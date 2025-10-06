* ==========================================================
* CCS Network Model 
* ==========================================================

$eolCom #
$inlineCom (* *)

* ---------- 0) Global options ----------
option mip = cplex, optcr = 0.02, limrow = 1000;

* ---------- 1) Sets ----------
Set
    t       "years"
    n       "regions"
    p       "coordinates (Latitude/Longitude)"
    c       "capture technologies" / coalPostCombustion, gasPostCombustion, coalOxyCombustion, preCombustion /
    q       "flowrate discretisation" / q1*q4 /
    ;
Alias (n,nn);   # needed for arcs (n,nn)

* ---------- 2) Parameters (loaded from GDX) ----------
* Core data
Parameter
    cn(n,nn)      "1 if regions n and nn are connected; 0 otherwise [-]"
    coastal_n(n)  "Coastal regions indicator [-]"
    mar_n(n)      "Maritime regions indicator [-]"
    inland_n(n)   "Inland regions indicator [-]"
    E(n,t)        "Emissions [kt/y]"
    coord(n,p)    "Coordinates [deg]"
    S(n)          "Storage capacity/inject limit [kt/y]"
    depth(n)      "Storage depth [m]"
    n_size(n)     "Region size [km]"
    d(n,nn)       "Distance between regions [km]"
    ;

* Cost inputs
Parameter
    aco(n,nn)       "Offshore pipeline factor (1 onshore; >1 offshore) [-]"
    USC_offshore(n) "Offshore sequestration multiplier (2.5 offshore, 1 onshore) [-]"
    UCC(c)          "Unitary Capture Cost [€/t]"
    UTC(q)          "Unit pipeline transport cost per distance [€/t·km]";
Scalar
    USC             "Unitary Sequestration Cost [€/t]";

* Capture feasibility & efficiency (from GDX)
Parameter
    cf(n,c)         "Tech feasibility share in region n for tech c [-]"
    ce(c)           "Capture efficiency of technology c in [0,1] [-]";
Scalar
    alpha           "Overall CO2 reduction target share in [0,1] [-]";

* ---------- Load from GDX ----------
$gdxin "C:\Users\josea\OneDrive - UNIVERSIDAD ALICANTE\Escritorio\OneDrive - UNIVERSIDAD ALICANTE\Doctorado\Tesis\Modelo\ALR_JACOBI_PROX\Matlab Files_ Jacobi Prox\Model_Data.gdx"
$load t, n, p
$load cn, coastal_n, mar_n, inland_n
$load E, coord, S, depth, n_size, d
$load aco, USC_offshore, UCC, UTC, USC
$load ce, cf
$gdxin

* Default target (kept here so it's easily tweakable)
alpha = 0.5;

* ---------- 3) Flowrate discretisation (kt/y) ----------
Parameter
    Qmin(q)  "minimum throughput when level q is used [kt/y]"
    Qmax(q)  "maximum throughput when level q is used [kt/y]";
Qmin('q1') =     1 ;   Qmax('q1') =   1500 ;
Qmin('q2') =  1500 ;   Qmax('q2') =   2500 ;
Qmin('q3') =  2500 ;   Qmax('q3') =   7500 ;
Qmin('q4') =  7500 ;   Qmax('q4') = 100000 ;

* ---------- 4) Feasible arcs for the model (exclude self-loops) ----------
Set A(n,nn) "feasible transport arcs";
A(n,nn) = yes$((cn(n,nn) > 0.5) and (ord(n) <> ord(nn)));

* ---------- 5) Decision variables ----------
Positive Variable
    pCO2(c,n,t)        "CO2 processed [kt/y]"
    cCO2(c,n,t)        "CO2 captured  [kt/y]"
    qCO2(q,n,nn,t)     "CO2 transported at level q on arc (n,nn) [kt/y]"
    Seq(n,t)           "CO2 sequestered [kt/y]";
Binary Variable
    y(q,n,nn,t)        "1 if flowrate level q is active on arc (n,nn) at time t";

* ---------- 6) Transport: linking constraints & level selection ----------
Equation
    flow_min(q,n,nn,t)  "qCO2 >= y * Qmin"
    flow_max(q,n,nn,t)  "qCO2 <= y * Qmax"
    choose_one(n,nn,t)  "at most one discrete level per arc and time"
    anti_parallel(n,nn,t) "no simultaneous activation in both directions";
flow_min(q,n,nn,t)$A(n,nn)..  qCO2(q,n,nn,t) =g= y(q,n,nn,t) * Qmin(q);
flow_max(q,n,nn,t)$A(n,nn)..  qCO2(q,n,nn,t) =l= y(q,n,nn,t) * Qmax(q);
choose_one(n,nn,t)$A(n,nn)..  sum(q, y(q,n,nn,t)) =l= 1;

* NEW: Anti-parallel activation (at most one direction per unordered pair & time)
anti_parallel(n,nn,t)$((ord(n) < ord(nn)) and (A(n,nn) or A(nn,n)))..
    sum(q, y(q,n,nn,t)) + sum(q, y(q,nn,n,t)) =l= 1;

* ---------- 7) Capture feasibility & effectiveness ----------
Equation
    reg_cap_limit(n,t)     "sum_c pCO2(c,n,t) <= E(n,t)"
    tech_feas_limit(c,n,t) "pCO2(c,n,t) <= cf(n,c)*E(n,t)"
    global_target          "sum_{c,n,t} pCO2 >= alpha * sum_{n,t} E"
    capture_eff(c,n,t)     "cCO2 = ce * pCO2";
reg_cap_limit(n,t)..       sum(c, pCO2(c,n,t)) =l= E(n,t);
tech_feas_limit(c,n,t)..   pCO2(c,n,t) =l= cf(n,c) * E(n,t);
global_target..            sum((c,n,t), pCO2(c,n,t)) =g= alpha * sum((n,t), E(n,t));
capture_eff(c,n,t)..       cCO2(c,n,t) =e= ce(c) * pCO2(c,n,t);

* ---------- 8) Mass balances ----------
Equation mass_balance(n,t) "inflow + local capture = outflow + sequestration";
mass_balance(n,t)..
      (  sum((q,nn)$A(nn,n), qCO2(q,nn,n,t))   # inflow to n
       + sum(c,             cCO2(c,n,t)) )     # capture at n
   =e= ( sum((q,nn)$A(n,nn), qCO2(q,n,nn,t))   # outflow from n
       + Seq(n,t) );                           # sequestered at n

* ---------- 9) Sequestration capacity ----------
Equation seq_cap(n,t) "injection/sequestration capacity limit";
seq_cap(n,t)..  Seq(n,t) =l= S(n);

* ---------- 10) Cost model ----------
* UNITS (keep consistent):
*   Flows pCO2, cCO2, Seq, qCO2 : [kt/y]
*   Distances d(n,nn)           : [km]
*   UCC(c)                      : [€/t]
*   UTC(q)                      : [€/t·km]  (pipeline onshore; offshore via aco(n,nn))
*   USC                         : [€/t]
*   aco(n,nn), USC_offshore(n)  : [-] multipliers
Scalar kt2t "conversion kt to t"; kt2t = 1e3;

Variable
    TCC(t)   "Total Capture Cost at time t [€/y]"
    TTC(t)   "Total Transport Cost at time t [€/y]"
    TSC(t)   "Total Sequestration Cost at time t [€/y]"
    TC(t)    "Total Cost at time t [€/y]"
    Z        "Objective: sum_t TC(t) [€]";
Equation
    def_TCC(t)  "TCC(t) = sum_{n,c} UCC(c) * cCO2(c,n,t)"
    def_TTC(t)  "TTC(t) = sum_{q,n,nn} qCO2 * d * UTC(q) * aco"
    def_TSC(t)  "TSC(t) = sum_{n} Seq(n,t) * USC * USC_offshore(n)"
    def_TC(t)   "TC(t)  = TCC(t) + TTC(t) + TSC(t)"
    obj         "Minimize Z = sum_t TC(t)";
def_TCC(t)..  TCC(t) =e= kt2t * sum((n,c), UCC(c) * cCO2(c,n,t));
def_TTC(t)..  TTC(t) =e= kt2t * sum((q,n,nn)$A(n,nn), qCO2(q,n,nn,t) * d(n,nn) * UTC(q) * aco(n,nn));
def_TSC(t)..  TSC(t) =e= kt2t * sum(n, Seq(n,t) * USC * USC_offshore(n));
def_TC(t)..   TC(t) =e= TCC(t) + TTC(t) + TSC(t);
obj..         Z =e= sum(t, TC(t));

* ---------- 11) Housekeeping for infeasible arcs ----------
qCO2.up(q,n,nn,t)$(not A(n,nn)) = 0;   # forbid flow on non-feasible arcs
y.fx(q,n,nn,t)$(not A(n,nn))     = 0;   # and forbid activation

* ---------- 12) Solve ----------
Model CCS / all /;
solve CCS using MIP minimizing Z;

* ---------- 13) Quick report ----------
display Z.l, TCC.l, TTC.l, TSC.l, TC.l, UCC, UTC, aco, USC, USC_offshore, ce, cf;

* ---------- Export results to GDX ----------
* We export levels (.l) of vars to keep post-processing simple
Parameter
    pCO2_l(c,n,t)    "CO2 processed [kt/y] (level)"
    cCO2_l(c,n,t)    "CO2 captured [kt/y] (level)"
    qCO2_l(q,n,nn,t) "CO2 transported [kt/y] (level)"
    Seq_l(n,t)       "CO2 sequestered [kt/y] (level)"
    y_l(q,n,nn,t)    "Transport activation (binary level)"
    TCC_l(t)         "Total capture cost [€/y]"
    TTC_l(t)         "Total transport cost [€/y]"
    TSC_l(t)         "Total sequestration cost [€/y]"
    TC_l(t)          "Total cost [€/y]"
    Z_l              "Objective value [€]"
    ;

pCO2_l(c,n,t)    = pCO2.l(c,n,t);
cCO2_l(c,n,t)    = cCO2.l(c,n,t);
qCO2_l(q,n,nn,t) = qCO2.l(q,n,nn,t);
Seq_l(n,t)       = Seq.l(n,t);
y_l(q,n,nn,t)    = y.l(q,n,nn,t);
TCC_l(t)         = TCC.l(t);
TTC_l(t)         = TTC.l(t);
TSC_l(t)         = TSC.l(t);
TC_l(t)          = TC.l(t);
Z_l              = Z.l;

execute_unload 'Results.gdx'
    n_size, inland_n, coastal_n, mar_n 
    t, n, p, c, q
    E, d, coord, aco, USC_offshore, UCC, UTC, USC
    pCO2_l, cCO2_l, qCO2_l, Seq_l, y_l,
    TCC_l, TTC_l, TSC_l, TC_l, Z_l
;
