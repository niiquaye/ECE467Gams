*----------------------------------------------------------
* Notations Used
*----------------
*  A generator emissions function is represented as a quadratic polynomial
*  as follows: D*P**2 + E*P + F
*  D:     Quadratic component of emission function, in kg per MW**2-h
*  E:     Linear component of emission function, in kg per MWh
*  F:     Constant component of emission function, in kg per hour
*  PMax:  Generating unit maximum capacity in MW
*  PMin:  Generating unit minimum capacity in MW
*  PD:    Total system load demand in MW
*  SME:   System marginal emissions in kg per MWh
*----------------------------------------------------------


Option solprint = ON;
Option sysout   = ON;
Option limrow   =  6;
Option limcol   =  6;
Option decimals =  8 ;

Set
    g   'set of generating units' /g1*g4/
    k   'set of hours'  /h1*h24/;

Scalar
    Base    'base MVA'                  /100/;
*    Load    'Total System Load (MW)'    /800/;

**********  DATA SECTION ****************************************
Table gendata(g,*) 'generator cost characteristics and limits'
        'Pmin'  'PMax'      'D'         'E'         'F'
*       MW      MW        kg/h*(MW*MW)  kg/MWh      kg/h
g1      150     600         0.01265     1.3552      22.983
g2      100     500         0.01378     1.2489      137.37
g3      50      300         0.00767     0.8051      363.705
g4      50      300         0.0905      0.7560      198.5;

Parameter D(k)
/
h1  800
h2  900
h3  1000
h4  1100
h5  1250
h6  1100
h7  1150
h8  1200
h9  1350
h10 1450
h11 1500
h12 1450
h13 1400
h14 1100
h15 1200
h16 1200
h17 1250
h18 1300
h19 1400
h20 1450
h21 1200
h22 1100
h23 950
h24 750
/;

*------- Convert generator data to per unit quantities ----------
Parameter PMx(g,k), PMn(g,k), Dc(g), Ec(g), Fc(g), PD(k);

PMx(g,k) = gendata(g,"PMax")/Base;
PMn(g,k) = gendata(g,"PMin")/Base;

Dc(g)  = gendata(g,"D")*Base*Base;
Ec(g)  = gendata(g,"E")*Base;
Fc(g)  = gendata(g,"F");

PD(k)  = D(k)/Base;
********** DATA INPUT SECTION ENDS HERE ******************************


********* MODEL DEFINITION, SOLVE AND OUTPUT RESULTS SECTION *********
Variables
Emit                 'Total system emissions in kg'
Pg(g,k)              'Active power generated by unit g at time t (p.u)'
Loss(k)              'Losses in MW/MW';

Equations
EmitFn               'Total system emissions in kg'
LoadEq(k)            'Demand supply balance'
GenUp(g,k)           'Generation upper limit'
GenLo(g,k)           'Generation lower limit'
LossFn(k)            'Loss per hour';

EmitFn..   Emit =e= sum((k), sum((g),((Pg(g,k)*Pg(g,k)*Dc(g))+(Pg(g,k)*Ec(g))+(Fc(g)))));

LossFn(k)..   Loss(k) =e= 0.00003*Pg("g1",k)*Pg("g1",k)+0.00009*Pg("g2",k)*Pg("g2",k)+0.00012*Pg("g3",k)*Pg("g3",k)+0.00007*Pg("g4",k)*Pg("g4",k);

LoadEq(k)..    sum(g, Pg(g,k)) - PD(k) - Loss(k) =e= 0;

GenUp(g,k)..   Pg(g,k) =l= PMx(g,k);
GenLo(g,k)..   Pg(g,k) =g= Pmn(g,k);

Model ELD
/
EmitFn
LoadEq
GenUp
GenLo
LossFn
/
;

solve ELD Minimizing Emit using NLP;

parameter P(g,k),dualSME(k),dualUP(g,k),dualLO(g,k);
P(g,k)        = Pg.l(g,k)*100;
dualSME(k)    = LoadEq.m(k)/100;
dualUP(g,k)   = GenUp.m(g,k)/100;
dualLO(g,k)   = GenLo.m(g,k)/100;

display P, Emit.l, Loss.l, dualSME, dualUP, dualLO