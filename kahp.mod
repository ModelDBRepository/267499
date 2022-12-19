TITLE Calcium-dependent K (afterhyperpolarization) current

NEURON { 
	SUFFIX kahp
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gkahp, ik, m
}

UNITS { 
	(molar) = (1/liter)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(mM) = (millimolar)
    FARADAY	= 96485.309 (coul/mole)
	PI = (pi) (1) 
}

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER { 
	gkahp = 0.00005 	(mho/cm2) 
}
 
ASSIGNED { 
	ik      (mA/cm2) 
	alpha   (/ms)
	beta	(/ms)
	v	(mV)	
	ek 	(mV)
	cai	(mM)
}
 
STATE {	m }

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik = gkahp*m*(v-ek)
}

INITIAL { 
	rates(cai)
	m = alpha/(alpha+beta)
	ik = gkahp*m*(v-ek)
}
 
DERIVATIVE states { 
	rates(cai)
	m' = alpha*(1-m)-beta*m
}

PROCEDURE rates(chi (mM)) { 

if(cai<5e-5) {cai = 5e-5} :lower cai limit

if(cai<=5.5e-5) {
		alpha = (cai-5e-5)*1e-2*2e5
	} else {
		alpha = 1e-2
	}
	beta = 1e-2
}
