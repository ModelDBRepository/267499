TITLE Fast Ca-dependent K current

NEURON { 
	SUFFIX kc
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE  gkc, ik, alpha, beta
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
	gkc = 196 	(mho/cm2)
} 

ASSIGNED { 
	ik 	(mA/cm2) 
	alpha   (/ms)
    beta	(/ms)
	ek 	(mV)
	v       (mV)
	cai	(mM)
}
 
STATE {	m }

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik = gkc*min(cai/250(mM),1)*m*(v-ek)
}

INITIAL { 
	settables(v) 
	m = alpha/(alpha+beta)
	ik = gkc*min(cai/250(mM),1)*m*(v-ek)
}
 
DERIVATIVE states { 
	settables(v) 
	m' = alpha*(1-m)-beta*m
}

PROCEDURE settables(v (mV)) { 
	TABLE alpha, beta FROM -120 TO 40 WITH 641

	if(v<=-10.0) {
		alpha = 1/18.975*(exp((v-10+60)/11-(v-6.5+60)/27))
		beta  = 2*exp((-v+6.5-60)/27)
	} else {
		alpha = 2*exp((-v+6.5-60)/27)
		beta  = 0
	}
}

FUNCTION min(x,y) {
	if(x<=y) {
		min = x
	} else {
		min = y
	}
}
