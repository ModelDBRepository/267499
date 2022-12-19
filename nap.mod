TITLE Persistent sodium current

NEURON {
        SUFFIX nap
        USEION na READ ena WRITE ina
        RANGE  gna, ina
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
    (molar) = (1/liter)
	(mM) = (millimolar)
	PI   = (pi) (1)
	FARADAY	= 96485.309 (coul/mole)
	R = (k-mole) (joule/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {
	gna = 0.0006 (mho/cm2)
}
 
ASSIGNED {
	ena	(mV)
    v 	(mV)
    ina	(mA/cm2)
    minf
	hinf	
	taum	(ms)
	tauh	(ms)
}

STATE { m h }
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        ina = gna*m*m*h*(v-ena)
}
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
    ina = gna*m*m*h*(v-ena)
}

DERIVATIVE state { 
    rates(v)
    m' = (minf-m)/taum
    h' = (hinf-h)/tauh
}

PROCEDURE rates(v (mV)) {  
	TABLE minf, hinf, taum, tauh
	FROM -200 TO 200 WITH 201

	minf = 1/(1+exp(-(v+48.7)/4.4))

	hinf = 1/(1+exp((v+48.8)/9.98))

	taum = 1/(0.091*(v+38)/(1-exp(-(v+38)/5))-0.062*(v+38)/(1-exp((v+38)/5)))

	:The following function has been fitted to the data in Magistretti et al., 1999
	
	if (v<=-60) {
		tauh = 3700+2000*1/(0.091*(v+22+38)/(1-exp(-(v+22+38)/5))-0.062*(v+22+38)/(1-exp((v+22+38)/5)))
	} else {
		tauh = 1200+8000*1/(0.091*(v+36+38)/(1-exp(-(v+36+38)/5))-0.062*(v+36+38)/(1-exp((v+36+38)/5)))
	}
}
