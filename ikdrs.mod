TITLE Delayed rectifier potassium current in the interneuron soma

NEURON {
	SUFFIX ikdrs
	USEION k READ ek WRITE ik
	RANGE gkdr, gk, ik
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
	gkdr = 0.027 (mho/cm2)
}

ASSIGNED { 
	ik	(mA/cm2)
	v	(mV)
	ek	(mV)
	gk	(mho/cm2)
}

STATE { n c }

BREAKPOINT {
	SOLVE kstate METHOD sparse
	gk = gkdr*n*n*n*n
	ik = gk*(v-ek)
}

INITIAL {
	n = n_inf(v)
	c = 1-n
	gk = gkdr*n*n*n*n
	ik = gk*(v-ek)
}

LOCAL a1, a2

KINETIC kstate {
	a1 = a_n(v)
	a2 = a_c(v)
	~ c <-> n	(a1, a2)
	CONSERVE n + c = 1
}

FUNCTION a_n(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	a_n = n_inf(v)/tau_act(v)
}

FUNCTION a_c(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	a_c = (1-n_inf(v))/tau_act(v)
}

FUNCTION n_inf(v(mV)) {
             n_inf = 1/(1+exp((-4.8-v-37)/13.6))
}

FUNCTION tau_act(v(mV)) {
	tau_act = 1.6/(0.0338338*exp(-v/40)+(0.016*exp(v/5)*(64.9+v))/(-0.00000230599+exp(v/5)))
}
