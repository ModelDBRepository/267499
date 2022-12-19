TITLE Leak current - passive electrical properties

NEURON {
	SUFFIX leak
	USEION k READ ek WRITE ik
	USEION na READ ena WRITE ina
	USEION cl READ ecl WRITE icl VALENCE -1
	RANGE gk, ik, gna, ina, gcl, icl
}

UNITS { 
	(mV) = (millivolt)
    (mA) = (milliamp)
    PI   = (pi) (1)
	FARADAY	= 96485.309 (coul/mole)
}

PARAMETER {
	gk	= 3e-5 (mho/cm2)	
	gna	= 1.5101316e-5 (mho/cm2)	
	gcl	= 1e-5 (mho/cm2)	
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	ek (mV)
	ina (mA/cm2)
	ena (mV)
	icl (mA/cm2)
	ecl (mV)
}

BREAKPOINT {
	ik = gk*(v-ek)
	ina = gna*(v-ena)
	icl = gcl*(v-ecl)
}

INITIAL {
	ik = gk*(v-ek)
	ina = gna*(v-ena)
	icl = gcl*(v-ecl)
}
