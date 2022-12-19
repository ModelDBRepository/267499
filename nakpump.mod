TITLE Na-K Pump

INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}

NEURON {
	SUFFIX nakpump
	USEION k READ ko WRITE ik
	USEION na READ nai WRITE ina
	RANGE ik, ina, km_k, km_na, imax
}

UNITS {
	(mV)	= (millivolt)
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(mol)	= (1)
	FARADAY	= 96485.309 (coul/mole)
	PI	    = (pi) (1)
	R 	    = (k-mole)(joule/degC)
}

PARAMETER {
	km_k = 2	(mM) 
	km_na = 10	(mM)
	imax = 0.01 (mA/cm2)
}

ASSIGNED {
	ik		(mA/cm2)
	ina		(mA/cm2)
	ko		(mM)
	nai		(mM)
}

BREAKPOINT {
	ik = -2*imax*flux(nai,ko)
	ina = ik*-3/2
}

INITIAL {
	ik = -2*imax*flux(nai,ko)
	ina = ik*-3/2
}

FUNCTION flux(na,k) {
	flux = (1/((1+km_k/k)*(1+km_k/k)))*(1/((1+km_na/na)*(1+km_na/na)*(1+km_na/na)))
}
