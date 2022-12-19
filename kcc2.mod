TITLE K-Cl cotransporter KCC2

INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}

NEURON {
	SUFFIX kcc2
	USEION k READ ko, ki WRITE ik
	USEION cl READ clo, cli WRITE icl VALENCE -1
	RANGE ik, icl
}

UNITS {
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	(mol) = (1)
}

PARAMETER {
  	U = 0.0020939418    (mA/cm2)
}

ASSIGNED {
	ik	(mA/cm2)
	icl	(mA/cm2)
	ko	(mM)
	ki	(mM)
  	clo     (mM)
  	cli     (mM)
}

BREAKPOINT {
	LOCAL rate
	rate = U*log((ki*cli)/(ko*clo))	
	ik =  rate
	icl = -rate
}
