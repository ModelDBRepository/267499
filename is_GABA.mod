COMMENT

Synaptic GABAergic mechanism 

Reversal potential Egaba is changing according to [Cl-]i change (due to Cl- influx, which we hypothesize to be significant). Bicarbonate (HCO3) flows through the GABAR too, and therefore Egaba is also [HCO3]i/[HCO3]o -dependent
igaba = icl + ihco3 (we assume icl and ihco3 to be mutually independent)

The GABAa model file is based on Jedlicka et al., 2010; https://doi.org/10.1002/hipo.20804
ENDCOMMENT

TITLE Inhibitory synapse in the pyramidal cells (GABA_A) with changing Cl- concentration

NEURON {	
	POINT_PROCESS is_GABA
	RANGE tau1, tau2, icl, ihco3, i, g, e, ehcl, ehco3, P
	USEION cl READ cli, clo, ecl WRITE icl VALENCE -1
	USEION hco3 READ hco3i, hco3o, ehco3 WRITE ihco3 VALENCE -1    		   
}

PARAMETER {
	tau1 = 2   (ms)
	tau2 = 6   (ms)
	P = 0.18  :HCO3/Cl relative permeability
}

UNITS	{
	FARADAY		= 96485.309 (coul/mole)
	R 		    = (k-mole) (joule/degC)
	(molar)		= (1/liter)
	(mM)		= (millimolar)
    (mV)        = (millivolt)
}

ASSIGNED {
	v    (mV)
	e    (mV)
	ecl    (mV)
	ehco3    (mV)
	i   (nanoamp)
	icl   (nanoamp)
	ihco3   (nanoamp)
  	g (microsiemens)
  	factor
	cli (mM)
	clo (mM)
	hco3i (mM)
	hco3o (mM)
}

STATE {
	A (microsiemens)
	B (microsiemens)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor	
	e = P*ehco3 + (1-P)*ecl
}

BREAKPOINT {
	SOLVE state METHOD cnexp	
	g = B - A
	icl = (1-P)*g*(v-ecl)
	ihco3 = P*g*(v-ehco3)
	i = icl + ihco3
	e = P*ehco3 + (1-P)*ecl	
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (microsiemens)) {
	A = A + weight*factor
	B = B + weight*factor
} 
