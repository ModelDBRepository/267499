TITLE calcium pump

NEURON {
	SUFFIX capump
	USEION ca READ cai WRITE ica	
	RANGE ica
	GLOBAL Vmax, Km, hill, Vconv, carest
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	FARADAY = (faraday) (coulombs)
}

PARAMETER {
	Vmax = 352 (uM/s)          :Units of Borgdorff
	vol_surf_ratio = 3.75 (um) :diam/4
	Km = .0069 (mM)	           :6.9 uM
	hill = 1	               :hill is 1.1, no significant diff from 1
	scale = 1e-4
}

ASSIGNED {
	ica (mA/cm2)
    cai (mM) 
	Vconv
	carest
}

INITIAL {
	VERBATIM
	cai = _ion_cai;
	carest = _ion_cai;
	ENDVERBATIM
	ica = pumprate(cai)*scale	
}

BREAKPOINT {
	ica = pumprate(cai)*scale
}

FUNCTION pumprate (ci) {
	Vconv = Vmax*vol_surf_ratio*FARADAY*2*(1e-4)
	if (fabs(ci-carest) < 1e-7) {
	  pumprate = (ci-carest)*Vconv/Km
	}
	else {
	  pumprate = Vconv/(1+Km/(ci-carest))
	}
}
