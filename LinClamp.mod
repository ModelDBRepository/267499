NEURON {
        POINT_PROCESS LinClamp
        RANGE del, tf, amp0, ampf
        ELECTRODE_CURRENT i
}

UNITS {
        (nA) = (nanoamp)
}

PARAMETER {
	del  = 60000   (ms)
	tf   = 100000  (ms)
	amp0 = 0.35    (nA)
	ampf = 0       (nA)
}

ASSIGNED {
        i (nA)
}

BREAKPOINT {
       IF (t < del) {
         i = 0   
	} ELSE { 
        IF (t > del && t < tf) {
          i = ((ampf-amp0)/(tf-del))*(t-del)+amp0
	 }
        }
}
