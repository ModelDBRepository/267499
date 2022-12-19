TITLE Transient sodium current in the pyram somas

NEURON {
	SUFFIX pnas
	USEION na READ ena WRITE ina
	RANGE gna, ina
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
        PI = (pi) (1)
	FARADAY	= 96485.309 (coul/mole)
}

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {
	gna = 0.014 	(mho/cm2)
}

ASSIGNED { 
	ina  (mA/cm2)
	v    (mV)
	ena  (mV)
}

STATE { ma mb ha hb }

BREAKPOINT {
	SOLVE nastate METHOD sparse
	ina = gna*ma*ma*ma*ha*(v-ena)
}

INITIAL {
	ma = m_inf(v)
	ha = h_inf(v)
	mb = 1-ma
	hb = 1-ha
	ina = gna*ma*ma*ma*ha*(v-ena)
}

LOCAL a1, a2, b1, b2

KINETIC nastate {	
	a1 = m_a(v)
	a2 = m_b(v)
	b1 = h_a(v)
	b2 = h_b(v)
	~ mb <-> ma (a1, a2)
	~ hb <-> ha (b1, b2)
	CONSERVE ma + mb = 1
	CONSERVE ha + hb = 1
}
	
FUNCTION m_a(v(mV)) {
	TABLE FROM -150 TO 150 WITH 301
    m_a = 0.8*(17.2-v-60+3)/(exp((17.2-v-60+3)/4)-1)
}

FUNCTION m_b(v(mV)) {
	TABLE FROM -150 TO 150 WITH 301
	m_b = 0.7*(v-3-42.2+60)/(exp((v-3-42.2+60)/5)-1)
}

FUNCTION h_a(v(mV)) {
	TABLE FROM -150 TO 150 WITH 301
	h_a = 0.32*exp((42-v-60+3)/18)
}

FUNCTION h_b(v(mV)) {
	TABLE FROM -150 TO 150 WITH 301
	h_b = 10/(1+exp((42-v-60+3)/5))
}

FUNCTION m_inf(v(mV)) {
	m_inf = m_a(v)/(m_a(v)+m_b(v))
}

FUNCTION h_inf(v(mV)) {
	h_inf = h_a(v)/(h_a(v)+h_b(v))
}
