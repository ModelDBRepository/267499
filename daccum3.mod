TITLE daccum3

NEURON {
	SUFFIX daccum3

	USEION na READ ina, nao, nai WRITE nai, nao, ina
	USEION k  READ ik, ko, ki WRITE ki, ko, ik
	USEION ca READ ica, cai WRITE cai, cao, ica VALENCE 2
	USEION cl READ icl, cli, clo WRITE cli, clo, icl VALENCE -1
	USEION hco3  READ hco3o, hco3i VALENCE -1
	USEION a  READ ia, ao, ai WRITE ai, ao VALENCE -1
	
	RANGE volin, volout, delta
	POINTER ko2, ko4, nao2, nao4, kis, kos, nais, naos, clis, clos, cais, caos
	GLOBAL TotalBuffer, TotalBuffer_calcium
}

UNITS	{
	(mV)		= (millivolt)
	(mA)		= (milliamp)
	FARADAY		= 96485.309 (coul/mole)
	(molar)		= (1/liter)
	(mM)		= (millimolar)
	PI	     	= (pi) (1)
	R 	    	= (k-mole) (joule/degC)
}

PARAMETER {
	Difna = 1.33	(um2/ms)
	Difk  = 1.96	(um2/ms)
	Difca = 0.6	    (um2/ms)
	Difcl = 2.03	(um2/ms)
	Difa = 0	    (um2/ms)

	k1buf = 0.0008 (/ms)
	TotalBuffer = 1100 (mM)
	x0 = 16 (mM)   
	q = -1.25 (mM)
    bf = 1    

	TotalBuffer_calcium = 1.562 (mM)
	Kd_calcium = .008 (mM)	            

	nab = 140 (mM) 
	kb  = 3.5 (mM)
	clb = 135 (mM)
	
	s = 4.4e4
	sk = 4.4e4
    bath = 1
	nadif = 1
    rad = 1 
    narad = 1 

	setvolin = 1
	setvolout = 1

	minvol = 0.04
	maxvol = 0.9
	tau = 250 (ms)
	swell = 1         

	L1 = 20 (um)	 
	L2 = 450 (um)	 
	L12 (um)	 
	S1 (um2)	 
	S2 (um2)
	S12 (um2)
	T12 (um2)
}

ASSIGNED {
	ina		(mA/cm2)
	ik		(mA/cm2)
	ica		(mA/cm2)
	icl		(mA/cm2)
	ia		(mA/cm2)

	diam		(um)

	naflux[2]
	kflux[2]
	caflux[2]
	clflux[2]
	aflux[2]

	dif[3]

	shell[4]

	klat[2]
	nalat[2]
	cllat[2]
	calat[2]
	
	nai ki cai cli hco3i ai
	nao ko cao clo hco3o ao

	volin     	
	volout		

	Kd (/mM)
	B0 (mM)

	B0_calcium (mM)

	delta (mM)

	ko2 (mM)
	ko4 (mM)
	nao2 (mM)
	nao4 (mM)

	kis (mM)
	kos (mM)
	nais (mM)
	naos (mM)
	clis (mM)
	clos (mM)
	cais (mM)
	caos (mM)
}

STATE { na[2] k[2] ca[2] cl[2] a[2] kbuf Buffer KBuffer catot Buffer_calcium CaBuffer (mM) vol[2] <1e-4> }

LOCAL b, c, d

BREAKPOINT {
	SOLVE state METHOD sparse
	L12 = (L1+L2)/2
	S1 = 15*15*PI/4
	S2 = 6.88*6.88*PI/4
	S12 = (S1+S2)/2
	T12 = setvolout*(S1+S2)/2
}

INITIAL {
	na[0]=nai
	na[1]=nao
	k[0]=ki
	k[1]=ko
	ca[0]=cai
	ca[1]=cao
	cl[0]=cli
	cl[1]=clo
	a[0]=ai
	a[1]=ao
	
	Kd = k2buf(k[1])/k1buf
    kbuf = 0
	B0 = TotalBuffer/(1+Kd*k[1])
	Buffer = B0
	KBuffer = TotalBuffer - B0

	B0_calcium = TotalBuffer_calcium/(1+Kd_calcium*cai)
	Buffer_calcium = B0_calcium
	CaBuffer = TotalBuffer_calcium - B0_calcium
	catot = cai*(1+(TotalBuffer_calcium/(cai+Kd_calcium)))
	
	vol[0]=setvolin
	vol[1]=setvolout

	volin=vol[0]
	volout=vol[1]
}

KINETIC state {

	delta = ((nai + ki + cli + cai + hco3i + ai) - (nao + ko + clo + cao + hco3o + ao))/tau

	IF (vol[1] <= minvol && delta > 0 ) {
	  delta = 0
	} ELSE {
	  IF (vol[0] <= maxvol && delta < 0 ) {
	    delta = 0
	  }
	~ vol[0] << (swell*delta/(diam*diam*PI/4))   :intracellular
	~ vol[1] << (-swell*delta/(diam*diam*PI/4))  :extracellular
	}

	COMPARTMENT i, vol[i]*PI/4*diam*diam { na k ca cl a }
	COMPARTMENT vol[0]*diam*diam*PI/4 { catot }

	naflux[0] = -delta*na[0] -(ina*diam)*PI*(1e4)/FARADAY
    kflux[0]  = -delta*k[0] -(ik*diam)*PI*(1e4)/FARADAY
	caflux[0] = -delta*catot -(ica*diam)*PI*(1e4)/(FARADAY*2)
	clflux[0] = -delta*cl[0]-(icl*diam)*PI*(1e4)/(FARADAY*-1)
    aflux[0]  = -delta*a[0]
	
	naflux[1] = delta*na[1] + (ina*diam)*PI*(1e4)/FARADAY
    kflux[1]  = delta*k[1]  + (ik*diam)*PI*(1e4)/FARADAY
	caflux[1] = delta*ca[1] + (ica*diam)*PI*(1e4)/(FARADAY*2)
	clflux[1] = delta*cl[1] + (icl*diam)*PI*(1e4)/(FARADAY*-1)
	aflux[1]  = delta*a[1]

	dif[0]    = nadif*bath*Difna*((geom(diam)/2)/2)*((nab-na[1])/s)*(PI*(diam+dr(diam)))
	dif[1]    = bath*Difk*((geom(diam)/2)/2)*((kb-k[1]-kbuf)/sk)*(PI*(diam+dr(diam)))
	dif[2]    = bath*Difcl*((geom(diam)/2)/2)*((clb-cl[1])/s)*(PI*(diam+dr(diam)))

	shell[0]  = rad*Difk*(geom(diam)/2)*surf(diam)*(ko2-k[1]-kbuf)
	shell[1]  = rad*Difk*(geom(diam)/2)*surf(diam)*(ko4-k[1]-kbuf)
	shell[2]  = narad*Difna*(geom(diam)/2)*surf(diam)*(nao2-na[1])
	shell[3]  = narad*Difna*(geom(diam)/2)*surf(diam)*(nao4-na[1])

	klat[0]    = -Difk*(k[0]-kis)*S12/(L12*L2)
	klat[1]    = -Difk*(k[1]+kbuf-kos)*T12/(L12*L2)
	nalat[0]    = -Difna*(na[0]-nais)*S12/(L12*L2)
	nalat[1]    = -Difna*(na[1]-naos)*T12/(L12*L2)
	cllat[0]    = -Difcl*(cl[0]-clis)*S12/(L12*L2)
	cllat[1]    = -Difcl*(cl[1]-clos)*T12/(L12*L2)
	calat[0]    = -Difca*(ca[0]-cais)*S12/(L12*L2)
	calat[1]    = -Difca*(ca[1]-caos)*T12/(L12*L2)
		
	IF (bf == 0 ) {
	  kbuf = 0
          Buffer = 0
          KBuffer = 0
	} ELSE {
	~  kbuf  << (-k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer+k1buf*KBuffer)
  	~  Buffer << (-k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer+k1buf*KBuffer)
  	~  KBuffer << (k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer-k1buf*KBuffer)
        }

	b = TotalBuffer_calcium*setvolin/vol[0]-catot+Kd_calcium
	c = -Kd_calcium*catot
	d = b*b-4*c
	cai = (-b+sqrt(d))/(2)
	CaBuffer = catot-cai
	Buffer_calcium = TotalBuffer_calcium*setvolin/vol[0] - CaBuffer

	~  na[0] << (naflux[0]+nalat[0])
        ~  k[0]  << (kflux[0]+klat[0])
	~  catot << (caflux[0]+calat[0])
	~  cl[0] << (clflux[0]+cllat[0])
	~  a[0] << (aflux[0])
	~  na[1] << (naflux[1]+dif[0]+shell[1]+nalat[1])
	~  k[1]  << (kflux[1]+dif[1]+shell[0]+shell[1]+klat[1])
	~  ca[1] << (caflux[1]+calat[1])
	~  cl[1] << (clflux[1]+dif[2]+cllat[1])
	~  a[1] << (aflux[1])

	volin=vol[0]
	volout=vol[1]
	nai=na[0] nao=na[1]
	ki=k[0]   ko=k[1]+kbuf
	cao=ca[1]
	cli=cl[0] clo=cl[1]
	ai=a[0]   ao=a[1]
}

FUNCTION k2buf(x(mM)) { 	        
	TABLE FROM 0 TO 199 WITH 200
	k2buf = k1buf/(1+exp((x-x0)/(q)))
}

FUNCTION dr(x(um)) { 			
	TABLE FROM 0 TO 15 WITH 10
	dr = x*(sqrt(1+setvolout)-1)
}

FUNCTION geom(x(um)) { 		
	TABLE FROM 0 TO 15 WITH 10
	geom = 1/dr(x)
}

FUNCTION surf(x(um)) {
	TABLE FROM 0 TO 15 WITH 10
	surf = (PI-2*atan2(x+2*dr(x/2),(x/2)))*((x+2*dr(x/2))/sin(atan2(x+2*dr(x/2),(x/2))))
}
