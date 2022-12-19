TITLE saccum1

NEURON {
	SUFFIX saccum1

	USEION na READ ina, nao, nai WRITE nai, nao, ina
	USEION k  READ ik, ko, ki WRITE ki, ko, ik
	USEION ca READ ica, cai WRITE cai, cao, ica VALENCE 2
	USEION cl READ icl, cli, clo WRITE cli, clo, icl VALENCE -1
	USEION hco3  READ hco3o, hco3i 	VALENCE -1
	USEION a  READ ia, ao, ai WRITE ai, ao VALENCE -1
	
	RANGE volin, volout, delta
	POINTER ko2, nao2, kid, kod, naid, naod, clid, clod, caid, caod
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
	:Ion diffusion constants
	Difna = 1.33	(um2/ms)   :Somjen 2008
	Difk  = 1.96	(um2/ms)   :Somjen 2008
	Difca = 0.6	    (um2/ms)   :Somjen 2008
	Difcl = 2.03	(um2/ms)   :Somjen 2008
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

	:Initial volume fractions
	setvolin = 1	:Pyram soma 1
	setvolout = 1   :Extracellular space

	minvol = 0.04      :Minimum feasible extracellular volume fraction
	maxvol = 0.9       :Maximum feasible intracellular shrinkage
	tau = 250 (ms)     :Timescale for cell swelling
	swell = 1          :Cell swelling decisional factor

	L1 = 20 (um)	   :soma1 length
	L2 = 450 (um)	   :dend1 length
	L12 (um)	   :distance between soma1-dend1 centers
	S1 (um2)	   :soma1 cross-section
	S2 (um2)	   :dend1 cross-section
	S12 (um2)	   :average soma-dend cross-section
	T12 (um2)	   :extracellular average soma-dend cross-section
}

ASSIGNED {
	:Ion currents
	ina		(mA/cm2)
	ik		(mA/cm2)
	ica		(mA/cm2)
	icl		(mA/cm2)
	ia		(mA/cm2)

	:Geometry	
	diam		(um)

	:Ionic fluxes (mM/ms)		
	naflux[2]
	kflux[2]
	caflux[2]
	clflux[2]
	aflux[2]

	:Cell-bath diffusion contributions (mM/ms)	
	dif[3]

	:Extracellular spaces contributions (mM/ms)
	shell[2]

	:Lateral diffusion contributions (mM/ms)
	klat[2]
	nalat[2]
	cllat[2]
	calat[2]
	
	nai ki cai cli hco3i ai		:Intracellular concentrations
	nao ko cao clo hco3o ao		:Extracellular concentrations

	volin     	:Interneuron volume
	volout		:Extracellular volume

	:K buffering parameters
	Kd (/mM)
	B0 (mM)

	:Ca buffering
	B0_calcium (mM)

	:Intra- and extracellular bulk concentrations gradient
	delta (mM)

	:Extracellular potassium and sodium exchange
	ko2 (mM)
	nao2 (mM)

	:Lateral diffusion
	kid (mM)
	kod (mM)
	naid (mM)
	naod (mM)
	clid (mM)
	clod (mM)
	caid (mM)
	caod (mM)
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
	
	:K buffering (implemented as in Somjen 2000)
	Kd = k2buf(k[1])/k1buf
    kbuf = 0
	B0 = TotalBuffer/(1+Kd*k[1])
	Buffer = B0
	KBuffer = TotalBuffer - B0

	:Ca buffering (Somjen)
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

	:CELLS SWELLING AND SHRINKING (Kager 2008)

	delta = ((nai + ki + cli + cai + hco3i + ai) - (nao + ko + clo + cao + hco3o + ao))/tau	:delta = 0 <--> no volume changes

	IF (vol[1] <= minvol && delta > 0 ) {
	  delta = 0
	} ELSE {
	  IF (vol[0] <= maxvol && delta < 0 ) {
	    delta = 0
	  }
	~  vol[0] << (swell*delta/(diam*diam*PI/4)) :intracellular
	~  vol[1] << (-swell*delta/(diam*diam*PI/4)) :extracellular
	}

	COMPARTMENT i, vol[i]*PI/4*diam*diam { na k ca cl a }
	COMPARTMENT vol[0]*diam*diam*PI/4 { catot }

	:INTRACELLULAR FLUXES
	naflux[0] = -delta*na[0] -(ina*diam)*PI*(1e4)/FARADAY
    kflux[0]  = -delta*k[0] -(ik*diam)*PI*(1e4)/FARADAY
	caflux[0] = -delta*catot -(ica*diam)*PI*(1e4)/(FARADAY*2)
	clflux[0] = -delta*cl[0]-(icl*diam)*PI*(1e4)/(FARADAY*-1)
    aflux[0]  = -delta*a[0]
    
	
	:EXTRACELLULAR FLUXES
	naflux[1] = delta*na[1] + (ina*diam)*PI*(1e4)/FARADAY
    kflux[1]  = delta*k[1] +  (ik*diam)*PI*(1e4)/FARADAY
	caflux[1] = delta*ca[1] + (ica*diam)*PI*(1e4)/(FARADAY*2)
	clflux[1] = delta*cl[1] + (icl*diam)*PI*(1e4)/(FARADAY*-1)
	aflux[1]  = delta*a[1]

	:DIFFUSION TO THE BATH
	dif[0]    = nadif*bath*Difna*((geom(diam)/2)/2)*((nab-na[1])/s)*(PI*(diam+dr(diam)))
	dif[1]    = bath*Difk*((geom(diam)/2)/2)*((kb-k[1]-kbuf)/sk)*(PI*(diam+dr(diam)))
	dif[2]    = bath*Difcl*((geom(diam)/2)/2)*((clb-cl[1])/s)*(PI*(diam+dr(diam)))

	:RADIAL DIFFUSION
	shell[0]  = rad*Difk*(geom(diam)/2)*surf(diam)*(ko2-k[1]-kbuf)
	shell[1]  = narad*Difna*(geom(diam)/2)*surf(diam)*(nao2-na[1])

	:LATERAL DIFFUSION
	klat[0]    = Difk*(kid-k[0])*S12/(L12*L1)
	klat[1]    = Difk*(kod-k[1]-kbuf)*T12/(L12*L1)
	nalat[0]    = Difna*(naid-na[0])*S12/(L12*L1)
	nalat[1]    = Difna*(naod-na[1])*T12/(L12*L1)
	cllat[0]    = Difcl*(clid-cl[0])*S12/(L12*L1)
	cllat[1]    = Difcl*(clod-cl[1])*T12/(L12*L1)
	calat[0]    = Difca*(caid-ca[0])*S12/(L12*L1)
	calat[1]    = Difca*(caod-ca[1])*T12/(L12*L1)
	
	:POTASSIUM BUFFERING
	IF (bf == 0 ) {
	  kbuf = 0
      Buffer = 0
      KBuffer = 0
	} ELSE {
	~  kbuf  << (-k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer+k1buf*KBuffer)
  	~  Buffer << (-k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer+k1buf*KBuffer)
  	~  KBuffer << (k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer-k1buf*KBuffer)
        }

	:CALCIUM BUFFERING
	b = TotalBuffer_calcium*setvolin/vol[0]-catot+Kd_calcium
	c = -Kd_calcium*catot
	d = b*b-4*c
	cai = (-b+sqrt(d))/(2)
	CaBuffer = catot-cai
	Buffer_calcium = TotalBuffer_calcium*setvolin/vol[0] - CaBuffer

	:DIFFERENTIAL EQUATIONS
	~  na[0] << (naflux[0]+nalat[0])
    ~  k[0]  << (kflux[0]+klat[0])
	~  catot << (caflux[0]+calat[0])
	~  cl[0] << (clflux[0]+cllat[0])
	~  a[0] << (aflux[0])
	~  na[1] << (naflux[1]+dif[0]+shell[1]+nalat[1])
	~  k[1]  << (kflux[1]+dif[1]+shell[0]+klat[1])
	~  ca[1] << (caflux[1]+calat[1])
	~  cl[1] << (clflux[1]+dif[2]+cllat[1])
	~  a[1] << (aflux[1])

	:VARIABLE IDENTIFICATION
	volin=vol[0]
	volout=vol[1]
	nai=na[0] nao=na[1]
	ki=k[0]   ko=k[1]+kbuf
	cao=ca[1]
	cli=cl[0] clo=cl[1]
	ai=a[0]   ao=a[1]
}



:Functions definition

FUNCTION k2buf(x(mM)) { 	         :Forward K buffering rate - Somjen 2000
	TABLE FROM 0 TO 199 WITH 200
	k2buf = k1buf/(1+exp((x-x0)/(q)))
}

FUNCTION dr(x(um)) { 			 :Interstitial space thickness
	TABLE FROM 0 TO 15 WITH 10
	dr = x*(sqrt(1+setvolout)-1)
}

FUNCTION geom(x(um)) {
	TABLE FROM 0 TO 15 WITH 10
	geom = 1/dr(x)
}

FUNCTION surf(x(um)) { 			 :Contact surface between shells
	TABLE FROM 0 TO 15 WITH 10
	surf = (PI-2*atan2(x+2*dr(x/2),(x/2)))*((x+2*dr(x/2))/sin(atan2(x+2*dr(x/2),(x/2))))
}
