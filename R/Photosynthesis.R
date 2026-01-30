#Implementation of the FvCB model described in Farquhar 1980, and adapted in Yin 2009 and Yin & Struik 2017

# Functions used to calculate effect of temperature on parameters
peaked25 = function(p25, Ep, Sp, Hdp, Tleaf) {
  return (p25*exp((Tleaf - 298.15)*Ep/(298.15*8.31*Tleaf))*(1 + exp((298.15*Sp - Hdp)/(8.31*298.15)))/(1 + exp((Tleaf*Sp - Hdp)/(8.31*Tleaf))))
}

arrhenius25 = function(p25, Ep, Tleaf) {
  return (p25*exp((Tleaf - 298.15)*Ep/(298.15*8.31*Tleaf)))
}

optimum25 = function(p25, sd, Tleaf) {
  T0 = 298.2	#reference temperature (K)	
  return (p25*exp(-(Tleaf-T0)**2/(2*sd**2)))
}

# Analytical solution of cubic equation due to coupling of A, gs and gm
CalcAnC3 = function(gm, gs0, fvpd, gb, x2, x1, gamma_star, Rd, Ca) {
  
  a = gs0*(x2 + gamma_star) + (gs0/gm + fvpd)*(x1 - Rd)
  b = Ca*(x1 - Rd) - gamma_star*x1 - Rd*x2
  c = Ca + x2 + (1/gm + 1/gb)*(x1 - Rd)
  d = x2 + gamma_star + (x1 - Rd)/gm
  m = 1/gm + (gs0/gm + fvpd)*(1/gm + 1/gb)
  p = -(d + (x1 - Rd)/gm + a*(1/gm + 1/gb) + (gs0/gm + fvpd)*c)/m
  q = (d*(x1 - Rd) + a*c + (gs0/gm + fvpd)*b)/m
  r = (-a*b)/m
  U = (2*p**3- 9*p*q + 27*r)/54.0
  Q = (p**2 - 3*q)/9.0
  psi = acos(U/sqrt(Q**3))
  A = -2*sqrt(Q)*cos(psi/3) - p/3.0
  
  return (A)
}

#solution of cubic equation when gs is known
CalcAnC3_gs = function(gm, gs, gb, x2, x1, gamma_star, Rd, Ca) {
  p = 1/gs + 1/gm + 1/gb
  q = -p * (x1 - Rd) - (Ca + x2)
  r = (Ca - gamma_star) * x1 - (Ca + x2) * Rd
  A = (-q - sqrt(q^2 - 4*p*r)) / (2*p*r)

  return(A)
}

# Calculation of internal CO2 concentration
CalcCi = function(gs0, An, Ca, Ci_star, Rd, fvpd) {
  
  a = gs0
  b = An - gs0*Ca - gs0*Ci_star + (An + Rd)*fvpd
  c = -An*Ci_star + gs0*Ca*Ci_star - (An + Rd)*Ca*fvpd
  Ci = (-b + sqrt(b**2 - 4*a*c))/(2*a)
  return (Ci)
}

# Calculate net photosynthesis without water limitation
calcA = function(PPFD, Ca, TleafC, VP, O2, LN, gs=NA) {
  #(PPFD = umol/m2/s Ca = ppm TleafC = dC VP = Pa O2 = ppm LN = gN/m2leaf, gs = mol/m2/s
  with(as.list(c(constants, parametersC3)),{
    
    # Calculate values of parameters at leafN
    k2ll = k2ll_a*LN + k2ll_b
    Vcmax25 = Vcmax25_a*LN + Vcmax25_b
    Jmax25 = Jmax25_a*LN + Jmax25_b
    TPU25 = TPU25_a*LN + TPU25_b
    
    # Need Kelvin for temperature responses
    Tleaf = TleafC + 273.15
    
    # Calculate values of parameters at Tleaf
    Rd = arrhenius25(Rd25, E_Rd, Tleaf)
    Sco = arrhenius25(Sco25, E_Sco, Tleaf)
    Vcmax = arrhenius25(Vcmax25, E_Vcmax, Tleaf)
    Kmc = arrhenius25(Kmc25, E_Kmc, Tleaf)
    Kmo = arrhenius25(Kmo25, E_Kmo, Tleaf)
    Jmax = peaked25(Jmax25, E_Jmax, S_Jmax, D_Jmax, Tleaf)
    TPU = peaked25(TPU25, E_TPU, S_TPU, D_TPU, Tleaf)
    gm = peaked25(gm25, E_gm, S_gm, D_gm, Tleaf)  
    
    #Parameters for the Pennman-Montheith equation
    #saturated vapour pressure as a function of temperature -> Tetens Equation (Pa/C)
    c1 = 610.78
    c2 = 17.269
    c3 = 237.3
    es = c1*c2*c3 * exp(c2*TleafC/(c3+TleafC))/(c3+TleafC)**2
    #saturation water vapour pressure (Pa)
    vps = c1*exp(c2*TleafC/(c3+TleafC))
    #Vapour Pressure Deficit (Pa)
    vpd = vps - VP
    #Latent heat of vaporisation of water in MJ/kg
    L = -0.0024*TleafC+2.5011
    #heat capacit of air in MJ/m3/k
    cp = (air_heat_capacity*1e-6)
    #phychometric constant (Pa K-1)
    y = cp*(Pair*1e3)/(L*0.622);
    
    # Parameters for CO2 diffusion
    fvpd = max(1/(1/(a1 - b1*(vpd/1000)) - 1),0.0)
    gamma_star = 0.5*O2/Sco
    Ci_star = gamma_star - Rd/gm
    
    # Limitation by Rubisco
    x1_c = Vcmax
    x2_c = Kmc*(1 + O2/Kmo)
    Ac = ifelse(gs==NA,
              CalcAnC3(gm, gs0, fvpd, gb, x2_c, x1_c, gamma_star, Rd, Ca),
              CalcAnC3_gs(gm, gs, gb, x2_c, x1_c, gamma_star, Rd, Ca))
    
    # Limitation by electron transport
    J = (k2ll*PPFD + Jmax - sqrt((k2ll*PPFD + Jmax)**2 - 4*theta*k2ll*Jmax*PPFD))/(2*theta)
    x1_j = (J/4.0)*(1 - fpseudo/(1 - fcyc))
    x2_j = 2*gamma_star
    Aj = ifelse(gs==NA,
              CalcAnC3(gm, gs0, fvpd, gb, x2_j, x1_j, gamma_star, Rd, Ca),
              CalcAnC3_gs(gm, gs, gb, x2_j, x1_j, gamma_star, Rd, Ca))
    
    # Limitation by TPU
    x1_p = 3.0*TPU
    x2_p = -gamma_star
    Ap = 3*TPU - Rd
    
    # Minimum of three potential rates
    An = min(Ac,min(Aj,Ap))
    
    # Transpiration
    Ag = An + Rd
    Ci = CalcCi(gs0, An, Ca, Ci_star, Rd, fvpd)
    Cc = Ci - An/gm
    gs = ifelse(gs==NA, gs0 + ((An + Rd)/(Ci - Ci_star))*fvpd, gs) #(mol/m2/s)
    gsw = gs*1.6 #stomatal conductance to water (mol/m2/s)
    gbw = gb*1.6 #boundary layer conductance conductance to water (mol/m2/s)
    
    # rbh = 100*sqrt(w/u) #boundary layer resistances to heat (rbh) and water (rbw) (in s m-1)
    # rbw = 0.93*rbh
    # rt = 0.74*(log((2-0.7*H)/(0.1*H))**2/(0.16*u))#Turbulent resistance (in s m-1)
    # rsw = (1/(1.6*gs))/mv               # stomatal resistance to water in s m-1 (/mv to convert from m2*s/mol)
    # PAR = PPFD/(1e6 * 0.55 * 4.55)      # in MJ/m2/s; 1e6: J to MJ, o.55: fraction PAR in Irradiance, 4.55: umol/J photons per J
    # Tr = (es*PAR + cp*vpd/(rbh+rt)) / (L*(es+y*(rbw+rsw+rt)/(rbh+rt))) #caluclate transpiration with Pennman-Monteith (1973)
    # WUE = An/Tr

    Tr = 1 / (1 / gsw + 1 / gbw) * vpd / (Pair*1e3) #Transpiration in (mol/m^2/s)
    Tr_L = Tr * 18.015 * 1e-3 # convert from mol to L
    
    return (list(An=An, Tr=Tr_L, gs=gs)) # if desired, we can also add Ag, gs, Ci, Ac, Aj, Ap
  })
}

##### Constants and parameters of the C3 reference paramerisation described in Yin 2009 and Yin & Struik 2017 ####
constants = c(
  pi = 3.14159265359,
  Pair = 101, # kPa
  emmisivity = 0.95,
  boltzmann = 5.67e-8, # J/m2/s/K4
  vap_heat_a = 56.897e3, # kJ/mol
  vap_heat_b = -0.043e3, # kJ/mol/K
  forcedA = 0.60,
  forcedn = 0.5,
  freeB = 0.5,
  freem = 0.25,
  diff = 20e-6, # m2/s
  visc = 15e-6, # m2/s
  Grcoef = 1.577e8, # 1/K/m3
  R = 8.31, # J/mol/K
  es0 = 0.61078, # kPa
  es_k = 17.269,
  Tzero = 273.15, # K
  es_Tref = 237.3, # K
  air_heat_capacity = 1200, # J/m3/K
  mv = 0.0245 #molair volume of air (m3/mol)
)

parametersC3 = c(
  # Rubisco CO2/O2 specificity
  Sco25 = 2800, # Sc/o parameter at 25 C
  E_Sco = -24460,	# Apparent activation energy of Sc/o (J/mol)
   # Michaelis-Menten constants carboxylation (J/mol/K)
  Kmc25 = 270, # Km for CO2 at 25 C (μmol/mol)
  E_Kmc = 80990,	# Activation energy of Kmc (J/mol)
  # Michaelis-Menten constants oxygenation
  Kmo25 = 165000, # Km for O2 at 25 C (umol/mol)
  E_Kmo = 23720,	# Activation energy of Kmo (J/mol)
  # Rubisco activity
  Vcmax25_a = 30.40,	# values for linear relation between Vcmax25 and leaf N (from Yin et al 2009 PCE)
  Vcmax25_b = 4.36,	#
  E_Vcmax = 65330,	# Activation energy of Vcmax (J/mol)
  # Electron transport
  theta = 0.7, # Curvature parameter
  k2ll_a = 0.044,	# values for linear relation between k2ll and leaf N (from Yin et al 2009 PCE)
  k2ll_b = 0.205,	#
  fpseudo = 0.1, # Fraction of electrons at PSI that are used by alternative electron sinks
  fcyc = 0.131, #  Fraction of electrons at PSI that follow cyclic transport around PSI
  Jmax25_a = 99.38,	# values for linear relation between Jmax25 and leaf N (from Yin et al 2009 PCE) Maximum rate of electron transport (μmol/m2/s)
  Jmax25_b = 5.75,	#
  E_Jmax = 30000,		# Activation energy Jmax (J/mol)
  D_Jmax = 200000,	# Deactivation energy of Jmax (J/mol)
  S_Jmax = 650,		# Entropy term for Jmax (J/K/mol)
  # Triose phosphate utilisation
  TPU25_a = 5.367,	# values for linear relation between TPU25 and leaf N (from Yin et al 2009 PCE) # Maximum rate of triose phosphate utilisation (μmol/m2/s)
  TPU25_b = 0.93,	#
  E_TPU = 53100, # Activation energy TPU (J/mol)
  D_TPU = 20180, # Deactivation energy of TPU (J/mol)
  S_TPU = 650, #  Entropy term for TPU (J/K/mol)
  # Respiration
  Rd25 = 1.2, # Respiration rate at 25 C (μmol/m2/s)
  E_Rd = 46390,	# Activation energy of Rd (J/mol)
  # Mesophyll conductance
  gm25 = 0.4, # Mesophyll conductance at 25 C (mol/m2/s)
  E_gm = 49600,	# Activation energy of gm (J/mol)
  D_gm = 437400,	# Deactivation energy of gm (J/mol)
  S_gm = 1400,	# Entropy term for gm (J/K/mol)
  # unit conversion
  mCO2 = 44.0095,		#molar mass of CO2 (g/mol)
  # Stomatal conductance
  a1 = 0.85, # Empirical parameter in gs formula
  b1 = 0.14, # Empirical parameter in gs formula (1/Pa)
  gs0 = 0.01, # Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)
  gb = 0.5, # boundary layer conductance for CO2 (mol/m2/s)
  #Vegetation stuff for Pennman-Monteith
  w = 0.1, #leaf width in m
  u = 2, #wind velocity in m s-1
  H = 2, #canopy height in m
  rH = 65 #relative humidity (%)
)

