from ThrustCalc import atmosphere
import math

def inst_params(R):
	##########################################################################################
										#Instantaneous Parameters
	##########################################################################################
	''' This function has an input of rocket altitude and returns the instantaneous thrust force (in Newtons) based on
	momentum thrust and pressure thrust from given engine parameters.
	F= mdot*v2 + (p2-p3)*a2
	Where
	mdot = mass flow rate
	v_2 = exhaust velocity
	p2 = Exit pressure
	p3 = atmospheric pressure

	'''


	#   P H Y S I C A L   C O N S T A N T S
	TZERO = 288.15  # sea-level temperature, kelvins
	PZERO = 101325.0  # sea-level pressure, N/sq.m
	RHOZERO = 1.225  # sea-level density, kg/cu.m
	k = R.k_2_phase # ratio of specific heats, 2-phase

	# E N G I N E  C O N S T A N T S
	a_2 = R.A_throat # Throat Diameter
	a_3 = R.A_exit # Exit Diameter
	p_3 = 101325 #Pascals # Exit Pressure - ideally expanded

	##############################################################
	#                      Chamber Pressure                      #
	##############################################################
	Ab = R.N_grain * (0.5 * math.pi * (R.D_grain** 2 - R.d** 2) + math.pi * R.L_grain * R.d)  # instantaneous grain burning surface area

	K_n = Ab/a_2 # ratio of burning surfaces to throat area
	r,a,n = R.burnRate()
	alpha = 1000000**n
	R.p_0 = (K_n*(a/alpha)*(R.grain_density/1000)*R.c_star)**(1/(1-n))

	# V A R I A B L E S
	Atmo = atmosphere(R.z/1000) # calculate atmospheric conditions
	p_a = Atmo[1]*PZERO # atmospheric pressure
	density_value = Atmo[0] * RHOZERO

	# S I M P L I F Y I N G   C A L C S
	aye = (2*k**2/(k-1))
	bee = (2/(k+1))**((k+1)/(k-1))
	cee = 1 - (p_3/R.p_0)**((k-1)/k)

	Total_Thrust = a_2*R.p_0*math.sqrt(aye*bee*cee) + (p_3-p_a)*a_3

	return (Total_Thrust, p_a, density_value)

def gravity(altitude, mass):
	G = 6.67408e-11
	M = 5.9722e24
	d = 6371000 + altitude

	F_g = (G*M*mass)/d**2

	# #Ellipsoid heights
	# Latitude =
	# a = 6378.137 # [km] Equatorial Radius
	# b = 6356.7523 # [km] Polar Radius
	return (F_g)


def dragCoeff(R):
	'''This Function models the drag coefficient of a rocket. It accounts for Body Drag, Fin Drag, Protuberance Drag,
	 Excrescencies Drag, Interference Drag, Base Drag, Transonic Wave Drag and Supersonic Wave Drag.
	 This approach is based off the DATCOM method.
	 '''
	# G A S   C O N S T A N T S
	gamma = 1.4 # Adiabatic constant of air
	R_const = 287.058 # Specific gas constant of air

	# A T M O S P H E R I C   C O N D I T I O N S
	Atmo = atmosphere(R.z/1000)
	TZERO = 288.15  # sea-level temperature, kelvins
	RHOZERO = 1.225  # sea-level density, kg/cu.m
	density = Atmo[0]*RHOZERO
	T = Atmo[2]*TZERO # Temperature

	# Constants for Sutherland's Forumla
	b = 1.458e-6
	S = 110.4

	Visc_dyn = (b*T**1.5)/(T+S)
	Visc_kin = Visc_dyn/density

	a = math.sqrt(gamma*R_const*T) # Local speed of sound
	M = abs(R.v/a) # Mach Number

	a_imperial = a*3.280839895 #1/0.3048
	Visc_kin_imperial = Visc_kin*0.09290304 # 0.3048^2


	L = R.L*39.37007874 # [inches] Length of rocket in inches
	D = R.D_outer*39.37007874 #[inches]  Maximum body diameter
	r = D/2
	L_cone = R.L_nc*39.37007874 # [inches] Length of Nose Cone
	S_b = math.pi*D*(L-L_cone)+ math.pi*r*(r+math.sqrt(L_cone**2+r**2)) # wetted surface of body and cone


	K = 0.00025
	'''	K = 0.0, for smooth surface
	= 0.00002 to 0.00008, for polished metal or wood
	= 0.00016, for natural sheet metal
	= 0.00025, for smooth matte paint, carefully applied
	= 0.0004 to 0.0012, for standard camouflage paint'''

	##############################################################
	#                          Body Drag                         #
	##############################################################

	#Compressible Renold's number
	R_comp = (a_imperial*M*L/(12*Visc_kin_imperial))*(1+0.0283*M-0.043*M**2 + 0.2107*M**3 - 0.03829*M**4 + 0.002709*M**5)


	cf_star = 0.037036*R_comp**-0.155079 # Incompressible skin friction coefficient
	cf = cf_star*(1+ 0.00798*M - 0.1813*M**2 + 0.0632*M**3 - 0.00933*M**4 + 0.000549*M**5 ) # Compressible skin friction coefficient

	cf_star_term = 1/((1.89+1.62*math.log10(L/K))**2.5)
	cf_term = cf_star_term/(1+0.02044*M**2)

	if cf >= cf_term:
		cf_final = cf
	else:
		cf_final = cf_term

	Cd_body = cf_final*(1 + (60/(L/D)**3) + 0.0025*(L/D))*(4*S_b/(math.pi*D**2))

	##############################################################
	#                          Fin Drag                          #
	##############################################################

	C_r = R.rootchord_fin*39.37007874  # [inches] root chord length of the fins
	C_t = R.tipchord_fin*39.37007874
	t = R.thickness_fin*39.37007874
	S_f = (R.span_fin*39.37007874/2)*(C_r + C_t)  # [inches^2] wetted area of each fin
	R_comp_fin = (a_imperial*M*C_r/(12*Visc_kin_imperial))*(1+0.0283*M-0.043*M**2 + 0.2107*M**3 - 0.03829*M**4 + 0.002709*M**5)

	cf_star_f = 0.037036 * R_comp_fin ** -0.155079  # Incompressible skin friction coefficient
	cf_f = cf_star * (1+0.00798*M - 0.1813*M**2 + 0.0632*M**3 - 0.00933*M**4 + 0.000549*M**5)  # Compressible skin friction coefficient

	cf_star_term_f = 1 / ((1.89 + 1.62 * math.log10(C_r / K)) ** 2.5)
	cf_term_f = cf_star_term_f / (1 + 0.02044 * M ** 2)

	if cf_f >= cf_term_f:
		cf_final_f = cf_f
	else:
		cf_final_f = cf_term_f

	Rn = (a_imperial*M*C_r/(12*Visc_kin_imperial)) # Incompressible Reynolds number
	lam = R.tiproot_ratio_fin

	if lam == 0:
		cf_lam = cf_final_f*(1+(0.5646/math.log10(Rn)))
	else:
		cf_lam = cf_final_f*((math.log10(Rn)**2.6)/(lam**2 - 1))*(((lam**2)/(math.log10(Rn*lam)**2.6))-(1/((math.log10(Rn))**2.6))+
																0.5646*(((lam**2)/(((math.log10(Rn*lam))**3.6)))-((1)/(((math.log10(Rn))**3.6)))))

	Cd_fin = cf_lam * (1 + (60 * ((t/ C_r) ** 4)) + 0.8*(1+5*R.x_t_fin**2)*(t / C_r)) * (4 *R.N_fin*S_f / (math.pi * D ** 2))

	##############################################################
	#                      Protuberance Drag                     #
	##############################################################
	# IGNORED

	##############################################################
	#                      Excrescencies Drag                    #
	##############################################################
	# Drag due to scratches, gouges, joints slots and holes. Assumed to be equally distributed over the wetted
	# surface of the rocket

	if M < 0.78:
		K_e = 0.00038 #coefficient for drag increment
	elif M >= 0.78 and M <= 1.04:
		K_e = -0.4501*M **4 +1.5954*M**3 - 2.1062*M** 2 +1.2288*M - 0.26717
	elif M > 1.04:
		K_e = 0.0002*M **2 - 0.0012*M + 0.0018
	else: return(False)

	Cd_E = K_e*(4*(S_b + S_f)/(math.pi*D**2))

	##############################################################
	#             Total Friction & Interference Drag             #
	##############################################################
	K_f = 1.04 # interference factor of fins and protuberances
	Cd_friction = Cd_body + K_f*Cd_fin + Cd_E

	##############################################################
	#                           Base Drag                        #
	##############################################################
	L_0 = L-L_cone
	K_b = 0.0274*math.atan((L_0/D)+ 0.0116)
	n = 3.6542*(L_0/D)**-0.2733

	D_b = D  # max diameter of rocket

	Cd_base = K_b * (((D_b / D) ** n) / (math.sqrt(Cd_friction)))
	f_b = 1
	if M > 0.6 and M< 1:
		f_b = 1 + 215.8*(M-0.6)**6
	elif M > 1 and M < 2:
		f_b = 2.0881*(M-1)**3 - 3.7938*(M-1)**2 + 1.4618*(M-1) + 1.883917
	elif M > 2:
		f_b = 0.297*(M-2)*3- 0.7937*(M- 2)**2 - 0.1115*(M- 2) + 1.64006

	Cd_base = Cd_base*f_b

	##############################################################
	#                    Transonic Wave Drag                     #
	##############################################################
	L_fairing = 0.3 # Doesn't exist

	if D_b > D:
		L_e = L_fairing
	else:
		L_e = L

	M_d =  -0.0156*(L_cone/D)**2 + 0.136*(L_cone/D) + 0.6817  # Transonic drag divergence Mach Number

	if (L_cone/L_e) < 0.2:
		a_tran = 2.4
		b_tran = -1.05
	elif (L_cone/L_e) >= 0.2:
		a_tran = -321.94*(L_cone/L_e)**2 + 264.07*(L_cone/L_e) - 36.348
		b_tran = 19.634*(L_cone/L_e)**2 + 18.369*(L_cone/L_e) - 1.7434


	M_f = a_tran*(L_e/D)**b_tran + 1.0275# Final Mach Number of Transonic Region

	#Calculate the maximum rise over transonic region:
	c = 50.676*(L_cone/L)**2 - 51.734*(L_cone/L) + 15.642
	g = -2.2538*(L_cone/L)**2 - 1.3108*(L_cone/L) - 1.7344

	if (L_e/D) >= 6:
		dCd_max_tran = c*(L_e/D)**g
	elif (L_e/D) < 6:
		dCd_max_tran = c*6**g

	#Transonic drag rise for a given Mach Number:
	x = ((M-M_d)/(M_f-M_d))
	F = -8.3474*x**5 + 24.543*x**4 - 24.946*x**3 + 8.6321*x**2 +1.1195*x
	if M >= M_d and M <= M_f:
		dCd_tran = dCd_max_tran*F
	elif M < M_d or M > M_f:
		dCd_tran = 0

	##############################################################
	#                   Supersonic Wave Drag                     #
	##############################################################

	if M >= M_f:
		dCd_super = dCd_max_tran
	else:
		dCd_super = 0

	##############################################################
	#                          Total Drag                        #
	##############################################################
	Cd_total = Cd_friction + Cd_base + dCd_tran + dCd_super


	# print(Cd_total)
	return(Cd_total)

if __name__ == "__main__":
	dragCoeff(R)