import math
'''
Design/Adaptation of engineA aka


This function initialises this rocket motor to be used in a delta-v calculator
'''

#   P H Y S I C A L   C O N S T A N T S
TZERO = 288.15  # sea-level temperature, kelvins
PZERO = 101325.0  # sea-level pressure, N/sq.m
RHOZERO = 1.225  # sea-level density, kg/cu.m
g0 = 9.81


''' 

Oxidiser: Potassium nitrate 65%
Fuel: Dextrose (sugar) 35%
Doped with red iron oxide ( a further 1%)
'''
def massFlowRate (R, dt):
    # P R O P E L L A N T  C O N S T A N T S
    # R=8314 # [J/kmol-K] Universal Gas Constant
    # M = 42.39 # [kg/kmol] Effective molecular weight of products
    # R_s = R/M # [J/kg-K] Specific gas constant
    # k = 1.0435  # ratio of specific heats, 2-phase
    # comb_eff = 0.95 # combustion efficiency
    # T_0act = 1625 # [K] Actual chamber temperature
    # c_star = 889 #[m/s] characteristic exhaust velocity
    # grain_den_act = 1785 # [kg/m^3] actual grain density

    # # E N G I N E  C O N S T A N T S
    # p_0 = 6.89e6  # Chamber Pressure - taken as 68 atmospheres
    #
    # a_2 = (0.01 ** 2) * (pi * 0.25)  # Throat Diameter
    # a_3 = (0.02 ** 2) * (pi * 0.25)  # Exit Diameter
    # p_3 = 72326  # Pascals # Exit Pressure



    # BATES grain
    # N = 1 # Number of BATES grains http://www.nakka-rocketry.net/design1.html
    # D =0.075 # initial outer diameter


    # Mass flow rate
    r = R.burnRate()[0]  # = 0.146*1000**0.342 # [cm/s] burn rate http://www.nakka-rocketry.net/bntest.html
    t = 0.5*(R.D_grain-R.d) # web thickness of motor
    Ab = R.N_grain*(0.5*math.pi*(R.D_grain**2-R.d**2)+math.pi*R.L_grain*R.d) # instantaneous grain burning surface area

    mass_flow = Ab*(r*dt/100)*R.grain_density
    R.d += 2*(r*dt/100)
    R.L_grain -= 2*(r*dt/100)
    if R.d >= R.D_grain or R.L_grain <=0:
        return False
    else:
        return (mass_flow,R.d,R.L_grain)


