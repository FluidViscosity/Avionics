import math
import numpy as np
import matplotlib.pyplot as plt

'''version = "1.0 (2004 Oct 04)"
greeting = "Tables - A Python program to compute atmosphere tables"
author = "Ralph L. Carmichael, Public Domain Aeronautical Software"
modifier = "Tomer Simhony"

'''



#   P H Y S I C A L   C O N S T A N T S

FT2METERS = 0.3048  # mult. ft. to get meters (exact)
KELVIN2RANKINE = 1.8  # mult deg K to get deg R
PSF2NSM = 47.880258  # mult lb/sq.ft to get sq.m
SCF2KCM = 515.379  # mult slugs/cu.ft to get kg/cu.m
TZERO = 288.15  # sea-level temperature, kelvins
PZERO = 101325.0  # sea-level pressure, N/sq.m
RHOZERO = 1.225  # sea-level density, kg/cu.m
AZERO = 340.294  # speed of sound at S.L.  m/sec
BETAVISC = 1.458E-6  # viscosity constant
SUTH = 110.4  # Sutherland's constant, kelvins



def atmosphere(alt):
    """ Compute temperature, density, and pressure in standard atmosphere.
  Correct to 86 km.  Only approximate thereafter.
  Input:
  alt geometric altitude, km.
  Return: (sigma, delta, theta)
  sigma   density/sea-level standard density
  delta   pressure/sea-level standard pressure
  theta   temperature/sea-level std. temperature
  """

    REARTH = 6371.0  # radius of the Earth (km)
    GMR = 34.163195


    htab = [0.0, 11.0, 20.0, 32.0, 47.0,
            51.0, 71.0, 84.852]
    ttab = [288.15, 216.65, 216.65, 228.65, 270.65,
            270.65, 214.65, 186.946]
    ptab = [1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3, 1.0945601E-3,
            6.6063531E-4, 3.9046834E-5, 3.68501E-6]
    gtab = [-6.5, 0.0, 1.0, 2.8, 0, -2.8, -2.0, 0.0]

    h = alt * REARTH / (alt + REARTH)  # geometric to geopotential altitude

    i = 0
    j = len(htab)
    while (j > i + 1):
        k = (i + j) // 2  # this is integer division in Python 3
        if h < htab[k]:
            j = k
        else:
            i = k
    tgrad = gtab[i]  # temp. gradient of local layer
    tbase = ttab[i]  # base  temp. of local layer
    deltah = h - htab[i]  # height above local base
    tlocal = tbase + tgrad * deltah  # local temperature
    theta = tlocal / ttab[0]  # temperature ratio

    if 0.0 == tgrad:
        delta = ptab[i] * math.exp(-GMR * deltah / tbase)
    else:
        delta = ptab[i] * math.pow(tbase / tlocal, GMR / tgrad)
    sigma = delta / theta
    return (sigma, delta, theta)


def thrust():
    # F= mv2 + (p2-p3)*a2
    # Define Engine Constants
    m_dot = 1409/3
    v_2 = 3560
    a_2 = (2.4**2)*(math.pi*0.25)
    p_2 = 72326 #Pascals

    # alt = input("How high are you flying today? \n")
    alt = int(86)

    thrust_value = []
    atmospheric_value = []
    temp_value = []
    density_value = []
    geometric_heights = []

    # every half kilometer,
    for altitude in np.linspace(0, alt, alt * 2, endpoint=False):
        thrust_value.append((m_dot * v_2) +
                            ((p_2 - atmosphere(altitude)[1] * PZERO) * a_2))
        atmospheric_value.append(atmosphere(altitude)[1] * PZERO)
        temp_value.append(atmosphere(altitude)[2] * TZERO)
        density_value.append(atmosphere(altitude)[0] * RHOZERO)
        geometric_heights.append(altitude)

    return (thrust_value, atmospheric_value, temp_value, density_value, geometric_heights)


def plot_thrust():
    alt = 86
    fig = plt.figure()
    axt = fig.add_subplot(111)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    # Big plot
    axt.set_ylabel('Altitude')
    axt.spines['top'].set_color('none')
    axt.spines['bottom'].set_color('none')
    axt.spines['left'].set_color('none')
    axt.spines['right'].set_color('none')
    axt.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    # Plot 1
    x = thrust()[1]
    y = np.linspace(0, alt, alt * 2, endpoint=False)
    ax1.title.set_text('Atmospheric Pressure')
    ax1.set_xlabel('Pressure [Pa]')
    ax1.plot(x, y)

    # Plot 2
    x = thrust()[0]
    ax2.title.set_text('Combined Thrust')
    ax2.set_xlabel('Force [N]')
    ax2.plot(x, y)

    # Plot 3
    x = thrust()[2]
    ax3.title.set_text('Temperature')
    ax3.set_xlabel('Temperature [K]')
    ax3.plot(x, y)

    # Plot 4
    x = thrust()[3]
    ax4.title.set_text('Density')
    ax4.set_xlabel('Density [kg/m^3]')
    ax4.plot(x, y)

    # adjust plot margins
    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.10, right=0.9, hspace=0.4, wspace=0.2)
    plt.suptitle("Rocket Parameters Vs. Altitude")
    plt.show()

if __name__ == "__main__":
    thrust()
    plot_thrust()

