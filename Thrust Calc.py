import sys
import math
import numpy


'''version = "1.0 (2004 Oct 04)"
greeting = "Tables - A Python program to compute atmosphere tables"
author = "Ralph L. Carmichael, Public Domain Aeronautical Software"
modifier = "Tomer Simhony"

'''

#   P H Y S I C A L   C O N S T A N T S

FT2METERS = 0.3048    # mult. ft. to get meters (exact)
KELVIN2RANKINE = 1.8    # mult deg K to get deg R
PSF2NSM = 47.880258   # mult lb/sq.ft to get sq.m
SCF2KCM = 515.379   # mult slugs/cu.ft to get kg/cu.m
TZERO = 288.15    # sea-level temperature, kelvins
PZERO = 101325.0    # sea-level pressure, N/sq.m
RHOZERO = 1.225     # sea-level density, kg/cu.m
AZERO = 340.294   # speed of sound at S.L.  m/sec
BETAVISC = 1.458E-6   # viscosity constant
SUTH = 110.4     # Sutherland's constant, kelvins


def Atmosphere(alt):
  """ Compute temperature, density, and pressure in standard atmosphere.
  Correct to 86 km.  Only approximate thereafter.
  Input:
  alt geometric altitude, km.
  Return: (sigma, delta, theta)
  sigma   density/sea-level standard density
  delta   pressure/sea-level standard pressure
  theta   temperature/sea-level std. temperature
  """

  REARTH = 6369.0     # radius of the Earth (km)
  GMR = 34.163195
  NTAB = 8            # length of tables

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
    k = (i + j) // 2      # this is integer division in Python 3
    if h < htab[k]:
      j = k
    else:
      i = k
  tgrad = gtab[i]     # temp. gradient of local layer
  tbase = ttab[i]     # base  temp. of local layer
  deltah = h - htab[i]        # height above local base
  tlocal = tbase + tgrad * deltah   # local temperature
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
    m_dot = 24.875
    v_2 = 2430
    a_2 = 0.0574
    p_2 = 70000

    alt = input("How high are you flying today? \n")
    alt = int(alt)

    thrust_value = []
    atmospheric_value = []

    # every half kilometer,
    for altitude in numpy.linspace(0, alt, alt * 2, endpoint=False):

        thrust_value.append((m_dot * v_2) + \
                           ((p_2 - Atmosphere(altitude)[1] * PZERO) * a_2))
        atmospheric_value.append(Atmosphere(altitude)[1] * PZERO)
    return print(thrust_value, "\n \n \n", atmospheric_value)

thrust()
