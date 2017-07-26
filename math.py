from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math
from ThrustCalc import atmosphere

# Define constants
mdot =  273.6 #1409/3  # ~constant
exhaust_v = 2239 #3560  # ~constant
g = 9.81  # Actually changes with altitude
# density = 1.225  # Actually changes with altitude
mass = 541300  # Actually changes with altitude
A_ref = 3.66**2 * math.pi  # reference area
C_d = 0.3  # coefficient of drag
a_nozzle = (2.4 ** 2) * (math.pi * 0.25)
p_2 = 72326  # Pascals

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


class Rocket():
    def __init__(self, height, velocity, mass):
        self.height = height
        self.velocity = velocity
        self.mass = mass


shuttle = Rocket(0, 0, mass)


def local_density(height):
    return atmosphere(height)[0] * RHOZERO


def local_pressure(height):
    return atmosphere(height)[1] * PZERO


def current_mass():
    shuttle.mass = shuttle.mass - mdot
    return shuttle.mass


def d(z, t, mdot, exhaust_v, g, density, A_ref, C_d, mass):
    return  np.array((
        (mdot*exhaust_v)/shuttle.mass - g + ((p_2 - local_pressure(shuttle.height))*a_nozzle)/shuttle.mass -
        (local_density(shuttle.height)*A_ref*C_d*(z[1])**2)/(2*shuttle.mass), z[0]
    ))


x0 = np.array([0, 0])
tmax, dt = 1, 0.1
t = np.linspace(0, tmax, num=np.round(tmax/dt)+1)
total_time = 50

recorded_height = []
recorded_velocity = []
# sol = odeint(d,x0,t, args=(mdot, exhaust_v, g, density(shuttle.height), A_ref, C_d, mass)) #0 = height, 1 = velocity


for time in range(total_time):
    sol = odeint(d, x0, t, args=(mdot, exhaust_v, g, local_density(shuttle.height), A_ref, C_d, shuttle.mass))
    shuttle.height = sol[:, 0][-1]
    shuttle.velocity = sol[:, 1][-1]
    shuttle.mass = current_mass()
    recorded_height.append(shuttle.height)
    recorded_velocity.append(shuttle.velocity)
    plt.plot(time, recorded_height[-1])
    print(x0)

    x0 = np.array([shuttle.height, shuttle.velocity])

#Plotting
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

plt.subplots_adjust(top=0.9, bottom=0.1, left=0.10, right=0.9, hspace=0.4, wspace=0.2)

# Plot 1 Height Vs Time
x = np.linspace(0, total_time, total_time)
y = recorded_height
ax1.title.set_text('Height Vs Time')
ax1.set_ylabel('Height [m]')
ax1.plot(x, y)

# Plot 2 Velocity Vs Time
x = np.linspace(0, total_time, total_time)
y = recorded_velocity
ax2.title.set_text('Velocity Vs Time')
ax2.set_ylabel('Velocity [m/s]')
ax2.plot(x, y)
plt.subplots_adjust(top=0.9, bottom=0.1, left=0.10, right=0.9, hspace=0.4, wspace=0.2)
plt.show()



