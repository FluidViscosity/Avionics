from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import itertools
from ThrustCalc import thrust, plot_thrust, atmosphere
##############################################################
# Need a function to define all constants used by Atmosphere #
##############################################################

mass = 580
empty_mass = 70
mdot = 3
c = 2000
CD = 0.8  # Drag coeffiecient of a long blunt cylinder



fig = plt.figure()

#Define integration arrays and boundary conditions (integrate over 1 second)
ts = np.linspace(0, 1, 101)

ts2 = np.linspace(0,1,101*(mass-empty_mass)//mdot)
#only works when mass minus empty mass is a multiple of mdot....

''' Define a function which calculates the derivative i.e. acceleration.
This includes a momentum thrust and pressure thrust terms (positive) as well as
 gravity loses and a drag term
'''


def dv_dt(v, t):
    return mdot * c / mass - 9.81


def integrator():

    global mass
    v0 = 0.0
    final_v = []

    tanks_empty = False

    while tanks_empty is False:

        if mass > empty_mass:
            mass = mass - mdot
            vs = odeint(dv_dt, v0, ts)  # need to check on the output format of odeint function. returning a weird array
            vs = np.array(vs).flatten()

            v0 = vs[-1]
            final_v.append(vs.tolist())

        else:
            tanks_empty = True

    merged = list(itertools.chain.from_iterable(final_v))
    #print(merged)

    return merged


ax1 = fig.add_subplot(111)
x = ts2
y = integrator()[:]

ax1.title.set_text('Integrated Velocity')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Velocity [m/s]')
ax1.plot(x, y, 'g')

plt.show()
thrust()
#plot_thrust()

