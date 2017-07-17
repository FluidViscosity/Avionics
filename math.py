from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

#Define constants
mdot = 3 #~constant
exhaust_v = 2000 #~constant
g = 9.81 #Actually changes with altitude
density = 1.225 #Actually changes with altitude
mass = 500 #Actually changes with altitude
A_ref = 2 #reference area
C_d = 0.3 #coefficient of drag

def d(z, t, mdot, exhaust_v, g, density, A_ref, C_d, mass):
    return  np.array((
        (mdot*exhaust_v)/mass - g - (density*A_ref*C_d*(z[1])**2)/(2*mass), z[0]
    ))

x0 = np.array([100,0])
tmax, dt = 10, 0.01
t = np.linspace(0,tmax, num=np.round(tmax/dt)+1)
sol = odeint(d,x0,t, args=(mdot, exhaust_v, g, density, A_ref, C_d, mass)) #0 = height, 1 = velocity

#Plotting
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

# Plot 1 Height Vs Time
x = t
y = sol[...,0]
ax1.title.set_text('Height Vs Time')
ax1.set_ylabel('Height [m]')
ax1.plot(x, y)

# Plot 2 Velocity Vs Time
x = t
y = sol[...,1]
ax2.title.set_text('Velocity Vs Time')
ax2.set_ylabel('Velocity [m/s]')
ax2.plot(x, y)
plt.subplots_adjust(top=0.9, bottom=0.1, left=0.10, right=0.9, hspace=0.4, wspace=0.2)
plt.show()



# To solve a second-order ODE using scipy.integrate.odeint, you should write it as a system of first-order ODEs:
#
# I'll define z = [x', x], then z' = [x'', x'], and that's your system! Of course, you have to plug in your real relations:
#
# x'' = -(b*x'(t) + k*x(t) + a*(x(t))^3 + m*g) / m
# becomes:
#
# z[0]' = -1/m * (b*z[0] + k*z[1] + a*z[1]**3 + m*g)
# z[1]' = z[0]
# Or, just call it d(z):
#
# def d(z, t):
#     return np.array((
#                      -1/m * (b*z[0] + k*z[1] + a*z[1]**3 + m*g),  # this is z[0]'
#                      z[0]                                         # this is z[1]'
#                    ))
# Now you can feed it to the odeint as such:
#
# _, x = odeint(d, x0, t).T
# (The _ is a blank placeholder for the x' variable we made)
#
# In order to minimize b subject to the constraint that the maximum of x is always negative, you can use scipy.optimize.minimize. I'll implement it by actually maximizing the maximum of x, subject to the constraint that it remains negative, because I can't think of how to minimize a parameter without being able to invert the function.
# #
# from scipy.optimize import minimize
# from scipy.integrate import odeint
# import numpy as np
#
# m = 1220
# k = 35600
# g = 17.5
# a = 450000
# z0 = np.array([-.5, 0])
# t = 5
#
# def d(z, t, m, k, g, a, b):
#     return np.array([-1/m * (b*z[0] + k*z[1] + a*z[1]**3 + m*g), z[0]])
#
# def func(b, z0, *args):
#     _, x = odeint(d, z0, t, args=args+(b,)).T
#     return -x.max()  # minimize negative max
#
# cons = [{'type': 'ineq', 'fun': lambda b: b - 1000, 'jac': lambda b: 1},   # b > 1000
#         {'type': 'ineq', 'fun': lambda b: 10000 - b, 'jac': lambda b: -1}, # b < 10000
#         {'type': 'ineq', 'fun': lambda b: func(b, z0, m, k, g, a)}] # func(b) > 0 means x < 0
#
# b0 = 10000
# b_min = minimize(func, b0, args=(z0, m, k, g, a), constraints=cons)

