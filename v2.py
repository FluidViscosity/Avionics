from inst_acc import *
from engineA import massFlowRate
import math
import configparser
import os
import sys


# path = os.path
path = os.path.join(sys.path[0],'Engines\\engineA_001.ini')
config = configparser.ConfigParser()

config.read(path)


class Rocket:

	def __init__(self, config ):

		self.mass = config.getfloat('PHYSICAL','mass')
		self.D_outer = config.getfloat('PHYSICAL','Diameter_outer')
		self.L = config.getfloat('PHYSICAL','Length')
		self.L_nc = config.getfloat('PHYSICAL','Nosecone_length')

		self.N_fin = config.getfloat('FIN', 'Number_of_fins')
		self.rootchord_fin = config.getfloat('FIN', 'Root_chord')
		self.tiproot_ratio_fin = config.getfloat('FIN', 'Tip_root_ratio')
		self.tipchord_fin = self.rootchord_fin*self.tiproot_ratio_fin
		self.x_t_fin = config.getfloat('FIN', 'Pos_max_thickness')
		self.thickness_fin = config.getfloat('FIN', 'Max_thickness')
		self.span_fin = config.getfloat('FIN', 'Span')

		self.N_grain = config.getfloat('GRAIN', 'Number_of_grains')
		self.d = config.getfloat('GRAIN', 'Diameter_core')
		self.D_grain = config.getfloat('GRAIN', 'Diameter_grain_outer')
		self.L_grain = config.getfloat('GRAIN', 'Segment_Length')

		self.grain_density = config.getfloat('PROPELLANT', 'Prop_density')
		self.k_2_phase = config.getfloat('PROPELLANT', 'Specific_heat_ratio')
		self.c_star = config.getfloat('PROPELLANT', 'c_star')

		self.A_throat = config.getfloat('NOZZLE', 'Throat_diameter')**2 * (math.pi/4)
		self.A_exit = config.getfloat('NOZZLE', 'Exit_diameter')**2 * (math.pi/4)
		self.p_0 = config.getfloat('NOZZLE', 'Chamber_pressure')

		self.A_cs = (math.pi/4)*self.D_outer**2 # [m^2] Cross sectional area
		self.z = 0
		self.v = 0.0001
		self.a = 0

	def burnRate(self):
		if self.p_0 < 779000:
			a = 8.88
			n = 0.619
		if self.p_0 >= 779000 and self.p_0 <2570000:
			a = 7.55
			n = -0.009
		if self.p_0 >= 2570000 and self.p_0 <5930000:
			a = 3.84
			n = 0.688
		if self.p_0 >= 5930000 and self.p_0 < 8500000:
			a = 17.2
			n = -0.148
		if self.p_0 >= 8500000 and self.p_0 < 11200000:
			a = 4.78
			n = 0.442

		r = a*(self.p_0/1000000)**n
		return(r,a,n)


R = Rocket(config)
# R O C K E T   P A R A M E T E R S
cd = 0.6 # Constant drag coefficient
# mass = 5.65 # [kg] Wet mass of Rocket
diameter = 0.1 # [m] diameter for cross sectional area


#
# d = 0.01 # grain core diameter
# L = 0.15 # grain length

dt = 0.1 # time step
time = 0

burning = True

while burning is True:
	InstantParams = inst_params(R) # checks atmospheric conditions and motor thrust
	R.a = (-gravity(R.z, R.mass)/R.mass) + (InstantParams[0]/R.mass) -((0.5*InstantParams[2]*R.v**2 *dragCoeff(R)*R.A_cs)/R.mass)

	R.v += R.a*dt
	R.z += R.v*dt
	time += dt

	mass_flow_condition = massFlowRate(R, dt)
	if mass_flow_condition == False:
		burning = False
	else:
		R.mass-=mass_flow_condition[0]
	print (R.a,R.v,R.z)

print ('Burnout!')
print ('Time to burnout:', time)

while R.v>0:
	InstantParams = inst_params(R)
	R.a = (-gravity(R.z, R.mass)/R.mass)  - ((0.5*InstantParams[2]*R.v**2*dragCoeff(R)*R.A_cs) / R.mass)

	R.v += R.a * dt
	R.z += R.v * dt
	time += dt
	print(R.a, R.v, R.z)

print ('Apogee!')
print ('Time to apogee:', time)
print(R.p_0)

# print ('Gravity:', gravity(z/1000, mass)/mass)

# if __name__ == "__main__":
#     inst_thrust(86)
#     # plot_thrust()

if __name__ == "__main__":
	R.v = 332*1.5
	R.z = 2000
	# print(dragCoeff(R))
