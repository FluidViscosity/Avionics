import matplotlib.pyplot as plt
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import splrep, splev, UnivariateSpline

'''
This program reads in a CSV file obtained from Chris Monson at 
https://plot.ly/~chris.monson/6/?share_key=Q0cEqeD0iGYhILBF8IaK90#data . This was taken from the SpaceX CRS 8 Mission 
to the ISS.



'''
#  unpack the csv file to lists time, velocity and altitude
time, velocity, altitude = np.loadtxt('CRS_12.csv', delimiter=',', unpack=True, skiprows=1)

dt = time[0]

#  Using numpy's gradient function, the acceleration can be found
acceleration = np.gradient(velocity, dt, edge_order=2)
velocity_r = np.gradient(altitude, dt, edge_order=2)

# filter noise results, frac values from 0.0025 to 0.005 work well
accel_filtered = lowess(acceleration, time, is_sorted=True, frac=0.025, it=0)
velocity_r_filtered = lowess(velocity_r, time, is_sorted=True, frac=0.005, it=0)

poly_a_spline = splrep(time,altitude, k=3, s=3)
plt.plot(time, splev(time, poly_a_spline, der=1))
#plt.plot(time, np.gradient(poly_a_spline, dt, edge_order=2))

plt.show()
plt.plot(time, splev(time, poly_a_spline))
plt.show()


def plot_csv():
    fig = plt.figure()
    axt = fig.add_subplot(111)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    #
    # fig2 = plt.figure()
    # axt2 = fig2.add_subplot(111)
    # ax3 = fig2.add_subplot(211)
    # ax4 = fig2.add_subplot(212)



    # Big plot
    axt.spines['top'].set_color('none')
    axt.spines['bottom'].set_color('none')
    axt.spines['left'].set_color('none')
    axt.spines['right'].set_color('none')
    axt.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    # axt2.spines['top'].set_color('none')
    # axt2.spines['bottom'].set_color('none')
    # axt2.spines['left'].set_color('none')
    # axt2.spines['right'].set_color('none')
    # axt2.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    # Plot 1 - Height vs Time
    x = time
    y = altitude
    ax1.title.set_text('Height vs Time')
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Height [km]')
    ax1.plot(x, y, color='g')

    # Plot 2 Velocity vs Time
    x = time
    y = velocity
    ax2.title.set_text('Velocity vs Time')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Velocity [km/h]')
    ax2.plot(x, y, label='Total Velocity')
    ax2.plot(x, velocity_r, label='Radial Velocity')


    # # Plot 3
    # x = time
    # y = acceleration
    # ax3.title.set_text('Acceleration, raw')
    # ax3.set_xlabel('Time [s]')
    # ax3.plot(x, y)
    #
    # # # Plot 4
    # x = accel_filtered[:, 0]
    # y = accel_filtered[:, 1]
    # ax4.title.set_text('Acceleration, filtered')
    # ax4.set_xlabel('Time [s]')
    # ax4.plot(x, y)
    #
    # #Annotation
    # ax4.annotate('Fairing Separation', (x[8601], y[8601]), xytext = (0.5, 0.2), textcoords = 'axes fraction',
    #              arrowprops = dict(arrowstyle="->",facecolor='grey', color='grey'))
    # ax4.annotate('MECO', (x[4500], y[4500]), xytext = (0.3, 0.5), textcoords = 'axes fraction',
    #              arrowprops = dict(arrowstyle="->",facecolor='red', color='red'))
    # #ax4.annotate('SECO', (x[17711], y[17711]), xytext = (0.9, 0.9), textcoords = 'axes fraction',
    #            #  arrowprops = dict(arrowstyle="->",facecolor='blue', color='blue'))
    # ax4.annotate('Max Q ''Bucket''', (x[1957], y[1597]), xytext = (0.05, 0.9), textcoords = 'axes fraction',
    #              arrowprops = dict(arrowstyle="->",facecolor='green', color='green'))

    # adjust plot margins
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.10, right=0.9, hspace=0.4, wspace=0.2)
   # fig2.subplots_adjust(top=0.9, bottom=0.1, left=0.10, right=0.9, hspace=0.4, wspace=0.2)
    fig.suptitle("Falcon 9 Full Thrust, CRS 8")
   # fig2.suptitle("Falcon 9 Full Thrust, CRS 8")
    plt.show()

if __name__ == "__main__":
    plot_csv()

