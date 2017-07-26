import matplotlib.pyplot as plt
import numpy as np

'''
This program reads in a CSV file obtained from Chris Monson at 
https://plot.ly/~chris.monson/6/?share_key=Q0cEqeD0iGYhILBF8IaK90#data . This was taken from the SpaceX CRS 8 Mission 
to the ISS.



'''

time, velocity, altitude = np.loadtxt('CRS_8.csv', delimiter=',', unpack=True, skiprows=1)


def acceleration(list, time):
    #  Uses forward difference

    accel = []
    for i in range(len(list)):
        if i < len(list)-1:
            delta_t = time[i+1] - time[i]
            accel.append((list[i+1] - list[i])/delta_t)
        else:
            accel.append((list[i] - list[i-1]) / delta_t)

    return accel

accel = acceleration(velocity, time)
vertical_velocity = acceleration(altitude, time)


def plot_csv():
    fig = plt.figure()
    axt = fig.add_subplot(111)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    # Big plot
    axt.spines['top'].set_color('none')
    axt.spines['bottom'].set_color('none')
    axt.spines['left'].set_color('none')
    axt.spines['right'].set_color('none')
    axt.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

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
    ax2.plot(x, y)

    # Plot 3
    x = time
    y = accel
    ax3.title.set_text('Acceleration')
    ax3.set_xlabel('Time [s]')
    ax3.plot(x, y)

    # # Plot 4
    # x = time
    # y = vertical_velocity
    # ax4.title.set_text('Density')
    # ax4.set_xlabel('Density [kg/m^3]')
    # ax4.plot(x, y)

    # adjust plot margins
    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.10, right=0.9, hspace=0.4, wspace=0.2)
    plt.suptitle("Falcon 9 Full Thrust, CRS 8")
    plt.show()

if __name__ == "__main__":
    plot_csv()
