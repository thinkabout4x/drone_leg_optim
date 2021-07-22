import scipy
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import animation


class Four_bar:

    def __init__(self, theta_init):

        self.phi_my = math.radians(160)
        self.theta_2 = theta_init
        self.x_4_1 = 0.007
        self.y_4_1 = 0.023

        self.phi_0 = 2*math.pi - self.phi_my
        self.theta_1 = -math.atan2(self.x_4_1, self.y_4_1)

        self.x_1_0 = 0.1
        self.y_1_0 = 0.045
        self.r_1 = math.sqrt(self.x_4_1*self.x_4_1+self.y_4_1*self.y_4_1)
        self.r_2 = 0.07736
        self.r_3 = 0.015
        self.r_4 = 0.065
        self.l_1 = 0.0325
        self.l_2 = math.sqrt(self.l_1*self.l_1+self.r_3*self.r_3-2*self.l_1*self.r_3*math.cos(self.phi_my))

        self.points = self.mech_calc(self.theta_2)

    def mech_calc(self, t_2):
        points = np.zeros((5, 2))
        z = math.sqrt(self.r_1 * self.r_1 + self.r_2 * self.r_2 - 2 * self.r_1 * self.r_2 * math.cos(t_2))  # +-
        gamma = -math.acos((self.r_3 * self.r_3 + self.r_4 * self.r_4 - z * z) / (2 * self.r_3 * self.r_4))  # +-
        alpha = -math.acos((z * z - self.r_3 * self.r_3 + self.r_4 * self.r_4) / (2 * z * self.r_4))  # +-
        beta = math.acos((z * z - self.r_2 * self.r_2 + self.r_1 * self.r_1) / (2 * z * self.r_1))  # ?
        theta_3 = math.pi - (alpha+beta+gamma)
        theta_4 = math.pi - (alpha+beta)
        points[0] = [self.x_1_0, self.y_1_0]
        points[1] = [points[0, 0] + self.r_2 * math.cos(t_2 + self.theta_1),
                     points[0, 1] + self.r_2 * math.sin(t_2 + self.theta_1)]
        points[4] = [points[0, 0] + self.r_1 * math.cos(self.theta_1),
                     points[0, 1] + self.r_1 * math.sin(self.theta_1)]
        points[3] = [points[4, 0] + self.r_4 * math.cos(theta_4 + self.theta_1),
                     points[4, 1] + self.r_4 * math.sin(theta_4 + self.theta_1)]
        points[2] = [points[1, 0] + self.l_1*math.cos(self.phi_0+theta_3+self.theta_1), points[1, 1] +
                     self.l_1 * math.sin(self.phi_0 + theta_3 + self.theta_1)]
        return points

    def mech_draw(self, l_1, l_2, l_3):
        l_1.set_data(self.points[:, 0], self.points[:, 1])
        l_2.set_data([self.points[1, 0], self.points[3, 0]], [self.points[1, 1], self.points[3, 1]])
        l_3.set_data([self.points[0, 0], self.points[4, 0]], [self.points[0, 1], self.points[4, 1]])
        return [l_1, l_2, l_3]


def animate(i, m, t_2, l_1, l_2, l_3):
    m.points = m.mech_calc(math.radians(t_2[i]))
    return m.mech_draw(l_1, l_2, l_3)


if __name__ == '__main__':
    theta_2 = np.linspace(80, 20, num=100)
    mech = Four_bar(math.radians(theta_2[0]))
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(0.05, 0.2), ylim=(0, 0.18))
    ax.set_aspect('equal')
    ax.grid()
    line_1, = ax.plot([], [], 'o-', lw=2)
    line_2, = ax.plot([], [], 'o-', lw=2)
    line_3, = ax.plot([], [], 'o-', lw=2)
    # mech.mech_draw(line_1, line_2, line_3)
    ani = animation.FuncAnimation(fig, animate, len(theta_2), fargs=[mech, theta_2, line_1, line_2, line_3],
                                  interval=50)
    plt.show()
