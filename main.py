import scipy
import numpy as np
import math
import matplotlib.pyplot as plt


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

    def mech_calc(self, theta_2):
        points = np.zeros((5, 2))
        z = math.sqrt(self.r_1 * self.r_1 + self.r_2 * self.r_2 - 2 * self.r_1 * self.r_2 * math.cos(theta_2))  # +-
        gamma = math.acos((self.r_3 * self.r_3 + self.r_4 * self.r_4 - z * z) / (2 * self.r_3 * self.r_4))
        alpha = math.acos((z * z - self.r_3 * self.r_3 + self.r_4 * self.r_4) / (2 * z * self.r_4))
        beta = math.acos((z * z - self.r_2 * self.r_2 + self.r_1 * self.r_1) / (2 * z * self.r_1))
        theta_3 = math.pi - (alpha+beta+gamma)
        theta_4 = math.pi - (alpha+beta)
        points[0] = [self.x_1_0, self.y_1_0]
        points[1] = [points[0, 0] + self.r_2 * math.cos(theta_2 + self.theta_1),
                     points[0, 1] + self.r_2 * math.sin(theta_2 + self.theta_1)]
        points[3] = [points[0, 0] + self.r_1 * math.cos(self.theta_1),
                     points[0, 1] + self.r_1 * math.sin(self.theta_1)]
        points[2] = [points[3, 0] + self.r_4 * math.cos(theta_4 + self.theta_1),
                     points[3, 1] + self.r_4 * math.sin(theta_4 + self.theta_1)]
        points[4] = [points[1, 0] + self.l_1*math.cos(self.phi_0+theta_3+self.theta_1), points[1, 1] +
                     self.l_1 * math.sin(self.phi_0 + theta_3 + self.theta_1)]
        print(points[4])
        return points

    def mech_draw(self):
        print(self.points[0])
        print(self.points[1])
        plt.plot([self.points[0, 0], self.points[1, 0]], [self.points[0, 1], self.points[1, 1]], marker='o')
        plt.plot([self.points[0, 0], self.points[3, 0]], [self.points[0, 1], self.points[3, 1]], marker='o')
        plt.plot([self.points[1, 0], self.points[2, 0]], [self.points[1, 1], self.points[2, 1]], marker='o')
        plt.plot([self.points[2, 0], self.points[3, 0]], [self.points[2, 1], self.points[3, 1]], marker='o')
        plt.plot([self.points[1, 0], self.points[4, 0]], [self.points[1, 1], self.points[4, 1]], marker='o')
        plt.plot([self.points[2, 0], self.points[4, 0]], [self.points[2, 1], self.points[4, 1]], marker='o')
        plt.show()
        return


if __name__ == '__main__':
    mech = Four_bar(math.radians(58.5))
    mech.mech_draw()
