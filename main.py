from scipy.optimize import differential_evolution
from scipy.optimize import NonlinearConstraint
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import animation


class Four_bar:

    def __init__(self, params_init):

        self.phi_my = params_init[0]
        self.x_4_1 = params_init[1]
        self.y_4_1 = params_init[2]

        self.phi_0 = 2*math.pi - self.phi_my
        self.theta_1 = -math.atan2(self.x_4_1, self.y_4_1)

        self.x_1_0 = params_init[3]
        self.y_1_0 = params_init[4]
        self.r_1 = math.sqrt(self.x_4_1*self.x_4_1+self.y_4_1*self.y_4_1)
        self.r_2 = params_init[5]
        self.r_3 = params_init[6]
        self.r_4 = params_init[7]
        self.l_1 = params_init[8]
        self.l_2 = math.sqrt(self.l_1*self.l_1+self.r_3*self.r_3-2*self.l_1*self.r_3*math.cos(self.phi_my))
        self.theta_2 = params_init[9]

        self.points = self.mech_calc(self.theta_2)

        print('s:', self.r_3)
        print('l', self.r_2)
        print('p', self.r_1)
        print('q', self.r_4)

    def mech_calc(self, t_2):
        points = np.zeros((5, 2))
        z = math.sqrt(self.r_1 * self.r_1 + self.r_2 * self.r_2 - 2 * self.r_1 * self.r_2 * math.cos(t_2))  # +
        gamma = -math.acos((self.r_3 * self.r_3 + self.r_4 * self.r_4 - z * z) / (2 * self.r_3 * self.r_4))  # -
        alpha = -math.acos((z * z - self.r_3 * self.r_3 + self.r_4 * self.r_4) / (2 * z * self.r_4))  # -
        beta = math.acos((z * z - self.r_2 * self.r_2 + self.r_1 * self.r_1) / (2 * z * self.r_1))  # +
        theta_3 = math.pi - (alpha+beta+gamma)  # +-
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

    def mech_draw(self, l_1_init, l_2_init, l_3_init):
        l_1_init.set_data(self.points[:, 0], self.points[:, 1])
        l_2_init.set_data([self.points[1, 0], self.points[3, 0]], [self.points[1, 1], self.points[3, 1]])
        l_3_init.set_data([self.points[0, 0], self.points[4, 0]], [self.points[0, 1], self.points[4, 1]])
        return [l_1_init, l_2_init, l_3_init]


def animate(i, m, t_2, l_1_init, l_2_init, l_3_init):
    m.points = m.mech_calc(math.radians(t_2[i]))
    return m.mech_draw(l_1_init, l_2_init, l_3_init)


def de_mech_calc(params_init, *args):
    phi_my_init, x_4_1_init, y_4_1_init, x_1_0_init, y_1_0_init, r_2_init, r_3_init, r_4_init, l_1_init, theta_2_1_init, theta_2_2_init = params_init
    points_des = np.array(args)
    theta_2_vec = np.array([[theta_2_1_init], [theta_2_2_init]])
    points_result = np.zeros((2, 3))
    phi_0 = 2 * math.pi - phi_my_init
    theta_1 = -math.atan2(x_4_1_init, y_4_1_init)
    r_1 = math.sqrt(x_4_1_init**2 + y_4_1_init**2)
    for i in range(np.shape(points_des)[0]):
        z = math.sqrt(r_1**2 + r_2_init**2 - 2 * r_1 * r_2_init * math.cos(theta_2_vec[i]))  # +-
        try:
            gamma = -math.acos((r_3_init**2 + r_4_init**2 - z * z) / (2 * r_3_init * r_4_init))  # +-
        except ValueError:
            # print([r_1, r_2_init, r_3_init, r_4_init])
            # print(smallest_link(params_init))
            # print(longest_link(params_init))
            return 1
        alpha = -math.acos((z**2 - r_3_init**2 + r_4_init**2) / (2 * z * r_4_init))  # +-
        beta = math.acos((z**2 - r_2_init**2 + r_1**2) / (2 * z * r_1))  # ?
        theta_3 = math.pi - (alpha + beta + gamma)
        x_p = x_1_0_init + r_2_init * math.cos(theta_2_vec[i] + theta_1) + l_1_init*math.cos(phi_0+theta_3+theta_1)
        y_p = y_1_0_init + r_2_init * math.sin(theta_2_vec[i] + theta_1) + l_1_init*math.sin(phi_0+theta_3+theta_1)
        points_result[i, :] = [x_p, y_p, theta_3]
    obj_func = 0
    for i in range(np.shape(points_des)[0]):
        obj_func += (points_des[i, 0]-points_result[i, 0])**2 + (points_des[i, 1]-points_result[i, 1])**2+(points_des[i, 2]-points_result[i, 2])**2
    return obj_func


def theta_2_1_limit_up(params_init):
    _, x_4_1_init, y_4_1_init, _, _, r_2_init, r_3_init, r_4_init, _, theta_2_1_init, _ = params_init
    r_1 = math.sqrt(x_4_1_init**2 + y_4_1_init**2)
    try:
        r = math.acos(((r_3_init + r_4_init) ** 2 - r_1 ** 2 - r_2_init ** 2) / (-2 * r_1 * r_2_init)) - theta_2_1_init
    except ValueError:
        r = -1
    return r


def theta_2_2_limit_up(params_init):
    _, x_4_1_init, y_4_1_init, _, _, r_2_init, r_3_init, r_4_init, _, _, theta_2_2_init = params_init
    r_1 = math.sqrt(x_4_1_init**2 + y_4_1_init**2)
    try:
        r = math.acos(((r_3_init+r_4_init)**2-r_1**2-r_2_init**2)/(-2*r_1*r_2_init))-theta_2_2_init
    except ValueError:
        r = -1
    return r


def theta_comp(params_init):
    _, _, _, _, _, _, _, _, _, theta_2_1_init, theta_2_2_init = params_init
    return theta_2_1_init - theta_2_2_init


def non_grashof_crit(params_init):
    _, x_4_1_init, y_4_1_init, _, _, r_2_init, r_3_init, r_4_init, _, _, _ = params_init
    r_1 = math.sqrt(x_4_1_init**2 + y_4_1_init**2)
    return r_3_init + r_2_init - r_1 - r_4_init


def smallest_link(params_init):
    _, x_4_1_init, y_4_1_init, _, _, r_2_init, r_3_init, r_4_init, _, _, _ = params_init
    r_1 = math.sqrt(x_4_1_init**2 + y_4_1_init**2)
    r = min([r_1, r_2_init, r_4_init])-r_3_init
    return r


def longest_link(params_init):
    _, x_4_1_init, y_4_1_init, _, _, r_2_init, r_3_init, r_4_init, _, _, _ = params_init
    r_1 = math.sqrt(x_4_1_init**2 + y_4_1_init**2)
    r = r_2_init - max([r_1, r_3_init, r_4_init])
    return r


if __name__ == '__main__':
    phi_my = math.radians(160)
    x_4_1 = 0.007
    y_4_1 = 0.023
    x_1_0 = -0.1
    y_1_0 = 0.045
    r_2 = 0.07736
    r_3 = 0.015
    r_4 = 0.065
    l_1 = 0.0325

    bounds = [(phi_my*0.8, phi_my*1.2), (0, x_4_1*1.2), (0, y_4_1*1.2), (x_1_0*1.2, 0), (0, y_1_0*1.2),
              (r_2*0.8, r_2*1.2), (r_3*0.8, r_3*1.2), (r_4*0.8, r_4*1.2), (0, l_1*1.5),
              (math.radians(70), math.radians(90)), (0, math.radians(30))]
    points_desired = ([-0.072, 0.15, math.radians(268.83)], [-0.023, 0.018, math.radians(88.83)])

    de_params = [phi_my, x_4_1, y_4_1, x_1_0, y_1_0, r_2, r_3, r_4, l_1, math.radians(80), math.radians(20)]
    # mech = Four_bar(de_params[:-1])
    nlc_1 = NonlinearConstraint(theta_2_1_limit_up, 0, np.inf)
    nlc_2 = NonlinearConstraint(theta_2_2_limit_up, 0, np.inf)
    nlc_3 = NonlinearConstraint(theta_comp, 0, np.inf)
    nlc_4 = NonlinearConstraint(non_grashof_crit, 0, np.inf)
    nlc_5 = NonlinearConstraint(smallest_link, 0, np.inf)
    nlc_6 = NonlinearConstraint(longest_link, 0, np.inf)
    result = differential_evolution(de_mech_calc, bounds, args=points_desired, strategy="best2bin",
                                    constraints=(nlc_1, nlc_2, nlc_3, nlc_4, nlc_5, nlc_6))
    print(result)
    theta_2 = np.linspace(math.degrees(result.x[-2]), math.degrees(result.x[-1]), num=100)
    # theta_2 = np.linspace(21, 21, num=1)
    params = result.x[:-1]
    mech = Four_bar(params)
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.15, 0.15), ylim=(-0.1, 0.2))
    ax.set_aspect('equal')
    ax.grid()
    line_1, = ax.plot([], [], 'o-', lw=2)
    line_2, = ax.plot([], [], 'o-', lw=2)
    line_3, = ax.plot([], [], 'o-', lw=2)
    ax.plot(points_desired[0][0], points_desired[0][1], "or")
    ax.plot(points_desired[1][0], points_desired[1][1], "or")
    ani = animation.FuncAnimation(fig, animate, len(theta_2), fargs=[mech, theta_2, line_1, line_2, line_3],
                                  interval=50)
    plt.show()
