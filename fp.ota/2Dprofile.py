import numpy as np
from scipy.ndimage import map_coordinates
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn


def polar2cartesian(r, t, grid, x, y, order=3):

    X, Y = np.meshgrid(x, y)

    new_r = np.sqrt(X*X+Y*Y)
    new_t = np.arctan2(X, Y)

    ir = interp1d(r, np.arange(len(r)), bounds_error=False)
    it = interp1d(t, np.arange(len(t)), kind='cubic')

    new_ir = ir(new_r.ravel())
    new_it = it(new_t.ravel())

    new_ir[new_r.ravel() > r.max()] = len(r)-1
    new_ir[new_r.ravel() < r.min()] = 0

    return map_coordinates(grid, np.array([new_ir, new_it]),
                           order=order).reshape(new_r.shape)


class plot:
    def __init__(self, path, name, type, xname="R", yname="Z"):
        self.name = name
        self.val = np.loadtxt(path, delimiter=',')
        self.raxis = np.arange(
            0, self.val.shape[0])*1.0/self.val.shape[0]+0.50/self.val.shape[0]
        self.zaxis = np.arange(
            0, self.val.shape[1])*1.0/self.val.shape[1]+0.50/self.val.shape[1]
        self.type = type
        self.xname = xname
        self.yname = yname
        # type = 0 : 等値線図
        # type = 1 : 3次元wire
        # type = 2 : 3次元surface
        # type = 3 : 等値線図 color map


class polar(plot):
    def __init__(self, path, name, type, xname=r"$p_\parallel$", yname=r"$p_\perp$", mode=0):
        super().__init__(path, name, type, xname, yname)
        self.zaxis_org = self.zaxis
        self.raxis_org = self.raxis
        self.val_org = self.val

        if mode == 0:
            self.costheta = 1.0-2*self.zaxis
            self.sintheta = (1.0-self.costheta**2)**0.5
            self.theta = np.arccos(self.costheta)
        elif mode == 1:
            self.theta = self.zaxis_org

        self.zaxis = np.linspace(0, 1, 50)
        self.raxis = np.linspace(-1, -1, 100)

        self.val = polar2cartesian(
            self.raxis_org, self.theta, self.val_org, self.zaxis, self.raxis, order=3)


def plot_profile2D(z, option):
    # option[0] = 0 : len(z)枚のグラフに描写
    # option[0] = 1 : 1枚のグラフ内に描写

    seaborn.set()
    n = len(z)
    if option[0] == 0:
        j = 1
        while True:
            if n <= 2:
                m1 = 1
                m2 = n
                break
            elif j**2 < n and n <= (j+1)**2:
                m1 = j+1
                m2 = j+1
                break
            j += 1

        fig = plt.figure()
        ax = []

        for i in range(n):
            x, y = np.meshgrid(z[i].zaxis, z[i].raxis)
            if z[i].type == 0:
                ax.append(fig.add_subplot(m1, m2, i+1))
                ax[i].contour(x, y, z[i].val, levels=10)
                plt.xlabel(z[i].xname)
                plt.ylabel(z[i].yname)

            elif z[i].type == 1:
                ax.append(fig.add_subplot(m1, m2, i+1, projection="3d"))
                ax[i].plot_wireframe(x, y, z[i].val)
                plt.xlabel(z[i].xname)
                plt.ylabel(z[i].yname)

            elif z[i].type == 2:
                ax.append(fig.add_subplot(m1, m2, i+1, projection="3d"))
                ax[i].plot_surface(x, y, z[i].val)
                plt.xlabel(z[i].xname)
                plt.ylabel(z[i].yname)

            elif z[i].type == 3:
                ax.append(fig.add_subplot(m1, m2, i+1))
                contour = ax[i].contourf(x, y, z[i].val, cmap="bwr")
                plt.xlabel(z[i].xname)
                plt.ylabel(z[i].yname)
                fig.colorbar(contour)

            plt.title(z[i].name)

        plt.show()

    elif option[0] == 1:
        fig = plt.figure()
        a = 0
        for i in range(n):
            a += z[i].type
        if a == 0:
            ax = fig.add_subplot(1, 1, 1)
        else:
            ax = fig.add_subplot(1, 1, 1, projection="3d")

        for i in range(n):
            x, y = np.meshgrid(z[i].raxis, z[i].zaxis)
            if z[i].type == 0:
                ax.contour(x, y, z[i].val, levels=10)

            elif z[i].type == 1:
                ax.plot_wireframe(x, y, z[i].val)

            elif z[i].type == 2:
                ax.plot_surface(x, y, z[i].val)

            elif z[i].type == 3:
                ax.contourf(x, y, z[i].val)

        plt.xlabel(z[0].xname)
        plt.ylabel(z[0].yname)
        plt.show()


if __name__ == "__main__":

    CSVDIR = "/Users/ota/git/task/fp.ota/csv/"
    theta_m = r"$\theta_m$"
    psi_m = r"$\psi_m$"
    p = "p"
    xi = r"$\xi=\cos\theta_m$"
    gtypez1 = 2
    gtypez2 = 3

    trapped_ion = plot(CSVDIR+"check_orbit_trapped_ion.csv",
                       r"$trapped_ion$", gtypez1, xi, psi_m)
    forbitten_ion = plot(CSVDIR+"check_orbit_forbitten_ion.csv",
                         r"$forbitten_ion$", gtypez1, xi, psi_m)
    trapped_ele = plot(CSVDIR+"check_orbit_trapped_ele.csv",
                       r"$trapped_ele$", gtypez1, xi, psi_m)
    forbitten_ele = plot(CSVDIR+"check_orbit_forbitten_ele.csv",
                         r"$forbitten_ele$", gtypez1, xi, psi_m)

    # fu = plot(CSVDIR+"fu.csv", r"$F_u$", 3, "p", r"$\theta$")
    fI_center = plot(CSVDIR+"fI_center.csv",
                     r"$f_I (nr=1)$", gtypez2, p, xi)
    fI_edge = plot(CSVDIR+"fI_edge.csv", r"$f_I( nr=nrmax)$", gtypez2, p, xi)
    fI_quarter = plot(CSVDIR+"fI_quarter.csv",
                      r"$f_I( nr=nrmax/2)$", gtypez2, p, xi)
    J = plot(CSVDIR+"J.csv", r"$J_I$", 0, "p", r"$\theta$")

    # fI_center_p = polar(CSVDIR+"fI_center.csv",
    #                     r"$f_I (nr=1)$", gtypez2, p, xi)

    # print(fI_center_p.theta)

    z1 = [trapped_ele, forbitten_ele, trapped_ion, forbitten_ion]
    z2 = [fI_center, fI_edge, fI_quarter]
    # zp = [fI_center_p]
    option = [0]
    # plot_profile2D(z1, option)
    plot_profile2D(z2, option)
