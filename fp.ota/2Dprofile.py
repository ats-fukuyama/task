import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn


class plot:
    def __init__(self, val, name, type, val_contour=0, zmin=0, zmax=0):
        self.name = name
        self.val = val
        self.raxis = np.arange(0, self.val.shape[0])*1.0/self.val.shape[0]
        self.zaxis = np.arange(0, self.val.shape[1])*1.0/self.val.shape[1]
        self.type = type
        # type = 0 : conter plot in 2D-plane, R-Z
        # type = 1 : 3D wire
        # type = 2 : 3D surface
        # type = 3 : 3D contour for G(R,Z)=val_contour
        self.val_contour = val_contour
        self.zmin = zmin
        self.zmax = zmax


def plot_profile2D(z, option):
    # option[0] = 0 : len(z)枚のグラフに描写，グラフ別居
    # option[0] = 1 : 1枚のグラフ内に描写，グラフ同居

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
                m2 = j+1
                for k in range(j+2):
                    if k*m2 >= n:
                        m1 = k
                        break
                break
            j += 1

        fig = plt.figure()
        ax = []

        for i in range(n):
            x, y = np.meshgrid(z[i].raxis, z[i].zaxis)
            if z[i].type == 0:
                ax.append(fig.add_subplot(m1, m2, i+1))
                ax[i].contour(x, y, z[i].val, levels=10)
                ax[i].set_xlabel("R")
                ax[i].set_ylabel("Z")

            elif z[i].type == 1:
                ax.append(fig.add_subplot(m1, m2, i+1, projection="3d"))
                ax[i].plot_wireframe(x, y, z[i].val)
                ax[i].set_xlabel("R")
                ax[i].set_ylabel("Z")
                ax[i].set_zlabel(z[i].name)

            elif z[i].type == 2:
                ax.append(fig.add_subplot(m1, m2, i+1, projection="3d"))
                ax[i].plot_surface(x, y, z[i].val)
                ax[i].set_xlabel("R")
                ax[i].set_ylabel("Z")
                ax[i].set_zlabel(z[i].name)

            else:
                ax.append(fig.add_subplot(m1, m2, i+1, projection="3d"))
                ax[i].contour(x, y, z[i].val, [z[i].val_contour])
                ax[i].set_xlabel("R")
                ax[i].set_ylabel("Z")
                ax[i].set_zlabel(z[i].name)

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

            else:
                ax.contour(x, y, z[i].val, [z[i].val_contour])

        plt.xlabel("R")
        plt.ylabel("Z")
        plt.show()


CSVDIR = "/Users/ota/git/task/fp/csv/"

# B = plot(np.loadtxt(CSVDIR+"Brz.csv", delimiter=','), r"$B$",1)
# psirztmp = plot(np.loadtxt(CSVDIR+"psirztmp.csv", delimiter=','), r"$\psi rztmp$", 1)
# psirzfp = plot(np.loadtxt(CSVDIR+"psirzfp.csv", delimiter=','), r"$\psi rzfp$", 1)
# f = plot(np.loadtxt(CSVDIR+"f.csv", delimiter=','), r"$f$", 2)
# f2 = plot(np.loadtxt(CSVDIR+"f2.csv", delimiter=','), r"$f2$", 2)
Brz = plot(np.loadtxt(CSVDIR+"Brz.csv", delimiter=','), r"$B_{rz}$", 1)
psirz = plot(np.loadtxt(CSVDIR+"psirz.csv", delimiter=','), r"$\Psi_{RZ}$", 0)
Frz = plot(np.loadtxt(CSVDIR+"Frz.csv", delimiter=','), r"$F_{rz}$", 0)

psirz1 = plot(np.loadtxt(CSVDIR+"psirz.csv", delimiter=','), r"$\Psi_{RZ}$", 1)
Frz1 = plot(np.loadtxt(CSVDIR+"Frz.csv", delimiter=','), r"$F_{rz}$", 1)

z = [Brz, Frz, psirz, Frz1, psirz1]
option = [0]
plot_profile2D(z, option)
