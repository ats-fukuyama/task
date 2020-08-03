import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn


class plot:
    def __init__(self, path, name, type):
        self.name = name
        self.val = np.loadtxt(path, delimiter=',')
        self.raxis = np.arange(0, self.val.shape[0])*1.0/self.val.shape[0]
        self.zaxis = np.arange(0, self.val.shape[1])*1.0/self.val.shape[1]
        self.type = type
        # type = 0 : 等値線図
        # type = 1 : 3次元wire
        # type = 2 : 3次元surface


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
                m1 = j+1
                m2 = j+1
                break
            j += 1

        fig = plt.figure()
        ax = []

        for i in range(n):
            x, y = np.meshgrid(z[i].raxis, z[i].zaxis)
            if z[i].type == 0:
                ax.append(fig.add_subplot(m1, m2, i+1))
                ax[i].contour(x, y, z[i].val, levels=10)
                plt.xlabel("R")
                plt.ylabel("Z")

            elif z[i].type == 1:
                ax.append(fig.add_subplot(m1, m2, i+1, projection="3d"))
                ax[i].plot_wireframe(x, y, z[i].val)
                plt.xlabel("R")
                plt.ylabel("Z")

            else:
                ax.append(fig.add_subplot(m1, m2, i+1, projection="3d"))
                ax[i].plot_surface(x, y, z[i].val)
                plt.xlabel("R")
                plt.ylabel("Z")

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

            else:
                ax.plot_surface(x, y, z[i].val)

        plt.xlabel("R")
        plt.ylabel("Z")
        plt.show()


CSVDIR = "/Users/ota/git/task/fp.ota/csv/"

trapped = plot(CSVDIR+"check_orbit_trapped.csv", r"$trapped$", 0)
forbitten = plot(CSVDIR+"check_orbit_forbitten.csv", r"$forbitten$", 0)

z = [trapped, forbitten]
option = [0]
plot_profile2D(z, option)
