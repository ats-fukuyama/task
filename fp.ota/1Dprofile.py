import numpy as np
import matplotlib.pyplot as plt
import seaborn


class plot:
    def __init__(self, val, name, ymin=0, ymax=0):
        self.name = name
        self.val = val
        self.raxis = np.arange(0, self.val.shape[0])*1.0/self.val.shape[0]
        self.ymin = ymin
        self.ymax = ymax


def plot_profile1D(y, option):
    # option[0] = 0 : 別のグラフに描写
    # option[0] = 1 : 同じグラフ内に描写
    seaborn.set()
    n = len(y)
    if option[0] == 0:
        plt.figure()
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

        for i in range(n):
            plt.subplot(m1, m2, i+1)
            plt.xlabel('r')
            plt.ylabel(y[i].name)
            plt.plot(y[i].raxis, y[i].val)
            if y[i].ymin != y[i].ymax:
                plt.ylim(y[i].ymin, y[i].ymax)

        plt.show()

    elif option[0] == 1:
        plt.subplot(1, 1, 1)
        plt.xlabel('r')
        for i in range(n):
            plt.plot(y[i].raxis, y[i].val, linewidth=2, label=y[i].name)
            if y[i].ymin != y[i].ymax:
                plt.ylim(y[i].ymin, y[i].ymax)

        plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0)
        plt.show()


CSVDIR = "/Users/ota/git/task/fp/csv/"

psi = plot(np.loadtxt(CSVDIR+"psi.csv", delimiter=','), r"$\psi$")
# psitmp = plot(np.loadtxt(CSVDIR+"psitmp.csv", delimiter=','), r"$\psi_{tmp}$")
F = plot(np.loadtxt(CSVDIR+"F.csv", delimiter=','), r"$F$", 0, 35)
# Ftmp = plot(np.loadtxt(CSVDIR+"Ftmp.csv", delimiter=','), r"$F_{tmp}$")
# f = plot(np.loadtxt(CSVDIR+"f.csv", delimiter=','), r"$f$")
# f1 = plot(np.loadtxt(CSVDIR+"f1.csv", delimiter=','), r"$f1$")
# f2 = plot(np.loadtxt(CSVDIR+"f2.csv", delimiter=','), r"$f2$")
# ftest = plot(np.loadtxt(CSVDIR+"ftest.csv", delimiter=','), r"$f$")
# gtest = plot(np.loadtxt(CSVDIR+"gtest.csv", delimiter=','), r"$g$")
Bout = plot(np.loadtxt(CSVDIR+"Bout.csv", delimiter=','), r"$B_{out}$")
Bin = plot(np.loadtxt(CSVDIR+"Bin.csv", delimiter=','), r"$B_{in}$")
psieq = plot(np.loadtxt(CSVDIR+"psieq.csv", delimiter=','), r"$\psi_{eq}$")
Beq = plot(np.loadtxt(CSVDIR+"Beq.csv", delimiter=','), r"$B_{eq}$")


y = [Bin, Bout, psieq, Beq, F, psi]

# plot_profile1D(y, [1])
plot_profile1D(y, [0])
