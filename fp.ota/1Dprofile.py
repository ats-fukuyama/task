import numpy as np
import matplotlib.pyplot as plt
import seaborn


class plot:
    def __init__(self, path, name, ymin=0, ymax=0):
        self.name = name
        self.val = np.loadtxt(path, delimiter=',')
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


CSVDIR = "/Users/ota/git/task/fp.ota/csv/"

psim = plot(CSVDIR+"psim.csv", "psim")
Bin = plot(CSVDIR+"Bin.csv", "Bin")
Bout = plot(CSVDIR+"Bout.csv", "Bout")
Fpsi = plot(CSVDIR+"Fpsi.csv", "Fpsi")
psimg = plot(CSVDIR+"psimg.csv", "psimg")
Bing = plot(CSVDIR+"Bing.csv", "Bing")
Boutg = plot(CSVDIR+"Boutg.csv", "Boutg")
Fpsig = plot(CSVDIR+"Fpsig.csv", "Fpsig")


y = [psim, Bin, Bout, Fpsi, psimg, Bing, Boutg, Fpsig]

# plot_profile1D(y, [1])
plot_profile1D(y, [0])
