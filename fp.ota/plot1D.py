import numpy as np
import matplotlib.pyplot as plt
import seaborn


class plot1D:
    def __init__(self, path, yname, xpath=None, xname='r/a', xmin=None, xmax=None, ymin=None, ymax=None):
        self.yname = yname
        self.val = np.loadtxt(path, delimiter=',')
        self.xname = xname
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        if xpath == None:
            self.xaxis = np.arange(0, self.val.shape[0])*1.0/self.val.shape[0]
        else:
            self.xaxis = np.loadtxt(xpath, delimiter=',')

    def polar2cartesian(self, mode=0):
        if mode == 0:
            th = self.xaxis
            costh = np.cos(th)
            sinth = np.sin(th)
        else:
            costh = self.xaxis
            sinth = np.sqrt(1-self.xaxis**2)

        x = self.val*costh
        y = self.val*sinth

        self.xaxis = x
        self.val = y
        self.xname = r"$v_\parallel/c$"

        return self


def plot_profile1D(y, option):
    # option[0] = 0 : 別のグラフに描写
    # option[0] = 1 : 同じグラフ内に描写, if option[0] == 1, then option[1] = y-axis name
    # seaborn.set()
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
            plt.xlabel(y[i].xname)
            plt.ylabel(y[i].yname)

            plt.plot(y[i].xaxis, y[i].val, linewidth=0.9)

            plt.grid(b=True, which='major', color='#666666', linestyle='-')
            plt.minorticks_on()
            plt.grid(b=True, which='minor', color='#999999',
                     linestyle='-', alpha=0.2)
            plt.xlim(y[i].xmin, y[i].xmax)
            plt.ylim(y[i].ymin, y[i].ymax)

        plt.tight_layout()
        plt.show()

    elif option[0] == 1:
        plt.subplot(1, 1, 1)
        plt.xlabel(y[0].xname)
        plt.ylabel(option[1])
        for i in range(n):
            if i == 2 or i == 4:
                c = 'lime'
            elif i == 3 or i == 5:
                c = 'black'
            elif i == 0:
                c = 'blue'
            elif i == 1:
                c = 'red'
            plt.plot(y[i].xaxis, y[i].val, linewidth=1.4,
                     label=y[i].yname, color=c)

        plt.xlim(y[0].xmin, y[0].xmax)
        plt.ylim(y[0].ymin, y[0].ymax)
        plt.grid(b=True, which='major', color='#666666', linestyle='-')
        plt.minorticks_on()
        plt.grid(b=True, which='minor', color='#999999',
                 linestyle='-', alpha=0.2)
        # plt.legend(bbox_to_anchor=(1, 1.18),
        #            loc='upper right', borderaxespad=0)
        plt.legend(bbox_to_anchor=(1, 1.1), ncol=2,
                   loc='upper right', borderaxespad=0)
        plt.rcParams['mathtext.fontset'] = 'stix'
        plt.fill_between(y[0].xaxis, 0, y[0].val, facecolor='aqua', alpha=0.3)
        plt.fill_between(y[2].xaxis, 0, y[2].val,
                         facecolor='slategray', alpha=0.3)
        plt.fill_between(y[3].xaxis, 0, y[3].val,
                         facecolor='slategray', alpha=0.3)
        plt.fill_between(y[4].xaxis, 0, y[4].val,
                         facecolor='purple', alpha=0.3)
        plt.fill_between(y[5].xaxis, 0, y[5].val,
                         facecolor='purple', alpha=0.3)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":

    CSVDIR = "/Users/ota/git/task/fp.ota/csv/"

    # psimg = plot1D(CSVDIR+"psimg.csv", "psimg")
    # Bing = plot1D(CSVDIR+"Bing.csv", "Bing")
    # Boutg = plot1D(CSVDIR+"Boutg.csv", "Boutg")
    # Fpsig = plot1D(CSVDIR+"Fpsig.csv", "Fpsig")

    psim = plot1D(CSVDIR+"psim.csv", r"$\psi_m$")
    Bin = plot1D(CSVDIR+"Bin.csv", r"$B_{in}$", CSVDIR+"psim.csv", "psim")
    Bout = plot1D(CSVDIR+"Bout.csv", r"$B_{out}$", CSVDIR+"psim.csv", "psim")
    Fpsi = plot1D(CSVDIR+"Fpsi.csv", r"$F_m$", CSVDIR+"psim.csv", "psim")

    dfdpsi = plot1D(CSVDIR+"dFdpsi.csv", r"$df/d\psi_m$",
                    CSVDIR+"psim.csv", "psim")
    d2fdpsi = plot1D(CSVDIR+"d2Fdpsi.csv", r"$d^2f/d\psi_m^2$",
                     CSVDIR+"psim.csv", "psim")

    dbindpsi = plot1D(CSVDIR+"dBmdpsi_in.csv", r"$dB_{in}/d\psi_m$",
                      CSVDIR+"psim.csv", "psim")

    dboutdpsi = plot1D(CSVDIR+"dBmdpsi_out.csv",
                       r"$dB_{out}/d\psi_m$", CSVDIR+"psim.csv", r"$\psi_{m}$", ymin=-0.1, ymax=0)
    d2boutdpsi = plot1D(CSVDIR+"d2Bmdpsi_out.csv",
                        r"$d^2B_{out}/d\psi_m^2$", CSVDIR+"psim.csv", r"$\psi_{m}$", ymin=0, ymax=0.1)

    dboutgdpsi = plot1D(CSVDIR+"dBoutgdpsi.csv",
                        r"$dB_{out}/d\psi_m$", CSVDIR+"psimg.csv", r"$\psi_{m}$", ymin=-0.1, ymax=0)
    d2boutgdpsi = plot1D(CSVDIR+"d2Boutgdpsi.csv",
                         r"$d^2B_{out}/d\psi_m^2$", CSVDIR+"psimg.csv", r"$\psi_{m}$", ymin=0, ymax=0.1)

    y = []

    trapped = plot1D(CSVDIR+"trapped_boundary_quarter.csv",
                     "D orbit", CSVDIR+"xi.csv", r"$\xi=\cos\theta_m$")
    forbitten = plot1D(CSVDIR+"forbitten_boundary_quarter.csv",
                       "stagnation orbit", CSVDIR+"xi.csv", r"$\xi$")

    x_co = plot1D(CSVDIR+"beta_x_stagnation_co.csv",
                  "X-type stagnation orbit", CSVDIR+"xi_x_stagnation_co.csv", r"$\xi$")
    o_co = plot1D(CSVDIR+"beta_o_stagnation_co.csv",
                  "O-type stagnation orbit", CSVDIR+"xi_o_stagnation_co.csv", r"$\xi$")
    x_cnt = plot1D(CSVDIR+"beta_x_stagnation_cnt.csv",
                   None, CSVDIR+"xi_x_stagnation_cnt.csv", r"$\xi$")
    o_cnt = plot1D(CSVDIR+"beta_o_stagnation_cnt.csv",
                   None, CSVDIR+"xi_o_stagnation_cnt.csv", r"$\xi$")
    pinch = plot1D(CSVDIR+"beta_pinch.csv",
                   "Pinch orbit", CSVDIR+"xi_pinch.csv", r"$\xi$")

    y.append(trapped)
    y.append(pinch)
    y.append(x_co)
    y.append(o_co)
    y.append(x_cnt)
    y.append(o_cnt)
    plot_profile1D(y, [1, r"$v/c$"])

    # trapped.polar2cartesian()
    # forbitten.polar2cartesian()
    # y.append(trapped)
    # y.append(forbitten)
    # plot_profile1D(y, [1, r"$v_\perp/c$"])
    # y = []
    # y.append(dboutgdpsi)
    # y.append(d2boutgdpsi)
    # plot_profile1D(y, [0])
