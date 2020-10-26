import seaborn
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import plot1D as p1D


class plot2D_org:
    def __init__(self, path, zname, type, xname="R", yname="Z", xmin=None, xmax=None, ymin=None, ymax=None):
        self.name = "plot2D_org"
        self.zname = zname
        self.val = np.loadtxt(path, delimiter=',')
        self.xaxis = np.arange(
            0, self.val.shape[0])*1.0/self.val.shape[0]+0.50/self.val.shape[0]
        self.yaxis = np.arange(
            0, self.val.shape[1])*1.0/self.val.shape[1]+0.50/self.val.shape[1]
        self.type = type
        self.xname = xname
        self.yname = yname
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        # type = 0 : 等値線図
        # type = 1 : 3次元wire
        # type = 2 : 3次元surface
        # type = 3 : 等値線図 color map


class plot2D(plot2D_org):
    def __init__(self, path, zname, type, xaxis_path, yaxis_path, list_auxiliary_line_path=[], xname=r"$\xi=\cos\theta_m$", yname=r"$\bar{p}$", xmin=None, xmax=None, ymin=None, ymax=None):

        super().__init__(path, zname, type, xname, yname, xmin, xmax, ymin, ymax)
        self.name = "plot2D"
        self.xaxis = np.loadtxt(xaxis_path, delimiter=',')
        self.yaxis = np.loadtxt(yaxis_path, delimiter=',')
        self.auxiliary_line = []

        if list_auxiliary_line_path != []:
            i = 1
            for A in list_auxiliary_line_path:
                self.auxiliary_line.append(
                    p1D.plot1D(A, "auxiliary_line"+str(i), xaxis_path))
                i += 1

    def polar2cartesian(self, mode=0):

        if self.auxiliary_line != []:
            for A in self.auxiliary_line:
                A.polar2cartesian(mode)

        size_r = self.yaxis.shape[0]
        x = np.linspace(-self.yaxis[-1], self.yaxis[-1], 2*size_r)
        y = np.linspace(0.0, self.yaxis[-1], size_r)
        func_znew = interpolate.interp2d(
            self.xaxis, self.yaxis, self.val.T, kind='cubic')
        z = np.empty((2*size_r, size_r))

        if mode == 0:
            for ix in range(0, 2*size_r):
                for iy in range(0, size_r):
                    r = (x[ix]**2+y[iy]**2)**0.5
                    if r >= self.yaxis[-1]:
                        z[ix, iy] = 0.0
                    else:
                        costh = x[ix]/(x[ix]**2+y[iy]**2)**0.5
                        th = np.arccos(costh)
                        z[ix, iy] = func_znew(th, r)

        else:
            for ix in range(0, 2*size_r):
                for iy in range(0, size_r):
                    r = (x[ix]**2+y[iy]**2)**0.5
                    if r >= self.yaxis[-1]:
                        z[ix, iy] = 0.0
                    else:
                        costh = x[ix]/(x[ix]**2+y[iy]**2)**0.5
                        z[ix, iy] = func_znew(costh, r)

        self.xaxis = x
        self.yaxis = y
        self.val = z
        self.xname = r"$p_\parallel$"
        self.yname = r"$p_{\perp}$"

        return


class plot2D_x_depende_on_y(plot2D):
    def __init__(self, path, zname, type, xaxis_path, yaxis_path, list_auxiliary_line_path2=[], xname2=r"$\theta_m$", yname2=r"$\bar{p}$", xmin2=None, xmax2=None, ymin2=None, ymax2=None):
        super().__init__(path, zname, type, xaxis_path, yaxis_path, list_auxiliary_line_path=list_auxiliary_line_path2,
                         xname=xname2, yname=yname2, xmin=xmin2, xmax=xmax2, ymin=ymin2, ymax=ymax2)
        self.name = "plot2D_x_depende_on_y"


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
            if z[i].name == "plot2D_x_depende_on_y":
                size_ = z[i].xaxis.shape[0]
                dummy = np.linspace(0, 1, size_)
                xx, y = np.meshgrid(dummy, z[i].yaxis, indexing='ij')
                x = z[i].xaxis
            else:
                x, y = np.meshgrid(z[i].xaxis, z[i].yaxis, indexing='ij')

            if z[i].type == 0:
                ax.append(fig.add_subplot(m1, m2, i+1))
                ax[i].contour(x, y, z[i].val, levels=10, cmap="bwr")
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
                # contour = ax[i].contourf(
                #     x1, y1, z[i].val.T, cmap="bwr", norm=Normalize(vmin=0, vmax=0.032))
                plt.xlabel(z[i].xname)
                plt.ylabel(z[i].yname)
                fig.colorbar(contour)
                # contour.set_clim(0, 0.032)
                # plt.axvline(x=0, color='black')

            ax[i].set_xlim(left=z[i].xmin, right=z[i].xmax)
            ax[i].set_ylim(bottom=z[i].ymin, top=z[i].ymax)

            if z[i].name == 'plot2D':
                if z[i].auxiliary_line != []:
                    j = 0
                    for A in z[i].auxiliary_line:
                        if j % 3 == 0:
                            l_style = 'dashed'
                        elif j % 3 == 1:
                            l_style = 'dashdot'
                        else:
                            l_style = 'dotted'
                        j += 1

                        plt.plot(A.xaxis, A.val, color='black',
                                 linestyle=l_style)

            plt.title(z[i].zname)

        plt.rcParams['mathtext.fontset'] = 'stix'
        plt.tight_layout()
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
            x, y = np.meshgrid(z[i].xaxis, z[i].yaxis, indexing='ij')
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
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":

    CSVDIR = "/Users/ota/git/task/fp.ota/csv/"
    # theta_m = r"$\theta_m$"
    # psi_m = r"$\psi_m$"
    # p = "p"
    # xi = r"$\xi=\cos\theta_m$"

    # # fu = plot2D_org(CSVDIR+"fu.csv", r"$F_u$", 3, "p", r"$\theta$")
    # # fI_center = plot2D(CSVDIR+"fI_center.csv", r"$f_I (nr=2)$",
    # #                   gtypez2, CSVDIR+"trapped_boundary_center.csv", CSVDIR+"forbitten_boundary_center.csv", CSVDIR+"xi.csv", CSVDIR+"pm_ion.csv")
    # # fI_edge = plot2D(CSVDIR+"fI_edge.csv", r"$f_I( nr=nrmax)$",
    # #                 gtypez2, CSVDIR+"trapped_boundary_edge.csv", CSVDIR+"forbitten_boundary_edge.csv", CSVDIR+"xi.csv", CSVDIR+"pm_ion.csv")
    # # fI_quarter = plot2D(CSVDIR+"fI_quarter.csv", r"$f_I( nr=nrmax/2)$",
    # #                    gtypez2, CSVDIR+"trapped_boundary_quarter.csv", CSVDIR+"forbitten_boundary_quarter.csv", CSVDIR+"xi.csv", CSVDIR+"pm_ion.csv")

    aux_c = [CSVDIR+"trapped_boundary_center.csv",
             CSVDIR+"forbitten_boundary_center.csv"]
    aux_q = [CSVDIR+"trapped_boundary_quarter.csv",
             CSVDIR+"forbitten_boundary_quarter.csv"]
    aux_e = [CSVDIR+"trapped_boundary_edge.csv",
             CSVDIR+"forbitten_boundary_edge.csv"]
    gtype = 3

    # mean_r_center = plot2D(CSVDIR+"mean_r_center.csv",
    #                        r"${r/a_{max}-\langle r/a\rangle}\ \ (nr=2)$", gtype,  CSVDIR+"xi.csv", CSVDIR+"pm_ion.csv", list_auxiliary_line_path=aux_c)
    # mean_r_quarter = plot2D(CSVDIR+"mean_r_quarter.csv",
    #                         r"${r/a_{max}-\langle r/a\rangle}\ \ (nr=nrmax/2)$", gtype, CSVDIR+"xi.csv", CSVDIR+"pm_ion.csv", list_auxiliary_line_path=aux_q)
    # mean_r_edge = plot2D(CSVDIR+"mean_r_edge.csv",
    #                      r"${r/a_{max}-\langle r/a\rangle}\ \ (nr=nrmax)$", gtype, CSVDIR+"xi.csv", CSVDIR+"pm_ion.csv", list_auxiliary_line_path=aux_e)
    # fI = plot2D_x_depende_on_y(CSVDIR+"fI.csv", r"$f_I$", gtype,
    #                            CSVDIR+"thetam.csv", CSVDIR+"pm_ion.csv")
    # fu = plot2D_x_depende_on_y(CSVDIR+"fu.csv", r"$f_u$", gtype,
    #                            CSVDIR+"thetam.csv", CSVDIR+"pm_ion.csv", xname2=r"$\theta$")
    # plot_profile2D([fI, fu], [0])
    # mean_r_quarter.polar2cartesian()
    # # J = plot2D_org(CSVDIR+"J.csv", r"$J_I$", 0, "p", r"$\theta$")

    # # z1 = [trapped_ele, forbitten_ele, trapped_ion, forbitten_ion]
    # # z2 = [fI_center, fI_edge, fI_quarter]
    # # z3 = [mean_r_center, mean_r_quarter, mean_r_edge]
    # z3 = [mean_r_quarter]
    # # zp = [fI_center_p]
    # option = [0]
    # # plot_profile2D(z1, option)
    # # plot_profile2D(z2, option)
    # plot_profile2D(z3, option)
