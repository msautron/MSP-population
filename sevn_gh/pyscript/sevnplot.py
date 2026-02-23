import numpy as np
from argparse import ArgumentParser
from plotdue import VeuzBaseDoc, rainbowlist
import os
from amuse.units import units, constants, nbody_system
import itertools

timelab = "Worldtime_0"
m1lab = "Mass_0"
m2lab = "Mass_1"
r1lab = "Radius_0"
r2lab = "Radius_1"
alab = "Semimajor"
elab = "Eccentricity"
ph1lab = "Phase_0"
ph2lab = "Phase_1"
spin1lab = "nd"
spin2lab = "nd"
inert1lab = "Inertia_0"
inert2lab = "Inertia_1"
lum1lab = "Luminosity_0"
lum2lab = "Luminosity_1"
temp1lab = "Temperature_0"
temp2lab = "Temperature_1"

class PlotTrack():
    def __init__(self, fname, idbin, task, it_Myr, ft_Myr):
        self.fullpath = fname
        self.task = task
        self.simname = os.path.basename(fname)
        self.foldname = os.path.dirname(fname)
        self.idbin = idbin

        idcolumn = np.loadtxt(self.fullpath, usecols=(0,))
        totrows = len(idcolumn)
        mycol = np.argwhere(idcolumn == self.idbin).flatten()
        idstart, idend = mycol[0], mycol[-1]

        self.tracktab = np.genfromtxt(self.fullpath, names=True, skip_header=idstart, skip_footer=totrows-idend-1)

        time = self.tracktab[timelab]
        tott = len(time)
        ft_Myr = time[-1] if ft_Myr == -1 else ft_Myr
        iind = np.argwhere(time > it_Myr).flatten()
        find = np.argwhere(time > ft_Myr).flatten()

        if len(iind) == 0:
            iind = 0
        else:
            iind = iind[0]
        if len(find) == 0:
            find = tott - 1
        else:
            find = find[0]

        if iind == find == 0: find = 1
        if iind == find == tott: iind = tott - 1

        print("it, ft = ", iind, find)
        self.iind = iind
        self.find = find

        # filter
        self.tracktab = self.tracktab[self.iind:self.find]

    def plot_track(self, plotter, nbin=0):
        t = self.tracktab[timelab]
        m1 = self.tracktab[m1lab]
        m2 = self.tracktab[m2lab]
        r1 = self.tracktab[r1lab]
        r2 = self.tracktab[r2lab]
        a = self.tracktab[alab]
        e = self.tracktab[elab]
        ph1 = self.tracktab[ph1lab]
        ph2 = self.tracktab[ph2lab]
        spin1 = np.zeros(len(t))  #self.tracktab[spin1lab]
        spin2 = np.zeros(len(t))  #self.tracktab[spin2lab]
        inert1 = self.tracktab[inert1lab]
        inert2 = self.tracktab[inert2lab]

        orb = (constants.G*(m1+m2 | units.MSun)/(a | units.RSun)**3).sqrt()
        orb = orb.value_in(units.yr**-1)

        e2 = e*e
        orbsync = orb * (1 + 7.5*e2 + 5.625*e2*e2 + 0.3125*e2*e2*e2)
        orbsync = orbsync / (1 + 3*e2 + 0.375*e2*e2)
        orbsync = orbsync / (1 - e2)**1.5

        cl = rainbowlist   #["blue", "red", "green", "yellow", "black"]
        Nc = len(cl)
        c = [cl[nbin % Nc]] * 2
        ls = ['solid', "dashed"]

        ig=0
        plotter.lines([t, t], [a, a*(1-e)], size=u'2pt', line=ls, color=c, igrid=0, igraph=ig)

        ig=1
        plotter.lines([t, t], [r1, r2], size=u'2pt', line=ls, color=c, igrid=0, igraph=ig)

        ig=2
        plotter.lines([t, t], [m1, m2], size=u'2pt', line=ls, color=c, igrid=0, igraph=ig)

        ig=3
        plotter.lines([t, t], [ph1, ph2], size=u'2pt', line=ls, color=c, igrid=0, igraph=ig)

        ig=4
        plotter.lines([t, t], [inert1, inert2], size=u'2pt', line=ls, color=c, igrid=0, igraph=ig)

        ig=5
        plotter.lines([t, t], [spin1, spin2], size=u'2pt', line=ls, color=c, igrid=0, igraph=ig)

    def plot_track_finalize(plotter):
        ig=0
        plotter.xaxis(u'', 0, igrid=0, igraph=ig, grid=True)
        plotter.yaxis(u'a, a(1-e) [R_{\odot}]', igrid=0, igraph=ig, log=True, grid=True)
        plotter.get(igrid=0, igraph=ig).bottomMargin.val = '0cm'
        plotter.get(igrid=0, igraph=ig).leftMargin.val = '2.5cm'

        ig=1
        plotter.xaxis(u'', 0, igrid=0, igraph=ig, grid=True)
        plotter.yaxis(u'Radius [R_{\odot}]', igrid=0, igraph=ig, grid=True, log=True)
        plotter.get(igrid=0, igraph=ig).bottomMargin.val = '0cm'
        plotter.get(igrid=0, igraph=ig).leftMargin.val = '2.5cm'

        ig=2
        plotter.xaxis(u'', 0, igrid=0, igraph=ig, grid=True)
        plotter.yaxis(u'Mass [M_{\odot}]', igrid=0, igraph=ig, grid=True)
        plotter.get(igrid=0, igraph=ig).topMargin.val = '0cm'
        plotter.get(igrid=0, igraph=ig).leftMargin.val = '2.5cm'
        plotter.get(igrid=0, igraph=ig).bottomMargin.val = '0cm'

        ig=3
        plotter.xaxis(u'', 0, igrid=0, igraph=ig, grid=True)
        plotter.yaxis(u'Stellar phase', igrid=0, igraph=ig, grid=True)
        plotter.get(igrid=0, igraph=ig).topMargin.val = '0cm'
        plotter.get(igrid=0, igraph=ig).leftMargin.val = '2.5cm'
        plotter.get(igrid=0, igraph=ig).bottomMargin.val = '0cm'

        ig=4
        plotter.xaxis(u'Time [Myr]', 0, igrid=0, igraph=ig, grid=True)
        plotter.yaxis(u'Inertia [R_{\odot}^{2} M_{\odot}^{2}]', igrid=0, igraph=ig, grid=True, log=True)
        plotter.get(igrid=0, igraph=ig).topMargin.val = '0cm'
        plotter.get(igrid=0, igraph=ig).leftMargin.val = '2.5cm'

        ig=5
        plotter.xaxis(u'Time [Myr]', 0, igrid=0, igraph=ig, grid=True)
        plotter.yaxis(u'Spin [yr^{-1}]', igrid=0, igraph=ig, grid=True)
        plotter.get(igrid=0, igraph=ig).topMargin.val = '0cm'
        plotter.get(igrid=0, igraph=ig).leftMargin.val = '2.5cm'

    def plot_HR(self, plotter, nbin=0):
        t = self.tracktab[timelab]
        temp1 = self.tracktab[temp1lab]
        temp2 = self.tracktab[temp2lab]
        lum1 = self.tracktab[lum1lab]
        lum2 = self.tracktab[lum2lab]

        cl = rainbowlist   #["blue", "red", "green", "yellow", "black"]
        Nc = len(cl)
        c = [cl[nbin % Nc]] * 2
        ls = ['solid', "dashed"]

        ig=0
        plotter.lines([temp1, temp2], [lum1, lum2], size=u'2pt', line=ls, color=c, igrid=0, igraph=ig)

    def plot_HR_finalize(plotter):
        ig=0
        plotter.xaxis(u'Temperature [K]', igrid=0, igraph=ig, grid=True, log=True)
        plotter.yaxis(u'Luminosity [L_ {\odot}]', igrid=0, igraph=ig, log=True, grid=True)


def plotter_choice(task, plotter, close=True):
    if 'show' in task:
        plotter.show()
    if 'png' in task:
        # plotter.save()
        plotter.export_png(".")
    if 'pdf' in task:
        # plotter.save()
        plotter.export_pdf(".")
    if 'vsz' in task:
        plotter.save(".")
        plotter.close()

def new_argument_parser():
    result = ArgumentParser()
    result.add_argument("fname", nargs="+",
                        help="trackfile", type=str)
    result.add_argument("-task", "--t", dest="task", default=['show'], nargs='+', type=str,
                        help="task for output values")
    result.add_argument("-id", dest="idbin", default=[0], nargs="+",
                        type=int, help="binary id")
    result.add_argument("-k", dest="key", default='aei',
                        type=str, help="plot type")
    result.add_argument("-it", dest="it_Myr", default=0,
                        type=float, help="initial time in Myr")
    result.add_argument("-ft", dest="ft_Myr", default=-1,
                        type=float, help="final time in Myr")
    return result

if __name__ in ('__main__'):
    args = new_argument_parser().parse_args()
    key = args.key

    outname = "binplot"
    plotter = VeuzBaseDoc(outname)
    plotter.add_page(height='35cm', width='35cm')
    plotter.add_grid(rows=3, columns=2)
    plotter.figure.ResizeWindow(780, 700)

    itrack = 0
    for fin in args.fname:
        outname = outname + "_" + os.path.basename(fin).replace(".dat", "")
        for id in args.idbin:
            outname = outname + "_" + str(id)
            sbin = PlotTrack(fin, args.idbin,
                          args.task,
                          args.it_Myr, args.ft_Myr)
            plotter.rename(outname)

            sbin.plot_track(plotter, itrack)
            itrack += 1

    PlotTrack.plot_track_finalize(plotter)
    plotter_choice(args.task, plotter)
