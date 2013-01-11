import pickle

class Cut:
    '''Cut to be applied to data.

    Cuts are (min_value, max_value) tuples.

    :param r: Energy cut (MeV)
    :param r: Radius cut (mm)
    :param t: PMT hit time cut (ns)
    '''
    def __init__(self, e=(0.0,10.0), r=(0.0, 10000.0), t=(-500.0, 500.0)):
        assert(len(e) == 2)
        assert(len(r) == 2)
        assert(len(t) == 2)

        self.e = e
        self.r = r
        self.t = t

    def as_suffix(self):
        '''A short string representation.

        E.g. ``_0.0-10.0_0.0-10000.0_-500.0-500.0``

        :returns: Cut as a string
        '''
        dj = lambda x: '-'.join(map(str, x))
        return '_' + '_'.join(map(dj, [self.e, self.r, self.t]))


def load_pdf(filename):
    '''Load a PDF from disk.

    :param filename: Path to PDF pickle
    :returns: PDF as numpy array
    '''
    pdf = None
    with open(filename) as f:
        pdf = pickle.load(f)
    return pdf


def flatten(l):
    '''Flatten a listed list.

    E.g. [[1,2],[3,4]] -> [1,2,3,4]

    http://stackoverflow.com/questions/2158395

    :param l: The nested list
    :returns: Iterable over flattened list
    '''
    import collections
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def plot_pdf_overlay(tlpath, dbpath, suffix='', rebin=None, show=False):
    '''Plot two PDFs overlaid.

    Produces "figures/pdfs[suffix].pdf".

    :param tlpath: Path to 208Tl PDF pickle
    :param dbpath: Path to 0vbb PDF pickle
    :param suffix: Appended to output files
    :param show: If True, show the plot
    '''
    import pdf
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    with open(tlpath) as f:
        tl = pickle.load(f)

    with open(dbpath) as f:
        db = pickle.load(f)

    if len(tl[0]) == 0 or len(db[0]) == 0:
        print 'warning: cannot plot empty pdfs'
        return

    if rebin is not None:
        tl_x = pdf.rebin(tl[:,0], rebin, trim=True)
        tl_y = pdf.rebin(tl[:,1], rebin, trim=True)
        tl = np.array([tl_x, tl_y])

        db_x = pdf.rebin(db[:,0], rebin, trim=True)
        db_y = pdf.rebin(db[:,1], rebin, trim=True)
        db = np.array([db_x, db_y])

    f = plt.figure(1, facecolor='white')

    if np.isnan(tl).any() or np.isnan(db).any():
        print 'warning: cannot plot pdfs containing nans'
        return

    # log scale
    plt.semilogy(tl[0], tl[1], color='blue', linewidth=2, label='$^{208}$Tl MC')
    plt.semilogy(db[0], db[1], color='red', linewidth=2, label='$0\\nu\\beta\\beta$ MC')

    bin_width = np.diff(tl[0])[0]
    plt.xlabel('Time residual (ns)')
    plt.ylabel('Normalized counts per %f ns bin' % bin_width)
    plt.legend(loc='upper left')
    plt.axis((-50,50,1e-6,1e-1))
    f.savefig('figures/pdfs_linear%s.pdf' % suffix)

    if show:
        plt.show()
        print 'press enter to continue'
        raw_input()

    plt.clf()

    # linear scale
    plt.plot(tl[0], tl[1], color='blue', linewidth=2, label='$^{208}$Tl MC')
    plt.plot(db[0], db[1], color='red', linewidth=2, label='$0\\nu\\beta\\beta$ MC')

    bin_width = np.diff(tl[0])[0]
    plt.xlabel('Time residual (ns)')
    plt.ylabel('Normalized counts per %f ns bin' % bin_width)
    plt.legend(loc='upper left')
    plt.axis((-50,50,1e-6,8e-2))
    f.savefig('figures/pdfs%s.pdf' % suffix)

    if show:
        plt.show()
        print 'press enter to continue'
        raw_input()

    plt.clf()


def calculate_fiducial_cut(cuts, years=1.0, binsize=100, background='./b8.txt', show=False):
    '''Compare leaked events to the irreducible distributed background.

    The number of events per bin passing all cuts is determined from any
    existing likelihood ratio pickles. We fit an exponential to these, allowing
    extrapolation to small radii where statistics are poor. 

    The point of intersection between the integrated exponential and the
    irreducible background provides a fiducial radius cut.

    Produces "figures/fiducial_radius.pdf".

    :param cuts: Cuts with which ratios are calculated
    :param years: Years of 208Tl data used to calculate L ratios
    :param binsize: Radial bin size in mm
    :param background: File containing points for integrated background dist.
    :param show: If True, show the plot
    '''
    import numpy as np
    import matplotlib
    if not show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    expo = lambda x, a, b: a * np.exp(b * x)

    # load rates
    nevents = {}
    for cut in cuts:
        try:
            with open('lratios/lratio%s.pickle' % cut.as_suffix()) as f:
                lratios = pickle.load(f)
                nevents[cut.r[1]] = len(lratios['tl208'])
        except IOError:
            print 'warning: unable to open lratios/lratio%s.pickle' % cut.as_suffix()
            continue

    for k, v in nevents.items():
        nevents[k] = v / years

    x, y = zip(*sorted(nevents.items()))
    yy = np.array([np.sum(y[:i]) for i in range(1,len(y)+1)])
    print x
    print y
    print yy

    popt, pcov = curve_fit(expo, np.array(x), np.array(yy), p0=(1e-8, 5e-3))

    xx = np.arange(2500, 5501, binsize)
    yyint = expo(xx, *popt)

    b8 = 15.0 * 4/3 * np.pi * np.power(xx, 3) / np.power(6005, 3)

    xx_oversampled = np.linspace(2500, 6500, 4000)
    yyint_oversampled = np.interp(xx_oversampled, xx, yyint)
    b8_oversampled = np.interp(xx_oversampled, xx, b8)

    intersection = xx_oversampled[np.abs(yyint_oversampled - b8_oversampled).argmin()]

    f = plt.figure(1, facecolor='white')

    plt.errorbar(x, yy, yerr=np.sqrt(1.0/yy), label='$^{208}$Tl MC')
    plt.semilogy(xx, yyint, label='$^{208}$Tl fit')
    plt.semilogy(xx, b8, label='$^{8}$B')
    plt.axvline(x=intersection, color='black', label='%1.0f mm' % intersection)

    plt.xlim((2500, 5500))
    plt.ylim((1, 1e4))
    plt.xlabel('Fiducial radius (mm)')
    plt.ylabel('Events inside fiducial radius')
    plt.legend(loc='upper left')
    plt.grid(True)
    f.savefig('figures/fiducial_radius.pdf')

    if show:
        plt.show()
        print 'press enter to continue'
        raw_input()

    plt.clf()

    return x, y

def cfrtest():
    e_mean = 2.555
    e_1s = e_mean + 0.1927 # mean + 1 sigma
    cut_list = [
        Cut(e=(e_mean, e_1s), r=(5400,5500), t=(-50,50)),
        Cut(e=(e_mean, e_1s), r=(5300,5400), t=(-50,50)),
        Cut(e=(e_mean, e_1s), r=(5200,5300), t=(-50,50)),
        Cut(e=(e_mean, e_1s), r=(5100,5200), t=(-50,50)),
        Cut(e=(e_mean, e_1s), r=(5000,5100), t=(-50,50)),
        Cut(e=(e_mean, e_1s), r=(4000,5000), t=(-50,50)),
        Cut(e=(e_mean, e_1s), r=(3000,4000), t=(-50,50)),
        Cut(e=(e_mean, e_1s), r=(2000,3000), t=(-50,50)),
        Cut(e=(e_mean, e_1s), r=(1000,2000), t=(-50,50)),
        Cut(e=(e_mean, e_1s), r=(0,1000), t=(-50,50))
    ]

    return calculate_fiducial_cut(cut_list, 4.339, show=True)

