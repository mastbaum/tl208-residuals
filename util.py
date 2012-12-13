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


def plot_pdf_overlay(tlpath, dbpath, suffix='', show=False):
    '''Plot two PDFs overlaid.

    Produces "figures/pdfs[suffix].pdf".

    :param tlpath: Path to 208Tl PDF pickle
    :param dbpath: Path to 0vbb PDF pickle
    :param suffix: Appended to output files
    :param show: If True, show the plot
    '''
    import pickle
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    with open(tlpath) as f:
        tl = pickle.load(f)

    with open(dbpath) as f:
        db = pickle.load(f)

    f = plt.figure(1, facecolor='white')

    plt.semilogy(tl[:,0], tl[:,1], color='blue', linewidth=2, label='$^{208}$Tl MC')
    plt.semilogy(db[:,0], db[:,1], color='red', linewidth=2, label='$0\\nu\\beta\\beta$ MC')

    bin_width = np.diff(tl[:,0])[0]
    plt.xlabel('Time residual (ns)')
    plt.ylabel('Normalized counts per %f ns bin' % bin_width)
    plt.legend(loc='upper left')
    plt.axis((-50,50,1e-6,1e-1))
    f.savefig('figures/pdfs%s.pdf' % suffix)

    if show:
        plt.show()

