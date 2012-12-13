import sys
import multiprocessing
import numpy as np

import util

# db parameters
scintvg = 190.0
avvg = 1.93109181500140664e2
watervg = 2.17554021555098529e2
offset = 0.6
av_inner_r = 6005.3
av_outer_r = 6060.4


# global counters incremented in extract_residuals
# FIXME no globals!
total_events = multiprocessing.Value('i', 0)
events_reconstructed = multiprocessing.Value('i', 0)
events_passing_cuts = multiprocessing.Value('i', 0)


def extract_residuals(filename, cut=None):
    '''Loop over events and pull out time residuals for those passing cuts.

    For best results, call with a multiprocessing.Pool map.

    :param filename: Path to the RAT ROOT file
    :param cut: Cut instance with cuts to apply to data
    :returns: A list of residuals
    '''
    from rat import ROOT
    print '.',
    sys.stdout.flush()

    f = ROOT.TFile(filename)
    if f.IsZombie():
        return np.array([])

    if cut is None:
        cut = util.Cut()

    t = f.Get("T")
    ds = ROOT.RAT.DS.Root()
    t.SetBranchAddress('ds', ds)

    runt = f.Get("runT")
    run = ROOT.RAT.DS.Run()
    runt.SetBranchAddress('run', run)
    runt.GetEntry(0)

    slp = run.GetStraightLinePath()
    slp.Initialise(av_inner_r, av_outer_r)
    elv = run.GetEffectiveVelocityTime()
    elv.Initialise(scintvg, avvg, watervg, offset)
    pmtprop = run.GetPMTProp()

    res = []
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        total_events.value += 1
        if ds.GetEVCount() > 0:
            try:
                fitter_valid = ds.GetEV(0).GetFitResult('scintFitter').GetValid()
            except Exception:
                fitter_valid = False

            if not fitter_valid:
                continue

            events_reconstructed.value += 1
            v = ds.GetEV(0).GetFitResult('scintFitter').GetVertex(0)
            fit_t = v.GetTime()
            fit_e = v.GetEnergy()
            fit_p = v.GetPosition()
            fit_r = fit_p.Mag()

            if fit_r < cut.r[0] or fit_r > cut.r[1] or fit_e < cut.e[0] or fit_e > cut.e[1]:
                continue

            events_passing_cuts.value += 1

            for ipmt in range(ds.GetEV(0).GetPMTCalCount()):
                pmt = ds.GetEV(0).GetPMTCal(ipmt)
                pmtpos = pmtprop.GetPos(pmt.id)

                dscint, dav, dwater = ROOT.Double(0), ROOT.Double(0), ROOT.Double(0)
                slp.CalcByPosition(fit_p, pmtpos, dscint, dav, dwater)
 
                tof = elv.CalcByDistance(dscint, dav, dwater)

                r = ds.GetEV(0).GetPMTCal(ipmt).GetsPMTt() - fit_t - tof
                if r < cut.t[0] or r > cut.t[1]:
                    continue
                res.append(r)

    return np.array(res)


def plot(pdf, nevents=None, show=False, cut=None, suffix=''):
    '''Plot the distributions of the likelihood ratio parameter.

    Creates "figures/pdf[suffix].pdf".

    :param h: Bin values
    :param e: Bin locations
    :param nevents: Number of events (for display only)
    :param show: If True, show the plot interactively
    :param suffix: Appended to the written PDF filename
    '''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if cut is None:
        cut = util.Cut()

    e, h = zip(*pdf)

    f = plt.figure(1, facecolor='white')
    label = ('MC (%i events)' % nevents) if nevents is not None else 'MC'
    plt.semilogy(e, h, color='black', label=label)

    bin_width = np.diff(e)[0]
    plt.xlabel('Fit time (ns)')
    plt.ylabel('Normalized counts per %f ns bin' % bin_width)
    plt.legend(loc='upper left')
    plt.axis((cut.t[0],cut.t[1],2e-6,0.3))

    f.savefig('figures/pdf%s.pdf' % suffix)

    if show:
        plt.show()

    plt.clf()


def phat(x, data):
    '''P hat(x) is a Gaussian kernel density estimate of the PDF from which
    data is drawn.

    :param x: Value of x at which to evaluate the estimated PDF
    :param data: The data set (a 1d array of samples).
    :returns: The value of the estimated PDF (float)
    '''
    data = np.array(data)
    h = np.power(4.0/(3*len(data)), 1.0/5) * 50
    return 1.0/(h*np.sqrt(2.0*np.pi)*len(data)) * np.sum(np.exp(-0.5*np.square((data-x)/h)))

    # generating a smoothed pdf
    #x = np.linspace(min(e), max(e), 500)
    #p = [phat(xx, res) for xx in x]


def rebin(a, factor, trim=False):
    '''Combine adjacent array elements.

    :param a: The array
    :param factor: Number of adjacent elements to combine
    :param trim: If True, trim elements off the end to make length divide by
                 factor evenly.
    :returns: array, shorter by a factor or factor
    '''
    if trim:
        a = a[:len(a)-len(a)%factor]

    if not len(a) % factor == 0:
            raise Exception('Array length must be evenly divisible by factor')

    rebinned = [1.0/factor * np.sum(a[i:i+factor]) for i in range(0, len(a)-1, factor)]

    return np.array(rebinned)


class ERF:
    def __init__(self, cut=None):
        self.cut = (cut if cut is not None else util.Cut())

    def __call__(self, filename):
        import pdf
        return pdf.extract_residuals(filename, cut=self.cut)

def make(filenames, nprocs, cut):
    '''Create time residual PDF for a set of data files.

    Note: you may wish to use a smaller number of nprocs than you have CPUs;
    this function will almost certainly be I/O-bound.

    :param filenames: list of RAT ROOT files containing data
    :param cut: A Cut instance with cuts to apply to data
    :param nprocs: number of parallel jobs to run
    '''
    p = multiprocessing.Pool(nprocs)
    erf = ERF(cut=cut)

    res = np.array(list(util.flatten(p.map(erf, filenames))))

    print
    print len(res), 'entries'
    h, e = np.histogram(res, bins=750, range=(cut.t[0],cut.t[1]), normed=True)

    pdf = np.array(zip(e,h))

    print 'total events:', total_events.value
    print 'events reconstructed:', events_reconstructed.value
    print 'events passing cuts:', events_passing_cuts.value

    return pdf

