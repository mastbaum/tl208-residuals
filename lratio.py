import numpy as np
import util

# db constants
scintvg = 190.0
avvg = 1.93109181500140664e2
watervg = 2.17554021555098529e2
offset = 0.6
av_inner_r = 6005.3
av_outer_r = 6060.4


def evaluate_pdf(samples, pdf):
    edges = pdf[:,0]
    bins = pdf[:,1]
    widths = np.diff(edges)
    d = np.digitize(samples, edges)
    d = d[d < len(widths)]

    return np.sum(bins[d] * widths[d])


def extract(filename, pdf_tl, pdf_db, cut=None):
    '''Evaluate likelihood ratios for all events in a file that pass cuts.

    :param filename: Path to RAT ROOT file
    :param pdf_tl: 208Tl PDF
    :param pdf_db: 0vbb PDF
    :param cut: Cut instance with cuts to apply to data
    :returns: List of likelihood ratios
    '''
    print filename
    from rat import ROOT
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

    lratios = []
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if ds.GetEVCount() > 0 and ds.GetEV(0).GetFitResult('scintFitter').GetValid():
            v = ds.GetEV(0).GetFitResult('scintFitter').GetVertex(0)
            fit_t = v.GetTime()
            fit_e = v.GetEnergy()
            fit_p = v.GetPosition()
            fit_r = fit_p.Mag()

            if fit_r < cut.r[0] or fit_r > cut.r[1] or fit_e < cut.e[0] or fit_e > cut.e[1]:
                continue

            res = []
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

            l_tl = evaluate_pdf(res, pdf_tl)
            l_db = evaluate_pdf(res, pdf_db)

            l_ratio = np.log(l_db) - np.log(l_tl)

            lratios.append(l_ratio)

    return np.array(lratios)


def plot(tl_ratios, db_ratios, suffix='', show=False):
    '''Plot the distributions of the likelihood ratios for Tl and DBD events.

    Creates "figures/lratios[suffix].pdf".

    :param tl_ratios: List of L ratios for Tl events
    :param db_ratios: List of L ratios for 0vbb events
    :param suffix: String appended to output filename
    :param show: If True, show the plot interactively
    '''
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    f = plt.figure(1, facecolor='white')

    b_tl, e_tl = np.histogram(tl_ratios, bins=100, normed=True)
    b_dbd, e_dbd = np.histogram(db_ratios, bins=100, normed=True)

    plt.bar(e_tl[:-1], b_tl, alpha=0.5, color='blue', linewidth=0, width=np.diff(e_tl), label='$^{208}$Tl MC')
    plt.bar(e_dbd[:-1], b_dbd, alpha=0.5, color='red', linewidth=0, width=np.diff(e_tl), label='$0\\nu\\beta\\beta$ MC')

    plt.xlabel('$\Delta$')
    plt.ylabel('Probability')
    plt.legend(loc='upper left')
    f.savefig('figures/lratio%s.pdf' % suffix)

    if show:
        plt.show()

