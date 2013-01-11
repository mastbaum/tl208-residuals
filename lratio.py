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
    edges = pdf[0]
    bins = pdf[1]
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
        if ds.GetEVCount() > 0:
            try:
                fitter_valid = ds.GetEV(0).GetFitResult('scintFitter').GetValid()
            except Exception:
                fitter_valid = False

            if not fitter_valid:
                continue

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


def plot(tl_ratios, db_ratios, suffix='', bins=100, plot_sl=True, show=False):
    '''Plot the distributions of the likelihood ratios for Tl and DBD events.

    Creates "figures/lratios[suffix].pdf".

    :param tl_ratios: List of L ratios for Tl events
    :param db_ratios: List of L ratios for 0vbb events
    :param suffix: String appended to output filename
    :param bins: Number of bins in histogram
    :params plot_sl: Also plot sacrifice and leakage
    :param show: If True, show the plot interactively
    '''
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    f = plt.figure(1, facecolor='white')

    if len(tl_ratios) == 0 or len(db_ratios) == 0:
        print 'warning: cannot plot with no likelihood ratios'
        return

    r = (min(min(tl_ratios), min(db_ratios)), max(max(tl_ratios), max(db_ratios)))

    b_tl, e_tl = np.histogram(tl_ratios, bins=bins, range=r, normed=True)
    b_dbd, e_dbd = np.histogram(db_ratios, bins=bins, range=r, normed=True)

    plt.bar(e_tl[:-1], b_tl, alpha=0.5, color='blue', linewidth=0, width=np.diff(e_tl), label='$^{208}$Tl MC')
    plt.bar(e_dbd[:-1], b_dbd, alpha=0.5, color='red', linewidth=0, width=np.diff(e_dbd), label='$0\\nu\\beta\\beta$ MC')

    plt.xlabel('$\Delta$')
    plt.ylabel('Normalized counts per bin')
    plt.legend(loc='upper left')
    f.savefig('figures/lratio%s.pdf' % suffix)

    if show:
        plt.show()

    plt.clf()

    if plot_sl:
        tl = np.empty(shape=(len(b_tl), 2))
        tl[:,0] = e_tl[:-1]
        tl[:,1] = b_tl
        db = np.empty(shape=(len(b_dbd), 2))
        db[:,0] = e_dbd[:-1]
        db[:,1] = b_dbd

        return plot_sacrifice_leakage(tl, db, suffix, show)

def plot_sacrifice_leakage(tl_delta, db_delta, suffix='', show=False):
    '''Plot sacrifice and leakage.

    Creates "figures/sacrifice_leakage[suffix].pdf".

    :param tl_delta: L ratio distribution for Tl events
    :param db_delta: L ratio distribution for 0vbb events
    :param suffix: String appended to output filename
    :param show: If True, show the plot interactively
    :returns: ([cut values], [sacrifice], [leakage]) tuple of lists
    '''
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    sacrifice = []
    leakage = []

    for cut in tl_delta[:,0]:
        e_cut_bin = np.abs(tl_delta[:,0]-cut).argmin()
        d_cut_bin = np.abs(db_delta[:,0]-cut).argmin()
        leakage.append(1.0-np.sum(tl_delta[:,1][:e_cut_bin]*np.diff(tl_delta[:,0])[:e_cut_bin]) / np.sum(tl_delta[:,1][:-1]*np.diff(tl_delta[:,0])))
        sacrifice.append(np.sum(db_delta[:,1][:d_cut_bin]*np.diff(db_delta[:,0])[:d_cut_bin]) / np.sum(db_delta[:,1][:-1]*np.diff(db_delta[:,0])))

    f = plt.figure(1, facecolor='white')

    sacrifice = np.array(sacrifice)
    leakage = np.array(leakage)

    plt.plot(tl_delta[:,0], sacrifice, linewidth=2, linestyle='--',color='red', label='Sacrifice')
    plt.plot(tl_delta[:,0], leakage, linewidth=2, linestyle='-.', color='blue', label='Leakage')

    # 0vbb sensitivity
    sr = [(1.0-sacrifice[i])/np.sqrt(39.8+1000*(leakage[i]))/(1.0/np.sqrt(39.8+1000)) for i in range(len(tl_delta[:,0]))]
    delta_c = tl_delta[:,0][np.argmax(sr)]
    plt.plot(tl_delta[:,0], sr, linewidth=3, color='black', label='Relative $T^{1/2}_{0\\nu\\beta\\beta}$')
    plt.axvline(x=delta_c, color='black', label='$\Delta_c=%1.4f$' % delta_c)
    print 'delta_c:', delta_c
    print 'sacrifice:', sacrifice[np.argmax(sr)]
    print 'leakage:', leakage[np.argmax(sr)]

    plt.xlabel('$\Delta_c$')
    plt.ylabel('Ratio')
    plt.legend(loc='lower left')
    plt.grid(True)
    f.savefig('figures/sacrifice_leakage%s.pdf' % suffix)

    if show:
        plt.show()
        raw_input()

    plt.clf()

    return tl_delta[:,0], sacrifice, leakage


def plot_sensitivity(cut_list, years=4.339):
    r, nevents = util.calculate_fiducial_cut(cut_list, years)

    for cut in cut_list:
        try:
            with open('lratios/lratio%s' % cut.as_suffix()) as f:
                lr = pickle.load(f)
        except IOError:
            print 'warning: could not load l ratios from lrations/lratio%s' % cut.as_suffix()
            continue

        tl_ratios = lr['tl208']
        db_ratios = lr['dbd']

        if len(tl_ratios) == 0 or len(db_ratios) == 0:
            print 'warning: cannot plot with no likelihood ratios'
            continue

        r = (min(min(tl_ratios), min(db_ratios)), max(max(tl_ratios), max(db_ratios)))

        b_tl, e_tl = np.histogram(tl_ratios, bins=bins, range=r, normed=True)
        b_dbd, e_dbd = np.histogram(db_ratios, bins=bins, range=r, normed=True)

        tl = np.empty(shape=(len(b_tl), 2))
        tl[:,0] = e_tl[:-1]
        tl[:,1] = b_tl
        db = np.empty(shape=(len(b_dbd), 2))
        db[:,0] = e_dbd[:-1]
        db[:,1] = b_dbd

        x, sacrifice, leakage = plot_sacrifice_leakage(tl, db, cut.as_suffix())

        sr = [(1.0-sacrifice[i])/np.sqrt(39.8+268*(1.0-leakage[i]))/(1.0/np.sqrt(39.8+268)) for i in range(len(tl_delta[:,0]))]
        delta_c = tl_delta[:,0][np.argmax(sr)]

