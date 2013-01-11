import sys
import pickle
import multiprocessing
from glob import glob
import numpy as np

import util
import pdf
import lratio

te130 = glob('/mnt/dropbox/jail/home/dropbox/mastbaum/te130_0vbb/te130_0vbb-*.root') + \
        glob('/mnt/dropbox/jail/home/dropbox/te130_0vbb-*.root') + \
        glob('/home/mastbaum/snoplus/tl208/data/pdf/dbd/dbd*.root') + \
        glob('/home/mastbaum/snoplus/tl208/data/fakedata/dbd_*.root')

tl208 = glob('/mnt/dropbox/jail/home/dropbox/av_tl208*.root') + \
        glob('/mnt/dropbox/jail/home/dropbox/av_tl208/av_tl208*.root') + \
        glob('/mnt/dropbox/jail/home/dropbox/mastbaum/av_tl208/run*/av_tl208*.root') + \
        glob('/home/mastbaum/snoplus/tl208/data/pdf/tl208/run*/av_tl208*.root') + \
        glob('/home/mastbaum/snoplus/tl208/data/fakedata/av_tl208*.root')

print len(te130), 'te130 files'
print len(tl208), 'tl208 files'

pdf_mc = {
    'tl208': tl208[:len(tl208)/2],
    'dbd': te130[:len(te130)/2]
}

data_mc = {
    'tl208': tl208[len(tl208)/2:],
    'dbd': te130[len(te130)/2:]
}

class LRF:
    '''Wrapper to allow lratio.extract to be called in a multiprocessing map.

    :param pdf_tl: 208Tl PDF
    :param pdf_db: 0vbb PDF
    :param cut: Cut instance with cuts to apply to data
    '''
    def __init__(self, pdf_tl, pdf_db, cut=None):
        self.pdf_tl = pdf_tl
        self.pdf_db = pdf_db
        self.cut = cut if cut is not None else util.Cut()

    def __call__(self, filename):
        '''Call lratio.extract.

        :param filename: Path to file from which to extract L ratios
        :returns: List of L ratios
        '''
        import lratio
        return lratio.extract(filename, pdf_tl=self.pdf_tl, pdf_db=self.pdf_db, cut=self.cut)


def main(cut_list, nprocs=2, rebin_factor=6):
    '''Perform likelihood ratio analysis for several cuts.

    :param cuts: List of Cut instances with cuts for whih to perform analysis
    :param nprocs: Number of parallel processes to launch
    :param rebin_factor: Factor by which to rebin PDFs
    '''
    for cut in cut_list:
        suffix = cut.as_suffix()
        print suffix

        # load or build pdfs
        try:
            print 'loading tl pdf'
            pdf_tl = util.load_pdf('pdfs/pdf_tl208%s.pickle' % suffix)
        except Exception:
            print 'building tl208 pdf'
            pdf_tl = pdf.make(pdf_mc['tl208'], nprocs=nprocs, cut=cut)
            with open('pdfs/pdf_tl208%s.pickle' % suffix, 'w') as f:
                pickle.dump(pdf_tl, f)

        try:
            print 'loading 0vbb pdf'
            pdf_db = util.load_pdf('pdfs/pdf_dbd%s.pickle' % suffix)
        except Exception:
            print 'building 0vbb pdf'
            pdf_db = pdf.make(pdf_mc['dbd'], nprocs=nprocs, cut=cut)
            with open('pdfs/pdf_dbd%s.pickle' % suffix, 'w') as f:
                pickle.dump(pdf_db, f)

        # rebin
        pdf_tl_x = pdf.rebin(pdf_tl[:,0], rebin_factor)
        pdf_tl_y = pdf.rebin(pdf_tl[:,1], rebin_factor)
        pdf_tl = np.array((pdf_tl_x, pdf_tl_y))
        pdf.plot(pdf_tl, cut=cut, suffix='_tl208'+suffix)

        pdf_db_x = pdf.rebin(pdf_db[:,0], rebin_factor)
        pdf_db_y = pdf.rebin(pdf_db[:,1], rebin_factor)
        pdf_db = np.array((pdf_db_x, pdf_db_y))
        pdf.plot(pdf_db, cut=cut, suffix='_dbd'+suffix)

        util.plot_pdf_overlay('pdfs/pdf_tl208%s.pickle' % suffix, 'pdfs/pdf_dbd%s.pickle' % suffix, rebin=rebin_factor, suffix=suffix)

        # load or evaluate likelihood ratios
        try:
            print 'loading likelihoods'
            with open('lratios/lratio%s.pickle' % suffix) as f:
                lr = pickle.load(f)
            lr_tl = lr['tl208']
            lr_db = lr['dbd']
        except Exception:
            p = multiprocessing.Pool(nprocs)
            lrf = LRF(pdf_tl, pdf_db, cut)

            print 'calculating likelihoods for tl208 data'
            lr_tl = np.array(list(util.flatten(p.map(lrf, data_mc['tl208']))))

            print 'calculating likelihoods for dbd data'
            lr_db = np.array(list(util.flatten(p.map(lrf, data_mc['dbd']))))

            with open('lratios/lratio%s.pickle' % suffix, 'w') as f:
                d = {'tl208': lr_tl, 'dbd': lr_db}
                pickle.dump(d, f)

        lratio.plot(lr_tl, lr_db, suffix=suffix)

if __name__ == '__main__':
    e_mean = 2.555
    e_1s = e_mean + 0.1927 # mean + 1 sigma
    cut_list = [
        util.Cut(e=(e_mean, e_1s), r=(5400,5500), t=(-50,50)),
        util.Cut(e=(e_mean, e_1s), r=(5300,5400), t=(-50,50)),
        util.Cut(e=(e_mean, e_1s), r=(5200,5300), t=(-50,50)),
        util.Cut(e=(e_mean, e_1s), r=(5100,5200), t=(-50,50)),
        util.Cut(e=(e_mean, e_1s), r=(5000,5100), t=(-50,50)),
        util.Cut(e=(e_mean, e_1s), r=(4000,5000), t=(-50,50)),
        util.Cut(e=(e_mean, e_1s), r=(3000,4000), t=(-50,50)),
        util.Cut(e=(e_mean, e_1s), r=(2000,3000), t=(-50,50)),
        util.Cut(e=(e_mean, e_1s), r=(1000,2000), t=(-50,50)),
        util.Cut(e=(e_mean, e_1s), r=(0,1000), t=(-50,50))
    ]

    nprocs = int(sys.argv[1] if len(sys.argv)>1 else 2)

    main(cut_list, nprocs=nprocs)

#for cut in e_tl:
#    e_cut_bin = np.abs(e_tl-cut).argmin()
#    d_cut_bin = np.abs(e_dbd-cut).argmin()
#    sacrifice = np.sum(b_dbd[:d_cut_bin]*np.diff(e_dbd)[:d_cut_bin]) / np.sum(b_dbd*np.diff(e_dbd))
#    rejection = np.sum(b_tl[:e_cut_bin]*np.diff(e_tl)[:e_cut_bin]) / np.sum(b_tl*np.diff(e_tl))
#
#    print cut, sacrifice, rejection

