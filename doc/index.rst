208Tl Rejection Analysis
========================
This code uses data (or simulation) to evaluate the performance of a likelihood
ratio-based cut of 208Tl beta decays occuring in the bulk AV acrylic. The
gammas from the 208Tl decay will produce light more spread out in time,
providing a handle for discrimination.

It does this in two steps:

1. Use one data set to build 1D (PMT hit time residual) PDFs for signal and
   background hypotheses

2. Loop through events of known type in a second data set and compute the ratio
   of likelihoods for the hypotheses

Likelihood ratio distributions for cuts at several radii may then be analyzed
to choose a fiducial volume cut that balances the leaked backgrounds with the
irreducible 8B neutrino scatters.

.. note::

   For a detailed explanation of the theory, analysis, and results, see
   :download:`the full report <../report/tl208.pdf>`.

Usage
-----
To perform the full analysis, run::

    $ python tl208_analysis.py [ncpus]

You will need to edit the data file paths at the top of ``tl208_analysis.py``
to point to your data sets.

Any useful analysis will require very large event samples for PDF generation
and for fake data, on the order of 5 million 208Tl decays each. The analysis
may take a long time (several hours). It can be accelerated by running parallel
processes with the ``ncpus`` argument (defaults to 2). However, parallel
processes quickly become I/O bound.

If appropriate PDFs are found in ``pdfs/``, they will be loaded. Otherwise PDFs
will be generated from data.

The code is highly modular, and individual parts of the analysis may be used
in other code or within the Python interpreter. ``main`` in
``tl208_analysis.py`` provides a usage example.

File Locations
--------------
The scripts output various files in a few places:

* ``pdfs``: PDFs generated from data. These take a very long time to produce, and
  so are cached for later reanalysis. ``tl208_analysis.py`` will search this
  location for appropriate PDFs before computing them.
* ``figures``: Portable Document Format PDF plots of probability distribution
  PDFs, PDF overlays, likelihood ratios, etc.
* ``lratios``: Lists of per-event likelihood ratios for different event classes.
  These are also expensive to compute, so are saved for later re-plotting,
  evaluating sacrifice, etc.

API
---

.. toctree::
   :maxdepth: 2

   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

