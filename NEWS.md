# pclmbayes 0.1.0

* Initial release. Implements the method of Lambert and Eilers (2009)
  "Bayesian density estimation from grouped continuous data"
  (*Computational Statistics and Data Analysis* 53(4), 1388-1399).
* Frequentist fit `pclm()`: penalised scoring with BIC/AIC selection of
  the smoothing parameter.
* Bayesian fit `bpclm()`: modified Langevin-Hastings sampler with
  rotation by the Cholesky of a frequentist warm-start variance-
  covariance, adaptive tuning of the step size, and a Gibbs step on
  the smoothing precision.
* Optional shape-constraint priors: unimodality, log-concavity and
  monotonicity (Eq. 7 of the paper).
* Example datasets `bloodlead` (verbatim from Hasselblad et al. 1980 /
  the paper) and `tbdeaths1907` (illustrative reconstruction; see help
  page for caveats).
* Vignette reproducing the worked examples of the paper.
