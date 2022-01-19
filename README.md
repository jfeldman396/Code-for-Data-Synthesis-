# BFGMSD

The following scripts are a demonstration of utilizing the methods proposed in a “Bayesian Data Synthesis and the Utility-Risk Trade-Off for Mixed Epidemeological Data” to create a fully synthetic version of a simulated dataset. The simulated dataset is comprised of an unordered categorical, binary, continuous and count variable. The count variable in this case is simulated as a complex function of the other variables, and as a result, is left out of the initial copula model for targeted univariate synthesis.

The scripts are meant to be run sequentially, with the ordering speciﬁed by the numeric preﬁx in the naming convention. Each script will create variables stored in the global environment that may be accessed at subsequent points in the framework. In addition, the script ‘sampleMGP.R’ should be run or sourced prior to running the sampler, as it contains the functions necessary for sampling the parameters of the multiplicative gamma process prior utilized in our sampling algorithm.

Version 2.2 or later of the TruncatedNormal package is required to run ‘2_Synthesizer_demo.R’. At the bottom of the ‘3_NonlinearSynthesizer_demo.R’, the user can ﬁnd simple univariate and bivariate summaries that conﬁrm the utility of the synthetic data.
