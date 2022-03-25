The code is from http://kendrickkay.net/analyzePRF/

The population receptive field (pRF) technique described in Dumoulin and Wandell (2008) provides a means for estimating the aggregate receptive field of populations of neurons as measured with fMRI. On this web site, we provide:

Stimuli that can be used for pRF mapping
MATLAB code (analyzePRF) that can be used to estimate pRFs from fMRI data
Because the code is general, it is not necessary to use the stimuli provided here; the code can be used in a standalone fashion to analyze a variety of data.

The specific pRF model that we implement is the Compressive Spatial Summation (CSS) model, which is similar to the original model (2D isotropic Gaussian), except that a static power-law nonlinearity is added after summation across the visual field. The nonlinearity is typically compressive, hence the name of the model. The CSS model is described in the following paper:

Kay, K.N., Winawer, J., Mezer, A., & Wandell, B.A. Compressive spatial summation in human visual cortex. Journal of Neurophysiology (2013)

If you use analyzePRF, please cite the above paper as well as the original Dumoulin and Wandell (2008) paper.

Regarding the stimuli, they have not yet been used in a published paper. If you use the stimuli, please cite this web site. In addition, please cite Kriegeskorte et al. 2008 since the stimuli are based on the visual object stimuli used in that study.

Terms of use: The content provided by this web site is licensed under a Creative Commons Attribution 3.0 Unported License. You are free to share and adapt the content as you please, under the condition that you cite the appropriate manuscripts.

Acknowledgements: Thanks to An Vu for 7T data collection and advice on stimulus design.