# AudioRestorationFramework
 
Supplemental material to the paper Flexible framework for audio restoration by Ondřej Mokrý, Pavel Rajmic and Pavel Záviška, accepted to the 23rd International Conference on Digital Audio Effects (DAFx), Vienna, 2020.

The repository contains all the MATLAB codes needed to produce the results presented in the paper. The main file is multitest.m, the file generating all the figures is plot_multitest.m. Pregenerated set of figures is available in the subfolder figures.

Note that some of the codes have two variants, which differ by the suffix -g; this denotes the implementation using the version of Condat-Vu algorithm assuming the use of a tight frame as the linear transform.
