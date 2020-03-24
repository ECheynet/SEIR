# Generalized SEIR Epidemic Model (fitting and computation)

[![View Generalized SEIR Epidemic Model (fitting and computation) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation)






To learn how to solve systems of non-linear coupled differential equations, I ended up playing with compartmental models in epidemiology, in particular, compartmental SIR models and all their derived models [1]. I took the generalized SEIR model by ref. [2] because it deals with a contemporary topic. There exist other types of generalized SEIR model that can be explored, but here I only use a single one for the sake of simplicity. The numerical implementation is done from scratch except for the fitting, that relies on the function "lsqcurvfit".

In version 1.2 and above, I also implemented two of the parameters as explicit functions of the time: the death rate and recovery rate. The idea behind this time-dependency as that the death rate (due to the disease) should be zero after an infinite time. If the death rate is kept constant, the number of death may become overestimated. As the same time, the recovery rate is also increasing toward a threshold value. Births and natural death are not modelled here. This means that the total population, including the deads, is kept constant.

Note that ref. [2] is a preprint that is not peer-reviewed and I am not qualified enough to judge the quality of the paper.

The present submission contains:
- A function SEIQRDP.m that is used to simulate the time histories of the infectious, recovered and dead cases (among others)
- A function fit_SEIQRDP.m that estimates the eight parameters used in SEIQRDP.m in the least square sense.
- One example file Example1.mlx that use only simulated data
- One example file Example2.mlx that use data collected by the Johns Hopkins University for the COVID-19 epidemy [3] for Hubei province (China).
- One example file Example3.mlx that use data collected by the Johns Hopkins University for the COVID-19 epidemy [3] for Italy.
- One file "ItalianRegions.mlx" written by Matteo Secli (https://github.com/matteosecli) that I have modified for a slightly more robust fitting.
- One example file Example4.mlx that illustrates how the function fit_SEIQRDP.m is used in a for loop to be fitted to the data [3] from the different Chinese provinces.
- One example "uncertaintiesIssues.mlx" illustrating the danger of fitting limited data sets.
- One function getDataCOVID, which read from [3] the data collected by Johns Hopkins University.
- One function getDataCOVID_ITA written by Matteo Secli (https://github.com/matteosecli), that  collects the updated data of the COVID-19 pandemic in Italy from the Italian governement [4]

That is the first version of the submission. It's probably full of typos that will be gradually corrected.

Any question, comment or suggestion is welcome.

References:

[1] https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#Bio-mathematical_deterministic_treatment_of_the_SIR_model

[2] Peng, L., Yang, W., Zhang, D., Zhuge, C., & Hong, L. (2020). Epidemic analysis of COVID-19 in China by dynamical modeling. arXiv preprint arXiv:2002.06563.

[3] https://github.com/CSSEGISandData/COVID-19

[4] https://github.com/pcm-dpc/COVID-19
