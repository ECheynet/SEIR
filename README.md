# Generalized SEIR Epidemic Model (fitting and computation)

[![View Generalized SEIR Epidemic Model (fitting and computation) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation)

## Description
A generalized SEIR model with seven states, as proposed by ref. [2]  is numerically implemented. There exist other types of generalized SEIR model that can be explored, but here I only use a single one for the sake of simplicity. The numerical implementation is done from scratch except for the fitting, that relies on the function "lsqcurvfit".

One major difference with respect to ref. [2] is the expression of the death rate and recovery rate, which are here analytical and empirical functions of the time. The idea behind this time-dependency as that the death rate (due to the disease) should converge toward zero as the time increases. If the death rate is kept constant, the number of death may be overestimated. As the same time, the recovery rate is also converging toward a constant value. Births and natural death are not modelled here. This means that the total population, including the number of deceased cases, is kept constant.

Note that ref. [2] is a preprint that is not peer-reviewed and I am not qualified enough to judge the quality of the paper.

## Content
The present submission contains:
- A function SEIQRDP.m that is used to simulate the time histories of the infectious, recovered and dead cases (among others)
- A function fit_SEIQRDP.m that estimates the eight parameters used in SEIQRDP.m in the least square sense.
- One example file Example1.mlx, which presents the numerical implementation.
- One example file Example2.mlx, which uses data collected by the Johns Hopkins University for the COVID-19 epidemy [3] for Hubei province (China).
- One example file Example3.mlx, which uses data collected by the Johns Hopkins University for the COVID-19 epidemy [3] for South Korea.
- One file "ItalianRegions.mlx" written by Matteo Secli (https://github.com/matteosecli) that I have modified for a slightly more robust fitting.
- One file "FrenchRegions.mlx", which gives another example for Data collected in France. The data quality is not as good as expected, so the fitting is unlikely to provide reliable parameter estimates.
- One example file ChineseProvinces.mlx, which illustrates how the function fit_SEIQRDP.m is used in a for loop to be fitted to the data [3] from the different Chinese provinces.
- One example "uncertaintiesIssues.mlx", which illustrates the danger of fitting limited data sets.
- One example "Example_US_cities.mlx" that illustrates the fitting when "recovered" data are not available.
- One function getDataCOVID, which read from [3] the data collected by Johns Hopkins University.
- One function getDataCOVID_ITA written by Matteo Secli (https://github.com/matteosecli), that collects the updated data of the COVID-19 pandemic in Italy from the Italian government [4]
- One function getDataCOVID_FRA that collects the updated data in France from [5]
- One function getDataCOVID_US that collects the updated data in the USA from [3]

Any question, comment or suggestion is welcome.

## References

[1] https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#Bio-mathematical_deterministic_treatment_of_the_SIR_model

[2] Peng, L., Yang, W., Zhang, D., Zhuge, C., & Hong, L. (2020). Epidemic analysis of COVID-19 in China by dynamical modeling. arXiv preprint arXiv:2002.06563.

[3] https://github.com/CSSEGISandData/COVID-19

[4] https://github.com/pcm-dpc/COVID-19

[5] https://github.com/cedricguadalupe/FRANCE-COVID-19



## Example (case of COVID-19 in Italy) 

The fitting of the extended SEIR model to real data provides the following results:

![Active, recoverd and dceased cases in italy](Italy.png)



