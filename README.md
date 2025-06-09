# Code for: "When correlation matter: on uncertainty propagation in the case of data disaggregation"

This is the code for reproducing the results for our study on “When correlation matter: on uncertainty propagation in the case of data disaggregation” 
submitted to the Journal of Industrial Ecology
(JIE).

## Data required

To reproduce the results and run the scripts you need to download the
following data and put them into the `./data` folder: 

-   EXIOBASE V3.8.2 `IOT_2015_ixi.zip` and `MRSUT_2015.zip` from here:
    <https://zenodo.org/records/5589597>
-   Intermediate results from the article "“Estimating the uncertainty
    of the greenhouse gas emission accounts in Global Multi-Regional
    Input-Output analysis” available here (for details how those data
    was generated, see the paper and its code): 10.5281/zenodo.13806019
    . Download them and save them in the `./data` folder

## How to run the scripts

1.  [Clone the repository to create a local copy on you
    computer](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).

2.  Open the project by clicking on `uncertainty_disaggregation.Rproj`.

3.  [renv](https://cran.r-project.org/web/packages/renv/vignettes/renv.html)
    will automatically bootstrap itself, downloading and installing the
    appropriate version of renv. It will also ask you if you want to
    download and install all the packages it needs by running
    `renv::restore()`.
    
4. You might need to install the `MaxentDisaggregation` package manually: 
```
renv::install("simschul/MaxentDisaggregation")
```
