# dWit: an R-library for algorithms around the notion of distributional witnesses

## Overview & Examples

### Two sample total variation lower bound

A simulation file canbe found in `tests/test_twoSampleTvLb.R`.

### Total Variation Lower Bound for time ordered mixtures

A simulation file canbe found in `tests/test_dwlb.R`.

## Installation

To install the package from github (private repository), you first need to go to [this setting page](https://github.com/settings/tokens) to create an access token, say `TOKEN`, with all **repo** scopes, and then run the following commands in R (using your newly generated token, `TOKEN`):

``` r
install.packages("devtools")
devtools::install_github("lorismichel/dWit", auth_token = "TOKEN")
```

## Issues

To report an issue, please use the [issue tracker](http://github.com/lorismichel/dWit/issues) on github.com.
