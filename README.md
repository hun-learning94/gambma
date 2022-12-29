
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gambms

<!-- badges: start -->
<!-- badges: end -->

The goal of the package `gambms` is to provide RCpp implementations of
estimation methods for generalized additive models based on Bayesian
model selection with mixtures of g-priors.

## Installation

You can install the development version of gambms from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hun-learning94/gambms")
```

## Examples

We illustrate the use case of `gambms` for a simulated dataset and a
real data. For more details, please see the paper in `docs` folder.

### Simulated data (Poisson, FREE-knot)

``` r
library(gambms)
set.seed(2021311165)
f_list = list(f1 = function(x) 0.5 * (2*x^5 + 3*x^2 + cos(3*pi*x) - 1),
              f2 = function(x) x,
              f3 = function(x) 0.75*(0.0035 * (x*3 + 1.5)^3 + (x > -0.5 & x < 0.85) *
                                       0.07 *sin(1.7*pi*(x*3 + 1.5)^2 / 3.2)*(x*3 -2.5)^2 * exp(x*3 + 1.5)))
n = 300
dat = simmat(f_list, -1, 1, n = n, family = "poisson")
maxk = 20
mf = y~ncs(x1, nk=maxk)+ ncs(x2, nk = maxk)  + ncs(x3, nk = maxk)
fit_sim = tryCatch(
  gambms(mf, dat,
         knotConfig = "FREE",
         prior = "Robust",
         family = "poisson",
         printIter=500,
         freeCtrl=list(mcmcIter = 4000, thin = maxk)),
  error = function(cnd)cnd
)
plot(fit_sim)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

### Real data 1) Boston Housing Prices (Gaussian, VS-knot)

``` r
data("Boston")
maxk=15; lambda = 1
mf = log(medv) ~ chas +
  ncs(crim, nk = maxk, lambda=lambda) +
  ncs(zn, nk = maxk, lambda=lambda) +
  ncs(indus, nk = maxk, lambda=lambda) +
  ncs(nox, nk = maxk, lambda=lambda) +
  ncs(rm, nk = maxk, lambda=lambda) +
  ncs(age, nk = maxk, lambda=lambda) +
  ncs(dis, nk = maxk, lambda=lambda) +
  ncs(rad, nk = maxk, lambda=lambda) +
  ncs(tax, nk = maxk, lambda=lambda) +
  ncs(ptratio, nk = maxk, lambda=lambda) +
  ncs(black, nk = maxk, lambda=lambda) +
  ncs(lstat, nk = maxk, lambda=lambda)
fit =  gambms(mf, Boston,
                 knotConfig = "VS",
                 prior = "Robust",
                 family = "gaussian")
summary(fit)
plot(fit)
plotnumknot(fit)
plotresiduals(fit)
```

    #> fm: 
    #> log(medv) ~ chas + ncs(crim, nk = maxk, lambda = lambda) + ncs(zn, 
    #>     nk = maxk, lambda = lambda) + ncs(indus, nk = maxk, lambda = lambda) + 
    #>     ncs(nox, nk = maxk, lambda = lambda) + ncs(rm, nk = maxk, 
    #>     lambda = lambda) + ncs(age, nk = maxk, lambda = lambda) + 
    #>     ncs(dis, nk = maxk, lambda = lambda) + ncs(rad, nk = maxk, 
    #>     lambda = lambda) + ncs(tax, nk = maxk, lambda = lambda) + 
    #>     ncs(ptratio, nk = maxk, lambda = lambda) + ncs(black, nk = maxk, 
    #>     lambda = lambda) + ncs(lstat, nk = maxk, lambda = lambda)
    #> 
    #> Family:                                           gaussian 
    #> Link function:                                    identity 
    #> g-prior:                                          Robust 
    #> ------------------------------------------------------------ 
    #> Sample size:                                      506 
    #> Number of smooth covariates:                      12 
    #> Number of linear terms (including intercept):     2 
    #> Knot configuration:                               VS 
    #> Maximum number of knots:                          15 15 15 15 15 15 15 8 15 15 15 15 
    #> Prior mean on number of knots:                    1 1 1 1 1 1 1 1 1 1 1 1 
    #> ------------------------------------------------------------- 
    #> Marginal probability that sm is linear: 
    #>    crim      zn   indus     nox      rm     age     dis     rad     tax ptratio 
    #>  0.0015  0.8600  0.8000  0.0000  0.0000  0.9000  0.0000  0.7800  0.6100  0.7300 
    #>   black   lstat 
    #>  0.7400  0.1500 
    #> ------------------------------------------------------------- 
    #> Linear coefficients:
    #>                 mean   median       sd  z_score  lower95  upper95
    #> (Intercept)   3.0363   3.0362   0.0065 464.3725   3.0234   3.0493
    #> chas          0.0301   0.0295   0.0286   1.0521  -0.0260   0.0860
    #> ------------------------------------------------------------ 
    #> Sampled g values:
    #>       mean   median       sd  lower95  upper95
    #> g 334.3914 297.4666 155.6324 154.4896 756.6889
    #> ------------------------------------------------------------- 
    #> Sampled sigma values:
    #>       mean median     sd lower95 upper95
    #> sig 0.1419 0.1417 0.0049  0.1326  0.1521
    #> -------------------------------------------------------------

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-5-2.png" width="100%" />
<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

### Real data 2) Pima Indian Diabetes (Bernoulli, EVEN-knot)

``` r
data("Pima")
maxk=15; lambda = 1
mf = diabetes ~
  ncs(pregnant, nk = maxk, lambda = lambda) +
  ncs(glucose, nk = maxk, lambda = lambda) +
  ncs(pressure, nk = maxk, lambda = lambda) +
  ncs(triceps, nk = maxk, lambda = lambda) +
  ncs(mass, nk = maxk, lambda = lambda) +
  ncs(pedigree, nk = maxk, lambda = lambda) +
  ncs(age, nk = maxk, lambda = lambda)
fit =  gambms(mf, Pima,
              knotConfig = "EVEN",
              prior = "Robust",
              family = "bernoulli",
              printIter = 500)

summary(fit)
plot(fit)
plotnumknot(fit)
plotresiduals(fit)
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
