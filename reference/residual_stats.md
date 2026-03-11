# Calculate statistics that describe the residuals of a fit

Calculates several key statistics from the residuals of of a fit: the
residual sum of squares (`RSS`), the mean squared error (`MSE`), the
root mean squared error (`RMSE`), the residual standard error (`RSE`),
and the Akaike information criterion (`AIC`). This function is used
internally by all fitting functions in the `PhotoGEA` package, such as
[`fit_ball_berry`](https://eloch216.github.io/PhotoGEA/reference/fit_ball_berry.md)
and
[`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md).

## Usage

``` r
residual_stats(fit_residuals, units, nparam)
```

## Arguments

- fit_residuals:

  A numeric vector representing the residuals from a fit, i.e., the
  differences between the measured and fitted values.

- units:

  A string expressing the units of the residuals.

- nparam:

  The number of free parameters that were varied when performing the
  fit.

## Details

This function calculates several model-independent measures of the
quality of a fit. The basis for these statistics are the `residuals`
(also known as the `errors`). If the measured values of a quantity `y`
are given by `y_measured` and the fitted values are `y_fitted`, then the
residuals are defined to be `residual = y_measured - y_fitted`. The key
statistics that can be calculated from the residuals are as follows:

- The residual sum of squares (`RSS`) is also known as the sum of
  squared errors (`SSE`). As its name implies, it is simply the sum of
  all the squared residuals: `RSS = sum(residuals^2)`.

- The mean squared error (`MSE`) is the mean value of the squared
  residuals: `MSE = sum(residuals^2) / n = RSS / n`, where `n` is the
  number of residuals.

- The root mean squared error (`RMSE`) is the square root of the mean
  squared error: `RMSE = sqrt(MSE) = sqrt(RSS / n)`.

- The residual standard error `RSE` is given by `RSE = sqrt(RSS / dof)`,
  where `dof = n - nparam` is the number of degrees of freedom involved
  in the fit.

- The Akaike information criterion `AIC` is given by
  `AIC = npts * (log(2 * pi) + 1) + npts * log(MSE) + 2 * (nparam + 1)`.

For a given model, the `RMSE` is usually a good way to compare the
quality of different fits. When trying to decide which model best fits
the measured data, the `AIC` may be a more appropriate metric since it
controls for the number of parameters in the model.

The AIC definition used here is appropriate for the results of maximum
likelihood fitting with equal variance, or minimum least squares
fitting. For more details about the AIC equation above and its relation
to the more general definition of AIC, see Section 2 of Banks & Joyner
(2017).

**References**:

Banks, H. T. & Joyner, M. L. "AIC under the framework of least squares
estimation." Applied Mathematics Letters 74, 33–45 (2017)
\[[doi:10.1016/j.aml.2017.05.005](https://doi.org/10.1016/j.aml.2017.05.005)
\].

## Value

An `exdf` object with one row and the following columns: `npts` (the
number of residual values), `nparam`, `dof`, `RSS`, `MSE`, `RMSE`,
`RSE`, `AIC`.

## Examples

``` r
# Generate some random residuals
residuals <- runif(10, -1, 1)

# Calculate residual stats as if these values had units of `kg` and were related
# to a model with 3 free parameters
residual_stats(residuals, 'kg', 3)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   npts [residual_stats] (NA) nparam [residual_stats] (NA)
#> 1                         10                            3
#>   dof [residual_stats] (NA) RSS [residual_stats] ((kg)^2)
#> 1                         7                      3.774498
#>   MSE [residual_stats] ((kg)^2) RMSE [residual_stats] (kg)
#> 1                     0.3774498                  0.6143694
#>   RSE [residual_stats] (kg) AIC [residual_stats] ()
#> 1                 0.7343119                26.63559
```
