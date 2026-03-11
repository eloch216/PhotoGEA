# Optimizers

These functions return optimizers that meet requirements for the
`optim_fun` input argument of
[`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md),
[`fit_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_variable_j.md),
[`fit_c4_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci.md),
and
[`fit_c4_aci_hyperbola`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci_hyperbola.md).
Essentially, they are wrappers for optimizers from other libraries that
serve to standardize their inputs and outputs.

## Usage

``` r
optimizer_deoptim(itermax, VTR = -Inf)

  optimizer_hjkb(tol, maxfeval = Inf, target = Inf)

  optimizer_nlminb(rel.tol, eval.max = 200, iter.max = 200, abs.tol = 0)

  optimizer_nmkb(tol, maxfeval = 2000, restarts.max = 10)

  optimizer_null()
```

## Arguments

- tol:

  A convergence tolerance value; to be passed to
  [`nmkb`](https://rdrr.io/pkg/dfoptim/man/nmkb.html) or
  [`hjkb`](https://rdrr.io/pkg/dfoptim/man/hookejeeves.html) via their
  `control` input arguments. A typical value is `1e-7`.

- maxfeval:

  A maximum value for the number of function evaluations to allow during
  optimization; to be passed to
  [`nmkb`](https://rdrr.io/pkg/dfoptim/man/nmkb.html) or
  [`hjkb`](https://rdrr.io/pkg/dfoptim/man/hookejeeves.html) via their
  `control` input arguments.

- target:

  A real number restricting the absolute function value; to be passed to
  [`hjkb`](https://rdrr.io/pkg/dfoptim/man/hookejeeves.html) via its
  `control` input argument.

- rel.tol:

  A relative convergence tolerance value; to be passed to
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html) via its `control`
  input argument. A typical value is `1e-10`.

- eval.max:

  A maximum value for the number of function evaluations; to be passed
  to [`nlminb`](https://rdrr.io/r/stats/nlminb.html) via its `control`
  input argument.

- iter.max:

  A maximum value for the number of iterations; to be passed to
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html) via its `control`
  input argument.

- abs.tol:

  An absolute convergence tolerance value; to be passed to
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html) via its `control`
  input argument.

- restarts.max:

  A maximum value for the number of restarts allowed during
  optimization; to be passed to
  [`nmkb`](https://rdrr.io/pkg/dfoptim/man/nmkb.html) via its `control`
  input argument.

- itermax:

  The maximum number of generations to be used; to be passed to
  [`DEoptim`](https://rdrr.io/pkg/DEoptim/man/DEoptim.html) via its
  `control` input argument. Note that when `VTR` is `-Inf`, the
  optimizer will always use the maximum number of generations. A typical
  value is `200`.

- VTR:

  The value to be reached; to be passed to
  [`DEoptim`](https://rdrr.io/pkg/DEoptim/man/DEoptim.html) via its
  `control` input argument.

## Details

`optimizer_deoptim` is a wrapper for
[`DEoptim`](https://rdrr.io/pkg/DEoptim/man/DEoptim.html).

`optimizer_hjkb` is a wrapper for
[`hjkb`](https://rdrr.io/pkg/dfoptim/man/hookejeeves.html).

`optimizer_nlminb` is a wrapper for
[`nlminb`](https://rdrr.io/r/stats/nlminb.html).

`optimizer_nmkb` is a wrapper for
[`nmkb`](https://rdrr.io/pkg/dfoptim/man/nmkb.html).

`optimizer_null` simply returns the initial guess without doing any
optimization; it can be useful for viewing initial guesses.

See the documentation for those functions for more information about how
the optimizers work.

## Value

Each of these functions returns an optimizer function `optim_fun`. The
returned `optim_fun` function has four input arguments: an initial guess
(`guess`), an error function (`fun`), lower bounds (`lower`), and upper
bounds (`upper`). It returns a list with four named elements: `par`,
`convergence`, `feval`, and `convergence_msg`.

## Examples

``` r
# Here we just show examples of the optim_fun results. Other examples using the
# optimizers can be found throughout PhotoGEA, such as in the user guides and
# the documentation for fit_c3_aci, fit_c4_aci, etc.

optimizer_deoptim(200)
#> function (guess, fun, lower, upper) 
#> {
#>     NP <- 10 * length(guess)
#>     initialpop <- matrix(nrow = NP, ncol = length(guess))
#>     initialpop[1, ] <- guess
#>     varsize <- 0.25
#>     for (i in seq(2, NP)) {
#>         tmp <- guess * (1 + stats::runif(length(guess), -varsize, 
#>             varsize))
#>         initialpop[i, ] <- constrain_guess(tmp, lower, upper, 
#>             0)
#>     }
#>     res <- DEoptim::DEoptim(fun, lower, upper, control = list(VTR = VTR, 
#>         itermax = itermax, NP = NP, initialpop = initialpop, 
#>         trace = FALSE))
#>     list(convergence = NA, convergence_msg = NA, feval = res[["optim"]][["nfeval"]], 
#>         par = res[["optim"]][["bestmem"]], optimizer = "optimizer_deoptim")
#> }
#> <bytecode: 0x5561652676b8>
#> <environment: 0x5561623d48b0>

optimizer_hjkb(1e-7)
#> function (guess, fun, lower, upper) 
#> {
#>     guess <- constrain_guess(guess, lower, upper, 0.01)
#>     res <- dfoptim::hjkb(guess, fun, lower, upper, control = list(tol = tol, 
#>         maxfeval = maxfeval, target = target))
#>     list(convergence = res[["convergence"]], convergence_msg = NA, 
#>         feval = res[["feval"]], par = res[["par"]], optimizer = "optimizer_hjkb")
#> }
#> <bytecode: 0x5561689de020>
#> <environment: 0x5561689dd140>

optimizer_nlminb(1e-7)
#> function (guess, fun, lower, upper) 
#> {
#>     guess <- constrain_guess(guess, lower, upper, 0.01)
#>     res <- stats::nlminb(guess, fun, lower = lower, upper = upper, 
#>         control = list(rel.tol = rel.tol, eval.max = eval.max, 
#>             iter.max = iter.max, abs.tol = abs.tol))
#>     list(convergence = res[["convergence"]], convergence_msg = res[["message"]], 
#>         feval = res[["evaluations"]][1], par = res[["par"]], 
#>         optimizer = "optimizer_nlminb")
#> }
#> <bytecode: 0x55616377aef8>
#> <environment: 0x556164aa96c8>

optimizer_nmkb(1e-7)
#> function (guess, fun, lower, upper) 
#> {
#>     guess <- constrain_guess(guess, lower, upper, 0.01)
#>     res <- dfoptim::nmkb(guess, fun, lower, upper, control = list(tol = tol, 
#>         maxfeval = maxfeval, restarts.max = restarts.max))
#>     list(convergence = res[["convergence"]], convergence_msg = res[["message"]], 
#>         feval = res[["feval"]], par = res[["par"]], optimizer = "optimizer_nmkb")
#> }
#> <bytecode: 0x55616493d700>
#> <environment: 0x556167af58a0>

optimizer_null()
#> function (guess, fun, lower, upper) 
#> {
#>     list(convergence = NA, convergence_msg = NA, feval = 0, par = guess, 
#>         optimizer = "optimizer_null")
#> }
#> <bytecode: 0x556163f25af0>
#> <environment: 0x556166927e58>
```
