optimizer_nmkb <- function(tol = 1e-7, maxfeval = 2000, restarts.max = 10) {
    function(guess, fun, lower, upper) {
        dfoptim::nmkb(guess, fun, lower, upper, control = list(
            tol = tol,
            maxfeval = maxfeval,
            restarts.max = restarts.max
        ))
    }
}

optimizer_deoptim <- function(VTR = -Inf, itermax = 200) {
  function(guess, fun, lower, upper) {
    # Set population size following DEoptim default behavior
    NP <- 10 * length(guess)

    # Initialize the population with parameter values close to guess
    initialpop <- matrix(nrow = NP, ncol = length(guess))

    initialpop[1, ] <- guess

    varsize <- 0.25
    for (i in seq(2, NP)) {
      initialpop[i, ] <-
        guess * (1 + stats::runif(length(guess), -varsize, varsize))
    }

    # Set control options and run the optimizer
    res <- DEoptim::DEoptim(fun, lower, upper, control = list(
      VTR = VTR,
      itermax = itermax,
      NP = NP,
      initialpop = initialpop,
      trace = FALSE
    ))

    # Return results as a list with the proper names
    list(
      par = res$optim$bestmem,
      convergence = NA,
      value = res$optim$bestval,
      feval = res$optim$nfeval
    )
  }
}
