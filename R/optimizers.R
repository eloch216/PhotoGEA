optimizer_nmkb <- function(tol, maxfeval = 2000, restarts.max = 10) {
    function(guess, fun, lower, upper) {
        guess <- constrain_guess(guess, lower, upper, 0.01)
        dfoptim::nmkb(guess, fun, lower, upper, control = list(
            tol = tol,
            maxfeval = maxfeval,
            restarts.max = restarts.max
        ))
    }
}

optimizer_deoptim <- function(itermax, VTR = -Inf) {
  function(guess, fun, lower, upper) {
    # Set population size following DEoptim default behavior
    NP <- 10 * length(guess)

    # Initialize the population with parameter values close to guess
    initialpop <- matrix(nrow = NP, ncol = length(guess))

    initialpop[1, ] <- guess

    varsize <- 0.25
    for (i in seq(2, NP)) {
      tmp <- guess * (1 + stats::runif(length(guess), -varsize, varsize))
      initialpop[i, ] <- constrain_guess(tmp, lower, upper, 0)
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

# A helping function to ensure the guess lies within (or possibly on) the bounds
constrain_guess <- function(guess, lower, upper, pad) {
    lower_temp <- lower + pad * (upper - lower)
    upper_temp <- upper - pad * (upper - lower)

    guess <- pmax(guess, lower_temp)
    pmin(guess, upper_temp)
}
