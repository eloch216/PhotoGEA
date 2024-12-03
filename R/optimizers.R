optimizer_nmkb <- function(tol, maxfeval = 2000, restarts.max = 10) {
    function(guess, fun, lower, upper) {
        guess <- constrain_guess(guess, lower, upper, 0.01)
        res <- dfoptim::nmkb(guess, fun, lower, upper, control = list(
            tol = tol,
            maxfeval = maxfeval,
            restarts.max = restarts.max
        ))

        # Return the results as a list with the proper names
        list(
            convergence = res[['convergence']],
            convergence_msg = res[['message']],
            feval = res[['feval']],
            par = res[['par']],
            optimizer = 'optimizer_nmkb'
        )
    }
}

optimizer_hjkb <- function(tol, maxfeval = Inf, target = Inf) {
    function(guess, fun, lower, upper) {
        guess <- constrain_guess(guess, lower, upper, 0.01)
        res <- dfoptim::hjkb(guess, fun, lower, upper, control = list(
            tol = tol,
            maxfeval = maxfeval,
            target = target
        ))

        # Return the results as a list with the proper names
        list(
            convergence = res[['convergence']],
            convergence_msg = NA, # this optimizer does not report messages
            feval = res[['feval']],
            par = res[['par']],
            optimizer = 'optimizer_hjkb'
        )
    }
}

optimizer_nlminb <- function(rel.tol, eval.max = 200, iter.max = 200, abs.tol = 0) {
    function(guess, fun, lower, upper) {
        guess <- constrain_guess(guess, lower, upper, 0.01)
        res <- stats::nlminb(guess, fun, lower = lower, upper = upper, control = list(
            rel.tol = rel.tol,
            eval.max = eval.max,
            iter.max = iter.max,
            abs.tol = abs.tol
        ))

        # Return results as a list with the proper names
        list(
            convergence = res[['convergence']],
            convergence_msg = res[['message']],
            feval = res[['evaluations']][1],
            par = res[['par']],
            optimizer = 'optimizer_nlminb'
        )
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
            convergence = NA,     # this optimizer does not check for convergence
            convergence_msg = NA, # this optimizer does not report messages
            feval = res[['optim']][['nfeval']],
            par = res[['optim']][['bestmem']],
            optimizer = 'optimizer_deoptim'
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

# A helping function to confirm that the results from an optimizer have the
# correct elements
check_optim_result <- function(optim_res) {
    required_length_one <- c(
        'convergence',
        'convergence_msg',
        'feval',
        'optimizer'
    )

    required <- c(required_length_one, 'par')

    if (!all(required %in% names(optim_res))) {
        msg <- paste0(
            'The optimizer result must include the following elements: ',
            paste(required, collapse = ', '),
            '. Found the following elements: ',
            paste(names(optim_res), collapse = ', ')
        )
        stop(msg)
    }

    length_okay <- sapply(required_length_one, function(x) {length(optim_res[[x]]) == 1})

    if (!all(length_okay)) {
        bad_length <- names(optim_res)[!length_okay]
        msg <- paste0(
            'The following optimizer outputs must have length 1, but do not: ',
            paste(bad_length, collapse = ', ')
        )
        stop(msg)
    }
}
