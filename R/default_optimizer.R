default_optimizer <- function(tol = 1e-7, maxfeval = 2000, restarts.max = 10) {
    function(guess, fun, lower, upper) {
        dfoptim::nmkb(guess, fun, lower, upper, control = list(
            tol = tol,
            maxfeval = maxfeval,
            restarts.max = restarts.max
        ))
    }
}
