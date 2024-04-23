# Calculate the AIC: https://en.wikipedia.org/wiki/Akaike_information_criterion
akaike_information_criterion <- function(log_likelihood, nparam) {
    2 * nparam - 2 * log_likelihood
}
