# This function creates a function that calculates the temperature-dependent
# values of parameters whose Arrhenius parameters are defined in `trf_param`.
photosynthesis_TRF <- function(trf_param) {
    function(leaf_temperature) {
        lapply(
            trf_param,
            function(x) {arrhenius(x$c, x$Ea, leaf_temperature)}
        )
    }
}
