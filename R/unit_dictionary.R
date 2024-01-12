# Specify units for some important parameters
unit_dictionary <- list(
    alpha       = 'dimensionless',
    J_at_25     = 'micromol m^(-2) s^(-1)',
    Rd_at_25    = 'micromol m^(-2) s^(-1)',
    tau         = 'dimensionless',
    TPU         = 'micromol m^(-2) s^(-1)',
    Vcmax_at_25 = 'micromol m^(-2) s^(-1)',
    Vpmax_at_25 = 'micromol m^(-2) s^(-1)',
    Vpr         = 'micromol m^(-2) s^(-1)'
)

# A helping function for adding a parameters (and their units) to a list of
# required parameters.
add_required_parameters <- function(required_variables, param_names) {
    for (pn in param_names) {
         required_variables[[pn]] <- unit_dictionary[[pn]]
    }
    required_variables
}
