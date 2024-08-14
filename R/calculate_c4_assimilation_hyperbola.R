calculate_c4_assimilation_hyperbola <- function(
    exdf_obj,
    c4_curvature,   # dimensionless      (typically this value is being fitted)
    c4_slope,       # mol / m^2 / s      (typically this value is being fitted)
    rL,             # micromol / m^2 / s (typically this value is being fitted)
    Vmax,           # micromol / m^2 / s (typically this value is being fitted)
    ci_column_name = 'Ci',
    hard_constraints = 0,
    perform_checks = TRUE,
    return_exdf = TRUE
)
{
    if (perform_checks) {
        if (!is.exdf(exdf_obj)) {
            stop('calculate_c4_assimilation_hyperbola requires an exdf object')
        }

        # Make sure the required variables are defined and have the correct units
        required_variables <- list()
        required_variables[[ci_column_name]] <- unit_dictionary[['Ci']]

        flexible_param <- list(
            c4_curvature = c4_curvature,
            rL = rL,
            c4_slope = c4_slope,
            Vmax = Vmax
        )

        required_variables <-
            require_flexible_param(required_variables, flexible_param)

        check_required_variables(exdf_obj, required_variables)
    }

    # Retrieve values of flexible parameters as necessary
    if (!value_set(c4_curvature))   {c4_curvature   <- exdf_obj[, 'c4_curvature']}
    if (!value_set(c4_slope))       {c4_slope       <- exdf_obj[, 'c4_slope']}
    if (!value_set(rL))             {rL             <- exdf_obj[, 'rL']}
    if (!value_set(Vmax))           {Vmax           <- exdf_obj[, 'Vmax']}

    # Extract a few columns from the exdf object to make the equations easier to
    # read, converting units as necessary
    Ci <- exdf_obj[, ci_column_name]        # micromol / mol

    # Make sure key inputs have reasonable values
    msg <- character()

    # Optionally check whether Ci is reasonable
    if (hard_constraints >= 1) {
        if (any(Ci < 0, na.rm = TRUE)) {msg <- append(msg, 'Ci must be >= 0')}
    }

    # Optionally check reasonableness of parameters that can be fit
    if (hard_constraints >= 2) {
        if (any(c4_curvature < 0 | c4_curvature > 1, na.rm = TRUE)) {msg <- append(msg, 'c4_curvature must be >= 0 and <= 1')}
        if (any(c4_slope < 0, na.rm = TRUE))                        {msg <- append(msg, 'c4_slope must be >= 0')}
        if (any(rL < 0, na.rm = TRUE))                              {msg <- append(msg, 'rL must be >= 0')}
        if (any(Vmax < 0, na.rm = TRUE))                            {msg <- append(msg, 'Vmax must be >= 0')}
    }

    msg <- paste(msg, collapse = '. ')

    # We only bypass these checks if !perform_checks && return_exdf
    if (perform_checks || !return_exdf) {
        if (msg != '') {
            stop(msg)
        }
    }

    # Make sure c4_curvature is a vector with the same length as Ci
    c4_curvature <- rep_len(c4_curvature, length(Ci))

    # Calculate gross assimilation, which is the negative root of a quadratic
    # equation
    Vinitial <- Ci * c4_slope # micromol / m^2 / s

    qa <- c4_curvature       # dimensionless
    qb <- -(Vinitial + Vmax) # micromol / m^2 / s
    qc <- Vinitial * Vmax    # micromol / m^2 / s

    Ag <- sapply(seq_along(qb), function(i) {
        quadratic_root_minus(qa[i], qb[i], qc[i]) # micromol / m^2 / s
    })

    # Calculate net assimilation rates
    Ainitial <- Vinitial - rL # micromol / m^2 / s
    Amax     <- Vmax - rL     # micromol / m^2 / s
    An       <- Ag - rL       # micromol / m^2 / s

    if (return_exdf) {
        # Make a new exdf object from the calculated variables and make sure units
        # are included
        output <- exdf(data.frame(
            Ag = Ag,
            Ainitial = Ainitial,
            Amax = Amax,
            An = An,
            c4_curvature = c4_curvature,
            c4_slope = c4_slope,
            rL = rL,
            Vinitial = Vinitial,
            Vmax = Vmax,
            c4_assimilation_hyperbola_msg = msg,
            stringsAsFactors = FALSE
        ))

        document_variables(
            output,
            c('calculate_c4_assimilation_hyperbola', 'Ag',                            unit_dictionary$Vmax),
            c('calculate_c4_assimilation_hyperbola', 'Ainitial',                      unit_dictionary$Ainitial),
            c('calculate_c4_assimilation_hyperbola', 'Amax',                          unit_dictionary$Vmax),
            c('calculate_c4_assimilation_hyperbola', 'An',                            unit_dictionary$Vmax),
            c('calculate_c4_assimilation_hyperbola', 'c4_curvature',                  unit_dictionary$c4_curvature),
            c('calculate_c4_assimilation_hyperbola', 'rL',                            unit_dictionary$rL),
            c('calculate_c4_assimilation_hyperbola', 'Vinitial',                      unit_dictionary$Vmax),
            c('calculate_c4_assimilation_hyperbola', 'Vmax',                          unit_dictionary$Vmax),
            c('calculate_c4_assimilation_hyperbola', 'c4_assimilation_hyperbola_msg', '')
        )
    } else {
        return(An)
    }
}
