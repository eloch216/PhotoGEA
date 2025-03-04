# A helping function for calculating Jmax
jmax_from_j <- function(J, I2, theta) {
    J * (I2 - theta * J) / (I2 - J)
}

# A helping function for creating error messages
jmax_error_msg <- function(J, I2) {
    error_cond <- I2 < J

    sapply(error_cond, function(ec) {
        if (ec) {
            'I2 < J'
        } else {
            ''
        }
    })
}

calculate_jmax <- function(
    data_table,
    alpha_j_at_25, # dimensionless
    theta_j_at_25, # dimensionless
    alpha_j_norm_column_name = 'alpha_j_norm',
    j_at_25_column_name = 'J_at_25',
    j_tl_column_name = 'J_tl_avg',
    qin_column_name = 'Qin_avg',
    theta_j_norm_column_name = 'theta_j_norm',
    tleaf_column_name = 'TleafCnd_avg',
    ...
)
{
    # Get optional arguments
    optional_args <- list(...)

    potential_optional_args <- c(
        'ignore_restriction'
    )

    check_optional_arguments(optional_args, potential_optional_args)

    ignore_restriction <- get_optional_argument(optional_args, 'ignore_restriction', FALSE)

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[alpha_j_norm_column_name]] <- unit_dictionary('alpha_j_norm')
    required_variables[[j_at_25_column_name]]      <- unit_dictionary('J_at_25')
    required_variables[[j_tl_column_name]]         <- unit_dictionary('J_tl_avg')
    required_variables[[qin_column_name]]          <- unit_dictionary('Qin')
    required_variables[[theta_j_norm_column_name]] <- unit_dictionary('theta_j_norm')
    required_variables[[tleaf_column_name]]        <- unit_dictionary('Tleaf_avg')

    flexible_param <- list(
        alpha_j_at_25 = alpha_j_at_25,
        theta_j_at_25 = theta_j_at_25
    )

    required_variables <-
        require_flexible_param(required_variables, flexible_param)

    check_required_variables(data_table, required_variables)

    # Retrieve values of flexible parameters as necessary
    if (!value_set(alpha_j_at_25)) {alpha_j_at_25 <- data_table[, 'alpha_j_at_25']}
    if (!value_set(theta_j_at_25)) {theta_j_at_25 <- data_table[, 'theta_j_at_25']}

    # Extract a few columns from the exdf object to make the equations easier to
    # read
    alpha_j_tl <- alpha_j_at_25 * data_table[, alpha_j_norm_column_name] # dimensionless
    theta_j_tl <- theta_j_at_25 * data_table[, theta_j_norm_column_name] # dimensionless

    J_at_25 <- data_table[, j_at_25_column_name] # micromol / m^2 / s
    J_tl    <- data_table[, j_tl_column_name]    # micromol / m^2 / s

    Qin <- data_table[, qin_column_name] # micromol / m^2 / s

    # Make sure key inputs have reasonable values
    msg <- character()

    if (any(alpha_j_at_25 < 0, na.rm = TRUE))                     {msg <- append(msg, 'alpha_j_at_25 must be >= 0')}
    if (any(alpha_j_tl < 0, na.rm = TRUE))                        {msg <- append(msg, 'alpha_j_tl must be >= 0')}
    if (any(theta_j_at_25 < 0 | theta_j_at_25 > 1, na.rm = TRUE)) {msg <- append(msg, 'theta_j_at_25 must be >= 0 and <= 1')}
    if (any(theta_j_tl < 0 | theta_j_tl > 1, na.rm = TRUE))       {msg <- append(msg, 'theta_j_tl must be >= 0 and <= 1')}

    msg <- paste(msg, collapse = '. ')

    if (msg != '') {
        stop(msg)
    }

    # Calculate Jmax at 25 degrees C
    I2_at_25   <- Qin * alpha_j_at_25                           # micromol / m^2 / s
    Jmax_at_25 <- jmax_from_j(J_at_25, I2_at_25, theta_j_at_25) # micromol / m^2 / s

    # Calculate Jmax at leaf temperature
    I2_tl   <- Qin * alpha_j_tl                     # micromol / m^2 / s
    Jmax_tl <- jmax_from_j(J_tl, I2_tl, theta_j_tl) # micromol / m^2 / s

    # Check for restrictions
    Jmax_at_25_msg <- jmax_error_msg(J_at_25, I2_at_25)
    Jmax_tl_msg    <- jmax_error_msg(J_tl,    I2_tl)

    if (!ignore_restriction) {
        Jmax_at_25[Jmax_at_25_msg != ''] <- NA
        Jmax_tl[Jmax_tl_msg != '']       <- NA
    }

    # Add calculated variables and units to the table, and return it
    data_table[, 'I2_at_25']       <- I2_at_25
    data_table[, 'I2_tl']          <- I2_tl
    data_table[, 'Jmax_at_25']     <- Jmax_at_25
    data_table[, 'Jmax_at_25_msg'] <- Jmax_at_25_msg
    data_table[, 'Jmax_tl']        <- Jmax_tl
    data_table[, 'Jmax_tl_msg']    <- Jmax_tl_msg

    document_variables(
        data_table,
        c('calculate_jmax', 'I2_at_25',       unit_dictionary('I2_at_25')),
        c('calculate_jmax', 'I2_tl',          unit_dictionary('I2_tl')),
        c('calculate_jmax', 'Jmax_at_25',     unit_dictionary('Jmax_at_25')),
        c('calculate_jmax', 'Jmax_at_25_msg', ''),
        c('calculate_jmax', 'Jmax_tl',        unit_dictionary('Jmax_tl')),
        c('calculate_jmax', 'Jmax_tl_msg',    '')
    )

}
