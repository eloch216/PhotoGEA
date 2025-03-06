# A helping function for calculating Jmax
jmax_from_j <- function(J, I2, theta) {
    J * (I2 - theta * J) / (I2 - J)
}

# A helping function for creating error messages
jmax_error_msg <- function(J, I2) {
    if (length(J) != length(I2)) {
        stop('Length mismatch between J and I2')
    }

    msg <- rep_len('', length(J))

    msg[is.na(J)] <- 'J is NA'
    msg[J >= I2]  <- 'J >= I2'

    msg
}

calculate_jmax <- function(
    data_table,
    alpha_j_at_25 = 0.293, # dimensionless
    theta_j_at_25 = 0.979, # dimensionless
    alpha_j_norm_column_name = 'alpha_j_norm',
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

    required_variables[['J_at_25']]                <- unit_dictionary('J_at_25')
    required_variables[['J_tl_avg']]               <- unit_dictionary('J_tl_avg')
    required_variables[[alpha_j_norm_column_name]] <- unit_dictionary('alpha_j_norm')
    required_variables[[qin_column_name]]          <- unit_dictionary('Qin')
    required_variables[[theta_j_norm_column_name]] <- unit_dictionary('theta_j_norm')
    required_variables[[tleaf_column_name]]        <- unit_dictionary('Tleaf_avg')

    optional_param <- list(
        J_at_25_lower  = 'J_at_25_lower',
        J_at_25_upper  = 'J_at_25_upper',
        J_tl_avg_lower = 'J_tl_avg_lower',
        J_tl_avg_upper = 'J_tl_avg_upper'
    )

    required_variables <-
        require_optional_param(required_variables, optional_param, data_table)

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

    Qin <- data_table[, qin_column_name] # micromol / m^2 / s

    # Make sure key inputs have reasonable values
    msg <- character()

    if (any(alpha_j_at_25 < 0, na.rm = TRUE))                      {msg <- append(msg, 'alpha_j_at_25 must be >= 0')}
    if (any(alpha_j_tl < 0, na.rm = TRUE))                         {msg <- append(msg, 'alpha_j_tl must be >= 0')}
    if (any(theta_j_at_25 <= 0 | theta_j_at_25 > 1, na.rm = TRUE)) {msg <- append(msg, 'theta_j_at_25 must be > 0 and <= 1')}
    if (any(theta_j_tl <= 0 | theta_j_tl > 1, na.rm = TRUE))       {msg <- append(msg, 'theta_j_tl must be > 0 and <= 1')}

    msg <- paste(msg, collapse = '. ')

    if (msg != '') {
        stop(msg)
    }

    # Specify lists of Jmax types to calculate
    to_calculate <- list(
        at_25 = list(J = data_table[, 'J_at_25'],  alpha_j = alpha_j_at_25, theta_j = theta_j_at_25, include_I2 = TRUE),
        tl    = list(J = data_table[, 'J_tl_avg'], alpha_j = alpha_j_tl,    theta_j = theta_j_tl,    include_I2 = TRUE)
    )

    if ('J_at_25_lower' %in% colnames(data_table)) {
        to_calculate$at_25_lower <- list(
            J = data_table[, 'J_at_25_lower'],
            alpha_j = alpha_j_at_25,
            theta_j = theta_j_at_25,
            include_I2 = FALSE
        )
    }

    if ('J_at_25_upper' %in% colnames(data_table)) {
        to_calculate$at_25_upper <- list(
            J = data_table[, 'J_at_25_upper'],
            alpha_j = alpha_j_at_25,
            theta_j = theta_j_at_25,
            include_I2 = FALSE
        )
    }

    if ('J_tl_avg_lower' %in% colnames(data_table)) {
        to_calculate$tl_lower <- list(
            J = data_table[, 'J_tl_avg_lower'],
            alpha_j = alpha_j_tl,
            theta_j = theta_j_tl,
            include_I2 = FALSE
        )
    }

    if ('J_tl_avg_upper' %in% colnames(data_table)) {
        to_calculate$tl_upper <- list(
            J = data_table[, 'J_tl_avg_upper'],
            alpha_j = alpha_j_tl,
            theta_j = theta_j_tl,
            include_I2 = FALSE
        )
    }

    # Run calculations for each type of Jmax
    for (i in seq_along(to_calculate)) {
        # Get info
        x      <- to_calculate[[i]]
        x_name <- names(to_calculate)[[i]]

        # Get I2
        I2 <- Qin * x$alpha_j

        # Get Jmax
        Jmax <- jmax_from_j(x$J, I2, x$theta_j)

        # Check for restrictions
        Jmax_msg <- jmax_error_msg(x$J, I2)

        # Apply restrictions if necessary
        if (!ignore_restriction) {
            Jmax[Jmax_msg != ''] <- NA
        }

        # Add results to table
        if (x$include_I2) {
            data_table <- set_variable(
                data_table,
                paste0('I2_', x_name),
                unit_dictionary('I2'),
                'calculate_jmax',
                I2
            )
        }

        data_table <- set_variable(
            data_table,
            paste0('Jmax_', x_name),
            unit_dictionary('J'),
            'calculate_jmax',
            Jmax
        )

        data_table <- set_variable(
            data_table,
            paste0('Jmax_', x_name, '_msg'),
            '',
            'calculate_jmax',
            Jmax_msg
        )
    }

    # Return the results
    return(data_table)
}
