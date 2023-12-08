estimate_operating_point <- function(
    aci_exdf,
    ca_atmospheric,
    type = 'c3',
    a_column_name = 'A',
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    pcm_column_name = 'PCm',
    return_list = FALSE
)
{
    if (!is.exdf(aci_exdf)) {
        stop('estimate_operating_point requires an exdf object')
    }

    type <- tolower(type)

    if (!type %in% c('c3', 'c4')) {
        stop('`type` must be "c3" or "c4"')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]  <- 'micromol m^(-2) s^(-1)'
    required_variables[[ca_column_name]] <- 'micromol mol^(-1)'
    required_variables[[ci_column_name]] <- 'micromol mol^(-1)'

    if (type == 'c3') {
        required_variables[[cc_column_name]] <- 'micromol mol^(-1)'
    } else {
        required_variables[[pcm_column_name]] <- 'microbar'
    }

    check_required_variables(aci_exdf, required_variables)

    # Make sure the atmospheric Ca is included in the Ca range
    min_ca <- min(aci_exdf[, ca_column_name])
    max_ca <- max(aci_exdf[, ca_column_name])

    unreliable <- ca_atmospheric < min_ca || ca_atmospheric > max_ca
    if (unreliable) {
        warning(
            'The atmospheric CO2 concentration (', ca_atmospheric,
            ') is outside the measured range (', min_ca, ' - ', max_ca, ')',
            call. = FALSE
        )
    }

    # Use linear interpolation to estimate the operating point
    operating_An <- if (!unreliable) {
        stats::approx(
            aci_exdf[, ca_column_name],
            aci_exdf[, a_column_name],
            ca_atmospheric
        )$y
    } else {
        NA
    }

    operating_Ci <- if (!unreliable) {
        stats::approx(
            aci_exdf[, ca_column_name],
            aci_exdf[, ci_column_name],
            ca_atmospheric
        )$y
    } else {
        NA
    }

    operating_Cc <- if (!unreliable && type == 'c3') {
        stats::approx(
            aci_exdf[, ca_column_name],
            aci_exdf[, cc_column_name],
            ca_atmospheric
        )$y
    } else {
        NA
    }

    operating_PCm <- if (!unreliable && type == 'c4') {
        stats::approx(
            aci_exdf[, ca_column_name],
            aci_exdf[, pcm_column_name],
            ca_atmospheric
        )$y
    } else {
        NA
    }

    if (return_list) {
        # Prepare an exdf that can be used for calculating the operating An
        # using fit parameters
        if (type == 'c3') {
            cc_seq <- aci_exdf[, cc_column_name]

            operating_row <- if (!unreliable) {
                which(abs(cc_seq - operating_Cc) == min(abs(cc_seq - operating_Cc)))
            } else {
                1
            }

            operating_exdf <- aci_exdf[operating_row, , TRUE]

            operating_exdf[, cc_column_name] <- operating_Cc
            operating_exdf[, ci_column_name] <- operating_Ci

            list(
                operating_An = operating_An,
                operating_Cc = operating_Cc,
                operating_Ci = operating_Ci,
                operating_exdf = operating_exdf
            )
        } else {
            pcm_seq <- aci_exdf[, pcm_column_name]

            operating_row <- if (!unreliable) {
                which(abs(pcm_seq - operating_PCm) == min(abs(pcm_seq - operating_PCm)))
            } else {
                1
            }

            operating_exdf <- aci_exdf[operating_row, , TRUE]

            operating_exdf[, pcm_column_name] <- operating_PCm
            operating_exdf[, ci_column_name] <- operating_Ci

            list(
                operating_An = operating_An,
                operating_PCm = operating_PCm,
                operating_Ci = operating_Ci,
                operating_exdf = operating_exdf
            )
        }
    } else {
        # Get the replicate identifier columns
        aci_identifiers <- identifier_columns(aci_exdf)

        # Store the results
        aci_identifiers[, 'operating_An'] <- operating_An
        aci_identifiers[, 'operating_Ci'] <- operating_Ci

        if (type == 'c3') {
            aci_identifiers[, 'operating_Cc'] <- operating_Cc
        } else {
            aci_identifiers[, 'operating_PCm'] <- operating_PCm
        }

        # Document the new columns that were added and return the exdf
        aci_identifiers <- document_variables(
            aci_identifiers,
            c('estimate_operating_point', 'operating_An', aci_exdf$units[[a_column_name]]),
            c('estimate_operating_point', 'operating_Ci', aci_exdf$units[[ci_column_name]])
        )

        if (type == 'c3') {
            document_variables(
                aci_identifiers,
                c('estimate_operating_point', 'operating_Cc', aci_exdf$units[[cc_column_name]])
            )
        } else {
            document_variables(
                aci_identifiers,
                c('estimate_operating_point', 'operating_PCm', aci_exdf$units[[pcm_column_name]])
            )
        }
    }
}
