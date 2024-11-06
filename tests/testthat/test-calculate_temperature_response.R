# Define a simple exdf object
test_exdf <- exdf(
    data.frame(
        TleafCnd = c(20, 30),
        not_used = NA
    ),
    units = data.frame(
        TleafCnd = 'degrees C',
        tmp = '',
        stringsAsFactors = FALSE
    )
)

# Choose test tolerance
TOLERANCE <- 1e-4

test_that('temperature response is calculated properly', {
    tr_res <- calculate_temperature_response(
        test_exdf,
        list(
          Kc = list(type = 'Arrhenius', c = 38.05, Ea = 79.43, units = 'micromol mol^(-1)'),
          Jmax = list(type = 'Gaussian', optimum_rate = 1, t_opt = 43, sigma = 16, units = 'micromol m^(-2) s^(-1)')
        )
    )

    expect_equal(
        as.numeric(tr_res[, 'Kc']),
        c(235.5536, 690.1576),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(tr_res[, 'Jmax']),
        c(0.1266401, 0.5167706),
        tolerance = TOLERANCE
    )
})

test_that('mistakes are caught', {
    expect_error(
        calculate_temperature_response(
            test_exdf,
            list(Kc = list(c = 38.05))
        ),
        'Temperature response parameter set named `Kc` does not specify a `type` value'
    )

    expect_error(
        calculate_temperature_response(
            test_exdf,
            list(Kc = list(type = 'bad_type', units = 'micromol mol^(-1)'))
        ),
        'Temperature response parameter set named `Kc` specifies an unsupported `type` value: `bad_type`. The available options are: arrhenius, gaussian.'
    )

    expect_error(
        calculate_temperature_response(
            test_exdf,
            list(Kc = list(type = 'Arrhenius', c = 38.05, Ea = 79.43))
        ),
        'Temperature response parameter set named `Kc` does not specify a `units` value'
    )

    expect_error(
        calculate_temperature_response(
            test_exdf,
            list(Kc = list(type = 'Arrhenius', c = 38.05, units = 'micromol mol^(-1)'))
        ),
        'Arrhenius parameter named `Kc` has the following elements: type, c, units; elements named `c`, `Ea`, and `units` are required.'
    )

    expect_error(
        calculate_temperature_response(
            test_exdf,
            list(Jmax = list(type = 'Gaussian', t_opt = 43, units = 'micromol m^(-2) s^(-1)'))
        ),
        'Gaussian parameter named `Jmax` has the following elements: type, t_opt, units; elements named `optimum_rate`, `t_opt`, `sigma`, and `units` are required.'
    )
})
