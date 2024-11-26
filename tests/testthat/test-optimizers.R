# Get test curves to use
source('one_curve_c3_aci.R')

# Specify an infinite mesophyll conductance (so `Cc` = `Ci`)
one_curve <- set_variable(
  one_curve,
  'gmc', 'mol m^(-2) s^(-1) bar^(-1)', value = Inf
)

# Calculate Cc
one_curve <- apply_gm(one_curve)

test_that('the nmkb optimizer works', {
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c3_aci(
            one_curve,
            Ca_atmospheric = 420,
            optim_fun = optimizer_nmkb(1e-5)
        )
    )

    expect_true(is.finite(fit_res$parameters[1, 'RL_at_25']))
})

test_that('the hjkb optimizer works', {
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c3_aci(
            one_curve,
            Ca_atmospheric = 420,
            optim_fun = optimizer_hjkb(1e-5)
        )
    )

    expect_true(is.finite(fit_res$parameters[1, 'RL_at_25']))
})

test_that('the deoptim optimizer works', {
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c3_aci(
            one_curve,
            Ca_atmospheric = 420,
            optim_fun = optimizer_deoptim(20)
        )
    )

    expect_true(is.finite(fit_res$parameters[1, 'RL_at_25']))
})
