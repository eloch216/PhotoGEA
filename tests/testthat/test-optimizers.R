# Get test curves to use
source('one_curve_c3_aci.R')

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

test_that('the nlminb optimizer works', {
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c3_aci(
            one_curve,
            Ca_atmospheric = 420,
            optim_fun = optimizer_nlminb(1e-5)
        )
    )

    expect_true(is.finite(fit_res$parameters[1, 'RL_at_25']))
})

test_that('the null optimizer works', {
    fit_res <- expect_silent(
        fit_c3_aci(
            one_curve,
            Ca_atmospheric = 420,
            optim_fun = optimizer_null()
        )
    )

    expect_true(is.finite(fit_res$parameters[1, 'RL_at_25']))
})

test_that('optimizer outputs are checked', {
    expect_error(
        fit_c3_aci(
            one_curve,
            Ca_atmospheric = 420,
            optim_fun = function(...) {list(a = 1)}
        ),
        'The optimizer result must include the following elements: convergence, convergence_msg, feval, optimizer, par. Found the following elements: a'
    )

    expect_error(
        fit_c3_aci(
            one_curve,
            Ca_atmospheric = 420,
            optim_fun = function(...) {list(convergence = 1, convergence_msg = 1, feval = c(1, 2), optimizer = 'fake', par = c(1, 2))}
        ),
        'The following optimizer outputs must have length 1, but do not: feval'
    )
})
