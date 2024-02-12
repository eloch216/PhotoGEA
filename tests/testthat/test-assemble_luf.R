# Set some parameter names to use for testing
test_param_names <- c('alpha_g', 'Rd', 'Vcmax', 'TPU')

test_that('reasonable inputs work properly', {
    res <- PhotoGEA:::assemble_luf(
        test_param_names,
        default_lower = list(alpha_g = 0, Rd = 0, Vcmax = 0, TPU = 0),
        default_upper = list(TPU = 100, alpha_g = 1, Rd = 10, Vcmax = 1000),
        default_fit_options = list(alpha_g = 0, Rd = 'fit', Vcmax = 'column', TPU = 40),
        user_lower = list(),
        user_upper = list(Vcmax = 400),
        user_fit_options = list(Vcmax = 'FIT', alpha_g = 0, Rd = 'column')
    )

    expect_true('lower' %in% names(res))
    expect_true('upper' %in% names(res))
    expect_true('fit_options' %in% names(res))
    expect_true('fit_options_vec' %in% names(res))
    expect_true('param_to_fit' %in% names(res))

    expect_true(is.numeric(res$lower))
    expect_true(is.numeric(res$upper))
    expect_true(is.list(res$fit_options))
    expect_true(is.numeric(res$fit_options_vec))
    expect_true(is.logical(res$param_to_fit))

    expect_equal(names(res$fit_options), test_param_names)

    expect_equal(res$lower, c(0, 0, 0, 0))
    expect_equal(res$upper, c(1, 10, 400, 100))
    expect_equal(res$fit_options_vec, c(0, NA, NA, 40))

    expect_equal(res$fit_options$alpha_g, 0)
    expect_equal(res$fit_options$Rd, 'column')
    expect_equal(res$fit_options$Vcmax, 'fit')
    expect_equal(res$fit_options$TPU, 40)
})

test_that('luf should be lists', {
    expect_error(
        PhotoGEA:::assemble_luf(
            test_param_names,
            default_lower = list(),
            default_upper = list(),
            default_fit_options = list(),
            user_lower = list(alpha_g = 0, Rd = 0, Vcmax = 0),
            user_upper = c(alpha_g = 1, Rd = 10, Vcmax = 1000),
            user_fit_options = list(alpha_g = 'fit', Rd = 'column', Vcmax = 100)
        ),
        '`lower`, `upper`, and `fit_options` must be lists'
    )
})

test_that('luf should be have names', {
    expect_error(
        PhotoGEA:::assemble_luf(
            test_param_names,
            default_lower = list(),
            default_upper = list(),
            default_fit_options = list(),
            user_lower = list(alpha_g = 0, Rd = 0, Vcmax = 0),
            user_upper = list(1, 10, 1000),
            user_fit_options = list(alpha_g = 'fit', Rd = 'column', Vcmax = 100)
        ),
        '`lower`, `upper`, and `fit_options` must have named elements if they are not empty lists'
    )
})

test_that('luf names cannot be empty', {
    expect_error(
        PhotoGEA:::assemble_luf(
            test_param_names,
            default_lower = list(),
            default_upper = list(),
            default_fit_options = list(),
            user_lower = list(alpha_g = 0, Rd = 0, Vcmax = 0),
            user_upper = list(alpha_g = 1, Rd = 10, Vcmax = 1000),
            user_fit_options = list(alpha_g = 'fit', 'column', Vcmax = 100)
        ),
        '`lower`, `upper`, and `fit_options` must not have any empty names'
    )
})

test_that('luf names must be in the parameter set', {
    expect_error(
        PhotoGEA:::assemble_luf(
            test_param_names,
            default_lower = list(),
            default_upper = list(),
            default_fit_options = list(),
            user_lower = list(alpha_g = 0, Rd = 0, Vcmax = 0),
            user_upper = list(alpha_g = 1, Rd = 10, Vcmax = 1000),
            user_fit_options = list(alpha_g = 'fit', J = 'column', Vcmax = 100)
        ),
        '`lower`, `upper`, and `fit_options` must only provide settings for `alpha_g`, `Rd`, `Vcmax`, `TPU`'
    )
})


test_that('fit options must be `fit`, `column`, or numeric', {
    expect_error(
        PhotoGEA:::assemble_luf(
            test_param_names,
            default_lower = list(),
            default_upper = list(),
            default_fit_options = list(),
            user_lower = list(alpha_g = 0, Rd = 0, Vcmax = 0),
            user_upper = list(alpha_g = 1, Rd = 10, Vcmax = 1000),
            user_fit_options = list(alpha_g = 'WRONG VALUE', Rd = 'column', Vcmax = 100)
        ),
        'Each element of `fit_options` must be `fit`, `column`, or a numeric value'
    )
})

test_that('at least one parameter must be fit', {
    expect_error(
            PhotoGEA:::assemble_luf(
            test_param_names,
            default_lower = list(alpha_g = 0, Rd = 0, Vcmax = 0, TPU = 0),
            default_upper = list(TPU = 100, alpha_g = 1, Rd = 10, Vcmax = 1000),
            default_fit_options = list(alpha_g = 0, Rd = 'fit', Vcmax = 'column', TPU = 40),
            user_lower = list(),
            user_upper = list(Vcmax = 400),
            user_fit_options = list(alpha_g = 0, Rd = 0, Vcmax = 0, TPU = 0)
        ),
        'No entries in `fit_options` are set to `fit`'
    )
})
