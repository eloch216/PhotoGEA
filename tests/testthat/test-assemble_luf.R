# Set some parameter names to use for testing
test_param_names <- c('alpha', 'Rd', 'Vcmax', 'TPU')

test_that('reasonable inputs work properly', {
    res <- PhotoGEA:::assemble_luf(
        test_param_names,
        default_lower = list(alpha = 0, Rd = 0, Vcmax = 0, TPU = 0),
        default_upper = list(TPU = 100, alpha = 1, Rd = 10, Vcmax = 1000),
        default_fit_options = list(alpha = 0, Rd = 'fit', Vcmax = 'column', TPU = 40),
        user_lower = list(),
        user_upper = list(Vcmax = 400),
        user_fit_options = list(Vcmax = 'FIT', alpha = 0, Rd = 'column')
    )

    expect_true('lower' %in% names(res))
    expect_true('upper' %in% names(res))
    expect_true('fit_options' %in% names(res))

    expect_equal(names(res$lower), test_param_names)
    expect_equal(names(res$upper), test_param_names)
    expect_equal(names(res$fit_options), test_param_names)

    expect_equal(res$lower$alpha, 0)
    expect_equal(res$lower$Rd, 0)
    expect_equal(res$lower$Vcmax, 0)
    expect_equal(res$lower$TPU, 0)

    expect_equal(res$upper$alpha, 1)
    expect_equal(res$upper$Rd, 10)
    expect_equal(res$upper$Vcmax, 400)
    expect_equal(res$upper$TPU, 100)

    expect_equal(res$fit_options$alpha, 0)
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
            user_lower = list(alpha = 0, Rd = 0, Vcmax = 0),
            user_upper = c(alpha = 1, Rd = 10, Vcmax = 1000),
            user_fit_options = list(alpha = 'fit', Rd = 'column', Vcmax = 100)
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
            user_lower = list(alpha = 0, Rd = 0, Vcmax = 0),
            user_upper = list(1, 10, 1000),
            user_fit_options = list(alpha = 'fit', Rd = 'column', Vcmax = 100)
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
            user_lower = list(alpha = 0, Rd = 0, Vcmax = 0),
            user_upper = list(alpha = 1, Rd = 10, Vcmax = 1000),
            user_fit_options = list(alpha = 'fit', 'column', Vcmax = 100)
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
            user_lower = list(alpha = 0, Rd = 0, Vcmax = 0),
            user_upper = list(alpha = 1, Rd = 10, Vcmax = 1000),
            user_fit_options = list(alpha = 'fit', J = 'column', Vcmax = 100)
        ),
        '`lower`, `upper`, and `fit_options` must only provide settings for `alpha`, `Rd`, `Vcmax`, `TPU`'
    )
})


test_that('fit options must be `fit`, `column`, or numeric', {
    expect_error(
        PhotoGEA:::assemble_luf(
            test_param_names,
            default_lower = list(),
            default_upper = list(),
            default_fit_options = list(),
            user_lower = list(alpha = 0, Rd = 0, Vcmax = 0),
            user_upper = list(alpha = 1, Rd = 10, Vcmax = 1000),
            user_fit_options = list(alpha = 'WRONG VALUE', Rd = 'column', Vcmax = 100)
        ),
        'Each element of `fit_options` must be `fit`, `column`, or a numeric value'
    )
})

test_that('at least one parameter must be fit', {
    expect_error(
            PhotoGEA:::assemble_luf(
            test_param_names,
            default_lower = list(alpha = 0, Rd = 0, Vcmax = 0, TPU = 0),
            default_upper = list(TPU = 100, alpha = 1, Rd = 10, Vcmax = 1000),
            default_fit_options = list(alpha = 0, Rd = 'fit', Vcmax = 'column', TPU = 40),
            user_lower = list(),
            user_upper = list(Vcmax = 400),
            user_fit_options = list(alpha = 0, Rd = 0, Vcmax = 0, TPU = 0)
        ),
        'No entries in `fit_options` are set to `fit`'
    )
})
