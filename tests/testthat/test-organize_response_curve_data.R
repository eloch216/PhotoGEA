licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

test_that('Curves are only checked when points are removed', {
    expect_error(
        organize_response_curve_data(
            licor_file,
            c('species', 'plot'),
            c(3, 4),
            'Qin',
            0.1,
            print_information = FALSE
        ),
        'The curves do not all follow the same sequence of the driving variable.'
    )

    expect_silent(
        organize_response_curve_data(
            licor_file,
            c('species', 'plot'),
            c(),
            'Qin',
            0.1
        )
    )
})
