test_that('Useful messages occur for bad example file names', {
    expect_error(
        PhotoGEA_example_file_path('fake_file.xlsx'),
        'The PhotoGEA package does not include an example file called `fake_file.xlsx`. If you are loading your own data files, do not use the `PhotoGEA_example_file_path` function. Type `?PhotoGEA_example_file_path` for more information.',
        fixed = TRUE
    )
})

test_that('Example files are found', {
    fn <- expect_no_error(
        PhotoGEA_example_file_path('ball_berry_1.xlsx')
    )

    expect_true(
        file.exists(fn)
    )
})
