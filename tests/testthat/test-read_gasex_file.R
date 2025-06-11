test_that('file_name cannot be an empty string', {
    expect_error(
        read_gasex_file(''),
        'The `file_name` input argument is an empty string. If `file_name` was generated using `system.file`, this means that the desired file could not be found.'
    )
})

test_that('file_name must exist', {
    expect_error(
        read_gasex_file('fake_file.xlsx'),
        '`fake_file.xlsx` does not exist'
    )
})


test_that('standardizations are applied', {
    fpath <- PhotoGEA_example_file_path('c3_aci_1.xlsx')

    licor_file_standard <- read_gasex_file(fpath)

    licor_file <- read_gasex_file(fpath, standardize_columns = FALSE)

    expect_equal(licor_file_standard$units$PhiPS2, 'dimensionless')
    expect_equal(licor_file$units$PhiPS2, 'NA')
})

test_that('NA rows are removed (plaintext)', {
    fpath <- PhotoGEA_example_file_path('plaintext_licor_file_v2')

    licor_file_with_NA_row <- expect_silent(
        read_gasex_file(fpath, remove_NA_rows = FALSE, include_user_remark_column = FALSE)
    )

    licor_file_without_NA_row <- expect_silent(
        read_gasex_file(fpath, remove_NA_rows = TRUE, include_user_remark_column = FALSE)
    )

    expect_equal(
        nrow(licor_file_without_NA_row),
        nrow(licor_file_with_NA_row) - 1
    )
})
