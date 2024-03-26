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
    fpath <- system.file('extdata', 'c3_aci_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)

    licor_file_standard <- read_gasex_file(fpath)

    licor_file <- read_gasex_file(fpath, standardize_columns = FALSE)

    expect_equal(licor_file_standard$units$PhiPS2, 'dimensionless')
    expect_equal(licor_file$units$PhiPS2, 'NA')
})
