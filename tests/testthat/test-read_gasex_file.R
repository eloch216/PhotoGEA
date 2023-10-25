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
