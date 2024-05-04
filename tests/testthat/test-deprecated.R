test_that('Deprecated functions produce errors when they are called', {
    expect_error(read_tdl_file())
    expect_error(read_licor_file())
    expect_error(check_licor_data())
})
