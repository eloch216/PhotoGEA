# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Get a subset of columns
id_col <- c('species', 'K')
other_col <- c('A', 'gsw')
licor_file <- licor_file[ , c(id_col, other_col), TRUE]

# Set one row to NA
for (cn in other_col) {
    licor_file[3, cn] <- NA
}

# Choose test tolerance
TOLERANCE <- 1e-4

test_that('basic stats have not changed (removing NA)', {
    res <- basic_stats(
      licor_file[ , c('species', 'K', 'A', 'gsw'), TRUE],
      'species'
    )

    expect_false('species_avg' %in% colnames(res))
    expect_false('species_stderr' %in% colnames(res))

    expect_equal(
        res[, 'A_avg'],
        c(23.74270, 17.95311),
        tolerance = TOLERANCE
    )

    expect_equal(
        res[, 'A_stderr'],
        c(2.595648, 3.026751),
        tolerance = TOLERANCE
    )
})

test_that('basic stats have not changed (not removing NA)', {
    res <- basic_stats(
      licor_file[ , c('species', 'K', 'A', 'gsw'), TRUE],
      'species',
      na.rm = FALSE
    )

    expect_false('species_avg' %in% colnames(res))
    expect_false('species_stderr' %in% colnames(res))

    expect_equal(
        res[, 'A_avg'],
        c(NA, 17.95311),
        tolerance = TOLERANCE
    )

    expect_equal(
        res[, 'A_stderr'],
        c(NA, 3.026751),
        tolerance = TOLERANCE
    )
})
