test_that('failures to find TDL cycles are identified', {
    tdl_file <- read_gasex_file(
        PhotoGEA_example_file_path('tdl_sampling_1.dat'),
        'TIMESTAMP'
    )

    expect_error(
        identify_tdl_cycles(
            tdl_file,
            valve_column_name = 'valve_number',
            cycle_start_valve = 20,
            expected_cycle_length_minutes = 2.7,
            expected_cycle_num_valves = 11, # this should be 9
            timestamp_colname = 'TIMESTAMP'
        ),
        'No TDL cycles were found; the values of `cycle_start_valve`, `expected_cycle_length_minutes`, `expected_cycle_num_valves`, or other input arguments may be incorrect for this data set'
    )
})
