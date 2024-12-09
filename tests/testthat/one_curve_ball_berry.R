# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'])

# Organize the data
licor_file <- organize_response_curve_data(
    licor_file,
    'species_plot',
    c(),
    'Qin',
    1.0
)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate additional gas properties
licor_file <- calculate_gas_properties(licor_file)

# Calculate the Ball-Berry index
licor_file <- calculate_ball_berry_index(licor_file)

# Get just one curve
one_curve <- licor_file[licor_file[, 'species_plot'] == 'soybean - 1a', , TRUE]
