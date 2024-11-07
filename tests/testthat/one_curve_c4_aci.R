# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c4_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'] )

# Organize the data
licor_file <- organize_response_curve_data(
    licor_file,
    'species_plot',
    c(9, 10, 16),
    'CO2_r_sp'
)

# Calculate temperature-dependent values of C4 photosynthetic parameters. Here
# we modify the scaling of Jmax_norm so it equals 1 at the optimal temperature.
licor_file <- calculate_temperature_response(
    licor_file,
     within(c4_temperature_param_vc, {Jmax_norm$optimum_rate = 1})
)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Get just one curve
one_curve <- licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE]

# Purposely introduce negative Ci values
one_curve_bad <- one_curve
one_curve_bad[one_curve_bad[, 'CO2_r_sp'] == 20, 'Ci'] <- -5
