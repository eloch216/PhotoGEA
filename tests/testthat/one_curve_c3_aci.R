# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c3_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'] )

# Organize the data
licor_file <- organize_response_curve_data(
    licor_file,
    'species_plot',
    c(9, 10, 16),
    'CO2_r_sp',
    columns_to_average = c('Qin', 'TleafCnd')
)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters. Here
# we use the "Bernacchi" option, but override the Tp and gmc responses with a
# flat one.
licor_file <- calculate_temperature_response(
    licor_file,
    within(c3_temperature_param_bernacchi, {
        Tp_norm = c3_temperature_param_flat$Tp_norm
        gmc_norm = c3_temperature_param_flat$gmc_norm
    })
)

# Get just one curve
one_curve <- licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE]

# Purposely introduce negative Ci values
one_curve_bad <- one_curve
one_curve_bad[one_curve_bad[, 'CO2_r_sp'] == 20, 'Ci'] <- -5
