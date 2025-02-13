\name{calculate_jmax}

\alias{calculate_jmax}

\title{TITLE}

\description{
  DESCRIPTION
}

\usage{
  calculate_jmax(
    data_table,
    alpha_j_at_25,
    theta_j_at_25,
    alpha_j_norm_column_name = 'alpha_j_norm',
    j_at_25_column_name = 'J_at_25',
    j_tl_column_name = 'J_tl_avg',
    qin_column_name = 'Qin_avg',
    theta_j_norm_column_name = 'theta_j_norm',
    tleaf_column_name = 'Tleaf_avg',
    \dots
  )
}

\arguments{
  \item{data_table}{
    A table-like R object such as a data frame or an \code{exdf}.
  }

  \item{\dots}{
    Optional arguments; see below.
  }
}

\details{
  DETAILS
}

\value{
  VALUE
}

\examples{
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
    'CO2_r_sp'
)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_bernacchi)

# For these examples, we will use a faster (but sometimes less reliable)
# optimizer so they run faster
optimizer <- optimizer_nmkb(1e-7)

# Fit all curves in the data set (it is more common to do this)
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c3_aci,
  Ca_atmospheric = 420,
  optim_fun = optimizer
))

# Calculate temperature-dependent values of Jmax-related parameters
aci_results$parameters <- calculate_temperature_response(
    aci_results$parameters,
    jmax_temperature_param_bernacchi,
    'Tleaf_avg'
)

# Calculate Jmax
aci_results$parameters <- calculate_jmax(
    aci_results$parameters,
    0.6895,
    0.97875
)

# Print a few columns
col_to_view <- c('species_plot', 'J_at_25', 'J_tl_avg', 'Jmax_at_25', 'Jmax_tl')

print(aci_results$parameters[, col_to_view, TRUE])
}

\concept{exdf}
