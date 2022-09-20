\name{initial_guess_c3_aci}

\alias{initial_guess_c3_aci}

\title{Make an initial guess of FvCB model parameter values for one curve}

\description{
  Creates a function that makes an initial guess of FvCB model parameter values
  for one curve. Values estimated in this way should be considered inaccurate,
  and should always be improved upon by an optimizer.
}

\usage{
  initial_guess_c3_aci(
    cc_threshold_rd = 100,
    Oc = 210000,
    a_column_name = 'A',
    cc_column_name = 'Cc',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    gamma_star_column_name = 'Gamma_star',
    vcmax_norm_column_name = 'Vcmax_norm',
    rd_norm_column_name = 'Rd_norm',
    j_norm_column_name = 'J_norm'
  )
}

\arguments{
  \item{cc_threshold_rd}{
    An upper cutoff value for the chloroplast CO2 concentration in
    \code{micromol mol^(-1)} to be used when estimating \code{Rd}.
  }

  \item{Oc}{
    The concentration of O2 in the chloroplast, expressed in
    \code{micromol mol^(-1)}. This is often assumed to be the ambient value. For
    air measurements, this would be 21 percent O2, which is 210000 ppm. For low
    oxygen measurements, this would be 2 percent O2, which is about 20000 ppm.
  }

  \item{a_column_name}{
    The name of the column in \code{licor_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{cc_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    chloroplastic CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{kc_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    Michaelis-Menten constant for rubisco carboxylation in
    \code{micromol mol^(-1)}.
  }

  \item{ko_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    Michaelis-Menten constant for rubisco oxygenation in
    \code{mmol mol^(-1)}.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    \code{Gamma_star} values in \code{micromol mol^(-1)}.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{rd_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Rd} values (with units of \code{normalized to Rd at 25 degrees C}).
  }

  \item{j_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{J} values (with units of \code{normalized to J at 25 degrees C}).
  }
}

\details{
  Here we estimate values of \code{Rd}, \code{Vcmax}, \code{J}, and \code{TPU}
  from a measured C3 CO2 response curve. For more information about these
  parameters, see the documentation for \code{\link{calculate_c3_assimilation}}.
  To estimate these parameter values, we use several equations from S. von
  Caemmerer, "Biochemical Models of Leaf Photosynthesis" (CSIRO Publishing,
  2000) [\doi{10.1071/9780643103405}]. Any equation numbers referenced below are
  from this book.

  \itemize{
    \item \strong{Estimating Rd}: An estimate for \code{Rd} can be obtained using
          Equation 2.20. To start, note that the term \code{Kc * (1 + Oc / Ko)}
          from this equation typically has a value of around \code{1000 ppm}.
          So, when the concentration of CO2 in the chloroplast is low,
          \code{Cc + Kc * (1 + Oc / Ko)} can be approximated as
          \code{Kc * (1 + Oc / Ko)}. Under these conditions, Equation 2.20
          simplifies to a linear relationship between \code{Ac} and \code{Cc}:
          \code{Ac = (Cc - Gamma_star) * Vcmax / (Kc * (1.0 + Oc / Ko)) - Rd}.
          Additionally, \code{Ac = -Rd} when \code{Cc = Gamma_star}. So, to
          estimate \code{Rd}, we make a linear fit of \code{A} vs. \code{Cc} in
          the low \code{Cc} range where assimilation is expected to be rubisco
          limited and this linear equation is expected to be valid. Then, we
          evaluate the linear model at \code{Cc = Gamma_star}, estimating
          \code{Rd} as the negative of the predicted assimilation value at that
          point.
    \item \strong{Estimating Vcmax}: An estimate for \code{Vcmax} can also be
          obtained from Equation 2.20. In this case, we simply solve the
          equation for \code{Vcmax} and then use it to calculate a value of
          \code{Vcmax} at each point in the curve from the measured values of
          \code{A} and \code{Cc}, and the value of \code{Rd} estimated as above.
          In the rubisco-limited part of the curve, \code{A = Ac} and the
          estimated values of \code{Vcmax} should be reasonable. In other parts
          of the curve, \code{Ac} is not the limiting rate, so \code{A < Ac}.
          Consequently, the estimated values of \code{Vcmax} in these parts of
          the curve will be smaller. So, to make an overall estimate, we choose
          the the largest estimated \code{Vcmax} value.
    \item \strong{Estimating J}: An estimate for \code{J} can be obtained by
          solving Equation 2.23 for \code{J}, similar to the method used to
          estimate \code{Vcmax}.
    \item \strong{Estimating TPU}: An estimate for \code{TPU} can be obtained by
          solving Equation 2.26 for \code{TPU}, similar to the method used to
          estimate \code{Vcmax}.
  }
}

\value{
  A function with one input argument \code{rc_exdf}, which should be an
  \code{exdf} object representing one C3 CO2 response curve. The return value of
  this function will be a numeric vector with four elements, representing the
  values of \code{TPU}, \code{J}, \code{Rd}, and \code{Vcmax} (in that order) at
  25 degrees C.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_licor_file(
  system.file('extdata', 'c3_aci_1.xlsx', package = 'PhotoGEA'),
  c(3, 5, 7, 9, 11, 13), 14, 15, 16, 17, 'time'
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

# Specify an infinite mesophyll conductance (so `Cc` = `Ci`)
licor_file <- set_variable(
  licor_file,
  'gmc', 'mol m^(-2) s^(-1) bar^(-1)', value = Inf
)

# Calculate Cc
licor_file <- calculate_cc(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_arrhenius(licor_file, c3_arrhenius_bernacchi)

# Create the guessing function
guessing_func <- initial_guess_c3_aci()

# Apply it and see the initial guesses for each curve
str(by(licor_file, licor_file[, 'species_plot'], guessing_func))

# Calculate simulated A-Ci curves based on the guesses and compare them to the
# actual data
calculated_aci <- do.call(rbind, by(
  licor_file,
  licor_file[, 'species_plot'],
  function(x) {
    param <- guessing_func(x)
    calculate_c3_assimilation(x, param[1], param[2], param[3], param[4])
  }
))

lattice::xyplot(
  Ac + Aj + Ap + An + A ~ Cc | species_plot,
  data = cbind(licor_file, calculated_aci)$main_data,
  type = 'b',
  auto = TRUE,
  grid = TRUE
)
}

\concept{exdf}
