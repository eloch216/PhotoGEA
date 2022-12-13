\name{initial_guess_c4_aci}

\alias{initial_guess_c4_aci}

\title{Make an initial guess of C4 photosynthesis parameter values for one curve}

\description{
  Creates a function that makes an initial guess of C4 photosynthesis model
  parameter values for one curve. The returned function meets the requirements
  for the \code{initial_guess_fun} input argument of \code{\link{fit_c4_aci}}.

  Values estimated by this guessing function should be considered inaccurate,
  and should always be improved upon by an optimizer.
}

\usage{
  initial_guess_c4_aci(
    pcm_threshold_rm = 100,
    gbs = 0.003,
    a_column_name = 'A',
    pcm_column_name = 'PCm',
    kp_column_name = 'Kp',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm'
  )
}

\arguments{
  \item{pcm_threshold_rm}{
    An upper cutoff value for the partial pressure of CO2 in the mesophyll (in
    \code{microbar}) to be used when estimating \code{Rm}.
  }

  \item{gbs}{
    The bundle sheath conductance to CO2 in \code{mol m^(-2) s^(-1) bar^(-1)}.
  }

  \item{a_column_name}{
    The name of the column in \code{licor_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{pcm_column_name}{
    The name of the column in \code{exdf_obj} that contains the partial pressure
    of CO2 in the mesophyll, expressed in \code{microbar}.
  }

  \item{kp_column_name}{
    The name of the column in \code{exdf_obj} that contains the Michaelis-Menten
    constant for PEP carboxylase carboxylation in \code{microbar}.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{vpmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vpmax} values (with units of
    \code{normalized to Vpmax at 25 degrees C}).
  }
}

\details{
  Here we estimate values of \code{Rm}, \code{Vcmax}, \code{Vpmax}, and
  \code{Vpr} from a measured C4 CO2 response curve. For more information about
  these parameters, see the documentation for
  \code{\link{calculate_c4_assimilation}}. To estimate these parameter values,
  we use several equations from S. von Caemmerer, "Biochemical Models of Leaf
  Photosynthesis" (CSIRO Publishing, 2000) [\doi{10.1071/9780643103405}]. Any
  equation numbers referenced below are from this book.

  \itemize{
    \item \strong{Estimating Rm}: An estimate for \code{Rm} can be obtained
          using Equation 4.26, which applies for low values of \code{PCm}. In
          this situation, \code{PCm + Kp} can be approximated by \code{Kp}, and
          Equation 4.26 simplifies to a linear relationship between the net
          assimilation \code{An} and \code{PCm}: \code{An =
          (gbs + Vpmax / kP) * PCm - Rm}. So, to estimate \code{Rm}, we make a
          linear fit of \code{An} vs. \code{PCm} in the low \code{PCm} range
          where this equation is expected to be valid. Then \code{Rm} is given
          by the negative of the intercept from the fit. In the C4 assimilation
          model, we assume that \code{Rm = 0.5 * Rd}, so we can also estimate
          \code{Rd = 2 * Rm} from this value.
    \item \strong{Estimating Vpmax}: An estimate for \code{Vpmax} can also be
          obtained from Equation 4.26. In this case, we simply solve the
          equation for \code{Vpmax} and use it to calculate a value of
          \code{Vpmax} at each point in the curve from the measured values of
          \code{An} and \code{PCm}, the input value of \code{gbs}, and the value
          of \code{Rm} estimated above. In the PEP-carboxylation-limited range,
          the estimated values of \code{Vpmax} should be reasonable. In other
          parts of the curve, the assimilation rate is limited by other factors,
          so \code{An} will be smaller than the PEP-carboxylation-limited
          values, causing the estimated values of \code{Vpmax} to be smaller.
          So, to make an overall estimate, we choose the largest estimated
          \code{Vpmax} value.
    \item \strong{Estimating Vcmax}: An estimate for \code{Vcmax} can be
          obtained by solving \code{An = Vcmax - Rd} for \code{Vcmax}, similar
          to the method used to estimate \code{Vpmax}.
    \item \strong{Estimating Vpr}: An estimate for \code{Vpr} can be
          obtained by solving \code{An = Vpr + gbs * PCm - Rm} for \code{Vpr},
          similar to the method used to estimate \code{Vpmax}.
  }

  Note that a key assumption underlying this approach is that the net
  assimilation can be reasonably approximated by \code{An = min(Apc, Apr, Ar)}
  (Equation 4.25). While this approximation seems to work well for low values of
  \code{PCm}, it tends to deviate significantly from the more accurate version
  (Equation 4.21) at higher values of \code{PCm}, predicting values that are
  noticably smaller. Thus, the values of \code{Vcmax} and \code{Vpr} estimated
  using this procedure are unlikely to be accurate. (See the example below.)
  This is not a problem, but instead it simply highlights the importance of
  improving this initial guess using an optimizer, which can be accomplished via
  \code{\link{fit_c4_aci}}.
}

\value{
  A function with one input argument \code{rc_exdf}, which should be an
  \code{exdf} object representing one C4 CO2 response curve. The return value of
  this function will be a numeric vector with three elements, representing the
  values of \code{Vcmax}, \code{Vpmax}, and \code{Vpr} (in that order) at 25
  degrees C.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  system.file('extdata', 'c4_aci_1.xlsx', package = 'PhotoGEA')
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

# Specify a value for mesophyll conductance
licor_file <- set_variable(
  licor_file,
  'gmc', 'mol m^(-2) s^(-1) bar^(-1)', value = 3.0
)

# Calculate PCm
licor_file <- apply_gm(licor_file, 'C4')

# Calculate temperature-dependent values of C4 photosynthetic parameters
licor_file <- calculate_arrhenius(licor_file, c4_arrhenius_von_caemmerer)

# Create the guessing function
guessing_func <- initial_guess_c4_aci()

# Apply it and see the initial guesses for each curve
str(by(licor_file, licor_file[, 'species_plot'], guessing_func))

# Calculate simulated A-Ci curves based on the guesses and compare them to the
# actual data
calculated_aci <- do.call(rbind, by(
  licor_file,
  licor_file[, 'species_plot'],
  function(x) {
    param <- guessing_func(x)
    calculate_c4_assimilation(x, param[1], param[2], param[3])
  }
))

lattice::xyplot(
  Apr + Apc + Ar + An + A ~ PCm | species_plot,
  data = cbind(licor_file, calculated_aci)$main_data,
  type = 'b',
  auto = TRUE,
  grid = TRUE
)
}

\concept{exdf}
