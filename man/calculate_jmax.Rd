\name{calculate_jmax}

\alias{calculate_jmax}

\title{Calculate maximum electron transport rate}

\description{
  Calculates maximum electron transport rates (\code{Jmax}) from estimates of
  the electron transport rate (\code{J}) at particular values of incident light
  (\code{Qin}).

  This function is typically used after \code{\link{fit_c3_aci}},
  \code{\link{fit_c3_variable_j}}, or \code{\link{fit_c4_aci}} is used to
  estimate values of \code{J}.
}

\usage{
  calculate_jmax(
    data_table,
    alpha_j_at_25 = 0.293,
    theta_j_at_25 = 0.979,
    alpha_j_norm_column_name = 'alpha_j_norm',
    qin_column_name = 'Qin_avg',
    theta_j_norm_column_name = 'theta_j_norm',
    tleaf_column_name = 'TleafCnd_avg',
    \dots
  )
}

\arguments{
  \item{data_table}{
    A table-like R object such as a data frame or an \code{exdf}.
  }

  \item{alpha_j_at_25}{
    The apparent quantum efficiency of electron transport \code{alpha_j} at 25
    degrees C (dimensionless). If \code{alpha_j_at_25} is not a number, then
    there must be a column in \code{data_table} called \code{alpha_j_at_25} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{alpha_j_at_25} column of \code{data_table} if it exists.
  }

  \item{theta_j_at_25}{
    The empirical curvature parameter \code{theta_j_at_25} at 25 degrees C
    (dimensionless). If \code{theta_j_at_25} is not a number, then there must be
    a column in \code{data_table} called \code{theta_j_at_25} with appropriate
    units. A numeric value supplied here will overwrite the values in the
    \code{theta_j_at_25} column of \code{data_table} if it exists.
  }

  \item{alpha_j_norm_column_name}{
    The name of the column in \code{data_table} that contains the normalized
    \code{alpha_j} values (with units of
    \code{normalized to alpha_j at 25 degrees C}).
  }

  \item{qin_column_name}{
    The name of the column in \code{data_table} that contains values of the
    incident photosynthetically active flux density in
    \code{micromol m^(-2) s^(-1)}.
  }

  \item{theta_j_norm_column_name}{
    The name of the column in \code{data_table} that contains the normalized
    \code{theta_j} values (with units of
    \code{normalized to theta_j at 25 degrees C}).
  }

  \item{tleaf_column_name}{
    The name of the column in \code{data_table} that contains the leaf
    temperature in units of \code{degrees C}.
  }

  \item{\dots}{
    Optional arguments; see below.
  }
}

\details{
  \strong{Basic Requirements:}

  This function requires that \code{data_table} contains columns called
  \code{J_at_25} and \code{J_tl_avg}, as would be included in the output from
  one of the PhotoGEA fitting functions (\code{\link{fit_c3_aci}},
  \code{\link{fit_c3_variable_j}}, and \code{\link{fit_c4_aci}}). These will be
  used to calculate values of \code{Jmax} at 25 degrees C and at leaf
  temperature.

  If any columns for the \code{J} confidence intervals are included in
  \code{data_table} (\code{J_at_25_upper}, \code{J_at_25_lower},
  \code{J_tl_avg_upper}, or \code{J_tl_avg_lower}), the corresponding confidence
  intervals for \code{Jmax} will also be calculated.

  If \code{data_table} is an \code{\link{exdf}} object, units will be checked
  for any columns used in the calculations.

  \strong{Overview of Jmax Calculations:}

  The potential electron transport rate going to support RuBP regeneration
  (\code{J}) depends on the available light energy. \code{J} quickly increases
  with the incident photosynthetically active photon flux density (\code{Qin})
  at low light levels, gradually reaching a plateau at high values of
  \code{Qin}. Although other mathematical representations have been used (Walker
  et al. 2021), this dependence is typically represented as a non-rectangular
  hyperbola:

  \code{J = (I2 + Jmax - sqrt[(I2 + Jmax)^2 - 4 * theta_j * I2 * Jmax]) /
  (2 * theta_j)}, (Eq. 1)

  where \code{Jmax} is the maximum value of \code{J} that would be achieved at
  infinitely large \code{Qin}, \code{0 < theta_j <= 1} is an empirical
  curvature parameter, and \code{I2} is the useful energy absorbed by
  photosystem II. In turn, \code{I2} is calculated by

  \code{I2 = alpha_j * Qin},

  where \code{alpha_j} is the apparent quantum efficiency of electron transport.
  \code{alpha_j} is often defined as

  \code{alpha_j = absorptance * phi_psii,max * beta_psii},

  where \code{absorptance} is the leaf absorptance, \code{phi_psii,max} is the
  maximum quantum yield of photosytem II, and \code{beta_psii} is the fraction
  of light energy partitioned to photosystem II.

  Equation 1 can be understood as a "smooth minimum" of two potential rates of
  electron transport: \code{I2} (which increases linearly with \code{Qin}) and
  \code{Jmax} (which is independent of \code{Qin}). For lower light levels,
  \code{I2} is the smaller rate, and \code{J} is approximately equal to
  \code{I2}; for very high light levels, \code{Jmax} is the smaller rate, and
  \code{J} is approximately equal to \code{Jmax}. For intermediate values of
  \code{Qin}, \code{J} smoothly transitions from \code{I2} to \code{Jmax}.

  One key property of Equation 1 is that the largest possible value of \code{J}
  at a given \code{Qin} is \code{I2}, which only occurs when \code{Jmax} is much
  larger than \code{I2}. In other words, the range of the function in Equation 1
  is \code{0 <= J <= I2}.

  Equation 1 can be solved for \code{Jmax}, enabling calculations of \code{Jmax}
  from estimates of \code{J}:

  \code{Jmax = J * (I2 - theta_j * J) / (I2 - J)} (Eq. 2)

  Because the range of the function in Equation 1 is \code{0 <= J <= I2}, the
  domain of its inverse function (defined in Equation 2) is also
  \code{0 <= J <= I2}. In other words, \code{Jmax} can only be calculated using
  Equation 2 when \code{J < I2}. Otherwise, there is no value of \code{Jmax}
  that can reproduce the value of \code{J} for the given value of
  \code{alpha_j}. This restriction can also be derived more rigorously; see the
  \strong{Detailed algebra} section below for more information.

  If \code{J >= I2}, the \code{calculate_jmax} function will return \code{NA}
  for the value of \code{Jmax}. This behavior can be bypassed by setting the
  optional input argument \code{ignore_restriction} to \code{TRUE}, but this is
  not recommended outside of pedagogical purposes. See Example 2 below for a
  demonstration of what goes wrong when Equation 2 is used for \code{J >= I2}.

  Note that this issue is more significant at lower light levels. For example,
  assuming a typical value of \code{alpha_j} (0.293), \code{I2} for
  \code{Qin} = 1800 micromol / m^2 / s would be 527.4 micromol / m^2 / s. Values
  of \code{J} are typically smaller than this, so an estimate of \code{Jmax} can
  almost always be made. But if a curve were measured at \code{Qin = 300},
  \code{I2} would only be 87.9 micromol / m^2 / s, placing a stronger
  restriction on the values of \code{J} where \code{Jmax} can be estimated. Say
  the best-fit value of \code{J} was 88.9 micromol / m^2 / s for a curve
  measured with \code{Qin} = 300 micromol / m^2 / s; in this case, it would not
  be possible to estimate \code{Jmax}, potentially indicating that the assumed
  value of \code{alpha_j} was not correct.

  \strong{Typical values:}

  According to von Caemmerer (2000), typical values of \code{absorptance},
  \code{phi_psii,max}, and \code{beta_psii} are 0.85, 1 - 0.15, and 0.5,
  respectively, leading to \code{alpha_j} = 0.36125, and the curvature parameter
  \code{theta_j} is typically 0.7.

  Bernacchi et al. (2003) reports that \code{phi_psii,max} is 0.6895 for
  light-adapted leaves at 25 degrees C, while \code{theta_j} at 25 degrees C is
  0.97875. Using this value of of \code{phi_psii,max} with typical values of
  absorptance and \code{beta_psii} results in an \code{alpha_j} estimate of
  0.2930375. These values, rounded to three decimal places, are used as the
  defaults when calculating \code{Jmax}.

  \strong{Absorbed light basis:}

  Values of \code{Jmax} can also be estimated from the absorbed
  photosynthetically active photon flux density (\code{Qabs}). In that case, we
  can regroup the terms in the definition of \code{I2} as follows:

  \code{I2 = (Qin * absorptance) * (phi_psii,max * beta_psii) = Qabs * alpha_j_abs},

  where \code{alpha_j_abs} is given by \code{phi_psii,max * beta_psii}. When
  working in this basis, the default value of \code{alpha_j} at 25 degrees C
  should be divided by the assumed absorptance (0.85), yielding an
  \code{alpha_j_abs} value of 0.345. Along with this change, it would also be
  necessary to change the name of the light column, likely to \code{Qabs_avg}.

  \strong{Why PhotoGEA Uses a Separate Function for Jmax:}

  In principle, values of \code{Jmax} could be estimated by the fitting
  functions that estimate \code{J}: \code{\link{fit_c3_aci}},
  \code{\link{fit_c3_variable_j}}, and \code{\link{fit_c4_aci}}. Instead,
  PhotoGEA requires users to use a separate function (\code{calculate_jmax}) to
  estimate \code{Jmax}. This serves several purposes:

  \itemize{
    \item It highlights that estimates of \code{Jmax} are made using the same
          equations for C3 and C4 leaves.
    \item It leaves open the possibility of other estimates of \code{Jmax}, such
          as those based on a rectangular hyperbola instead of the
          non-rectangular hyperbola used here.
    \item It emphasizes that sometimes it is not possible to provide an estimate
          for \code{Jmax}, depending on the values of \code{Qin},
          \code{alpha_j}, and \code{J}, because of the requirement that
          \code{J < I2 = alpha_j * Qin}.
  }

  The last point is especially important. If \code{Jmax} were varied during the
  fitting process, and \code{J} was estimated from \code{Jmax} using Equation 1,
  there would be a restriction on the possible values of \code{J} that could be
  obtained: \code{J < alpha_j * Qin}. This could potentially bias the fitting
  results, since it may be the case that the best fit would be found for
  \code{J} outside this range.

  In other words, keeping estimates of \code{Jmax} separate from the fitting
  process ensures that the values of \code{alpha_j} and \code{theta_j} have no
  influence on the fits or best-fit values of \code{J}. This is important since
  the true values of these parameters for a particular leaf are difficult or
  impossible to determine.

  \strong{Detailed algebra:}

  Here we will solve Equation 1 for \code{Jmax}, arriving at Equation 2. This
  algebra is reproduced here to highlight the important restriction that
  \code{J < I2}.

  First, multiply both sides of Equation 1 by \code{2 * theta_j}:

  \code{2 * theta_j * J = I2 + Jmax - sqrt[(I2 + Jmax)^2 - 4 * theta_j * I2 * Jmax]}. (Eq. 3)

  Next, isolate the square root term on one side:

  \code{I2 + Jmax - 2 * theta_j * J = sqrt[(I2 + Jmax)^2 - 4 * theta_j * I2 * Jmax]}. (Eq. 4)

  A key point here is that the right hand side cannot be negative, since the
  square root of a real number is never negative. Thus, the left hand side also
  cannot be negative. In other words,

  \code{I2 + Jmax - 2 * theta_j * J >= 0}. (Eq. 5)

  We will return to this restriction later. For now, we square both sides of
  Equation 4:

  \code{(I2 + Jmax)^2 - 4 * theta_j * J * (I2 + Jmax) + 4 * theta_j^2 * J^2 =
  (I2 + Jmax)^2 - 4 * theta_j * I2 * Jmax}. (Eq. 6)

  The term \code{(I2 + Jmax)^2} appears on both sides of Equation 6 and can
  therefore be cancelled out. Grouping the remaining terms that contain
  \code{Jmax} on one side, we have:

  \code{4 * theta_j * Jmax * (I2 - J) = 4 * theta_j * J * (I2 - theta_j * J)} (Eq. 7)

  Finally, provided that \code{I2 - J} is not zero (in other words, that
  \code{I2} is not equal to \code{J}), we can divide both sides of Equation 7 by
  \code{4 * theta_j * (I2 - J)} to obtain Equation 2 above.

  Now, we can use this expression (Equation 2) to replace \code{Jmax} in
  Equation 5:

  \code{I2 + J * (I2 - theta_j * J) / (I2 - J) - 2 * theta_j * J >= 0}. (Eq. 8)

  This can be converted to a single ratio as follows:

  \code{[(I2 - 2 * theta_j * J) * (I2 - J) + J * (I2 - theta_j * J)] / (I2 - J) >= 0}. (Eq. 9)

  Multiplying out the factors in the numerator and collecting like terms,
  Equation 9 becomes

  \code{[I2^2 - 2 * theta_j * I2 * J + theta_j * J^2] / (I2 - J) >= 0}. (Eq. 10)

  Because \code{theta_j} must lie between 0 and 1, \code{theta_j^2} is always
  less than or equal to \code{theta_j}. This allows us to place a lower bound on
  the value of the numerator of the left hand side of Equation 10:

  \code{I2^2 - 2 * theta_j * I2 * J + theta_j * J^2 >=
  I2^2 - 2 * theta_j * I2 * J + theta_j^2 * J^2}. (Eq. 11)

  The right hand side of Equation 11 can be refactored:

  \code{I2^2 - 2 * theta_j * I2 * J + theta_j * J^2 >=
  (I2 - theta_j J)^2}. (Eq. 12)

  The right hand side of Equation 12 can never be negative, so from this we can
  also conclude that the numerator of the left hand side of Equation 10 can also
  never be negative. Thus, the inequality in Equation 10 is satisfied whenever
  its denominator is positive. In other words, whenever \code{I2 - J > 0}, or,
  equivalently, \code{J < I2}.

  Thus, we have shown that Equation 2 holds whenever \code{J < I2}, since, when
  this inequality is satisfied, Equation 5 is also satisfied.

  Although we do not do so here, it can be shown that when \code{I2 < J}, the
  value of \code{Jmax} that would be calculated by Equation 2 is the inverse of

  \code{J = (I2 + Jmax + sqrt[(I2 + Jmax)^2 - 4 * theta_j * I2 * Jmax]) /
  (2 * theta_j)} (Eq. 13)

  rather than the inverse of Equation 1. Note the difference: in Equation 13,
  the square root term is added to \code{I2 + Jmax} rather than subtracted. This
  is a "smooth maximum" function, rather than a smooth minimum. In fact,
  whenever \code{I2 > Jmax}, Equation 13 would predict \code{J > Jmax}, clearly
  a nonsensical result. Likewise, the inverse of the function in Equation 13
  would predict some values of \code{Jmax} that are smaller than \code{J}.
  Example 2 below shows that it can even return negative values of \code{Jmax},
  which is clearly not reasonable from a biological perspective.

  \strong{References:}

  \itemize{
    \item von Caemmerer, S. "Biochemical Models of Leaf Photosynthesis" (CSIRO
          Publishing, 2000) [\doi{10.1071/9780643103405}].

    \item Walker, A. P. et al. "Multi-hypothesis comparison of Farquhar and
          Collatz photosynthesis models reveals the unexpected influence of
          empirical assumptions at leaf and global scales." Global Change
          Biology 27, 804–822 (2021) [\doi{10.1111/gcb.15366}].

    \item Bernacchi, C. J., Pimentel, C. & Long, S. P. "In vivo temperature
          response functions of parameters required to model RuBP-limited
          photosynthesis" Plant, Cell & Environment 26, 1419–1430 (2003)
          [\doi{10.1046/j.0016-8025.2003.01050.x}].
  }
}

\value{
  The return value is a table based on \code{data_table} that includes several
  new columns: \code{I2_at_25}, \code{Jmax_at_25}, \code{Jmax_at_25_msg},
  \code{I2_tl}, \code{Jmax_tl}, and \code{Jmax_tl_msg}. The \code{_msg} columns
  indicate when the error condition \code{J >= I2} has occurred.

  If \code{J} confidence intervals were provided in the inputs, then there will
  be correspoding columns for the related \code{Jmax}, and \code{msg} values;
  for example, \code{Jmax_at_25_lower} and \code{Jmax_at_25_lower_msg}.
}

\examples{
## Example 1: Estimating Jmax after fitting several C3 A-Ci curves

# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c3_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'] )

# Organize the data; we will need average values of leaf temperature and
# incident PPFD in order to calculate Jmax later
licor_file <- organize_response_curve_data(
    licor_file,
    'species_plot',
    c(9, 10, 16),
    'CO2_r_sp',
    columns_to_average = c('TleafCnd', 'Qin')
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
    'TleafCnd_avg'
)

# Calculate Jmax
aci_results$parameters <- calculate_jmax(aci_results$parameters)

# Print a few columns
col_to_view <- c('species_plot', 'J_at_25', 'J_tl_avg', 'Jmax_at_25', 'Jmax_tl')

print(aci_results$parameters[, col_to_view, TRUE])

## Example 2: Illustrating the importance of requiring I2 > J

# Define a data frame with input values
npts <- 200
J_seq <- seq_len(npts)

jmax_df <- data.frame(
  J_at_25 = J_seq,
  J_tl_avg = J_seq,
  alpha_j_norm = 1,
  Qin_avg = 300,
  theta_j_norm = 1,
  TleafCnd_avg = 25
)

# Calculate Jmax values, overriding the default behavior so that values of Jmax
# are returned even when I2 < J.
jmax_df <- calculate_jmax(jmax_df, ignore_restriction = TRUE)

# Plot the Jmax values, distinguishing between cases where J < I2 and where
# J > I2. Here we can see that when J > I2, values of Jmax are smaller than J,
# and can even be negative, which is clearly unreasonable from a biological
# perspective. To highlight these considerations, J = I2 is plotted as a dashed
# black line, Jmax = J is plotted as a black long-dashed line, and Jmax = 0 is
# plotted as a solid black line.
ymin <- -50
ymax <- 250

xmin <- min(J_seq)
xmax <- max(J_seq)

I2 <- jmax_df$I2_at_25[1]

jmax_df$Jmax_at_25_msg[jmax_df$Jmax_at_25_msg == ''] <- 'J < I2'

lattice::xyplot(
  Jmax_at_25 ~ J_at_25,
  group = Jmax_at_25_msg,
  data = jmax_df,
  auto = TRUE,
  type = 'l',
  xlim = c(xmin, xmax),
  ylim = c(ymin, ymax),
  xlab = 'J (micromol / m^2 / s)',
  ylab = 'Jmax (micromol / m^2 / s)',
  panel = function(x, y, ...) {
    lattice::panel.lines(c(0, 0) ~ c(xmin, xmax), lty = 1, col = 'black')
    lattice::panel.lines(c(ymin, ymax) ~ c(I2, I2), lty = 2, col = 'black')
    lattice::panel.lines(J_seq ~ J_seq, lty = 5, col = 'black')
    lattice::panel.xyplot(x, y, ...)
  }
)
}

\concept{exdf}
