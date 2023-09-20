\name{calculate_arrhenius}

\alias{calculate_arrhenius}

\title{Calculate temperature-dependent values using Arrhenius equations}

\description{
  Calculate leaf-temperature-dependent values of various parameters using
  Arrhenius equations. This function can accomodate alternative names for the
  leaf temperature column. It also checks the units of the leaf temperature
  column and will produce an error of the units are not correct.
}

\usage{
  calculate_arrhenius(
    exdf_obj,
    arrhenius_parameters,
    tleaf_column_name = 'TleafCnd'
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{arrhenius_parameters}{
    A list of named lists. Each list element should describe the Arrhenius
    scaling factor (\code{c}), activation energy in \code{kJ / mol} (\code{Ea}),
    and units (\code{units}) for a variable that follows an Arrhenius
    temperature dependence. The name of each list element should be the
    corresponding name of the variable.
  }

  \item{tleaf_column_name}{
    The name of the column in \code{exdf_obj} that contains the leaf temperature
    in units of \code{degrees C}.
  }
}

\details{
  The Arrhenius equation is often used to calculate the temperature dependence
  of the rate of a chemical reaction. It is often stated as follows:

  (1) \code{rate = A * exp(-Ea / (R * T))}

  where \code{A} is the "pre-exponential factor" that sets the overall scaling,
  \code{Ea} is the activation energy, \code{R} is the ideal gas constant, and
  \code{T} is the temperature in Kelvin. See, for example, the
  \href{https://en.wikipedia.org/wiki/Arrhenius_equation}{Wikipedia page for the equation}.

  In photosynthesis research, it is common to use an alternative form of the
  equation, where the pre-exponential factor \code{A} is rewritten as an
  exponent \code{A = exp(c)}, where \code{c} is a "scaling factor" whose value
  can be calculated from \code{A} according to \code{c = ln(A)}). In this
  formulation, the equation becomes:

  (2) \code{rate = exp(c) * exp(-Ea / (R * T)) = exp(c - Ea / (R * T))}

  The advantage of this version is that the natural logarithm of the rate is
  equal to \code{c - Ea / (R * T)}. This means that the Arrhenius paramerer
  values can be easily determined from a linear fit of \code{log(rate)} against
  \code{1 / (R * T)}; \code{c} is the y-intercept and \code{-Ea} is the slope.

  In \code{calculate_arrhenius}, the scaling factor (\code{c}), activation
  energy (\code{Ea}), and units (\code{units}) for a variable must be specified
  as elements of a list, which itself is a named element of
  \code{arrhenius_parameters}. For example, if a variable called \code{Kc} has
  \code{c = 38.05}, {Ea = 79.43}, and units of \code{micromol mol^(-1)}, the
  \code{arrhenius_parameters} argument could be specified as follows:
  \code{list(Kc = list(c = 38.05, Ea = 79.43, units = 'micromol mol^(-1)'))}. It
  is rare to directly specify the Arrhenius parameters; instead, it is more
  typical to use one of the pre-set values such as
  \code{\link{c3_arrhenius_sharkey}}.

  Sometimes a publication will specify the value of a variable at 25 degrees C
  instead of the Arrhenius scaling factor \code{c}. In this case, there is a
  "trick" for determining the value of \code{c}. For example, if the Arrhenius
  exponent should be \code{X} at 25 degrees C, then we have the following:
  \code{X = exp(c - Ea / (R * (25 + 273.15)))}, which we can solve algebraically
  for \code{c} as follows: \code{c = ln(X) + Ea / f}, where
  \code{f = R * (25 + 273.15)}. As a special case, for parameters normalized to
  1 at 25 degrees C, we have \code{c = Ea / f}. The value of \code{f} can be
  accessed as \code{PhotoGEA:::f}.

  Another common scenario is that we may wish to convert the units of a variable
  defined by Arrhenius exponents. For example, let's say \code{Y} is determined
  by an Arrhenius exponent, i.e., that \code{Y = exp(c - Ea / (R * T))}, and we
  want to convert \code{Y} to different units via a multiplicative conversion
  factor \code{cf}. Then, in the new units, \code{Y} becomes
  \code{Y_new = cf * Y = cf * exp(c - (R * T))}. Through algebra, it is possible
  to combine \code{cf} with the original value of \code{c} as
  \code{c_new = c + ln(cf)}. Then we can continue calculating \code{Y_new} using
  an Arrhenius factor as \code{Y_new = exp(c_new - Ea / (R * T))}.
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes one new column
  for each element of \code{arrhenius_parameters}, where the
  temperature-dependent values of these new columns are determined using the
  temperature values specified by the \code{tleaf_column_name} column. The
  category of each of these new columns is \code{calculate_arrhenius} to
  indicate that they were created using this function.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  system.file('extdata', 'ball_berry_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
)

licor_file <- calculate_arrhenius(licor_file, c3_arrhenius_sharkey)

licor_file$units$Kc      # View the units of the new `Kc` column
licor_file$categories$Kc # View the category of the new `Kc` column
licor_file[,'Kc']        # View the values of the new `Kc` column
}

\concept{exdf}
