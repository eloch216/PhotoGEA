\name{calculate_temperature_response_johnson}

\alias{calculate_temperature_response_johnson}

\title{Calculate temperature-dependent values using Johnson-Eyring-Williams equations}

\description{
  Calculate leaf-temperature-dependent values of various parameters using
  Johnson-Eyring-Williams equations. It is rare for users to call this function
  directly; instead, it is used internally by
  \code{\link{calculate_temperature_response}}.
}

\usage{
  calculate_temperature_response_johnson(
    exdf_obj,
    johnson_parameters,
    tleaf_column_name = 'TleafCnd'
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{johnson_parameters}{
    A list of named lists. Each list element should describe the scaling factor
    (\code{c}), enthalpy of activation in \code{kJ / mol} (\code{Ha}),
    enthalpy of deactivation in \code{kJ / mol} (\code{Hd}), entropy in
    \code{kJ / K / mol} (\code{S}), and units (\code{units}) for a variable that
    follows a Johnson-Eyring-Williams temperature dependence. The name of each
    list element should be the corresponding name of the variable.
  }

  \item{tleaf_column_name}{
    The name of the column in \code{exdf_obj} that contains the leaf temperature
    in units of \code{degrees C}.
  }
}

\details{
  The Johnson-Eyring-Williams equation is often used to calculate the
  temperature dependence of the rate of a chemical reaction. It can be stated
  as follows:

  \code{rate = exp(c - Ha / (R * T)) / (1 + exp(S / R - Hd / (R * T)))}

  where \code{c} is the scaling factor that sets the overall magnitude of the
  rate, \code{Ha} is the enthalpy of activation, \code{Hd} is the enthalpy of
  deactivation, \code{S} is the entropy, \code{R} is the ideal gas constant, and
  \code{T} is the temperature in Kelvin.

  This equation exhibits a peak; in other words, there is a particular
  temperature (the optimal temperature) where the rate is maximized. Thus, it is
  often used in place of an Arrhenius equation (see
  \code{\link{calculate_temperature_response_arrhenius}}) for photosynthetic
  parameters that exhibit a decrease at high temperatures.

  This equation was originally published in Johnson, Eyring, & Williams (1942)
  and has been used to model the temperature dependence of key photosynthetic
  parameters, as in Harley et al. (1992), Bernacchi et al. (2003), Sharkey et
  al. (2007), and others.

  In \code{calculate_temperature_response_johnson}, the scaling factor
  (\code{c}), enthalpy of activation (\code{Ha}), enthalpy of deactivation
  (\code{Hd}), entopy (\code{S}), and units (\code{units}) for a variable must
  be specified as elements of a list, which itself is a named element of
  \code{johnson_parameters}. For example, if a variable called
  \code{Tp} has \code{c = 21.46}, \code{Ha = 53.1}, \code{Hd = 201.8},
  \code{S = 0.65}, and units of \code{micromol mol^(-1)}, the
  \code{johnson_parameters} argument could be specified as follows:
  \code{list(Tp = list(c = 21.46, Ha = 53.1, Hd = 201.8, S = 0.65, units =
  'micromol mol^(-1)'))}.

  It is rare to directly specify these parameters; instead, it is more
  typical to use one of the pre-set values such as those included in
  \code{\link{c3_temperature_param_sharkey}}.

  \strong{References:}

  \itemize{
    \item Johnson, F. H., Eyring, H. & Williams, R. W. "The nature of enzyme
          inhibitions in bacterial luminescence: Sulfanilamide, urethane,
          temperature and pressure." Journal of Cellular and Comparative
          Physiology 20, 247–268 (1942) [\doi{10.1002/jcp.1030200302}].
    \item Harley, P. C., Thomas, R. B., Reynolds, J. F. & Strain, B. R.
          "Modelling photosynthesis of cotton grown in elevated CO2." Plant,
          Cell & Environment 15, 271–282 (1992)
          [\doi{10.1111/j.1365-3040.1992.tb00974.x}].
    \item Bernacchi, C. J., Pimentel, C. & Long, S. P. "In vivo temperature
          response functions of parameters required to model RuBP-limited
          photosynthesis." Plant, Cell & Environment 26, 1419–1430 (2003)
          [\doi{10.1046/j.0016-8025.2003.01050.x}].
    \item Sharkey, T. D., Bernacchi, C. J., Farquhar, G. D. & Singsaas, E. L.
          "Fitting photosynthetic carbon dioxide response curves for C3 leaves."
          Plant, Cell & Environment 30, 1035–1040 (2007)
          [\doi{10.1111/j.1365-3040.2007.01710.x}].
  }
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes one new column
  for each element of \code{johnson_parameters}, where the temperature-dependent
  values of these new columns are determined using the temperature values
  specified by the \code{tleaf_column_name} column. The category of each of
  these new columns is \code{calculate_temperature_response_johnson} to indicate
  that they were created using this function.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_temperature_response_johnson(
  licor_file,
  list(Tp_norm = c3_temperature_param_sharkey$Tp_norm)
)

licor_file$units$Tp_norm      # View the units of the new `Tp_norm` column
licor_file$categories$Tp_norm # View the category of the new `Tp_norm` column
licor_file[,'Tp_norm']        # View the values of the new `Tp_norm` column
}

\concept{exdf}
