library(PhotoGEA)

library(dfoptim)

INCLUDE_FLUORESCENCE <- TRUE

PERFORM_CALCULATIONS <- TRUE

CHOOSE_FILES_INTERACTIVELY <- TRUE

LICOR_FILES_TO_PROCESS <- c()

# Specify the filenames depending on the value of the CHOOSE_FILES_INTERACTIVELY
# and USE_GM_TABLE booleans
if (PERFORM_CALCULATIONS) {
    if (CHOOSE_FILES_INTERACTIVELY) {
        LICOR_FILES_TO_PROCESS <- choose_input_licor_files()
    }
    else {
        LICOR_FILES_TO_PROCESS <- c(
            "2021-04-07-site 11 vulcan cs 36627-1-17.xlsx",
            "20210407-pluto-site13-36627-WT-3.xlsx"
        )
    }
}                                                                   ###

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- "event"
REP_COLUMN_NAME <- "replicate"
MEASUREMENT_NUMBER_NAME <- "obs"
GM_COLUMN_NAME <- "gmc"
CI_COLUMN_NAME <- "Ci"
CC_COLUMN_NAME <- "Cc"
A_COLUMN_NAME <- "A"
GSW_COLUMN_NAME <- "gsw"
IWUE_COLUMN_NAME <- "iwue"
O2_COLUMN_NAME <- "O2"
F_PRIME_COLUMN_NAME <- "f_prime"
GAMMA_STAR_COLUMN_NAME <- "gamma_star"
KC_COLUMN_NAME <- "Kc"
KO_COLUMN_NAME <- "Ko"
TIME_COLUMN_NAME <- "time"

# Specify the variables to extract. Note that when the file is loaded, any
# Unicode characters such as Greek letters will be converted into `ASCII`
# versions, e.g. the character Δ will be become `Delta`. The conversion rules
# are defined in the `UNICODE_REPLACEMENTS` data frame (see `read_licor.R`).
VARIABLES_TO_EXTRACT <- c(
    "obs",
    TIME_COLUMN_NAME,
    "elapsed",
    "date",
    "hhmmss",
    EVENT_COLUMN_NAME,
    REP_COLUMN_NAME,
    A_COLUMN_NAME,
    CI_COLUMN_NAME,
    GSW_COLUMN_NAME,
    "CO2_r_sp",
    "Ca",
    "gbw",
    "Qin",
    "Qabs",
    "CO2_r",
    "TleafCnd"
)

if (INCLUDE_FLUORESCENCE) {
    VARIABLES_TO_EXTRACT <- c(VARIABLES_TO_EXTRACT, "PhiPS2", "ETR")
}


###                                                                   ###
### COMMANDS THAT ACTUALLY CALL THE FUNCTIONS WITH APPROPRIATE INPUTS ###
###                                                                   ###

# Load the data and calculate the stats, if required
if (PERFORM_CALCULATIONS) {
    multi_file_info <- batch_read_licor_file(
        LICOR_FILES_TO_PROCESS,
        preamble_data_rows = c(3, 5, 7, 9, 11, 13),
        variable_category_row = 14,
        variable_name_row = 15,
        variable_unit_row = 16,
        data_start_row = 17,
        timestamp_colname = TIME_COLUMN_NAME
    )

    extracted_multi_file_info <- batch_extract_variables(
        multi_file_info,
        VARIABLES_TO_EXTRACT
    )

    combined_info <- combine_exdf(extracted_multi_file_info)

    combined_info_subset <- combined_info[combined_info[,'event'] == 'WT',]

    two_reps <- combined_info_subset[combined_info_subset$replicate == '17-12' | combined_info_subset$replicate == '10-5',]

    one_rep <- combined_info_subset[combined_info_subset$replicate == '17-12',]

    combined_info_data <- combined_info[['main_data']]
    combined_info_data$event_replicate <- paste(combined_info_data$event, combined_info_data$replicate)

    set.seed(123)
}

tmp2 <- fit_variable_j_std(two_reps, 'replicate', dfoptim::hjkb)

xyplot(A + An_estimate ~ Ci | replicate, data = tmp2$fits, type = 'b', pch = 20, auto = TRUE, grid = TRUE)

xyplot(gm ~ Ci, group = replicate, data = tmp2$fits, type = 'b', pch = 20, auto = TRUE, grid = TRUE)



tmp3 <- fit_variable_j_std(combined_info_subset, 'replicate', dfoptim::hjkb)

tmp3 <- fit_variable_j_std(
    combined_info_subset,
    'replicate',
    function(guess, fun, lower, upper) {
        dfoptim::nmkb(guess, fun, lower, upper, control = list(
            tol = 1e-7,
            maxfeval = 2000,
            restarts.max = 10
        ))
    }
)


tmp3 <- fit_variable_j_std(combined_info_subset, 'replicate', dfoptim::mads, upper = c(500, 500, 10, 0.6, 500))

tmp3 <- fit_variable_j_std(combined_info_subset, 'replicate', dfoptim::nmkb)

View(tmp3$parameters)

xyplot(A + An_estimate ~ Ci | replicate, data = tmp3$fits, type = 'b', pch = 20, auto = TRUE, grid = TRUE)

xyplot(gm ~ Ci, group = replicate, data = tmp3$fits, type = 'b', pch = 20, auto = TRUE, grid = TRUE, ylim=c(-1,1))

xyplot(gm ~ Cc, group = replicate, data = tmp3$fits, type = 'b', pch = 20, auto = TRUE, grid = TRUE, ylim=c(-1,1))







tmp4 <- fit_variable_j_std(
    combined_info_data,
    'event_replicate',
    function(guess, fun, lower, upper) {
        dfoptim::nmkb(guess, fun, lower, upper, control = list(
            tol = 1e-7,
            maxfeval = 2000,
            restarts.max = 10
        ))
    }
)

View(tmp4$parameters)

xyplot(A + An_estimate ~ Ci | event_replicate, data = tmp4$fits, type = 'b', pch = 20, auto = TRUE, grid = TRUE)

xyplot(gm ~ Ci | event, group = event_replicate, data = tmp4$fits, type = 'b', pch = 20, auto.key = list(space = 'right'), grid = TRUE, ylim=c(-1,3))

xyplot(gm ~ Cc | event, group = event_replicate, data = tmp4$fits, type = 'b', pch = 20, auto.key = list(space = 'right'),, grid = TRUE, ylim=c(-1,3))





tmp5 <- fit_variable_j_gjrttv(
    combined_info_data,
    'event_replicate',
    function(guess, fun, lower, upper) {
        dfoptim::nmkb(guess, fun, lower, upper, control = list(
            tol = 1e-7,
            maxfeval = 2000,
            restarts.max = 10
        ))
    },
    upper = c(500, 500, 10, 0.6, 20, 500)
)

View(tmp5$parameters)

a_cols <- c(
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#000000"
)

xyplot(
    Ac + Aj + Ap + An_estimate + A ~ Ci | event_replicate,
    data = tmp5$fits,
    type = 'b',
    pch = 20,
    auto = TRUE,
    grid = TRUE,
    ylim = c(-10, 60),
    par.settings=list(
        superpose.line=list(col=a_cols),
        superpose.symbol=list(col=a_cols)
    )
)

xyplot(
    Ac + Aj + Ap + An_estimate + A ~ Cc | event_replicate,
    data = tmp5$fits,
    type = 'b',
    pch = 20,
    auto = TRUE,
    grid = TRUE,
    ylim = c(-10, 50),
    xlim = c(0, 300),
    par.settings=list(
        superpose.line=list(col=a_cols),
        superpose.symbol=list(col=a_cols)
    )
)


xyplot(gm ~ Ci | event, group = event_replicate, data = tmp5$fits, type = 'b', pch = 20, auto.key = list(space = 'right'), grid = TRUE, ylim=c(-1,3))

xyplot(gm ~ Cc | event, group = event_replicate, data = tmp5$fits, type = 'b', pch = 20, auto.key = list(space = 'right'),, grid = TRUE, ylim=c(-1,3))






tmp <- list(event = "WT", replicate = "10-5", Tleaf2 = 1000, Gamma_star = 104, J_high = 235, Rd = 0.783, tau = 0.47, Vcmax = 158, convergence = 0, convergence_msg = NULL, optimum_val = 95.7)

 error      : num  3.699 0.127 6.322 8.855 7.338 ...
List of 11
 $ event          : chr "WT"
 $ replicate      : chr "10-5"
 $ Tleaf2         : num 1000
 $ Gamma_star     : num 104
 $ J_high         : num 235
 $ Rd             : num 0.783
 $ tau            : num 0.47
 $ Vcmax          : num 158
 $ convergence    : num 0
 $ convergence_msg: NULL
 $ optimum_val    : num 95.7
Error in (function (..., row.names






    if (FALSE) {
        A_COLUMN_NAME = 'A'
        CI_COLUMN_NAME = 'Ci'
        PHIPS2_COLUMN_NAME = 'PhiPS2'
        QIN_COLUMN_NAME = 'Qin'

        Kc = 404    # micromol / mol
        Ko = 278    # mmol / mol
        O = 210     # mmol / mol
        TPU = 1000  # micromol / m^2 / s

        initial_guess = c(Gamma_star = 37, J_high = 110, Rd = 1, tau = 0.47, Vcmax = 100)
        lower = c(0, 0, 0, 0.2, 0)
        upper = c(Inf, Inf, Inf, 0.6, Inf)

        ERROR_FUN <- dpmn_error_std(Kc, Ko, O, TPU)

        tmp <- fit_variable_j_replicate(
            one_rep,
            A_COLUMN_NAME,
            CI_COLUMN_NAME,
            PHIPS2_COLUMN_NAME,
            QIN_COLUMN_NAME,
            ERROR_FUN,
            dfoptim::hjkb,
            initial_guess,
            lower,
            upper
        )

        View(tmp$fits)

        str(tmp$parameters)

        xyplot(A + An_estimate ~ Ci, data = tmp$fits, type = 'b', pch = 20, auto = TRUE, grid = TRUE)
    }
