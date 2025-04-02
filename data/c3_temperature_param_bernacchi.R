c3_temperature_param_bernacchi <- list(
    Gamma_star_norm  = list(type = 'Arrhenius',  c = 37.83 / PhotoGEA:::f, Ea = 37.83,       units = 'normalized to Gamma_star at 25 degrees C'),
    J_norm           = list(type = 'Arrhenius',  c = 17.57,                Ea = 43.5,        units = 'normalized to J at 25 degrees C'),
    Kc_norm          = list(type = 'Arrhenius',  c = 79.43 / PhotoGEA:::f, Ea = 79.43,       units = 'normalized to Kc at 25 degrees C'),
    Ko_norm          = list(type = 'Arrhenius',  c = 36.38 / PhotoGEA:::f, Ea = 36.38,       units = 'normalized to Ko at 25 degrees C'),
    RL_norm          = list(type = 'Arrhenius',  c = 18.72,                Ea = 46.39,       units = 'normalized to RL at 25 degrees C'),
    Vcmax_norm       = list(type = 'Arrhenius',  c = 26.35,                Ea = 65.33,       units = 'normalized to Vcmax at 25 degrees C'),
    Vomax_norm       = list(type = 'Arrhenius',  c = 22.98,                Ea = 60.11,       units = 'normalized to Vcmax at 25 degrees C'),
    gmc_norm         = list(type = 'Johnson',    c = 20.01, Ha = 49.6, Hd = 437.4, S = 1.4,  units = 'normalized to gmc at 25 degrees C'),
    Tp_norm          = list(type = 'Johnson',    c = 21.46, Ha = 53.1, Hd = 201.8, S = 0.65, units = 'normalized to Tp at 25 degrees C'),
    Gamma_star_at_25 = list(type = 'Polynomial', coef = 42.93205,                            units = 'micromol mol^(-1)'),
    Kc_at_25         = list(type = 'Polynomial', coef = 406.8494,                            units = 'micromol mol^(-1)'),
    Ko_at_25         = list(type = 'Polynomial', coef = 277.1446,                            units = 'mmol mol^(-1)')
)
