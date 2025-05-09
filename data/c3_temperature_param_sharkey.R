c3_temperature_param_sharkey <- list(
    Gamma_star_norm  = list(type = 'Arrhenius', c = 24.46 / PhotoGEA:::f, Ea = 24.46,       units = 'normalized to Gamma_star at 25 degrees C'),
    J_norm           = list(type = 'Arrhenius', c = 17.71,                Ea = 43.9,        units = 'normalized to J at 25 degrees C'),
    Kc_norm          = list(type = 'Arrhenius', c = 80.99 / PhotoGEA:::f, Ea = 80.99,       units = 'normalized to Kc at 25 degrees C'),
    Ko_norm          = list(type = 'Arrhenius', c = 23.72 / PhotoGEA:::f, Ea = 23.72,       units = 'normalized to Ko at 25 degrees C'),
    RL_norm          = list(type = 'Arrhenius', c = 18.7145,              Ea = 46.39,       units = 'normalized to RL at 25 degrees C'),
    Vcmax_norm       = list(type = 'Arrhenius', c = 26.355,               Ea = 65.33,       units = 'normalized to Vcmax at 25 degrees C'),
    gmc_norm         = list(type = 'Johnson',   c = 20.01, Ha = 49.6, Hd = 437.4, S = 1.4,  units = 'normalized to gmc at 25 degrees C'),
    Tp_norm          = list(type = 'Johnson',   c = 21.46, Ha = 53.1, Hd = 201.8, S = 0.65, units = 'normalized to Tp at 25 degrees C'),
    Gamma_star_at_25 = list(type = 'Polynomial', coef = 36.94438,                           units = 'micromol mol^(-1)'),
    Kc_at_25         = list(type = 'Polynomial', coef = 269.3391,                           units = 'micromol mol^(-1)'),
    Ko_at_25         = list(type = 'Polynomial', coef = 163.7146,                           units = 'mmol mol^(-1)')
)
