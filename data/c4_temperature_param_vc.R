c4_temperature_param_vc <- list(
    Vcmax_norm = list(type = 'Arrhenius', c = 78 / PhotoGEA:::f,                     Ea = 78,   units = 'normalized to Vcmax at 25 degrees C'),
    Vpmax_norm = list(type = 'Arrhenius', c = 50.1 / PhotoGEA:::f,                   Ea = 50.1, units = 'normalized to Vpmax at 25 degrees C'),
    RL_norm =    list(type = 'Arrhenius', c = 66.4 / PhotoGEA:::f,                   Ea = 66.4, units = 'normalized to RL at 25 degrees C'),
    Kc =         list(type = 'Arrhenius', c = log(1210) + 64.2 / PhotoGEA:::f,       Ea = 64.2, units = 'microbar'),
    Ko =         list(type = 'Arrhenius', c = log(292) + 10.5 / PhotoGEA:::f,        Ea = 10.5, units = 'mbar'),
    Kp =         list(type = 'Arrhenius', c = log(82) + 38.3 / PhotoGEA:::f,         Ea = 38.3, units = 'microbar'),
    gamma_star = list(type = 'Arrhenius', c = log(0.5 / 1310) + 31.1 / PhotoGEA:::f, Ea = 31.1, units = 'dimensionless'),
    ao =         list(type = 'Arrhenius', c = log(0.047) + 1.63 / PhotoGEA:::f,      Ea = 1.63, units = 'dimensionless'),
    gmc_norm =   list(type = 'Arrhenius', c = 49.8 / PhotoGEA:::f,                   Ea = 49.8, units = 'normalized to gmc at 25 degrees C'),
    J_norm =     list(type = 'Gaussian', optimum_rate = exp((25 - 43)^2 / 26^2), t_opt = 43, sigma = 26, units = 'normalized to J at 25 degrees C')
)
