c4_arrhenius_von_caemmerer <- list(
    Vcmax_norm = list(c = 78 / PhotoGEA:::f,                     Ea = 78,   units = 'normalized to Vcmax at 25 degrees C'),
    Vpmax_norm = list(c = 50.1 / PhotoGEA:::f,                   Ea = 50.1, units = 'normalized to Vpmax at 25 degrees C'),
    Rd_norm =    list(c = 66.4 / PhotoGEA:::f,                   Ea = 66.4, units = 'normalized to Rd at 25 degrees C'),
    Kc =         list(c = log(1210) + 64.2 / PhotoGEA:::f,       Ea = 64.2, units = 'microbar'),
    Ko =         list(c = log(292) + 10.5 / PhotoGEA:::f,        Ea = 10.5, units = 'mbar'),
    Kp =         list(c = log(82) + 38.3 / PhotoGEA:::f,         Ea = 38.3, units = 'microbar'),
    gamma_star = list(c = log(0.5 / 1310) + 31.1 / PhotoGEA:::f, Ea = 31.1, units = 'dimensionless'),
    ao =         list(c = log(0.047) + 1.63 / PhotoGEA:::f,      Ea = 1.63, units = 'dimensionless'),
    gmc=         list(c = log(1) + 49.8 / PhotoGEA:::f,          Ea = 49.8, units = 'mol m^(-2) s^(-1) bar^(-1)')
)
