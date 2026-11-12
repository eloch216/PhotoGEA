# PhotoGEA Publications

## Overview

This page contains a list of all publications about or using the
PhotoGEA R package, including links to publicly-available analysis
scripts whenever possible. If we are missing any publications, please
let us know!

## Publications About PhotoGEA

The PhotoGEA R package was first described in:

1.  Lochocki, E. B., Salesse-Smith, C. E. & McGrath, J. M. “PhotoGEA: An
    R Package for Closer Fitting of Photosynthetic Gas Exchange Data
    With Non-Gaussian Confidence Interval Estimation.” *Plant, Cell &
    Environment* **48**, 5104–5119 (2025).

    > - [DOI:10.1111/pce.15501](https://doi.org/10.1111/pce.15501)
    >
    > - PhotoGEA version: `v1.1.0`
    >
    > - Analysis scripts:
    >   <https://github.com/ripeproject/PhotoGEA-paper/>

If you use PhotoGEA for your own work, please cite this publication and
specify the version of the package you used.

## Publications Using PhotoGEA

2.  Salesse-Smith, C. E. *et al*. “Greater mesophyll conductance and
    leaf photosynthesis in the field through modified cell wall porosity
    and thickness via AtCGR3 expression in tobacco.” *Plant
    Biotechnology Journal* **22**, 2504–2517 (2024).

    > - [DOI:10.1111/pbi.14364](https://doi.org/10.1111/pbi.14364)
    >
    > - PhotoGEA used for: Calculating mesophyll conductance from
    >   tunable diode laser absorption measurements
    >   (`calculate_gm_busch`), C₃*A-C_(i)* fitting (`fit_c3_aci`),
    >   C₃*A-C_(i)* + CF Variable *J* fitting (`fit_c3_variable_j`)
    >
    > - PhotoGEA version: `v0.10.0`
    >
    > - Analysis scripts:
    >   <https://github.com/ripeproject/CGR3-tobacco-2024>

3.  Pelech, E. A., Stutz, S. S., Wang, Y., Lochocki, E. B. & Long, S. P.
    “Have We Selected for Higher Mesophyll Conductance in Domesticating
    Soybean?” *Plant, Cell & Environment* **48**, 1594–1607 (2025).

    > - [DOI:10.1111/pce.15206](https://doi.org/10.1111/pce.15206)
    >
    > - PhotoGEA used for: Calculating limitations to C₃ photosynthesis
    >   (`calculate_c3_limitations_grassi` and
    >   `calculate_c3_limitations_warren`)
    >
    > - PhotoGEA version: `v1.0.0`
    >
    > - Analysis scripts: <https://doi.org/10.13012/B2IDB-7809185_V2>

4.  Salesse-Smith, C. E. *et al*. “Adapting C₄ photosynthesis to
    atmospheric change and increasing productivity by elevating Rubisco
    content in sorghum and sugarcane.” *Proceedings of the National
    Academy of Sciences* **122**, e2419943122 (2025).

    > - [DOI:10.1073/pnas.2419943122](https://doi.org/10.1073/pnas.2419943122)
    >
    > - PhotoGEA used for: Fitting C₄*A-C_(i)* curves using the
    >   mechanistic model (`fit_c4_aci`), reading and plotting induction
    >   curves (`read_gasex_file`)
    >
    > - PhotoGEA version: `v0.11.0`
    >
    > - Analysis scripts:
    >   <https://github.com/cabbi-bio/sorghum-sugarcane-RBCS-RAF1-2024>

5.  Lochocki, E. B. & McGrath, J. M. “Widely Used Variants of the
    Farquhar-von-Caemmerer-Berry Model Can Cause Errors in Parameter
    Estimation.” *in silico Plants* diaf014 (2025)

    > - [DOI:10.1093/insilicoplants/diaf014](https://doi.org/10.1093/insilicoplants/diaf014)
    >
    > - PhotoGEA used for: Fitting C₃*A-C_(i)* curves (`fit_c3_aci`)
    >   using several variants of the FvCB model
    >
    > - PhotoGEA version: `v1.2.0`
    >
    > - Analysis scripts: <https://github.com/ripeproject/FvCB-min-A>
    >
    > - Preprint available from *bioRxiv*:
    >   [DOI:10.1101/2025.03.11.642611](https://doi.org/10.1101/2025.03.11.642611)

6.  Tamang, B. G., Bernard, G., Bernacchi, C. J., Diers, B. W. &
    Ainsworth, E. A. “Bigger is not always better: Optimizing leaf area
    index with narrow leaf shape in soybean.” *Plant Physiology* kiaf663
    (2025).

    > - [DOI:10.1093/plphys/kiaf663](https://doi.org/10.1093/plphys/kiaf663)
    >
    > - PhotoGEA used for: Fitting C₃*A-C_(i)* curves (`fit_c3_aci`)
    >
    > - PhotoGEA version: `v1.1.0`
    >
    > - Analysis scripts: Not yet available
    >
    > - Preprint available from *bioRxiv*:
    >   [DOI:10.1101/2025.07.07.663573](https://doi.org/10.1101/2025.07.07.663573)

7.  Garen, J. C. & Crous, K. Y. “Choose Wisely: Parameter Choice is Key
    for Ensuring Consistent Estimates of Photosynthetic Capacity From
    A–Ci Response Curves.” *Plant, Cell & Environment* **49**, 3251–3253
    (2026).

    > - [DOI:10.1111/pce.70462](https://doi.org/10.1111/pce.70462)
    >
    > - PhotoGEA used for: Fitting C₃*A-C_(i)* curves (`fit_c3_aci`)
    >
    > - PhotoGEA version: `v1.3.2`
    >
    > - Analysis scripts: Not yet available

8.  Mbaluto, C. M., Martínez-Goñi, X. S., Tripathi, A. & Singh, P.
    “Cereal grafting in rice and pearl millet preserves photosynthetic
    performance and stomatal dynamics, establishing a platform for
    root-shoot communication studies.” *bioRxiv* preprint (2026).

    > - [DOI:10.64898/2026.06.25.734459](https://doi.org/10.64898/2026.06.25.734459)
    >
    > - PhotoGEA used for: Fitting C₃ and C₄*A-C_(i)* curves
    >   (`fit_c3_aci` and `fit_c4_aci`)
    >
    > - PhotoGEA version: Not specified
    >
    > - Analysis scripts: Not yet available

9.  Jeong, J. *et al*. “Chloroplast expression of *Chlamydomonas*
    glycolate dehydrogenase en route to an improved photorespiratory
    bypass.” *bioRxiv* preprint (2026).

    > - [DOI:10.64898/2026.07.23.739952](https://doi.org/10.64898/2026.07.23.739952)
    >
    > - PhotoGEA used for: Fitting C₃*A-C_(i)* curves (`fit_c3_aci`),
    >   calculating mesophyll conductance from tunable diode laser
    >   absorption measurements (`calculate_gm_busch`)
    >
    > - PhotoGEA version: `v1.3.2`
    >
    > - Analysis scripts: Not yet available

10. Piao, X., Lochocki, E. B., McGrath, J. M. & Matthews, M. L.
    “Integrating carbon utilization and transport processes into a crop
    growth model enables the prediction of emergent soybean carbon
    allocation behavior. *bioRxiv* preprint (2026).

    > - [DOI:10.64898/2026.08.27.747615](https://doi.org/10.64898/2026.08.27.747615)
    >
    > - PhotoGEA used for: Fitting C₃*A-C_(i)* curves (`fit_c3_aci`)
    >
    > - PhotoGEA version: `v1.3.3`
    >
    > - Analysis scripts: Not yet available

11. Kwafo, T. S., Yerkes, K. & McGrath, J. M. “Nutrient Concentration
    Declines but Nutrient Yield Increases Under Elevated CO2
    Concentration in Soybean: Disentangling the relative effects of
    Dilution, Transpiration, and Uptake Activity.” *bioRxiv* preprint
    (2026).

    > - [DOI:10.64898/2026.08.28.747630](https://doi.org/10.64898/2026.08.28.747630)
    >
    > - PhotoGEA used for: Fitting Ball-Berry curves (`fit_ball_berry`)
    >
    > - PhotoGEA version: Not specified
    >
    > - Analysis scripts: Not yet available

12. de Oliveira, L. P. *et al*. “A century of soybean breeding increased
    photosynthetic capacity but not NPQ relaxation.” *Plant Physiology*
    kiag811 (2026).

    > - [DOI:10.1093/plphys/kiag811](https:doi.org/10.1093/plphys/kiag811)
    >
    > - PhotoGEA used for: Fitting C₃*A-C_(i)* curves (`fit_c3_aci`)
    >
    > - PhotoGEA version: Not specified
    >
    > - Analysis scripts: Not yet available
