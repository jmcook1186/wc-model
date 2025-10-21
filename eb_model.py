# -*- coding: utf-8 -*-
"""
ebmodel.py

This module implements the point-based glacier surface energy balance model
originally developed as a spreadsheet application by:

Brock, B. W., and Arnold, N. (2000) Technical Communication: A spreadsheet-
based (Microsoft Excel) point surface energy balance model for glacier and
snow melt studies. Earth Surface Processes and Landforms, 25, p 649-658.

This code was modified from an earlier Python implementation by
Andrew Tedstone, May 2018
Which can be found here: https://github.com/atedstone/ebmodel/tree/master


"""

import numpy as np
import math
from typing import Tuple, Union

# Physical Constants
DEG_TO_RAD = math.pi / 180
RAD_TO_DEG = 180 / math.pi
STEFAN_BOLTZMANN = 5.7e-8  # W m-2 K-4
LATENT_HEAT_VAPORIZATION = 2500000  # J kg-1
SPECIFIC_GAS_CONSTANT_DRY_AIR = 0.622
VON_KARMAN_CONSTANT = 0.4
AIR_DENSITY = 1.225  # kg m-3
KINEMATIC_VISCOSITY = 0.00001461  # m2 s-1
GRAVITY = 9.81  # m s-2
MELT_CONVERSION_FACTOR = 0.0107784  # Conversion from W m-2 to mm w.e.
SOLAR_CONSTANT = 1368  # W m-2
ATMOSPHERIC_PRESSURE_SEA_LEVEL = 100000  # Pa


def deg_to_rad(degrees: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert degrees to radians."""
    return degrees * DEG_TO_RAD


def rad_to_deg(radians: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert radians to degrees."""
    return radians * RAD_TO_DEG


def calc_dayang(day: float) -> float:
    """Calculate day angle (degrees) from day of year."""
    return (day / 365.25) * 360


def calc_eqtim(dayang: float) -> float:
    """Calculate equation of time correction."""
    return (-0.128 * np.sin(deg_to_rad(dayang - 2.8)) -
            0.165 * np.sin(deg_to_rad(2 * dayang + 19.7)))


def calc_solhour(time: float, lon: float, lon_ref: float,
                 eqtim: float, summertime: float) -> float:
    """Calculate solar hour from time and location."""
    return (time / 100) + ((lon - lon_ref) / 15) + eqtim - summertime



def calc_soldec(day: float) -> float:
    """Calculate solar declination angle (degrees)."""
    i1 = -23.2559 * np.cos(((2 * math.pi * day) / 365) + 0.1582)
    i2 = -0.3915 * np.cos(((4 * math.pi * day) / 365) + 0.0934)
    i3 = -0.1764 * np.cos(((6 * math.pi * day) / 365) + 0.4539)
    return 0.3948 + (i1 + i2 + i3)


def calc_solhran(solhour: float) -> float:
    """Calculate solar hour angle (degrees)."""
    return 15 * (solhour - 12)


def calc_solaltr(lat: float, soldec: float, solhran: float) -> float:
    """Calculate solar altitude in radians."""
    lat_rad = deg_to_rad(lat)
    soldec_rad = deg_to_rad(soldec)
    solhran_rad = deg_to_rad(solhran)

    return np.arcsin(
        np.sin(lat_rad) * np.sin(soldec_rad) +
        np.cos(lat_rad) * np.cos(soldec_rad) * np.cos(solhran_rad)
    )


def calc_solaltd(solaltr: float) -> float:
    """Calculate solar altitude in degrees."""
    return rad_to_deg(solaltr)


def calc_cossolaz(lat: float, solaltr: float, soldec: float) -> float:
    """Calculate cosine of solar azimuth."""
    lat_rad = deg_to_rad(lat)
    soldec_rad = deg_to_rad(soldec)

    return ((np.sin(lat_rad) * np.sin(solaltr) - np.sin(soldec_rad)) /
            (np.cos(lat_rad) * np.cos(solaltr)))


def calc_sinolaz(soldec: float, solhran: float, solaltr: float) -> float:
    """Calculate sine of solar azimuth."""
    soldec_rad = deg_to_rad(soldec)
    solhran_rad = deg_to_rad(solhran)

    return (np.cos(soldec_rad) * np.sin(solhran_rad) / np.cos(solaltr))


def calc_solaz(sinsolaz: float, cossolaz: float) -> float:
    """Calculate solar azimuth (radians)."""
    if sinsolaz < 0:
        return -np.arccos(cossolaz)
    else:
        return np.arccos(cossolaz)


def calc_solaz360(solaz: float) -> float:
    """Calculate solar azimuth in 0-360 degree format."""
    return rad_to_deg(solaz) + 180



def calc_cloudn(inswrad: float, elevation: float, solaltr: float) -> float:
    """Calculate cloud cover.

    B+A eqn. 3
    """
    elevation_factor = ATMOSPHERIC_PRESSURE_SEA_LEVEL / (
        ATMOSPHERIC_PRESSURE_SEA_LEVEL * (1 - elevation * 0.0001)
    )
    clear_sky_radiation = (SOLAR_CONSTANT * (0.75 ** elevation_factor) *
                          np.cos(1.57 - solaltr))
    return np.min([1 - (inswrad / clear_sky_radiation), 1])



def calc_diffuser(cloudn: float) -> float:
    """Calculate diffuse fraction of total incoming shortwave radiation.

    B+A eqn. 2
    """
    if cloudn > 0.8:
        return 0.8
    else:
        return 0.65 * np.max([cloudn, 0]) + 0.15


def calc_directr(diffuser: float) -> float:
    """Calculate direct (incident) fraction of total incoming shortwave radiation.

    Following B+A eqn. 2
    """
    return 1 - diffuser


def calc_Qn(directr: float, inswrad: float, solaltr: float) -> float:
    """Calculate radiation received at surface normal to sun's rays.

    B+A eqn. 4
    """
    return directr * inswrad / np.sin(solaltr)


def calc_Qi(Qn: float, solaltr: float, slope: float,
            solaz: float, aspect: float) -> float:
    """Calculate direct (incident) component of shortwave radiation.

    B+A eqn. 5
    """
    slope_rad = deg_to_rad(slope)
    aspect_rad = deg_to_rad(aspect)

    return Qn * (
        np.sin(solaltr) * np.cos(slope_rad) +
        np.cos(solaltr) * np.sin(slope_rad) * np.cos(solaz - aspect_rad)
    )


def calc_Q_Wm2(albedo: float, Qi: float, diffuser: float,
               inswrad: float, slope: float) -> float:
    """Calculate net shortwave radiation (W m-2).

    B+A eqn. 7
    """
    slope_rad = deg_to_rad(slope)
    half_slope = slope_rad / 2

    direct_term = (1 - albedo) * Qi
    diffuse_term = ((1 - albedo) * diffuser *
                   (inswrad * (np.cos(half_slope) ** 2) +
                    albedo * (np.sin(half_slope) ** 2)))

    return direct_term + diffuse_term


def calc_Q_melt(Q_Wm2: float) -> float:
    """Calculate melting by (net) shortwave radiation (mm w.e.)."""
    return Q_Wm2 * MELT_CONVERSION_FACTOR



def calc_eo(airtemp: float, elevation: float,
            met_elevation: float, lapse: float) -> float:
    """Calculate clear sky emissivity.

    B+A see after eqn. 10
    """
    adjusted_temp = airtemp - (lapse * (elevation - met_elevation)) + 273.16
    return 8.733 * (0.001 * (adjusted_temp ** 0.788))


def calc_e_star(cloudn: float, eo: float) -> float:
    """Calculate sky effective emissivity.

    B+A eqn. 10
    """
    return (1 + 0.26 * np.max([cloudn, 0])) * eo


def calc_lwin(e_star: float, airtemp: float, lapse: float,
              elevation: float, met_elevation: float) -> float:
    """Calculate downwelling longwave radiation (W m-2).

    B+A eqn. 9
    """
    adjusted_temp = airtemp - (lapse * (elevation - met_elevation)) + 273.16
    return e_star * STEFAN_BOLTZMANN * (adjusted_temp ** 4)


def calc_lnet_Wm2(lwin: float, lwout: float = 316.) -> float:
    """Calculate net longwave radiation (W m-2).

    B+A eqn. 8
    """
    return lwin - lwout


def calc_lnet_melt(lnet_Wm2: float) -> float:
    """Calculate melt due to (net) longwave radiation (mm w.e.)."""
    return lnet_Wm2 * MELT_CONVERSION_FACTOR



def calc_atmpres(elevation: float) -> float:
    """Calculate atmospheric pressure (Pa)."""
    return ATMOSPHERIC_PRESSURE_SEA_LEVEL * (1 - elevation * 0.0001)


def calc_de(avp: float, elevation: float,
            met_elevation: float, lapse: float) -> float:
    """Calculate vapor pressure deficit (Pa)."""
    return (avp - (45 * (elevation - met_elevation) * lapse)) - 610.8


def calc_spehum(avp: float, atmpres: float) -> float:
    """Calculate specific humidity of air."""
    return avp / atmpres


def calc_spehtair(spehum: float) -> float:
    """Calculate specific heat of air at constant pressure (J kg-1 K-1).

    B+A para following eqn. 12
    """
    return 1004.67 * (1 + 0.84 * spehum)


def calc_Re_star(U_star: float, roughness: float) -> float:
    """Estimate scaling length for roughness.

    B+A para following eqn. 12 and 16
    """
    return (U_star * roughness) / KINEMATIC_VISCOSITY


def calc_ln_zt(roughness: float, Re_star: float) -> float:
    """Estimate scaling length for temperature.

    B+A para following eqn. 12 and 16
    """
    ln_Re = np.log(Re_star)
    return np.log(roughness) + 0.317 - 0.565 * ln_Re - 0.183 * (ln_Re ** 2)


def calc_ln_ze(roughness: float, Re_star: float) -> float:
    """Estimate scaling length for humidity.

    B+A para following eqn. 12 and 16
    """
    ln_Re = np.log(Re_star)
    return np.log(roughness) + 0.396 - 0.512 * ln_Re - 0.18 * (ln_Re ** 2)



def calc_L(U_star: float, airtemp: float, lapse: float, elevation: float,
           met_elevation: float, spehtair: float, Qs_Wm2: float) -> float:
    """Estimate Monin-Obukhov length scale.

    B+A eqn. 13
    """
    adjusted_temp = airtemp - (lapse * (elevation - met_elevation)) + 273.16
    avg_temp = (adjusted_temp + 273.16) / 2

    numerator = AIR_DENSITY * (U_star ** 3) * avg_temp * spehtair
    denominator = VON_KARMAN_CONSTANT * GRAVITY * Qs_Wm2

    return numerator / denominator


def calc_U_star_z_L(windspd: float, roughness: float) -> float:
    """Initial estimate of friction velocity.

    B+A para following eqn. 14
    """
    return (VON_KARMAN_CONSTANT * windspd) / np.log(2 / roughness)


def calc_U_star_full(windspd: float, roughness: float, L: float) -> float:
    """Full estimate of friction velocity.

    B+A eqn. 14
    """
    return (VON_KARMAN_CONSTANT * windspd) / (np.log(2 / roughness) + 5 * (2 / L))


def calc_Qsz_L0(spehtair: float, windspd: float, airtemp: float, lapse: float,
                elevation: float, met_elevation: float, roughness: float,
                ln_zt: float) -> float:
    """Initial estimate of sensible heat flux (without L).

    B+A para following eqn. 14
    """
    adjusted_temp = airtemp - (lapse * (elevation - met_elevation))
    numerator = AIR_DENSITY * spehtair * 0.16 * windspd * adjusted_temp
    denominator = np.log(2 / roughness) * (np.log(2) - ln_zt)

    return numerator / denominator


def calc_Qlz_L0(windspd: float, de: float, atmpres: float,
                roughness: float, ln_ze: float) -> float:
    """Initial estimate for latent heat flux (without L).

    B+A para following eqn. 14
    """
    numerator = (AIR_DENSITY * SPECIFIC_GAS_CONSTANT_DRY_AIR *
                LATENT_HEAT_VAPORIZATION * 0.16 * windspd * (de / atmpres))
    denominator = np.log(2 / roughness) * (np.log(2) - ln_ze)

    return numerator / denominator


def calc_Qs_full(spehtair: float, windspd: float, airtemp: float, lapse: float,
                 elevation: float, met_elevation: float, roughness: float,
                 L: float, ln_zt: float) -> float:
    """Full estimate of Sensible Heat Flux W m-2.

    B+A eqn. 11
    """
    adjusted_temp = airtemp - (lapse * (elevation - met_elevation))
    stability_term = 5 * (2 / L)

    numerator = AIR_DENSITY * spehtair * 0.16 * windspd * adjusted_temp
    denominator = ((np.log(2 / roughness) + stability_term) *
                  (np.log(2) - ln_zt + stability_term))

    return numerator / denominator


def calc_QI_full(windspd: float, de: float, atmpres: float, roughness: float,
                 ln_ze: float, L: float) -> float:
    """Full estimate of Latent Heat Flux W m-2.

    B+A eqn. 12
    """
    stability_term = 5 * (2 / L)

    numerator = (AIR_DENSITY * SPECIFIC_GAS_CONSTANT_DRY_AIR *
                LATENT_HEAT_VAPORIZATION * 0.16 * windspd * (de / atmpres))
    denominator = ((np.log(2 / roughness) + stability_term) *
                  (np.log(2) - ln_ze + stability_term))

    return numerator / denominator



def calc_turbulent_fluxes(
        windspd: float, roughness: float, spehtair: float, airtemp: float,
        de: float, atmpres: float, lapse: float, elevation: float,
        met_elevation: float, max_n: int = 25, tol: float = 0.001,
        verbose: bool = False, return_steps: bool = False
) -> Union[Tuple[float, float], Tuple[float, float, float, float, float, float, float]]:
    """Iteratively solve the turbulent fluxes and the Monin-Obukhov scale.

    B+A page 653.

    Args:
        windspd: Wind speed (m s-1)
        roughness: Surface roughness length (m)
        spehtair: Specific heat of air (J kg-1 K-1)
        airtemp: Air temperature (°C)
        de: Vapor pressure deficit (Pa)
        atmpres: Atmospheric pressure (Pa)
        lapse: Temperature lapse rate (°C m-1)
        elevation: Surface elevation (m)
        met_elevation: Meteorological station elevation (m)
        max_n: Maximum number of iterations
        tol: Tolerance threshold for iteration
        verbose: If True then print information about iterations
        return_steps: If True then return values of intermediary calculations

    Returns:
        tuple (SHF, LHF) in units W m-2

        if return_steps=True then returns:
        tuple (U_star, Re_star, ln_zt, ln_ze, Qs_Wm2, QI_Wm2, L)
    """
    # Solve initial conditions
    U_star = calc_U_star_z_L(windspd, roughness)
    Re_star = calc_Re_star(U_star, roughness)

    if verbose:
        print(f'Initial conditions: U_star: {U_star}, Re_star: {Re_star}')

    # Iteration block 1
    if verbose:
        print('Iteration block 1 . . . ')
        print('ln_zt  ln_ze  Qs_Wm2  Ql_Wm2')

    for _ in range(max_n):
        ln_zt = calc_ln_zt(roughness, Re_star)
        ln_ze = calc_ln_ze(roughness, Re_star)
        Qs_Wm2 = calc_Qsz_L0(spehtair, windspd, airtemp, lapse, elevation,
                             met_elevation, roughness, ln_zt)
        Ql_Wm2 = calc_Qlz_L0(windspd, de, atmpres, roughness, ln_ze)

        if verbose:
            print(ln_zt, ln_ze, Qs_Wm2, Ql_Wm2)

    # Iteration block 2
    if verbose:
        print('Iteration block 2 . . . ')
        print('U_star  Re_star  ln_zt  ln_ze  Qs_Wm2  QI_Wm2  L')

    L = 0.
    L_old = 1.
    n = 0

    while (n < max_n) and (np.abs(L - L_old) > tol):
        L_old = L
        L = calc_L(U_star, airtemp, lapse, elevation, met_elevation, spehtair, Qs_Wm2)
        Re_star = calc_Re_star(U_star, roughness)
        ln_zt = calc_ln_zt(roughness, Re_star)
        ln_ze = calc_ln_ze(roughness, Re_star)
        Qs_Wm2 = calc_Qs_full(spehtair, windspd, airtemp, lapse, elevation,
                              met_elevation, roughness, L, ln_zt)
        QI_Wm2 = calc_QI_full(windspd, de, atmpres, roughness, ln_ze, L)
        U_star = calc_U_star_full(windspd, roughness, L)

        if verbose:
            print(U_star, Re_star, ln_zt, ln_ze, Qs_Wm2, QI_Wm2, L)

        n += 1

    if return_steps:
        return (U_star, Re_star, ln_zt, ln_ze, Qs_Wm2, QI_Wm2, L)
    else:
        return (Qs_Wm2, QI_Wm2)



def _should_suppress_turbulent_melt(windspd: float, airtemp: float) -> bool:
    """Check if turbulent flux melt should be suppressed based on conditions.

    Helper function for calc_shf_melt and calc_lhf_melt.
    """
    if windspd > 2:
        return False

    if (windspd / airtemp) < 0.3:
        return True

    if (-1.5 < airtemp < 1.5) and (windspd < 1.5):
        return True

    if (-2 < airtemp < 2) and (windspd < 1):
        return True

    return False


def calc_shf_melt(windspd: float, Qs_Wm2: float, airtemp: float) -> float:
    """Calculate melting (mm w.e.) by sensible heat flux."""
    if _should_suppress_turbulent_melt(windspd, airtemp):
        return 0
    return Qs_Wm2 * MELT_CONVERSION_FACTOR


def calc_lhf_melt(windspd: float, QI_Wm2: float, airtemp: float) -> float:
    """Calculate melting (mm w.e.) by latent heat flux."""
    if _should_suppress_turbulent_melt(windspd, airtemp):
        return 0
    return QI_Wm2 * MELT_CONVERSION_FACTOR



def calc_melt_total(SWR: float, LWR: float, SHF: float, LHF: float) -> float:
    """Sum melt components (mm w.e.)."""
    melt_total = SWR + LWR + SHF + LHF
    if melt_total < 0:
        return 0.
    else:
        return melt_total


def calculate_seb(
        lat: float, lon: float, lon_ref: float, day: float, time: float,
        summertime: float, slope: float, aspect: float, elevation: float,
        met_elevation: float, lapse: float, inswrad: float, avp: float,
        airtemp: float, windspd: float, albedo: float, roughness: float
) -> Tuple[float, float, float, float]:
    """Convenience function to solve energy balance for a single timestep.

    Args:
        lat: Latitude (degrees)
        lon: Longitude (degrees)
        lon_ref: Reference longitude (degrees)
        day: Day of year
        time: Time (HHMM format)
        summertime: Summertime adjustment (hours)
        slope: Surface slope (degrees)
        aspect: Surface aspect (degrees)
        elevation: Surface elevation (m)
        met_elevation: Meteorological station elevation (m)
        lapse: Temperature lapse rate (°C m-1)
        inswrad: Incoming shortwave radiation (W m-2)
        avp: Actual vapor pressure (Pa)
        airtemp: Air temperature (°C)
        windspd: Wind speed (m s-1)
        albedo: Surface albedo (0-1)
        roughness: Surface roughness length (m)

    Returns:
        Tuple of (net_shortwave, net_longwave, sensible_heat, latent_heat) in W m-2
    """
    # Solar characteristics/quantities
    dayang = calc_dayang(day)
    eqtim = calc_eqtim(dayang)
    solhour = calc_solhour(time, lon, lon_ref, eqtim, summertime)
    soldec = calc_soldec(day)
    solhran = calc_solhran(solhour)
    solaltr = calc_solaltr(lat, soldec, solhran)
    solaltd = calc_solaltd(solaltr)
    cossolaz = calc_cossolaz(lat, solaltr, soldec)
    sinolaz = calc_sinolaz(soldec, solhran, solaltr)
    solaz = calc_solaz(sinolaz, cossolaz)
    solaz360 = calc_solaz360(solaz)

    # Calculate shortwave radiation
    cloudn = calc_cloudn(inswrad, elevation, solaltr)
    diffuser = calc_diffuser(cloudn)
    directr = calc_directr(diffuser)
    Qn = calc_Qn(directr, inswrad, solaltr)
    Qi = calc_Qi(Qn, solaltr, slope, solaz, aspect)
    # Net shortwave radiation
    Q_Wm2 = calc_Q_Wm2(albedo, Qi, diffuser, inswrad, slope)

    # Calculate longwave radiation
    eo = calc_eo(airtemp, elevation, met_elevation, lapse)
    e_star = calc_e_star(cloudn, eo)
    lwin = calc_lwin(e_star, airtemp, lapse, elevation, met_elevation)
    # Net longwave radiation
    lnet_Wm2 = calc_lnet_Wm2(lwin)

    # Calculate turbulent fluxes
    atmpres = calc_atmpres(elevation)
    de = calc_de(avp, elevation, met_elevation, lapse)
    spehum = calc_spehum(avp, atmpres)
    spehtair = calc_spehtair(spehum)
    # Return individual sensible and latent heat fluxes
    Qs_Wm2, QI_Wm2 = calc_turbulent_fluxes(
        windspd, roughness, spehtair, airtemp, de, atmpres,
        lapse, elevation, met_elevation
    )

    return (Q_Wm2, lnet_Wm2, Qs_Wm2, QI_Wm2)


def calculate_melt(
        swnet: float, lwnet: float, shf: float, lhf: float,
        windspd: float, airtemp: float
) -> Tuple[float, float, float, float, float]:
    """Convert energy fluxes to melt quantities.

    Args:
        swnet: Net shortwave radiation (W m-2)
        lwnet: Net longwave radiation (W m-2)
        shf: Sensible heat flux (W m-2)
        lhf: Latent heat flux (W m-2)
        windspd: Wind speed (m s-1)
        airtemp: Air temperature (°C)

    Returns:
        Tuple of (swnet_melt, lwnet_melt, shf_melt, lhf_melt, melt_total) in mm w.e.
    """
    swnet_melt = calc_Q_melt(swnet)
    lwnet_melt = calc_lnet_melt(lwnet)
    SHF_melt = calc_shf_melt(windspd, shf, airtemp)
    LHF_melt = calc_lhf_melt(windspd, lhf, airtemp)
    melt_total = calc_melt_total(swnet_melt, lwnet_melt, SHF_melt, LHF_melt)

    return (swnet_melt, lwnet_melt, SHF_melt, LHF_melt, melt_total)
