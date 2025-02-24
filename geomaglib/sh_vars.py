import numpy as np 
def comp_sh_vars(lon: float, r: float, theta: float, sh_degrees: int) -> dict[str, list[float]]:
    """
    Computes Spherical Variabale and put the results in a dictionary
    variable computed are (a/r) ^ (n+2), cos_m(lon), and sin_m(lon)

    Parameters:
    lon (float or int): longitude in degrees
    r (float or int): distance from the center of ellipsoid computed by geod_to_geoc function
    theta (float or int): geocentric latitude computed by geod_to_geoc function
    sh_degrees(int): Total number of spherical harmonic degrees

    Returns:
    dictionary: with keys relative_radius_power, cos_mlon, and sin_mlon for the three different value arrays
    """
    ellip_re = 6371.2 # Mean radius of ellipsoid

    sh_degrees = int(sh_degrees)
    sh_vars = {}
    cos_lon = np.cos(np.radians(lon))
    sin_lon = np.sin(np.radians(lon))

    sh_vars["relative_radius_power"] = []

    ellip_re_over_r = ellip_re/r

    sh_vars["relative_radius_power"].append((ellip_re_over_r) ** 2)

    for n in range(1,sh_degrees+1):
        sh_vars["relative_radius_power"].append(sh_vars["relative_radius_power"][n-1] * ellip_re_over_r)

    sh_vars["cos_mlon"] = []
    sh_vars["sin_mlon"] = []

    sh_vars["cos_mlon"].append(1.0)
    sh_vars["sin_mlon"].append(0.0)

    sh_vars["cos_mlon"].append(cos_lon)
    sh_vars["sin_mlon"].append(sin_lon)

    for m in range(2, sh_degrees+1):
        sh_vars["cos_mlon"].append(sh_vars["cos_mlon"][m-1] * cos_lon -  sh_vars["sin_mlon"][m-1] * sin_lon)
        sh_vars["sin_mlon"].append(sh_vars["cos_mlon"][m-1] * sin_lon +  sh_vars["sin_mlon"][m-1] * cos_lon)

    return sh_vars


