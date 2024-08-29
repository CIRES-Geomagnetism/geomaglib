import math
import numpy as np

def rad2deg(rad: float) -> float:

    """
        Convert radius to degree
    """
    
    return rad*180.0/math.pi

def deg2rad(deg):

    """
        Convert degree to radius
    """

    return deg*math.pi/180.0



def calc_Bp_Pole(nmax, geoc_lat, sph, coef_dict):

    PcupS = [0.0]*(nmax+1)

    PcupS[0] = 1.0

    schmidtQuasiNorm1 = 1.0


    Bp = 0.0
    sin_phi = math.sin(deg2rad(geoc_lat))

    for n in range(1, nmax):
        idx = int(n * (n + 1) / 2 + 1)

        schmidtQuasiNorm2 = schmidtQuasiNorm1 * (2 * n - 1) / n
        schmidtQuasiNorm3 = schmidtQuasiNorm2 * math.sqrt((n * 2) / (n + 1))
        schmidtQuasiNorm1 = schmidtQuasiNorm2

        if n == 1:
            PcupS[1] = 1.0
        else:
            k = (((n - 1) * (n - 1)) - 1) / ((2 * n - 1) * (2 * n - 3))
            PcupS[n] = sin_phi * PcupS[n - 1] - k * PcupS[n - 2]

        Bp += sph["relative_radius_power"][n] * (coef_dict["g"][idx] * sph["sin_mlon"][n] - coef_dict["h"][idx] * sph["cos_mlon"][n]) * PcupS[n] * schmidtQuasiNorm3

    return Bp

def mag_SPH_summation(nmax, sph, coef_dict, Leg, geoc_lat)->tuple:


    Br, Bt, Bp = 0.0, 0.0, 0.0

    legP = np.array(Leg[0]).flatten()
    legdP = np.array(Leg[1]).flatten()

    pidx = 1

    for m in range(nmax + 1):
        # degree
        for n in range(m, nmax + 1):
            if n == 0:
                continue
            gidx = int(n * (n + 1) / 2 + m)
            #gidx =gidx - 1

            Bt -= sph["relative_radius_power"][n] * (
                        coef_dict["g"][gidx] * sph["cos_mlon"][m] + coef_dict["h"][gidx] * sph["sin_mlon"][m]) * legdP[
                      pidx]

            Bp += sph["relative_radius_power"][n] * (
                    coef_dict["g"][gidx] * sph["sin_mlon"][m] - coef_dict["h"][gidx] * sph["cos_mlon"][m]) * m * legP[pidx]

            Br -= sph["relative_radius_power"][n] * (
                        coef_dict["g"][gidx] * sph["cos_mlon"][m] + coef_dict["h"][gidx] * sph["sin_mlon"][m]) * (
                              n + 1) * legP[pidx]
            pidx += 1

    cos_phi = math.cos(deg2rad(geoc_lat))

    if math.fabs(cos_phi) < 1.0e-10:
        Bp += calc_Bp_Pole(nmax, geoc_lat, sph, coef_dict)
    else:
        Bp = Bp / cos_phi

    Bt = -Bt

    return Bt, Bp, Br



def mag_SPH_summation_alf(nmax, sph, coef_dict, legP, legdP, geoc_lat)->tuple:


    Br, Bt, Bp = 0.0, 0.0, 0.0



    for n in range(1, nmax + 1):
        # degree
        for m in range(n+1):

            gidx = int(n * (n + 1) / 2 + m)

            Bt -= sph["relative_radius_power"][n] * (
                        coef_dict["g"][gidx] * sph["cos_mlon"][m] + coef_dict["h"][gidx] * sph["sin_mlon"][m]) * legdP[
                      gidx]

            Bp += sph["relative_radius_power"][n] * (
                    coef_dict["g"][gidx] * sph["sin_mlon"][m] - coef_dict["h"][gidx] * sph["cos_mlon"][m]) * m * legP[gidx]

            Br -= sph["relative_radius_power"][n] * (
                        coef_dict["g"][gidx] * sph["cos_mlon"][m] + coef_dict["h"][gidx] * sph["sin_mlon"][m]) * (
                              n + 1) * legP[gidx]


    cos_phi = math.cos(deg2rad(geoc_lat))

    if math.fabs(cos_phi) < 1.0e-10:
        Bp += calc_Bp_Pole(nmax, geoc_lat, sph, coef_dict)
    else:
        Bp = Bp / cos_phi



    return Bt, Bp, Br



def rotate_magvec(Bt, Bp, Br, geoc_lat, geod_lat):
    """
            Convert magnetic vector from spherical to geodetic

            Parameters:
            ___________

            B: magnetic vector based on geocentric
            geoc_lat: geocentric latitude
            geod_lat: geeodetic latitude

            Returns:
            _________

            B:array the magnetic vector based on geodetic
            B = [Bx, By, Bz]
    """


    psi = (math.pi/180.0) * (geoc_lat - geod_lat)


    Bz = Bt * math.sin(psi) + Br * math.cos(psi)
    Bx = Bt * math.cos(psi) - Br * math.sin(psi)
    By = Bp



    return Bx, By, Bz

    
def get_magh(B):
    """
        Compute the magnetic horizontal

        Parameters:
        ____________

        B:array the gepdetic magnetic vector
        B = [Bx, By, Bz]
        Returns:
        ____________

        h:float the magneitc horizontal

    """

    return math.sqrt(B[0]**2 + B[1]**2)

def get_magf(B):
    """
        Get the total intensity

        Parameters:
        ____________
        B:array the magnetic vector

        Returns:
        f:float the total intensity value
        _________
    """

    f  = math.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    return f 


