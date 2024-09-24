import math
from typing import Optional, Tuple

import numpy as np


def rad2deg(rad: float) -> float:
    """
        Convert radius to degree
    """

    return rad * 180.0 / math.pi


def deg2rad(deg):
    """
        Convert degree to radius
    """

    return deg * math.pi / 180.0


def calc_Bp_Pole(nmax: int, geoc_lat: float, sph:dict[str, list[float]], g: list[float], h: list[float]) -> float:
    """
    Calculate the B_phi magnetic elements at pole
    Args:
        nmax: maximum degree
        geoc_lat: geocentric latitude in degree
        sph: the dict svaed with spherical harmonic varialbles like (a/r) ^ (n+2), cos_m(lon), and sin_m(lon)
        g: g coefficients
        h: h coefficients

    Returns:

    """
    PcupS = [0.0] * (nmax + 1)

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

        Bp += sph["relative_radius_power"][n] * (g[idx] * sph["sin_mlon"][n] - h[idx] * sph["cos_mlon"][n]) * PcupS[
            n] * schmidtQuasiNorm3

    return Bp


def mag_SPH_summation(nmax: int, sph: dict[str, list[float]], g: list[float], h: list[float], Leg: list[list[float]], geoc_lat: float) -> tuple:
    """
    Compute the magnetic eelements
    Args:
        nmax: max degree
        sph: the dict svaed with spherical harmonic varialbles like (a/r) ^ (n+2), cos_m(lon), and sin_m(lon)
        g: g coefficients
        h: h coefficients
        Leg: legendre function array. Leg[0] for Plm array; Leg[1] for dPlm array.
        geoc_lat: geocentric latitude in degree

    Returns:

    """
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

            Bt -= sph["relative_radius_power"][n] * (
                    g[gidx] * sph["cos_mlon"][m] + h[gidx] * sph["sin_mlon"][m]) * legdP[
                      pidx]

            Bp += sph["relative_radius_power"][n] * (
                    g[gidx] * sph["sin_mlon"][m] - h[gidx] * sph["cos_mlon"][m]) * m * legP[pidx]

            Br -= sph["relative_radius_power"][n] * (
                    g[gidx] * sph["cos_mlon"][m] + h[gidx] * sph["sin_mlon"][m]) * (
                          n + 1) * legP[pidx]
            pidx += 1

    cos_phi = math.cos(deg2rad(geoc_lat))

    if math.fabs(cos_phi) < 1.0e-10:
        Bp += calc_Bp_Pole(nmax, geoc_lat, sph, g, h)
    else:
        Bp = Bp / cos_phi

    Bt = -Bt

    return Bt, Bp, Br


def mag_SPH_summation_alf(nmax, sph, coef_dict, legP, legdP, geoc_lat) -> tuple:
    """
    Compute the magnetic elements based on Legendre function from geomag's team C library
    Args:
        nmax:
        sph:
        coef_dict:
        legP:
        legdP:
        geoc_lat:

    Returns:

    """
    Br, Bt, Bp = 0.0, 0.0, 0.0

    for n in range(1, nmax + 1):
        # degree
        for m in range(n + 1):
            gidx = int(n * (n + 1) / 2 + m)

            Bt -= sph["relative_radius_power"][n] * (
                    coef_dict["g"][gidx] * sph["cos_mlon"][m] + coef_dict["h"][gidx] * sph["sin_mlon"][m]) * legdP[
                      gidx]

            Bp += sph["relative_radius_power"][n] * (
                    coef_dict["g"][gidx] * sph["sin_mlon"][m] - coef_dict["h"][gidx] * sph["cos_mlon"][m]) * m * legP[
                      gidx]

            Br -= sph["relative_radius_power"][n] * (
                    coef_dict["g"][gidx] * sph["cos_mlon"][m] + coef_dict["h"][gidx] * sph["sin_mlon"][m]) * (
                          n + 1) * legP[gidx]

    cos_phi = math.cos(deg2rad(geoc_lat))

    if math.fabs(cos_phi) < 1.0e-10:
        Bp += calc_Bp_Pole(nmax, geoc_lat, sph, coef_dict["g"],coef_dict["h"])
    else:
        Bp = Bp / cos_phi

    return Bt, Bp, Br


def rotate_magvec(Bt, Bp, Br, geoc_lat, geod_lat) -> Tuple[float, float, float]:
    """
            Convert magnetic vector from spherical to geodetic

            Parameters:
            ___________

            Bt: magnetic elements theta
            Bp: magnetic elements phi
            Br: magnetic elements radius
            geoc_lat: geocentric latitude
            geod_lat: geeodetic latitude

            Returns:
            _________

            B:array the magnetic vector based on geodetic
            B = [Bx, By, Bz]
    """

    psi = (math.pi / 180.0) * (geoc_lat - geod_lat)

    Bz = Bt * math.sin(psi) + Br * math.cos(psi)
    Bx = Bt * math.cos(psi) - Br * math.sin(psi)
    By = Bp

    return Bx, By, Bz


class GeomagElements:

    def __init__(self, Bx: float, By: float, Bz: float, dBx: Optional[float] = None, dBy: Optional[float] = None, dBz: Optional[float] = None):
        """
        Compute magnetic elements
        Args:
            Bx: float type
            By: float type
            Bz: float type
            dBx: float type
            dBy: float type
            dBz: float type
        """
        self.Bx = float(Bx)
        self.By = float(By)
        self.Bz = float(Bz)

        self.dBx = dBx
        self.dBy = dBy
        self.dBz = dBz

        if isinstance(self.dBx, float): self.dBx = float(self.dBx)
        if isinstance(self.dBy, float): self.dBy = float(self.dBy)
        if isinstance(self.dBz, float): self.dBz = float(self.dBz)




    def get_Bh(self) -> float:
        """
            Compute the magnetic horizontal


            Returns:
            ____________

            h: magneitc horizontal elements

        """

        return math.sqrt(self.Bx ** 2 + self.By ** 2)

    def get_Bf(self) -> float:
        """
            Get the total intensity
            Returns:
            f: the total intensity value
            _________
        """

        f = math.sqrt(self.Bx ** 2 + self.By ** 2 + self.Bz ** 2)
        return f

    def get_Bdec(self) -> float:
        """
        Get the declination value
        """
        dec = rad2deg(math.atan2(self.By, self.Bx))

        return dec

    def get_Binc(self) -> float:
        """
        Get the inclination value
        Returns:

        """
        Bh = self.get_Bh()
        inc = rad2deg(math.atan2(self.Bz, Bh))

        return inc

    def get_all_base(self) -> dict[str, float]:
        """
        Get Bx, By, Bz, Bh, Bf, Bdec and Binc in dict

        """
        mag_map = {}

        mag_map["x"] = float(self.Bx)
        mag_map["y"] = float(self.By)
        mag_map["z"] = float(self.Bz)
        mag_map["h"] = float(self.get_Bh())
        mag_map["f"] = float(self.get_Bf())
        mag_map["dec"] = float(self.get_Bdec())
        mag_map["inc"] = float(self.get_Binc())

        return mag_map

    def get_all(self) -> dict[str, float]:
        """

        Returns: all of magnetic elements:
        Bx, By, Bz, Bh, Bf, Bdec, Binc,
        dBx, dBy, dBz, dBh, dBf, dBdec and dBinc in dict

        """
        mag_map = {}

        mag_map["x"] = float(self.Bx)
        mag_map["y"] = float(self.By)
        mag_map["z"] = float(self.Bz)
        h = float(self.get_Bh())
        f = float(self.get_Bf())

        mag_map["h"] = h
        mag_map["f"] = f
        mag_map["dec"] = self.get_Bdec()
        mag_map["inc"] = self.get_Binc()

        mag_map["dx"] = self.dBx
        mag_map["dy"] = self.dBy
        mag_map["dz"] = self.dBz
        mag_map["dh"] = (self.Bx * self.dBx + self.By * self.dBy) / h
        mag_map["df"] = (self.Bx * self.dBx + self.By * self.dBy + mag_map["z"] * self.dBz) / mag_map["f"]
        mag_map["ddec"] = 180 / math.pi * (self.Bx * self.dBy - self.By * self.dBx) / (h ** 2)
        mag_map["dinc"] = float(180 / math.pi * (h * self.dBz - self.Bz * mag_map["dh"])) / (f ** 2)

        return mag_map

    def get_dBh(self) -> float:
        """

        Returns: delta horizontal

        """
        h = self.get_Bh()
        return (self.Bx * self.dBx + self.By * self.dBy) / h

    def get_dBf(self) -> float:
        """
        Returns: delta total intensity
        """
        f = self.get_Bf()
        return (self.Bx * self.dBx + self.By * self.dBy + self.Bz * self.dBz) / f

    def get_dBdec(self) -> float:
        """
        Returns: delta declination value
        """
        h = self.get_Bh()
        return 180 / math.pi * (self.Bx * self.dBy - self.By * self.dBx) / (h ** 2)

    def get_dBinc(self) -> float:
        """
        Returns: delta inclination value

        """
        h = self.get_Bh()
        f = self.get_Bf()
        dh = (self.Bx * self.dBx + self.By * self.dBy) / h
        return 180 / math.pi * (h * self.dBz - self.Bz * dh) / (f ** 2)
