import math
from typing import Optional, Tuple, Union

import numpy as np


def rad2deg(rad: float) -> float:
    """
        Convert radius to degree
    """

    return rad * 180.0 / math.pi


def deg2rad(deg) ->float:
    """
        Convert degree to radius
    """

    return deg * math.pi / 180.0


def calc_Bp_Pole(nmax: int, geoc_lat: Union[list[float], np.ndarray], sph:dict[str, list[float]], g: list[float], h: list[float]) -> float:
    """
    Calculate the B_phi magnetic elements at pole
    Args:
        nmax: maximum degree
        geoc_lat: geocentric latitude in degree
        sph: the dict svaed with spherical harmonic varialbles like (a/r) ^ (n+2), cos_m(lon), and sin_m(lon)
        g: g coefficients
        h: h coefficients

    Returns:
        B_phi
    """

    if isinstance(geoc_lat, list):
        geoc_lat = np.array(geoc_lat)
    PcupS = [0.0] * (nmax + 1)

    PcupS[0] = 1.0

    schmidtQuasiNorm1 = 1.0

    Bp = 0.0
    sin_phi = np.sin(deg2rad(geoc_lat))

    for n in range(1, nmax):
        idx = int(n * (n + 1) / 2 + 1)

        schmidtQuasiNorm2 = schmidtQuasiNorm1 * (2 * n - 1) / n
        schmidtQuasiNorm3 = schmidtQuasiNorm2 * np.sqrt((n * 2) / (n + 1))
        schmidtQuasiNorm1 = schmidtQuasiNorm2

        if n == 1:
            PcupS[1] = 1.0
        else:
            k = (((n - 1) * (n - 1)) - 1) / ((2 * n - 1) * (2 * n - 3))
            PcupS[n] = sin_phi * PcupS[n - 1] - k * PcupS[n - 2]

        Bp += sph["relative_radius_power"][n] * (g[idx] * sph["sin_mlon"][n] - h[idx] * sph["cos_mlon"][n]) * PcupS[
            n] * schmidtQuasiNorm3

    return Bp


def mag_SPH_summation(nmax: int, sph: dict[str, list[float]], g: list[float], h: list[float], Leg: list[list[float]],geoc_lat: Union[list[float], np.ndarray]) -> tuple:
    """
    Compute the magnetic elements based on pole (B_theta, B_phi, B_radius)
    Args:
        nmax: max degree
        sph: the dict saved with spherical harmonic varialbles like (a/r) ^ (n+2), cos_m(lon), and sin_m(lon)
        g: g coefficients
        h: h coefficients
        Leg: legendre function array. Leg[0] for Plm array; Leg[1] for dPlm array.
        geoc_lat: geocentric latitude in degree

    Returns:
        B_theta, B_phi, B_radius
    """

    if isinstance(geoc_lat, list):
        geoc_lat = np.array(geoc_lat)

    Br, Bt, Bp = np.zeros(len(geoc_lat)),np.zeros(len(geoc_lat)),np.zeros(len(geoc_lat))

    
    legP = np.array(Leg[0])
    legdP = np.array(Leg[1])
    
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

    cos_phi = np.cos(deg2rad(geoc_lat))

    mask = np.abs(cos_phi) < 1.0e-10 
    # Apply calc_Bp_Pole where the mask is True, otherwise perform division
    Bp = np.where(mask, Bp + calc_Bp_Pole(nmax, geoc_lat, sph, g, h), Bp / cos_phi)

    Bt = -Bt

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

    Bz = Bt * np.sin(psi) + Br * np.cos(psi)
    Bx = Bt * np.cos(psi) - Br * np.sin(psi)
    By = Bp

    return Bx, By, Bz


class GeomagElements:
    def __init__(self, 
                 Bx: Union[float, np.ndarray, list],
                 By: Union[float, np.ndarray, list],
                 Bz: Union[float, np.ndarray, list],
                 dBx: Optional[Union[float, np.ndarray, list]] = None,
                 dBy: Optional[Union[float, np.ndarray, list]] = None,
                 dBz: Optional[Union[float, np.ndarray, list]] = None):
        """
        Compute magnetic elements.
        
        Args:
            Bx: float or np.ndarray
            By: float or np.ndarray
            Bz: float or np.ndarray
            dBx: Optional float or np.ndarray
            dBy: Optional float or np.ndarray
            dBz: Optional float or np.ndarray
        """
        self.Bx = np.asarray(Bx, dtype=np.float64)
        self.By = np.asarray(By, dtype=np.float64)
        self.Bz = np.asarray(Bz, dtype=np.float64)

        self.dBx = np.asarray(dBx, dtype=np.float64) if dBx is not None else None
        self.dBy = np.asarray(dBy, dtype=np.float64) if dBy is not None else None
        self.dBz = np.asarray(dBz, dtype=np.float64) if dBz is not None else None



    def get_Bh(self) -> float:
        """
            Compute the magnetic horizontal


            Returns:
            ____________

            h: magneitc horizontal elements

        """

        return np.sqrt(self.Bx ** 2 + self.By ** 2)

    def get_Bf(self) -> float:
        """
            Get the total intensity
            Returns:
            f: the total intensity value
            _________
        """

        f = np.sqrt(self.Bx ** 2 + self.By ** 2 + self.Bz ** 2)
        return f

    def get_Bdec(self) -> float:
        """
        Get the declination value
        """
        dec = rad2deg(np.arctan2(self.By, self.Bx))

        return dec

    def get_Binc(self) -> float:
        """
        Get the inclination value
        Returns:

        """
        Bh = self.get_Bh()
        inc = rad2deg(np.arctan2(self.Bz, Bh))

        return inc

    def get_all_base(self) -> dict[str, np.ndarray]:
        """
        Get Bx, By, Bz, Bh, Bf, Bdec and Binc in dict

        """
    
        mag_map = {}

        mag_map["x"] = np.asarray(self.Bx, dtype=np.float64)
        mag_map["y"] = np.asarray(self.By, dtype=np.float64)
        mag_map["z"] = np.asarray(self.Bz, dtype=np.float64)
        mag_map["h"] = np.asarray(self.get_Bh(), dtype=np.float64)
        mag_map["f"] = np.asarray(self.get_Bf(), dtype=np.float64)
        mag_map["dec"] = np.asarray(self.get_Bdec(), dtype=np.float64)
        mag_map["inc"] = np.asarray(self.get_Binc(), dtype=np.float64)
        return mag_map

    def get_all(self) -> dict[str, np.ndarray]:
        """

        Returns: all of magnetic elements:
        Bx, By, Bz, Bh, Bf, Bdec, Binc,
        dBx, dBy, dBz, dBh, dBf, dBdec and dBinc in dict

        """
        mag_map = {}

        mag_map["x"] = np.asarray(self.Bx, dtype=np.float64)
        mag_map["y"] = np.asarray(self.By, dtype=np.float64)
        mag_map["z"] = np.asarray(self.Bz, dtype=np.float64)
        h = np.asarray(self.get_Bh(), dtype=np.float64)
        f = np.asarray(self.get_Bf(), dtype=np.float64)

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
        mag_map["dinc"] = np.asarray((180 / math.pi * (h * self.dBz - self.Bz * mag_map["dh"])) / (f ** 2), dtype=np.float64)
        
        return mag_map

    def get_dBh(self) -> np.ndarray:
        """

        Returns: delta horizontal

        """
        h = self.get_Bh()
        return (self.Bx * self.dBx + self.By * self.dBy) / h

    def get_dBf(self) -> np.ndarray:
        """
        Returns: delta total intensity
        """
        f = self.get_Bf()
        return (self.Bx * self.dBx + self.By * self.dBy + self.Bz * self.dBz) / f

    def get_dBdec(self) -> np.ndarray:
        """
        Returns: delta declination value
        """
        h = self.get_Bh()
        return 180 / math.pi * (self.Bx * self.dBy - self.By * self.dBx) / (h ** 2)

    def get_dBinc(self) -> np.ndarray:
        """
        Returns: delta inclination value

        """
        h = self.get_Bh()
        f = self.get_Bf()
        dh = (self.Bx * self.dBx + self.By * self.dBy) / h
        return 180 / math.pi * (h * self.dBz - self.Bz * dh) / (f ** 2)
