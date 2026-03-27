from .magmath import rad2deg, deg2rad, mag_SPH_summation, GeomagElements, calc_Bp_Pole, rotate_magvec
from .util import (cart_to_sph_deg, sph_deg_to_cart, geod_to_geoc_lat, alt_to_ellipsoid_height,
                    calc_dec_year, calc_dec_year_array, decimalYearToDateTime, jd2000)
from .sh_vars import comp_sh_vars
from .legendre import Flattened_Chaos_Legendre1
from .dipole import Dipole