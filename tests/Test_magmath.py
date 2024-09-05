import unittest
from geomaglib import util, sh_loader, sh_vars, Leg_SHA_for_import, magmath



class Test_magmath(unittest.TestCase):

    def setUp(self):

        self.Bh = [1510.0, 1910.8, 2487.8, 24377.2, 21666.6]
        self.Bx = [-575.7, 1518.0, 1555.6, 24375.3, 21556.3]
        self.By = [-1396.0, -1160.5, 1941.4, 303.2, -2183.2]
        self.Bz = [56082.3, 55671.9, 56520.5, 49691.4, -52676.0]
        self.Bf = [56102.7,  55704.7, 56575.2, 55348.7, 56957.9]
        self.wmm_coeff = "WMM2020.cof"
    def test_magallB(self):

        N = len(self.Bx)

        for i in range(N):

            self.assertEqual(True, False)  # add assertion here



    def test_mag_SPH_summation(self):


        lat = -18
        lon = 138
        alt = 77

        dec_year = 2024.5

        # load g, h
        coef_dict = sh_loader.load_coef(self.wmm_coeff, skip_two_columns=True)
        timly_coef_dict = sh_loader.timely_modify_magnetic_model(coef_dict, dec_year)
        nmax = sh_loader.calc_num_elems_to_sh_degrees(len(coef_dict["g"]))

        # alt = util.alt_to_ellipsoid_height(alt, lat, lon)

        r, theta = util.geod_to_geoc_lat(lat, alt)
        sph_dict = sh_vars.comp_sh_vars(lon, r, theta, nmax)
        cotheta = 90 - theta
        colats = [cotheta]

        Leg = Leg_SHA_for_import.Flattened_Chaos_Legendre1(nmax, colats)

        Bt, Bp, Br = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict, Leg, theta)
        Bx, By, Bz = magmath.rotate_magvec(Bt, Bp, Br, theta, lat)

        self.assertAlmostEqual(round(Bx, 1), 31722.0, delta=0.01)
        self.assertAlmostEqual(round(By, 1), 2569.6, delta=0.01)
        self.assertAlmostEqual(round(Bz, 1), -34986.2, delta=0.01)

    def test_calc_Bp_Pole(self):
        lat = 90
        lon = 138
        alt = 77

        dec_year = 2024.5

        # load g, h
        coef_dict = sh_loader.load_coef(self.wmm_coeff, skip_two_columns=True)
        timly_coef_dict = sh_loader.timely_modify_magnetic_model(coef_dict, dec_year)
        nmax = sh_loader.calc_num_elems_to_sh_degrees(len(coef_dict["g"]))

        # Assume alt already in WGS, otherwise, use util.alt_to_ellipsoid_height(alt, lat, lon) to transform

        r, theta = util.geod_to_geoc_lat(lat, alt)
        sph_dict = sh_vars.comp_sh_vars(lon, r, theta, nmax)
        cotheta = 90 - theta
        colats = [cotheta]

        Leg = Leg_SHA_for_import.Flattened_Chaos_Legendre1(nmax, colats)

        Bt, Bp, Br = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict, Leg, theta)
        Bpole_t, Bpole_p, Bpole_r = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict, Leg, theta)

        self.assertEqual(Bt, Bpole_t)  # add assertion here





if __name__ == '__main__':
    unittest.main()
