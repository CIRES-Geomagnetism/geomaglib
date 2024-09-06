import unittest
import os
from geomaglib import util, sh_loader, sh_vars, magmath, legendre



class Test_magmath(unittest.TestCase):

    def setUp(self):

        self.Bh = [1510.0, 1910.8, 2487.8, 24377.2, 21666.6]
        self.Bx = [-575.7, 1518.0, 1555.6, 24375.3, 21556.3]
        self.By = [-1396.0, -1160.5, 1941.4, 303.2, -2183.2]
        self.Bz = [56082.3, 55671.9, 56520.5, 49691.4, -52676.0]
        self.Bf = [56102.7,  55704.7, 56575.2, 55348.7, 56957.9]
        self.Bdec = [-112.41,  -37.40, 51.30, 0.71, -5.78]
        self.Binc = [88.46, 88.03, 87.48, 63.87, -67.64]

        self.tol = 1e-1

        self.curr_dir = os.path.dirname(os.path.abspath(__file__))
        self.top_dir = os.path.dirname(self.curr_dir)
        self.wmm_coeff = os.path.join(self.curr_dir, "WMM2020.cof")
    def test_get_Bh(self):

        N = len(self.Bx)

        for i in range(N):
            h = magmath.get_Bh(self.Bx[i], self.By[i])

            self.assertAlmostEqual(h, self.Bh[i], delta=self.tol)  # add assertion here


    def test_get_Bf(self):

        N = len(self.Bf)

        for i in range(N):
            f = magmath.get_Bf(self.Bx[i], self.By[i], self.Bz[i])
            print(f)
            self.assertAlmostEqual(f, self.Bf[i], delta=self.tol)  # add assertion here

    def test_get_Bdec(self):

        N = len(self.Bdec)

        for i in range(N):
            dec = magmath.get_Bdec(self.Bx[i], self.By[i])

            self.assertAlmostEqual(round(dec, 2), self.Bdec[i], delta=0.01)  # add assertion here

    def test_get_Binc(self):

        N = len(self.Bdec)

        for i in range(N):
            inc = magmath.get_Binc(self.Bx[i], self.By[i], self.Bz[i])

            self.assertAlmostEqual(round(inc, 2), self.Binc[i], delta=0.01)  # add assertion here

    def test_get_allB(self):

        N = len(self.Bx)

        for i in range(N):
            map = magmath.get_allB(self.Bx[i], self.By[i], self.Bz[i])

            self.assertAlmostEqual(map["x"], self.Bx[i], delta=self.tol)
            self.assertAlmostEqual(map["y"], self.By[i], delta=self.tol)
            self.assertAlmostEqual(map["z"], self.Bz[i], delta=self.tol)
            self.assertAlmostEqual(map["h"], self.Bh[i], delta=self.tol)
            self.assertAlmostEqual(map["f"], self.Bf[i], delta=self.tol)
            self.assertAlmostEqual(map["dec"], self.Bdec[i], delta=0.01)
            self.assertAlmostEqual(map["inc"], self.Binc[i], delta=0.01)
    def test_mag_SPH_summation(self):


        lats = [89, 80, 82, 43, -33]
        lons = [-121, -96, 87, 93, 109]
        alts = [28, 48, 54, 65, 51]
        N = len(lats)

        dec_year = 2020.0

        # load g, h
        coef_dict = sh_loader.load_coef(self.wmm_coeff, skip_two_columns=True)
        timly_coef_dict = sh_loader.timely_modify_magnetic_model(coef_dict, dec_year)
        nmax = sh_loader.calc_num_elems_to_sh_degrees(len(coef_dict["g"]))

        # alt = util.alt_to_ellipsoid_height(alt, lat, lon)
        for i in range(N):
            r, theta = util.geod_to_geoc_lat(lats[i], alts[i])
            sph_dict = sh_vars.comp_sh_vars(lons[i], r, theta, nmax)
            cotheta = 90 - theta
            colats = [cotheta]

            Leg = legendre.Flattened_Chaos_Legendre1(nmax, colats)

            Bt, Bp, Br = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict, Leg, theta)
            x, y, z = magmath.rotate_magvec(Bt, Bp, Br, theta, lats[i])

            self.assertAlmostEqual(x, self.Bx[i], delta=self.tol)
            self.assertAlmostEqual(y, self.By[i], delta=self.tol)
            self.assertAlmostEqual(z, self.Bz[i], delta=self.tol)

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

        Leg = legendre.Flattened_Chaos_Legendre1(nmax, colats)


        theta_nopole = theta + 1e-3
        Bt, Bp, Br = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict, Leg, theta_nopole)
        print(f"Bt: {Bt}, Br: {Br}, Bp: {Bp}")
        Bpole_t, Bpole_p, Bpole_r = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict, Leg, theta)
        print(f"Bt: {Bpole_t}, Br: {Bpole_r}, Bp: {Bpole_p}")

        #self.assertEqual(Bt, Bpole_t)  # add assertion here
        self.assertEqual(Bp, Bpole_p)
        #self.assertEqual(Br, Bpole_r)





if __name__ == '__main__':
    unittest.main()
