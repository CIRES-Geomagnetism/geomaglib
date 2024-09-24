import unittest
import os
import math
from geomaglib import util, sh_loader, sh_vars, magmath, legendre

from geomaglib import GeomagElements



class Test_magmath(unittest.TestCase):

    def setUp(self):

        self.Bh = [1510.0, 1910.8, 2487.8, 24377.2, 21666.6]
        self.Bx = [-575.7, 1518.0, 1555.6, 24375.3, 21556.3]
        self.By = [-1396.0, -1160.5, 1941.4, 303.2, -2183.2]
        self.Bz = [56082.3, 55671.9, 56520.5, 49691.4, -52676.0]
        self.Bf = [56102.7,  55704.7, 56575.2, 55348.7, 56957.9]
        self.Bdec = [-112.41,  -37.40, 51.30, 0.71, -5.78]
        self.Binc = [88.46, 88.03, 87.48, 63.87, -67.64]

        self.dBx = [69.5, 72.5, -48.7, -33.7, 34.4]
        self.dBy = [-9.2, 26.3, -20.4, -32.4, 9.2]
        self.dBz = [20.1, -11.3, 39.7, 101.1, -16.8]

        self.dBdec = [2.6, 1.9, 0.6, -0.1, 0.0]
        self.dBinc = [0.0, -0.0, 0.0, 0.1, 0.0]
        self.dBh = [-18.0, 41.6, -46.4, -34.1, 33.3]
        self.dBf = [19.6, -9.8, 37.7, 75.7, 28.2]

        self.tol = 1e-1

        self.curr_dir = os.path.dirname(os.path.abspath(__file__))
        self.top_dir = os.path.dirname(self.curr_dir)
        self.wmm_coeff = os.path.join(self.curr_dir, "WMM2020.cof")
    def test_get_Bh(self):

        N = len(self.Bx)


        for i in range(N):
            results = magmath.GeomagElements(self.Bx[i], self.By[i], self.Bz[i])
            h = results.get_Bh()

            self.assertAlmostEqual(h, self.Bh[i], delta=self.tol)  # add assertion here


    def test_get_dBh(self):

        N = len(self.Bx)


        for i in range(N):
            results = magmath.GeomagElements(self.Bx[i], self.By[i], self.Bz[i], self.dBx[i], self.dBy[i], self.dBz[i])

            dh = results.get_dBh()

            self.assertAlmostEqual(dh, self.dBh[i], delta=self.tol)  # add assertion here
    def test_get_Bf(self):

        N = len(self.Bf)

        for i in range(N):
            results = magmath.GeomagElements(self.Bx[i], self.By[i], self.Bz[i])
            f = results.get_Bf()
            self.assertAlmostEqual(f, self.dBf[i], delta=self.tol)  # add assertion here

    def test_get_dBf(self):

        N = len(self.Bf)

        for i in range(N):
            results = magmath.GeomagElements(self.Bx[i], self.By[i], self.Bz[i], self.dBx[i], self.dBy[i], self.dBz[i])
            df = results.get_dBf()
            print(math.fabs(df - self.dBf[i]))
            self.assertAlmostEqual(df, self.dBf[i], delta=self.tol)  # add assertion here

    def test_get_Bdec(self):

        N = len(self.Bdec)

        for i in range(N):
            results = magmath.GeomagElements(self.Bx[i], self.By[i], self.Bz[i])
            dec = results.get_Bdec()

            self.assertAlmostEqual(round(dec, 2), self.Bdec[i], delta=0.01)  # add assertion here

    def test_get_dBdec(self):

        N = len(self.Bdec)

        for i in range(N):
            results = magmath.GeomagElements(self.Bx[i], self.By[i], self.Bz[i], self.dBx[i], self.dBy[i], self.dBz[i])
            ddec = results.get_dBdec()

            self.assertAlmostEqual(round(ddec, 1), self.dBdec[i], delta=0.01)  # add assertion here

    def test_get_Binc(self):

        N = len(self.Bdec)

        for i in range(N):
            results = magmath.GeomagElements(self.Bx[i], self.By[i], self.Bz[i])
            dinc = results.get_Binc()

            self.assertAlmostEqual(round(dinc, 2), self.Binc[i], delta=0.01)  # add assertion here

    def test_get_dBinc(self):

        N = len(self.Bdec)

        for i in range(N):
            results = magmath.GeomagElements(self.Bx[i], self.By[i], self.Bz[i], self.dBx[i], self.dBy[i], self.dBz[i])
            dinc = results.get_dBinc()

            self.assertAlmostEqual(round(dinc, 1), self.dBinc[i], delta=0.01)  # add assertion here

    def test_get_all_base(self):

        N = len(self.Bx)

        for i in range(N):
            results = GeomagElements(self.Bx[i], self.By[i], self.Bz[i])
            map = results.get_all_base()

            self.assertAlmostEqual(map["x"], self.Bx[i], delta=self.tol)
            self.assertAlmostEqual(map["y"], self.By[i], delta=self.tol)
            self.assertAlmostEqual(map["z"], self.Bz[i], delta=self.tol)
            self.assertAlmostEqual(map["h"], self.Bh[i], delta=self.tol)
            self.assertAlmostEqual(map["f"], self.Bf[i], delta=self.tol)
            self.assertAlmostEqual(map["dec"], self.Bdec[i], delta=0.01)
            self.assertAlmostEqual(map["inc"], self.Binc[i], delta=0.01)

    def test_get_all(self):

        for i in range(len(self.dBx)):
            results = GeomagElements(self.Bx[i], self.By[i], self.Bz[i], self.dBx[i], self.dBy[i], self.dBz[i])

            map = results.get_all()


            self.assertAlmostEqual(round(map["dh"], 1), self.dBh[i], delta=self.tol)

            self.assertAlmostEqual(round(map["df"], 1), self.dBf[i], delta=0.01)

            self.assertAlmostEqual(round(map["ddec"], 1), self.dBdec[i], delta=0.01)
            self.assertAlmostEqual(round(map["dinc"], 1), self.dBinc[i], delta=0.01)



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

        print(coef_dict.keys())

        # alt = util.alt_to_ellipsoid_height(alt, lat, lon)
        for i in range(N):
            r, theta = util.geod_to_geoc_lat(lats[i], alts[i])
            sph_dict = sh_vars.comp_sh_vars(lons[i], r, theta, nmax)
            cotheta = 90 - theta
            colats = [cotheta]

            Leg = legendre.Flattened_Chaos_Legendre1(nmax, colats)

            Bt, Bp, Br = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict["g"], timly_coef_dict["h"], Leg, theta)
            x, y, z = magmath.rotate_magvec(Bt, Bp, Br, theta, lats[i])

            dBt, dBp, dBr = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict["g_sv"], timly_coef_dict["h_sv"], Leg,
                                                   theta)
            dx, dy, dz = magmath.rotate_magvec(dBt, dBp, dBr, theta, lats[i])

            self.assertAlmostEqual(x, self.Bx[i], delta=self.tol)
            self.assertAlmostEqual(y, self.By[i], delta=self.tol)
            self.assertAlmostEqual(z, self.Bz[i], delta=self.tol)

            self.assertAlmostEqual(dx, self.dBx[i], delta=self.tol)
            self.assertAlmostEqual(dy, self.dBy[i], delta=self.tol)
            self.assertAlmostEqual(dz, self.dBz[i], delta=self.tol)

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
        Bt, Bp, Br = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict["g"], timly_coef_dict["h"], Leg, theta_nopole)
        print(f"Bt: {Bt}, Br: {Br}, Bp: {Bp}")
        Bpole_t, Bpole_p, Bpole_r = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict["g_sv"], timly_coef_dict["h_sv"], Leg, theta)
        print(f"Bt: {Bpole_t}, Br: {Bpole_r}, Bp: {Bpole_p}")

        #self.assertEqual(Bt, Bpole_t)  # add assertion here
        self.assertEqual(Bp, Bpole_p)
        #self.assertEqual(Br, Bpole_r)





if __name__ == '__main__':
    unittest.main()
