import unittest
import os
import math
import numpy as np
from geomaglib import util, sh_loader, sh_vars, magmath, legendre

from geomaglib import GeomagElements



class Test_magmath(unittest.TestCase):

    def setUp(self):


        self.tol = 1e-1

        self.curr_dir = os.path.dirname(os.path.abspath(__file__))
        self.top_dir = os.path.dirname(self.curr_dir)
        self.igrf_coeff = os.path.join(self.curr_dir, "coefs", "IGRF14_sv.COF")
        self.igrf_testval = os.path.join(self.curr_dir, "IGRF14_TEST_VALUES.txt")

        self.setup_IGRF_TestValues(self.igrf_testval)

    def setup_IGRF_TestValues(self, file_name):

        self.dyears, self.alts, self.lats, self.lons = [], [], [], []
        self.Bh, self.Bf, self.Bx, self.By, self.Bz, self.Bdec, self.Binc = [], [], [], [], [], [], []
        self.dBh, self.dBf, self.dBx, self.dBy, self.dBz, self.dBdec, self.dBinc  = [], [], [], [], [], [], []

        with open(file_name, "r") as file:
            for line in file:
                vals = line.split()
                if vals[0] == "#":
                    continue
                else:
                    for i in range(len(vals)):
                        vals[i] = float(vals[i])
                    dyear, alt, lat, lon = vals[0], vals[1], vals[2], vals[3]
                    dec, inc, h, f, x, y, z = vals[4], vals[5], vals[6], vals[7], vals[8], vals[9], vals[10]
                    ddec, dinc, dh, df, dx, dy, dz = vals[11], vals[12], vals[13], vals[14], vals[15], vals[16], vals[17]

                    self.dyears.append(dyear)
                    self.alts.append(alt)
                    self.lats.append(lat)
                    self.lons.append(lon)

                    self.Bdec.append(dec)
                    self.Binc.append(inc)
                    self.Bx.append(x)
                    self.By.append(y)
                    self.Bz.append(z)
                    self.Bh.append(h)
                    self.Bf.append(f)

                    self.dBdec.append(ddec)
                    self.dBinc.append(dinc)
                    self.dBx.append(dx)
                    self.dBy.append(dy)
                    self.dBz.append(dz)
                    self.dBh.append(dh)
                    self.dBf.append(df)

    def test_get_Bh(self):

        N = len(self.Bx)

        results = magmath.GeomagElements(self.Bx, self.By, self.Bz)
        h = results.get_Bh()

        self.assertEqual(type(h), np.ndarray)


        for i in range(N):

            self.assertAlmostEqual(h[i], self.Bh[i], delta=self.tol)  # add assertion here


    def test_get_dBh(self):

        N = len(self.Bx)
        results = magmath.GeomagElements(self.Bx, self.By, self.Bz, self.dBx, self.dBy, self.dBz)

        dh = results.get_dBh()
        self.assertEqual(type(dh), np.ndarray)
        for i in range(N):
            self.assertAlmostEqual(dh[i], self.dBh[i], delta=self.tol)  # add assertion here
    def test_get_Bf(self):

        N = len(self.Bf)
        results = magmath.GeomagElements(self.Bx, self.By, self.Bz)
        f = results.get_Bf()

        self.assertEqual(type(f), np.ndarray)
        for i in range(N):

            self.assertAlmostEqual(f[i], self.Bf[i], delta=self.tol)  # add assertion here

    def test_get_dBf(self):

        N = len(self.Bf)
        results = magmath.GeomagElements(self.Bx, self.By, self.Bz, self.dBx, self.dBy, self.dBz)


        df = results.get_dBf()
        self.assertEqual(type(df), np.ndarray)
        for i in range(N):

            self.assertAlmostEqual(df[i], self.dBf[i], delta=self.tol)  # add assertion here

    def test_get_Bdec(self):

        N = len(self.Bdec)
        results = magmath.GeomagElements(self.Bx, self.By, self.Bz)
        dec = results.get_Bdec()

        self.assertEqual(type(dec), np.ndarray)

        for i in range(N):


            self.assertAlmostEqual(round(dec[i], 2), self.Bdec[i], delta=0.01)  # add assertion here

    def test_get_dBdec(self):

        N = len(self.Bdec)
        results = magmath.GeomagElements(self.Bx, self.By, self.Bz, self.dBx, self.dBy, self.dBz)
        ddec = results.get_dBdec()

        self.assertEqual(type(ddec), np.ndarray)
        for i in range(N):
            self.assertAlmostEqual(ddec[i]*60, self.dBdec[i], delta=0.1)  # add assertion here

    def test_get_Binc(self):

        N = len(self.Bdec)
        results = magmath.GeomagElements(self.Bx, self.By, self.Bz)
        inc = results.get_Binc()

        self.assertEqual(type(inc), np.ndarray)

        for i in range(N):

            self.assertAlmostEqual(round(inc[i], 2), self.Binc[i], delta=self.tol)  # add assertion here

    def test_get_dBinc(self):

        N = len(self.Bdec)
        results = magmath.GeomagElements(self.Bx, self.By, self.Bz, self.dBx, self.dBy, self.dBz)
        dinc = results.get_dBinc()

        self.assertEqual(type(dinc), np.ndarray)

        for i in range(N):


            self.assertAlmostEqual(round(dinc[i]*60, 1), self.dBinc[i], delta=self.tol)  # add assertion here

    def test_get_all_base(self):

        N = len(self.Bx)
        results = GeomagElements(self.Bx, self.By, self.Bz)
        map = results.get_all_base()

        for key in map:
            self.assertEqual(type(map[key]), np.ndarray)

        for i in range(N):

            self.assertAlmostEqual(map["x"][i], self.Bx[i], delta=self.tol)
            self.assertAlmostEqual(map["y"][i], self.By[i], delta=self.tol)
            self.assertAlmostEqual(map["z"][i], self.Bz[i], delta=self.tol)
            self.assertAlmostEqual(map["h"][i], self.Bh[i], delta=self.tol)
            self.assertAlmostEqual(map["f"][i], self.Bf[i], delta=self.tol)
            self.assertAlmostEqual(map["dec"][i], self.Bdec[i], delta=0.01)
            self.assertAlmostEqual(map["inc"][i], self.Binc[i], delta=0.01)

    def test_get_all(self):

        Bx = np.array(self.Bx)
        By = np.array(self.By)
        Bz = np.array(self.Bz)

        dBx = np.array(self.dBx)
        dBy = np.array(self.dBy)
        dBz = np.array(self.dBz)

        results = GeomagElements(Bx, By, Bz, dBx, dBy, dBz)
        map = results.get_all()

        for key in map:
            self.assertEqual(type(map[key]), np.ndarray)

        for i in range(len(self.dBx)):
            self.assertAlmostEqual(map["dh"][i], self.dBh[i], delta=self.tol)

            df = results.get_dBf()
            self.assertAlmostEqual(df[i], map["df"][i], delta=1e-6)

            self.assertAlmostEqual(map["df"][i], self.dBf[i], delta=self.tol)


            self.assertAlmostEqual(map["ddec"][i]*60, self.dBdec[i], delta=0.1)
            self.assertAlmostEqual(map["dinc"][i]*60, self.dBinc[i], delta=0.1)


    def test_arrinputs_base(self):


        Bx = np.array(self.Bx)
        By = np.array(self.By)
        Bz = np.array(self.Bz)
        results = GeomagElements(Bx, By, Bz)

        H = results.get_Bh()
        F = results.get_Bf()
        dec = results.get_Bdec()
        inc = results.get_Binc()



        self.assertEqual(len(H), len(Bx))
        self.assertEqual(len(F), len(Bx))
        self.assertEqual(len(dec), len(Bx))
        self.assertEqual(len(inc), len(Bx))



    def test_arrinputs_sv(self):


        Bx = np.array(self.Bx)
        By = np.array(self.By)
        Bz = np.array(self.Bz)

        dBx = np.array(self.dBx)
        dBy = np.array(self.dBy)
        dBz = np.array(self.dBz)
        results = GeomagElements(Bx, By, Bz, dBx, dBy, dBz)

        dH = results.get_dBh()
        dF = results.get_dBf()
        ddec = results.get_dBdec()
        dinc = results.get_dBinc()

        self.assertEqual(len(dH), len(Bx))
        self.assertEqual(len(dF), len(Bx))
        self.assertEqual(len(ddec), len(Bx))
        self.assertEqual(len(dinc), len(Bx))



    def test_mag_SPH_summation(self):

        lats = np.array(self.lats)
        lons = np.array(self.lons)
        alts = np.array(self.alts)
        N = len(lats)




        for i in range(N):

            epoch = 1900 + 5*(float(self.dyears[i]) - 1900)//5



            coef_dict = sh_loader.load_coef(self.igrf_coeff, skip_two_columns=True, load_sv=True, load_year=int(epoch))



            timly_coef_dict = sh_loader.timely_modify_magnetic_model(coef_dict, float(self.dyears[i]))

            nmax = sh_loader.calc_num_elems_to_sh_degrees(len(timly_coef_dict["g"]))



            if float(self.dyears[i]) < 2000:
                nmax = 10
            else:
                nmax = 13
            # alt = util.alt_to_ellipsoid_height(alt, lat, lon)

            r, theta = util.geod_to_geoc_lat(lats[i], alts[i])
            r = np.array([r])
            #theta = np.array([theta])
            theta = [theta]
            sph_dict = sh_vars.comp_sh_vars(lons[i], r, theta, nmax)
            cotheta = 90 - theta[0]

            Leg = legendre.Flattened_Chaos_Legendre1(nmax, cotheta)

            Bt, Bp, Br = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict["g"], timly_coef_dict["h"], Leg,
                                                   theta)

            dBt, dBp, dBr = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict["g_sv"], timly_coef_dict["h_sv"],
                                                      Leg,
                                                      theta)



            x, y, z = magmath.rotate_magvec(Bt, Bp, Br, theta, lats[i])
            dx, dy, dz = magmath.rotate_magvec(dBt, dBp, dBr, theta, lats[i])
            results = GeomagElements(x, y, z, dx, dy, dz)
            h = results.get_Bh()
            f = results.get_Bf()
            dec = results.get_Bdec()
            inc = results.get_Binc()

            dh = results.get_dBh()
            df = results.get_dBf()
            ddec = results.get_dBdec()
            dinc = results.get_dBinc()



            self.assertAlmostEqual(x, self.Bx[i], delta=self.tol)
            self.assertAlmostEqual(y, self.By[i], delta=self.tol)
            self.assertAlmostEqual(z, self.Bz[i], delta=self.tol)
            self.assertAlmostEqual(h, self.Bh[i], delta=self.tol)
            self.assertAlmostEqual(f, self.Bf[i], delta=self.tol)
            self.assertAlmostEqual(dec, self.Bdec[i], delta=self.tol)
            self.assertAlmostEqual(inc, self.Binc[i], delta=self.tol)

            self.assertAlmostEqual(dx, self.dBx[i], delta=self.tol)
            self.assertAlmostEqual(dy, self.dBy[i], delta=self.tol)
            self.assertAlmostEqual(dz, self.dBz[i], delta=self.tol)
            self.assertAlmostEqual(dh, self.dBh[i], delta=self.tol)
            self.assertAlmostEqual(df, self.dBf[i], delta=self.tol)
            self.assertAlmostEqual(ddec*60, self.dBdec[i], delta=self.tol)
            self.assertAlmostEqual(dinc*60, self.dBinc[i], delta=self.tol)

    def test_calc_Bp_Pole(self):
        lat = 90
        lon = 138
        alt = 77

        dec_year = 2024.5

        # load g, h
        coef_dict = sh_loader.load_coef(self.igrf_coeff, skip_two_columns=True, load_year=2020)
        timly_coef_dict = sh_loader.timely_modify_magnetic_model(coef_dict, dec_year)
        nmax = sh_loader.calc_num_elems_to_sh_degrees(len(coef_dict["g"]))

        # Assume alt already in WGS, otherwise, use util.alt_to_ellipsoid_height(alt, lat, lon) to transform

        r, theta = util.geod_to_geoc_lat(lat, alt)
        sph_dict = sh_vars.comp_sh_vars(lon, r, theta, nmax)
        cotheta = 90 - theta
        colats = cotheta

        Leg = legendre.Flattened_Chaos_Legendre1(nmax, colats)


        theta_nopole = theta + 1e-9
        Bt, Bp, Br = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict["g"], timly_coef_dict["h"], Leg, np.array([theta_nopole]))


        Bpole_t, Bpole_p, Bpole_r = magmath.mag_SPH_summation(nmax, sph_dict, timly_coef_dict["g"], timly_coef_dict["h"], Leg, np.array([theta]))


        self.assertEqual(Bt, Bpole_t)  # add assertion here
        self.assertEqual(Bp, Bpole_p)
        self.assertEqual(Br, Bpole_r)





if __name__ == '__main__':
    unittest.main()
