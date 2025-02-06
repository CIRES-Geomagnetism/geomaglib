import sys
import os
import unittest
import numpy as np
import math

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from geomaglib import util

class TestUtil(unittest.TestCase):

    def _test_single_geod_to_geoc_conv(self, lat, alt, exp_r, exp_theta):
        act_r, act_theta = util.geod_to_geoc_lat(lat,alt)
        self.assertAlmostEqual(act_r, exp_r, places=10)
        self.assertAlmostEqual(act_theta, exp_theta,places=10)

    def test_geod_to_geoc_conv(self):
        self._test_single_geod_to_geoc_conv(0,0,6378.137,0)
        self._test_single_geod_to_geoc_conv(90,0,6356.75231424518,90)
        self._test_single_geod_to_geoc_conv(-90,2,6358.75231424518,-90)
        self._test_single_geod_to_geoc_conv(22.342,2.02, 6377.08899029245,22.2070535610)
        self._test_single_geod_to_geoc_conv(45.734,-5,6362.21562444842,45.5414722008)
        self._test_single_geod_to_geoc_conv(-36.254,3.45,6374.14924310438,-36.0707588822) 
        self._test_single_geod_to_geoc_conv(-70.89,-7,6352.06156995608,-70.7705043304) 
        self._test_single_geod_to_geoc_conv(50.5,0,6365.44775048439,50.3109909949)
        self._test_single_geod_to_geoc_conv(10.2,20,6397.47181557818,10.1333459251)
    
    def test_alt_to_ellipsoid(self):

        self.assertAlmostEqual(util.alt_to_ellipsoid_height(10,40,20),10.034187,places=10)
        self.assertAlmostEqual(util.alt_to_ellipsoid_height(5,20,60),4.957280,places=10)
        self.assertAlmostEqual(util.alt_to_ellipsoid_height(1,-90,0),0.970466,places=10)
        self.assertAlmostEqual(util.alt_to_ellipsoid_height(0,90,0),0.013606,places=10)
        self.assertAlmostEqual(util.alt_to_ellipsoid_height(-2,45,-50),-1.982954,places=10)
        self.assertAlmostEqual(util.alt_to_ellipsoid_height(5,-60,-50),5.019614,places=10)
        self.assertAlmostEqual(util.alt_to_ellipsoid_height(3.5,20,270),3.489277,places=10)
        self.assertAlmostEqual(util.alt_to_ellipsoid_height(2.005,5.5,40)[0],1.985737,places=10)
        self.assertAlmostEqual(util.alt_to_ellipsoid_height(7,-20,-20),6.999373,places=10)
    
    def test_calc_dec_year(self):
        self.assertAlmostEqual(util.calc_dec_year(2024,5,20),2024.3825136612,places=10)
        self.assertAlmostEqual(util.calc_dec_year(2023,5,20),2023.3808219178,places=10)
        self.assertAlmostEqual(util.calc_dec_year(2022,2,11),2022.11232876712,places=10)
        self.assertAlmostEqual(util.calc_dec_year(2021,8,30),2021.6602739726,places=10)
        self.assertAlmostEqual(util.calc_dec_year(2020,2,11),2020.11202185792,places=10)
        self.assertAlmostEqual(util.calc_dec_year(2015,10,1),2015.74794520547,places=10)
        self.assertAlmostEqual(util.calc_dec_year(2011,12,31),2011.99726027397,places=10)
        self.assertAlmostEqual(util.calc_dec_year(1900,3,13),1900.19452054794,places=10)

    def test_calc_dec_year_arr(self):

        years = np.array([2024, 2023, 2022, 2021, 2020, 2015, 2011, 1900])
        months = np.array([5, 5, 2, 8, 2, 10, 12, 3])
        days = np.array([20, 20, 11, 30, 11, 1, 31, 13])

        dyears = util.calc_dec_year_array(years, months, days)


        for i in range(len(years)):
            self.assertAlmostEqual(dyears[i], util.calc_dec_year(years[i], months[i], days[i]), places=10)


    def test_jd2000(self):

        year = 2024
        month = 12
        day = 31
        ut = 12
        minutes = 30
        jd = util.jd2000(year, month, day, ut, minutes)

        self.assertAlmostEqual(jd, 9131.52083333, places=6)











if __name__ == '__main__':
    unittest.main()
