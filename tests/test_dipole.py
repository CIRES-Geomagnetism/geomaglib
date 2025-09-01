import unittest
from numpy.testing import assert_allclose
from geomaglib.dipole import Dipole

class TestDipole(unittest.TestCase):
    
    def test_coord_accuracy(self):
        """Tests results against reference results from
        https://github.com/klaundal/dipole.
        We can only get this to .01 accuracy because 
        Laundal dipole does not use secular variation, but rather
        interpolates coefficients linearly"""
        dp = Dipole(2025.5)
        gclats = [-89,-45.,0.,45.,89.]
        glons = [-179.,-90.,0.,90.,179.]
        exp_mlats = [-81.03790935, -36.16814416, 2.70974606, 36.16814416, 80.45017686]
        exp_mlons = [353.8245409,  344.97297991,  72.99018428, 164.97297991, 185.73468103]
        exp_mlons180 = [mlon-360. if mlon>180. else mlon for mlon in exp_mlons]

        #unittest.TestCase.assertAlmostEqual doesn't work on arrays or lists
        with self.subTest('geo to mag'):
            test_mlats,test_mlons = dp.coords(gclats,glons)
            assert_allclose(exp_mlats,test_mlats,atol=.01,rtol=0)
            assert_allclose(exp_mlons180,test_mlons,atol=.01,rtol=0)
        with self.subTest('mag to geo'):
            test_gclats,test_glons = dp.coords(exp_mlats,exp_mlons180,inverse=True)
            assert_allclose(gclats,test_gclats,atol=.01,rtol=0)
            assert_allclose(glons,test_glons,atol=.01,rtol=0)
    
    def _coord_round_trip(self,epoch,gclati,gloni):
        dp = Dipole(epoch)
        mlat,mlon = dp.coords(gclati,gloni)
        gclato,glono = dp.coords(mlat,mlon,inverse=True)
        return gclato,glono
    
    def test_coord_type_handling(self):
        lat_tol=1.0e-10
        lon_tol=1.0e-10

        with self.subTest('floats in and out'):
            gclati = 45.
            gloni = 90.
            gclato,glono = self._coord_round_trip(2000.,gclati,gloni)
            self.assertAlmostEqual(gclati,gclato,places=10)
            self.assertAlmostEqual(gloni,glono,places=10)

        with self.subTest('lists in and out'):
            gclati = [-45.,45.]
            gloni = [-90.,90.]
            gclato,glono = self._coord_round_trip(1995.,gclati,gloni)
            assert_allclose(gclato,gclati,atol=lat_tol,rtol=0.)
            assert_allclose(glono,gloni,atol=lon_tol,rtol=0.)

        with self.subTest('list and float in and lists out'):
            gclati = [-45.,45.]
            gloni = 90.
            gclato,glono = self._coord_round_trip(1980.,gclati,gloni)
            assert_allclose(gclato,gclati,atol=lat_tol,rtol=0.)
            assert_allclose(glono,[gloni,gloni],atol=lon_tol,rtol=0.)
        