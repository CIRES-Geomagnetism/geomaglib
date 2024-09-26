import sys
import os
import unittest

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from geomaglib import util
from geomaglib import sh_vars

class TestUtil(unittest.TestCase):

    def _test_act_dict_to_exp_dict(self, act_dict, exp_dict):
        for key in exp_dict:
            act_arr = act_dict.get(key)
            exp_arr = exp_dict.get(key)
            for i in range(len(exp_arr)):
                self.assertAlmostEqual(act_arr[i], exp_arr[i], places=10)

    def test_geod_to_geoc_conv(self):
        r, theta = util.geod_to_geoc_lat(20,10) 
        sh_vars_dict = sh_vars.comp_sh_vars(50,r,theta, 5)
        exp_dict = {"relative_radius_power":[0.9954781233,0.9932248585,0.9909766939,0.9887336181,0.9864956195,0.9842626866],"cos_mlon":[1.0000000000,0.6427876097,-0.1736481777,-0.8660254038,-0.9396926208,-0.3420201433],"sin_mlon":[0.0000000000,0.7660444431,0.9848077530,0.5000000000,-0.3420201433,-0.9396926208]}
        self._test_act_dict_to_exp_dict(sh_vars_dict, exp_dict)

        r, theta = util.geod_to_geoc_lat(90,5)
        sh_vars_dict = sh_vars.comp_sh_vars(0,r,theta, 7)
        exp_dict = {"relative_radius_power":[1.0029723575,1.0044618477,1.0059535499,1.0074474673,1.0089436034,1.0104419614,1.0119425445,1.0134453561],"cos_mlon":[1.0000000000,1.0000000000,1.0000000000,1.0000000000,1.0000000000,1.0000000000,1.0000000000,1.0000000000],"sin_mlon":[0.0000000000,0.0000000000,0.0000000000,0.0000000000,0.0000000000,0.0000000000,0.0000000000,0.0000000000]}
        self._test_act_dict_to_exp_dict(sh_vars_dict, exp_dict)

        r, theta = util.geod_to_geoc_lat(-25,0)
        sh_vars_dict = sh_vars.comp_sh_vars(45,r,theta, 10)
        exp_dict = {"relative_radius_power":[0.9990138422,0.9985211281,0.9980286570,0.9975364288,0.9970444433,0.9965527005,0.9960612002,0.9955699423,0.9950789267,0.9945881533,0.9940976219],"cos_mlon":[1.0000000000,0.7071067812,0.0000000000,-0.7071067812,-1.0000000000,-0.7071067812,-0.0000000000,0.7071067812,1.0000000000,0.7071067812,0.0000000000],"sin_mlon":[0.0000000000,0.7071067812,1.0000000000,0.7071067812,0.0000000000,-0.7071067812,-1.0000000000,-0.7071067812,-0.0000000000,0.7071067812,1.0000000000]}
        self._test_act_dict_to_exp_dict(sh_vars_dict, exp_dict)

        r, theta = util.geod_to_geoc_lat(10.5,-2.6)
        sh_vars_dict = sh_vars.comp_sh_vars(-70,r,theta, 4)
        exp_dict = {"relative_radius_power":[0.9988606674,0.9982914879,0.9977226328,0.9971541019,0.9965858949],"cos_mlon":[1.0000000000,0.3420201433,-0.7660444431,-0.8660254038,0.1736481777],"sin_mlon":[0.0000000000,-0.9396926208,-0.6427876097,0.5000000000,0.9848077530]}
        self._test_act_dict_to_exp_dict(sh_vars_dict, exp_dict)



if __name__ == '__main__':
    unittest.main()
