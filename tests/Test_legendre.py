import unittest
import numpy as np
import os
import math 

from lib.Legendre_and_SHA import Leg_SHA_for_import
from load import sh_vars, util
class Test_legendre(unittest.TestCase):

    def setUp(self)->None:
        
        self.curr_dir = os.path.dirname(os.path.abspath(__file__))
        self.top_dir = os.path.dirname(self.curr_dir)    
        self.gsl_file = os.path.join(self.curr_dir, "gsl_results.csv")

    def read_gsl_results(self):
        
        Plm = [0]
        dPlm = [0]
        with open(self.gsl_file, "r") as file:

            for line in file:
                rows = line.split(",")
                
                Plm.append(float(rows[0]))
                dPlm.append(float(rows[1]))
        
        return Plm, dPlm

    def flatten(self, nested_list):
        for item in nested_list:
            if isinstance(item, list):
                yield from flatten(item)
            else:
                yield item

    def compare_gsl_flattenL(self, fP, fdP, gP, gdP, nmax, out_file):


        file = open(out_file, "w")
        file.write("degree,order,flatten_Plm,gsl_Plm,flatten_dPlm,gsl_dPlm\n")
        fidx = 1
        for m in range(nmax+1):
            for n in range(m, nmax+1):
                if n == 0:
                    continue

                gidx = int(n * (n + 1) / 2 + m)
                

                diff_P = math.fabs(fP[fidx] - gP[gidx])
                diff_dP = math.fabs(fdP[fidx] - gdP[gidx])

                if diff_P > 1e-6 or diff_dP > 1e-6: 
                    file.write(f"{n},{m},{round(fP[fidx],6)},{gP[gidx]},{round(fdP[fidx],6)},{gdP[gidx]}\n")        
                
                fidx += 1

        file.close()
        
    def compare_gsl_manoj(self, fP, fdP, gP, gdP, nmax, out_file):


        file = open(out_file, "w")
        file.write("degree,order,alf_Plm,gsl_Plm,alf_dPlm,gsl_dPlm\n")
        
        for n in range(1, nmax+1):
            for m in range(n+1):
                if n == 0:
                    continue

                gidx = int(n * (n + 1) / 2 + m)

                gdP[gidx] = -gdP[gidx]
                

                diff_P = math.fabs(fP[gidx] - gP[gidx])
                diff_dP = math.fabs(fdP[gidx] - gdP[gidx])

                if diff_P > 1e-6 or diff_dP > 1e-6: 
                    file.write(f"{n},{m},{round(fP[gidx],6)},{gP[gidx]},{round(fdP[gidx],6)},{gdP[gidx]}\n")        
                
            

        file.close()
         
    def test_Flattened_Chaos_Legendre1(self):
    
        lat = 20.1
        #colat = 90 - lat
        alt = 100
        r, theta = util.geod_to_geoc_lat(lat, alt)
        print(f"theta: {theta}")
        cotheta = 90 - theta
        colats = [cotheta]
        res_file = os.path.join(self.curr_dir, "compare_results_flattten.csv")
        
        nmax = 12
        Leg = Leg_SHA_for_import.Flattened_Chaos_Legendre1(nmax, colats)


        Plm = np.array(Leg[0])
        dPlm = np.array(Leg[1])

        fPlm = Plm.flatten()
        fdPlm = dPlm.flatten()

        print(fPlm)
    
        gPlm, gdPlm = self.read_gsl_results()
        
        self.assertEqual(len(fPlm),len(gPlm))        
        self.compare_gsl_flattenL(fPlm, fdPlm, gPlm, gdPlm, nmax, res_file)


    def test_legendre_manoj_low(self):

        lat = 20.1
        colat = 90 - lat 
        colats = [colat]
        alt = 100
        res_file = os.path.join(self.curr_dir, "compare_results_alf_low.csv")
        
        nmax = 12
        r, theta = util.geod_to_geoc_lat(lat, alt)
        Leg = Leg_SHA_for_import.legendre_manoj(theta, nmax)
        
        Plm = np.array(Leg[0])
        dPlm = np.array(Leg[1])

        fPlm = Plm.flatten()
        fdPlm = dPlm.flatten()
    
        gPlm, gdPlm = self.read_gsl_results()
        
        self.assertEqual(len(fPlm),len(gPlm))        
        self.compare_gsl_manoj(fPlm, fdPlm, gPlm, gdPlm, nmax, res_file)

    def test_legendre_manoj_high(self):

        lat = 20.1
        alt = 100

        r, theta = util.geod_to_geoc_lat(lat, alt)



        print(f"geoc_lat: {theta}")


        res_file = os.path.join(self.curr_dir, "compare_results_alf_high.csv")

        nmax = 12
        mPlm, mdPlm = Leg_SHA_for_import.PcupHigh(theta, nmax)


        gPlm, gdPlm = self.read_gsl_results()


        self.compare_gsl_manoj(mPlm, mdPlm, gPlm, gdPlm, nmax, res_file)
    def test_compare_flatten_manoj(self):

        lats = np.linspace(57, 61, 10)


        nmax = 790

        for lat in lats:

            colat = 90.0 - float(lat)
            colats = [colat]
            fLeg = Leg_SHA_for_import.Flattened_Chaos_Legendre1(nmax, colats)
            mPlm, mdPlm = Leg_SHA_for_import.PcupHigh(lat, nmax)

            fPlm = np.array(fLeg[0]).flatten()
            #mPlm = np.array(mLeg[0]).flatten()

            fdPlm = np.array(fLeg[1]).flatten()
            #mdPlm = np.array(mLeg[1]).flatten()

            tol = 1e-12

            fidx = 1
            for m in range(nmax+1):
                for n in range(m, nmax+1):
                    if n == 0:
                        continue

                    gidx = int(n * (n + 1) / 2 + m)

                    self.assertAlmostEqual(fPlm[fidx], mPlm[gidx], delta=tol)
                    self.assertAlmostEqual(fdPlm[fidx], -mdPlm[gidx], delta=tol)

                    fidx += 1






if __name__ == "__main__":
    unittest.main()
