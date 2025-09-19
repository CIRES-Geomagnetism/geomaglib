import pathlib
from typing import Union,List,Tuple
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike

from geomaglib import mag_SPH_summation
from geomaglib.sh_vars import comp_sh_vars
from geomaglib.sh_loader import load_coef,timely_modify_magnetic_model
from geomaglib.util import cart_to_sph_deg,sph_deg_to_cart

COEF_PATH = (pathlib.Path(__file__).parents[1] / 'tests' / 'coefs' / 'IGRF14_sv.COF').resolve()

class Dipole:
    """Calculate the geomagnetic dipole magnetic field
    Approach and math inspired by https://github.com/klaundal/dipole"""
    
    def __init__(self, epoch : float):
        self.base_year = int(epoch) - int(epoch % 5)
        self.base_coef = load_coef(COEF_PATH,
                               skip_two_columns=True,
                               load_sv=True,
                               load_year=self.base_year)
        #Note that Laundal just interpolates linearly whereas geomaglib
        #actually applies the secular variation
        self.coef = timely_modify_magnetic_model(self.base_coef,epoch)
        
        #Index math int(n * (n + 1) / 2 + m)
        g10 = self.coef['g'][1]
        g11 = self.coef['g'][2]
        h11 = self.coef['h'][2]
        self.g10,self.g11,self.h11 = g10,g11,h11

        #magnetic field strength constant
        self.B0 = np.sqrt(g10**2+g11**2+h11**2)

        #Unit vector in direction of geomagnetic north pole
        nhat_x = -1.*g11/self.B0
        nhat_y = -1.*h11/self.B0
        nhat_z = -1.*g10/self.B0
        nhat_r = np.sqrt(nhat_x**2+nhat_y**2+nhat_z**2)
        self.nhat = np.array([nhat_x,nhat_y,nhat_z])

        #Location of the northern dipole pole in geocentric spherical in degrees
        _,self.pole_colat,self.pole_lon = cart_to_sph_deg(nhat_x,nhat_y,nhat_z)
        if self.pole_lon>180: #put longitude in -180,180
            self.pole_lon-=360.

    def coords(self,
               lat: Union[float,List[float]], 
                lon: Union[float,List[float]],
                inverse=False) -> Tuple[Union[float,List[float]],Union[float,List[float]]]:
        """Calculate dipole magnetic latitude and longitude for any number of
        given geocentric latitudes and longitudes (or visa-versa if inverse==True)
        Takes and returns floats or lists of floats"""
        latarr, lonarr = np.broadcast_arrays(np.asarray(lat),
                                                np.asarray(lon))
        shape = latarr.shape
        latarr, lonarr = latarr.flatten(), lonarr.flatten()

        # make rotation matrix from geo to cd
        Zcd = self.nhat
        Zgeo_x_Zcd = np.cross(np.array([0, 0, 1]), Zcd)
        Ycd = Zgeo_x_Zcd / np.linalg.norm(Zgeo_x_Zcd)
        Xcd = np.cross(Ycd, Zcd)

        #geo to centered dipole rotation matrix
        rotmat = np.vstack((Xcd, Ycd, Zcd))

        # transpose rotation matrix to get inverse operation
        if inverse: 
            rotmat = rotmat.T

        #Omit 'r' from these spherical to cartesian 
        #formulae by assuming it's 1
        
        # convert input to cartesian:
        r_x_in,r_y_in,r_z_in = sph_deg_to_cart(np.ones_like(latarr),
                                               90.-latarr,
                                               lonarr)

        
        r_vec_in = np.vstack([r_x_in,r_y_in,r_z_in])
        
        # rotate
        r_vec_out = np.dot(rotmat,r_vec_in)
        
        # cart to spherical
        r_x_out,r_y_out,r_z_out = r_vec_out[0,:],r_vec_out[1,:],r_vec_out[2,:]
        _,colat_out,lon_out = cart_to_sph_deg(r_x_out,r_y_out,r_z_out)
        lat_out = 90.-colat_out
    
        if lat_out.size==1 and lon_out.size==1:
            
            return lat_out,(lon_out - 360 if lon_out > 180. else lon_out)
        else:
            lon_out[lon_out>180.]-=360.
            return lat_out.flatten().tolist(),lon_out.flatten().tolist()
        

