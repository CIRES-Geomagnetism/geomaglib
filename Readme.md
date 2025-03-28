# geomaglib

![PyPI - Version](https://img.shields.io/pypi/v/geomaglib)
![PyPI - License](https://img.shields.io/pypi/l/geomaglib)
![PyPI - Downloads](https://img.shields.io/pypi/dm/geomaglib)
[![PyPI Downloads](https://static.pepy.tech/badge/geomaglib)](https://pepy.tech/projects/geomaglib)

geomaglib Python API includes the core functions for building geomag model like WMM, WMMHR, HDGM-crust, HDGM-core and IGRF.
The API is developed by by the Cooperative Institute For Research in Environmental Sciences (CIRES), University of Colorado. 

## Installation

### Install the geomaglib deployed to PyPI
Install the geomaglib from PyPI
```commandline
pip install geomaglib
```

### Install the Public Git repository to your virtual environment

#### 1. Install specific branch
```commandline
pip install git+https://github.com/CIRES-Geomagnetism/geomaglib.git@<branch_name>
```
#### 2. Install specific tag
```commandline
pip install git+https://github.com/CIRES-Geomagnetism/geomaglib.git@<tag_number>
```
## API Reference

### geomaglib.legendre

- Compute the Legendre polynomials

<details>
<summary>See the detail of functions</summary>
<p><b>Flattened_Chaos_Legendre1()</b>: Outputs a 2 dimensional numpy array which contains the associated
    legendre polynomials (Pnm) and the respective derivatives (dPnm).</p>
<p><b>get_index()</b>: Get the index of Flattened_Chaos_Legendre1 </p>
</details>

### geomaglib.magmath

- Compute the Geomagnetic elements in geodetic: Bx, By, Bz, Bh(horizontal), Bf(total intensity), declination, inclination and the delta Bx, By, Bz, Bh(horizontal), Bf(total intensity), declination, inclination 
- Radius and degree conversion
- Convert magnetic vector from spherical to geodetic
<details>
<summary>See the detail of functions</summary>
<p>rad2deg(): Convert radius to degree</p>
<p>deg2rad(): Convert degree to radius</p>
<p>calc_Bp_Pole(): Calculate the B_phi magnetic elements at pole</p>
<p>mag_SPH_summation(): Compute the sperical magnetic vectors (B_theta, B_phi, B_radius)</p>
<p>rotate_magvec():   Convert magnetic vector from spherical to geodetic</p>
<p>GeomagElements(): The Python class for computing magnetic elements in geodetic coordinates </p>
    <ul>
    <li>Get all of the magnetic elements: get_all()</li>
    <li>Get the specific magnetic element: 
    <ul>
    <li>get_Bx()</li>
    <li>get_By()</li>
    <li>get_Bz()</li>
    <li>get_Bh()</li>
    <li>get_Bf()</li>
    <li>get_Bdec()</li>
    <li>get_Binc()</li>
    <li>get_dBx()</li>
    <li>get_dBy()</li>
    <li>get_dBz()</li>
    <li>get_dBh()</li>
    <li>get_dBf()</li>
    <li>get_dBdec()</li>
    <li>get_dBinc()</li>
    </ul>
    </li>
    </ul>

<p></p>
</details>




### geomaglib.util

- Compute decimal year
- Converting the coordinates based on gedetic to geocentric
- Converting the mean sea level altitude to ellipsoid height

<details>
<summary>See the detail of functions</summary>
<p> <b>geod_to_geoc_lat()</b>: Conversion from geodetic X,Z components to geocentric B_r, B_theta </p>
<p> <b>alt_to_ellipsoid_height()</b>: This function converts altitude MSL in kilometers to ellipsoid height in
    kilometers</p>
<p> <b>calc_dec_year()</b>: Takes year, month, and day and calculates the decimal year from those inputs </p>
<p> <b>calc_dec_year_array()</b>:Takes the array of year, month, and day and outputs the decimal year from those inputs </p>
<p> <b>decimalYearToDateTime()</b>: Covert the decimal year to year, month, day, hour and minute </p>
</details>

### geomaglib.sh_vars

- Computes Spherical Variabale
<details>
<summary>See the detail of functions</summary>

<p> <b>comp_sh_vars()</b>: Computes spherical variabale and put the results in a dictionary
    variable computed are (a/r) ^ (n+2), cos_m(lon), and sin_m(lon)</p>

</details>


### geomaglib.sh_loader

- Estimate the total elements from end degree ofr the model coefficients
- Funtions for loading WMM2025.COF, IGRF14.COF
- Compute time change the Model coefficients from the base year of the model(epoch) using secular variation coefficients.

<details>
<summary>See the detail of functions</summary>

<p> <b>calc_sh_degrees_to_num_elems()</b>: Gives you the number of elements from the number of 
    sphereical harmonic degree
</p>
<p> <b>calc_num_elems_to_sh_degrees()</b>: Gives you the number of spherical harmonic degrees from the total number
    of elements
</p>
<p> <b>load_wmm_coef()</b></p>
<p> <b>load_igrf_coef()</b></p>
<p> <b>timely_modify_magnetic_model_wmm()</b></p>
<p> <b>timely_modify_magnetic_model_HDGM()</b></p>

</details>

## Contacts and contributing to geomaglib:
If you have any questions, please email `geomag.models@noaa.gov`, submit issue or pull request at [https://github.com/CIRES-Geomagnetism/geomaglib](https://github.com/CIRES-Geomagnetism/geomaglib).