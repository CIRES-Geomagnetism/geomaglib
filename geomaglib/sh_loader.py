import math
from typing import Optional

import numpy as np
import copy

from geomaglib import util

def calc_sh_degrees_to_num_elems(sh_degrees):
    """
    Gives you the number of elements from the number of 
    sphereical harmonic degree

    Parameters:
    sh_degrees (int): The number of spherical harmonic degrees

    Returns:
    int: The number of elements
    """
    return int(sh_degrees * (sh_degrees+3)/2 ) + 1

def calc_num_elems_to_sh_degrees(num_elems):
    """
    Gives you the number of spherical harmonic degrees from the total number
    of elements

    Parameters:
    num_elems (int): The number of elements

    Returns:
    int: The number of spherical harmonic degrees
    """
    return int((-3+np.sqrt(8*num_elems + 1) )/2)

def load_wmm_coef(filename: str, skip_two_columns: bool =False, load_sv: bool = True, end_degree: Optional[int] = None, load_year: Optional[float] = None) -> dict:
    """
    Takes a coefficient filename of WMM or path and gives you back a dictionary with the
    components in arrays under the keys g,h,g_sv, and h_sv

    Parameters:
    filename (string): The relative path or just the name of the coefficient file

    skip_two_columns (boolean): Sometimes the coefficient file in the data lines
    lists the spherical harmonic degree numbers in the first two columns like in
    WMM and IGRF. If this parameter is true it lets the function know to skip those
    two columns

    load_sv (boolean): The coefficient file doesn't always have sv columns like
    in HDGM crust,set this false to let the function know to not load g_sv and h_sv

    end_degree (int): If you just want to load a specific amount of spherical harmonic
    degrees from the coefficient file set this to that number. Set it to None if
    you wish to ignore

    load_year (int): The coefficient file could be seperated by year components like
    in HDGM core, set this to a year value based on what year you want to be loaded.
    Set it to None if you wish to ignore this parameter

    Returns:
    dictionary: The dictionary loaded with g, h, g_sv, h_sv arrays under those keys
    """
    coef_dict = {"g": [], "h": []}
    if load_sv:
        coef_dict["g_sv"] = []
        coef_dict["h_sv"] = []

    skip_adder = 0
    if skip_two_columns:
        skip_adder = 2


    coef_file = open(filename, 'r')
    lines = coef_file.readlines()

    header_line = True
    load = True
    footer_line_split_len = 4
    if not load_sv:
        footer_line_split_len = 2
    if load_year is not None:
        load = False
    num_lines_load = None
    load_counter = 0
    if end_degree is not None:
        num_lines_load = calc_sh_degrees_to_num_elems(end_degree)
    for line in lines:
        split = line.split()
        # This will detect the footer and we can stop loading
        if not header_line:
            if len(split) < footer_line_split_len:
                break

        if load_year is not None and split[1] == (str(load_year) + ".0"):
            load = True
            header_line = False
            coef_dict["epoch"] = float(split[1])
            continue
        if load_year is not None and split[1] == (str(load_year + 1) + ".0"):
            break
        if header_line:
            header_line = False
            if load_sv:
                coef_dict["epoch"] = float(split[0])
            
            if('/' in split[2]): #modern WMM where second value is mm/dd/yyyy
                date_string = split[2]

                
                month, day, year = map(int, date_string.split('/'))
                month = np.array([month])
                day = np.array([day])
                year = np.array([year])
                coef_dict["min_year"] = util.calc_dec_year_array(year, month, day)
                coef_dict["min_date"] = str(f"{year}-{month}-{day}")
            elif('/' in split[3]):#Old WMM with second value being decimal year
                date_string = split[3]

                month, day, year = map(int, date_string.split('/'))
                month = np.array([month])
                day = np.array([day])
                year = np.array([year])
               

                coef_dict["min_year"] = float(split[2])
                year, month, day, hour, minute = util.decimalYearToDateTime(float(split[2]))

                coef_dict["min_date"] = str(f"{year}-{month}-{day} {hour}:{minute}")
            else:
                raise ValueError(f"Header line of WMM.COF file should have form: 2025.0           WMM              11/13/2024")
            continue
        if num_lines_load is not None and load_counter >= num_lines_load:
            break
        if load:
            load_counter = load_counter + 1
            coef_dict['g'].append(float(split[0 + skip_adder]))
            coef_dict['h'].append(float(split[1 + skip_adder]))
            if load_sv:
                coef_dict['g_sv'].append(float(split[2 + skip_adder]))
                coef_dict['h_sv'].append(float(split[3 + skip_adder]))

    coef_file.close()
    if len(coef_dict["g"]) > 0 and (coef_dict["g"][0] != 0 or coef_dict['h'][0] != 0):
        coef_dict["g"].insert(0, 0)
        coef_dict["h"].insert(0, 0)
        if load_sv:
            coef_dict["g_sv"].insert(0, 0)
            coef_dict["h_sv"].insert(0, 0)
    if end_degree is not None and len(coef_dict["g"]) > num_lines_load:
        coef_dict["g"].pop()
        coef_dict["h"].pop()
        if load_sv:
            coef_dict["g_sv"].pop()
            coef_dict["h_sv"].pop()
    return coef_dict

def load_igrf_coef(filename, skip_two_columns=False, load_sv=True, end_degree=None, load_year=None):
    """
    Takes a coefficient filename or path and gives you back a dictionary with the
    components in arrays under the keys g,h,g_sv, and h_sv

    Parameters:
    filename (string): The relative path or just the name of the coefficient file
    skip_two_columns (boolean): Sometimes the coefficient file in the data lines
    lists the spherical harmonic degree numbers in the first two columns like in
    WMM and IGRF. If this parameter is true it lets the function know to skip those
    two columns
    load_sv (boolean): The coefficient file doesn't always have sv columns like
    in HDGM crust,set this false to let the function know to not load g_sv and h_sv
    end_degree (int): If you just want to load a specific amount of spherical harmonic
    degrees from the coefficient file set this to that number. Set it to None if
    you wish to ignore
    load_year (int): The coefficient file could be seperated by year components like
    in HDGM core, set this to a year value based on what year you want to be loaded.
    Set it to None if you wish to ignore this parameter

    Returns:
    dictionary: The dictionary loaded with g, h, g_sv, h_sv arrays under those keys
    """
    coef_dict = {"g": [], "h": []}
    if load_sv:
        coef_dict["g_sv"] = []
        coef_dict["h_sv"] = []

    skip_adder = 0
    if skip_two_columns:
        skip_adder = 2

    coef_file = open(filename, 'r')
    lines = coef_file.readlines()
    header_line = True
    load = True
    footer_line_split_len = 4
    if not load_sv:
        footer_line_split_len = 2
    if load_year is not None:
        load = False
    num_lines_load = None
    load_counter = 0
    if end_degree is not None:
        num_lines_load = calc_sh_degrees_to_num_elems(end_degree)
    out_degree = 0
    for line in lines:
        split = line.split()
        # This will detect the footer and we can stop loading
        if len(split) < footer_line_split_len:
            break

        if load_year is not None and (split[1] == (str(load_year) + ".0") or split[1] == (str(load_year) + ".00")):
            load = True
            out_degree = int(split[2])
            header_line = False
            coef_dict["epoch"] = int(float(split[1]))
            continue
        if load_year is not None and (split[1] == (str(load_year) + ".0") or split[1] == (str(load_year+5) + ".00")):
            break
        if header_line:

            header_line = False
            if load_sv:
                coef_dict["epoch"] = int(float(split[1]))
            continue
        if num_lines_load is not None and load_counter >= num_lines_load:
            break
        if load:
            load_counter = load_counter + 1
            coef_dict['g'].append(float(split[0 + skip_adder]))
            coef_dict['h'].append(float(split[1 + skip_adder]))
            if load_sv:
                coef_dict['g_sv'].append(float(split[2 + skip_adder]))
                coef_dict['h_sv'].append(float(split[3 + skip_adder]))

    coef_file.close()

    if len(coef_dict["g"]) > 0 and (coef_dict["g"][0] != 0 or coef_dict['h'][0] != 0):
        coef_dict["g"].insert(0, 0)
        coef_dict["h"].insert(0, 0)
        if load_sv:
            coef_dict["g_sv"].insert(0, 0)
            coef_dict["h_sv"].insert(0, 0)

    if end_degree is not None and len(coef_dict["g"]) > num_lines_load:
        coef_dict["g"].pop()
        coef_dict["h"].pop()
        if load_sv:
            coef_dict["g_sv"].pop()
            coef_dict["h_sv"].pop()
    return coef_dict, out_degree

def load_coef(filename, skip_two_columns=False, load_sv=True, end_degree=None, load_year=None):
    """
    Takes a coefficient filename or path and gives you back a dictionary with the
    components in arrays under the keys g,h,g_sv, and h_sv

    Parameters:
    filename (string): The relative path or just the name of the coefficient file
    skip_two_columns (boolean): Sometimes the coefficient file in the data lines
    lists the spherical harmonic degree numbers in the first two columns like in
    WMM and IGRF. If this parameter is true it lets the function know to skip those
    two columns
    load_sv (boolean): The coefficient file doesn't always have sv columns like
    in HDGM crust,set this false to let the function know to not load g_sv and h_sv
    end_degree (int): If you just want to load a specific amount of spherical harmonic
    degrees from the coefficient file set this to that number. Set it to None if
    you wish to ignore
    load_year (int): The coefficient file could be seperated by year components like
    in HDGM core, set this to a year value based on what year you want to be loaded.
    Set it to None if you wish to ignore this parameter

    Returns:
    dictionary: The dictionary loaded with g, h, g_sv, h_sv arrays under those keys
    """
    coef_dict = {"g": [], "h": []}
    if load_sv:
        coef_dict["g_sv"] = []
        coef_dict["h_sv"] = []

    skip_adder = 0
    if skip_two_columns:
        skip_adder = 2

    coef_file = open(filename, 'r')
    lines = coef_file.readlines()

    header_line = True
    load = True
    footer_line_split_len = 4
    if not load_sv:
        footer_line_split_len = 2
    if load_year is not None:
        load = False
    num_lines_load = None
    load_counter = 0
    if end_degree is not None:
        num_lines_load = calc_sh_degrees_to_num_elems(end_degree)
    for line in lines:
        split = line.split()
        # This will detect the footer and we can stop loading
        if len(split) < footer_line_split_len:
            break

        if load_year is not None and (split[1] == (str(load_year) + ".0") or split[1] == (str(load_year) + ".00")):
            load = True
            header_line = False
            coef_dict["epoch"] = int(float(split[1]))
            continue
        if load_year is not None and (split[1] == (str(load_year) + ".0") or split[1] == (str(load_year+5) + ".00")):
            break
        if header_line:
            header_line = False
            if load_sv:
                coef_dict["epoch"] = int(float(split[1]))
            continue
        if num_lines_load is not None and load_counter >= num_lines_load:
            break
        if load:
            load_counter = load_counter + 1
            coef_dict['g'].append(float(split[0 + skip_adder]))
            coef_dict['h'].append(float(split[1 + skip_adder]))
            if load_sv:
                coef_dict['g_sv'].append(float(split[2 + skip_adder]))
                coef_dict['h_sv'].append(float(split[3 + skip_adder]))

    coef_file.close()

    if len(coef_dict["g"]) > 0 and (coef_dict["g"][0] != 0 or coef_dict['h'][0] != 0):
        coef_dict["g"].insert(0, 0)
        coef_dict["h"].insert(0, 0)
        if load_sv:
            coef_dict["g_sv"].insert(0, 0)
            coef_dict["h_sv"].insert(0, 0)

    if end_degree is not None and len(coef_dict["g"]) > num_lines_load:
        coef_dict["g"].pop()
        coef_dict["h"].pop()
        if load_sv:
            coef_dict["g_sv"].pop()
            coef_dict["h_sv"].pop()

    return coef_dict
def timely_modify_magnetic_model(sh_dict, dec_year):
    """
    Time change the Model coefficients from the base year of the model(epoch) using secular variation coefficients.
Store the coefficients of the static model with their values advanced from epoch t0 to epoch t.
Copy the SV coefficients.  If input "t�" is the same as "t0", then this is merely a copy operation.

    Parameters:
    sh_dict (dictionary): This is the input dictionary, you would get this dictionary from using the load_coef function
    dec_year(float or int): Decimal year input for calculating the time shift
    epoch (float or int): The base year of the model

`   Returns:
    dictionary: Copy of sh_dict with the elements timely shifted
    """

    sh_dict_time = copy.deepcopy(sh_dict)
    epoch = sh_dict.get("epoch", 0)
    #If the sh_dict doesn't have secular variations just return a copy
    #of the dictionary
    if  "g_sv" not in sh_dict or "h_sv" not in sh_dict:
        return sh_dict_time
    num_elems = len(sh_dict["g"])
    sh_degrees = calc_num_elems_to_sh_degrees(num_elems)
    date_diff = dec_year - epoch
    for n in range(1, (sh_degrees+1)):
        for m in range(n+1):
            index = int(n * (n + 1) / 2 + m)
            if index < num_elems:
                sh_dict_time["g"][index] = sh_dict["g"][index] + date_diff * sh_dict["g_sv"][index]
                sh_dict_time["h"][index] = sh_dict["h"][index] + date_diff * sh_dict["h_sv"][index]

    return sh_dict_time



def timely_modify_magnetic_model_wmm(sh_dict, dec_year, max_sv: Optional[int] = None):
    """
    Time change the Model coefficients from the base year of the model(epoch) using secular variation coefficients.
Store the coefficients of the static model with their values advanced from epoch t0 to epoch t.
Copy the SV coefficients.  If input "t�" is the same as "t0", then this is merely a copy operation.

    Parameters:
    sh_dict (dictionary): This is the input dictionary, you would get this dictionary from using the load_coef function
    dec_year(float or int): Decimal year input for calculating the time shift
    epoch (float or int): The base year of the model

`   Returns:
    dictionary: Copy of sh_dict with the elements timely shifted
    """

    sh_dict_time = copy.deepcopy(sh_dict)
    epoch = sh_dict.get("epoch", 0)
    #If the sh_dict doesn't have secular variations just return a copy
    #of the dictionary
    num_elems = len(sh_dict["g"])


    if max_sv is None:
        max_sv = calc_num_elems_to_sh_degrees(num_elems)
    if  "g_sv" not in sh_dict or "h_sv" not in sh_dict:
        return sh_dict_time

    date_diff = dec_year - epoch
    for n in range(1, (max_sv+1)):
        for m in range(n+1):
            index = int(n * (n + 1) / 2 + m)
            if index < num_elems:
                
                sh_dict_time["g"][index] = sh_dict["g"][index] + date_diff * sh_dict["g_sv"][index]
                sh_dict_time["h"][index] = sh_dict["h"][index] + date_diff * sh_dict["h_sv"][index]

    return sh_dict_time


def timely_modify_magnetic_model_HDGM(sh_dict, dec_year, max_sv: Optional[int] = None):
    """
    Time change the Model coefficients from the base year of the model(epoch) using secular variation coefficients.
Store the coefficients of the static model with their values advanced from epoch t0 to epoch t.
Copy the SV coefficients.  If input "t�" is the same as "t0", then this is merely a copy operation.

    Parameters:
    sh_dict (dictionary): This is the input dictionary, you would get this dictionary from using the load_coef function
    dec_year(float or int): Decimal year input for calculating the time shift
    epoch (float or int): The base year of the model

`   Returns:
    dictionary: Copy of sh_dict with the elements timely shifted
    """

    sh_dict_time = copy.deepcopy(sh_dict)
    epoch = sh_dict.get("epoch", 0)
    #If the sh_dict doesn't have secular variations just return a copy
    #of the dictionary
    num_elems = len(sh_dict["g"])


    if max_sv is None:
        max_sv = calc_num_elems_to_sh_degrees(num_elems)
    if  "g_sv" not in sh_dict or "h_sv" not in sh_dict:
        return sh_dict_time

    date_diff = dec_year - epoch
    for n in range(1, (max_sv+1)):
        for m in range(n+1):
            index = int(n * (n + 1) / 2 + m)
            if index < num_elems:
                
                sh_dict_time["g"][index] = sh_dict["g"][index] + date_diff * sh_dict["g_sv"][index]
                sh_dict_time["h"][index] = sh_dict["h"][index] + date_diff * sh_dict["h_sv"][index]

    return sh_dict_time

