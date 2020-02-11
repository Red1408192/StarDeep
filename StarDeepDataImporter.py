# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 16:42:05 2020

@author: Red
"""
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#Read Element and Isotopes in consideration from the csv file:
elem_iso = pd.read_csv('elements_and_isotopes.csv')

# Read the conversion file, and create for each metallicity percentage-element 
# combination a cubic interpolation function from the masses and values
conversion_table = pd.read_csv('conversion.csv')

def get_interpolation(group):
    masses = list(pd.to_numeric(group.columns))
    values = group.values.flatten()
    return interp1d(masses, values, kind='cubic', fill_value="extrapolate")

gby = conversion_table.groupby(['metallicity', 'element'])
interp_fns = gby.apply(get_interpolation)


def getStarYield(metallicity, starMass):
    """
    Get star yield, a function that take a metallicity precentage and the star mass of star and
    return a dictionary of element+renmant keys, for each key a the value of hydrogen converted
    in solar masses. If the metallicity is superior to any key, return the last one
    """
    
    metal_vals = conversion_table['metallicity'].unique()
    closest_metallicity = metal_vals[metal_vals < metallicity].max()
    
    solar_masses = interp_fns[closest_metallicity].apply(lambda x: x(starMass))
    return solar_masses.to_dict()
    
def getStarAllInterpolatedYield(metallicity):
    """
    ***
    """
    x = [i for i in range(100)]
    xi = [i for i in range(200)]

    metal_vals = conversion_table['metallicity'].unique()
    closest_metallicity = metal_vals[metal_vals < metallicity].max()
    
    for el in conversion_table['element'].unique():
        f = interp_fns[closest_metallicity, el]
        yi = f(xi)
        plt.plot(xi, yi, label = "intr.Data")
        plt.suptitle(el)
        y = f(x)
        plt.plot(x, y, "r--", label = "data")
        plt.legend(loc="best")
        plt.show()