from astropy.io import ascii as asctab
from astropy.coordinates import SkyCoord 
import numpy as np

def position_getter(location):
    '''
    Creates SkyCoord Catalogue from csv table.

    Parameters:
        location (str): path to csv table

    Returns:
        catalogue (SkyCoord Catalogue): SkyCoord Catalogue of positions
    '''
    
    ra=asctab.read(location)['RA_deg']
    dec=asctab.read(location)['DEC_deg']
    
    return SkyCoord(ra,dec,unit='deg')

def foreground_calculator(location,Map):
    '''
    Calculates the Milky Way foreground contribution for the Magellanic System. 
    Foreground equations were taken from:
        SMC: Mao+2008 ApJ.688.1029M (https://ui.adsabs.harvard.edu/abs/2008ApJ...688.1029M/abstract)
        Bridge: Kaczmarek+2017 MNRAS.467.1776K (https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.1776K/abstract)
        LMC: Mao+2012 ApJ.759.25M (https://ui.adsabs.harvard.edu/abs/2012ApJ...759...25M/abstract)

    Parameters:
        Positions (array, Nx2): An array of the positions of each source in ra and dec
        Map (str): Choice of 'SMC', 'Bridge', and 'LMC' for foreground model

    Returns:
        foreground (array): array of foreground values for each supplied position
    '''
    Positions=position_getter(location)
    
    foreground=[]
    if Map=='SMC': # Foreground map for the SMC
        for i in Positions.icrs:
            foreground.append(46.1-4.9*i.ra.deg*(np.cos(np.deg2rad(i.dec.deg))))
            
    if Map=='Bridge': # Foreground map for the Magellanic Bridge
        for i in Positions.galactic:
            foreground.append(-0.511*i.l.deg+1.28*i.b.deg+225)
            
    if Map=='LMC': # Foreground map for the LMC
        #Centre of LMC
        LMC_centre=SkyCoord(ra='5h16m00s',dec='-68d41m00s',frame='icrs')
        for i in Positions.icrs:
            foreground.append(28.7+0.68*(i.ra.deg-LMC_centre.ra.deg)-0.39*(i.dec.deg-LMC_centre.dec.deg))
        
    return np.array(foreground)
