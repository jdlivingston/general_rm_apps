from astropy.io import ascii as asctab
from astropy.coordinates import SkyCoord 
import numpy as np

def position_getter(location):
    '''
    Creates SkyCoord Catalogue from csv table.
    Args:
        location ([str]): path to csv table
    Returns:
        ret ([SkyCoord Catalogue]): SkyCoord Catalogue of positions
    '''
    if 'RA_deg' in asctab.read(location).colnames:
        ra = asctab.read(location)['RA_deg']
        dec = asctab.read(location)['DEC_deg']
        ret = SkyCoord(ra,dec,unit='deg')
    elif 'GLON_deg' in asctab.read(location).colnames:
        glon = asctab.read(location)['GLON_deg']
        glat = asctab.read(location)['GLAT_deg']
        ret = SkyCoord(glon,glat,frame='galactic',unit='deg')
    return ret.icrs

def foreground_calculator(positions,fore_map,table_read=False):
    '''
    Calculates the Milky Way foreground contribution for the Magellanic System. 
    Foreground equations were taken from:
        SMC: Mao+2008 ApJ.688.1029M (https://ui.adsabs.harvard.edu/abs/2008ApJ...688.1029M/abstract)
        Bridge: Kaczmarek+2017 MNRAS.467.1776K (https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.1776K/abstract)
        LMC: Mao+2012 ApJ.759.25M (https://ui.adsabs.harvard.edu/abs/2012ApJ...759...25M/abstract)
    Args:
        Positions ([array, Nx2]): An array of the positions of each source in ra and dec
        Map ([str]): Choice of 'SMC', 'Bridge', and 'LMC' for foreground model
    Returns:
        foreground ([array]): array of foreground values for each supplied position
    '''
    if table_read == True:
        pos=position_getter(positions)
    else:
        pos=positions
        
    foreground=np.zeros((len(pos)))
    if fore_map=='SMC': # Foreground map for the SMC
        for s,i in enumerate(pos.icrs):
            foreground[s]=46.1-4.9*i.ra.deg*(np.cos(np.deg2rad(i.dec.deg)))
            
    if fore_map=='Bridge': # Foreground map for the Magellanic Bridge
        for s,i in enumerate(pos.galactic):
            foreground[s]=-0.511*i.l.deg+1.28*i.b.deg+225
            
    if fore_map=='LMC': # Foreground map for the LMC
        #Centre of LMC
        LMC_centre=SkyCoord(ra='5h16m00s',dec='-68d41m00s',frame='icrs')
        for s,i in enumerate(pos.icrs):
            foreground[s]=28.7+0.68*(i.ra.deg-LMC_centre.ra.deg)-0.39*(i.dec.deg-LMC_centre.dec.deg)
        
    return foreground
