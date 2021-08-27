import numpy as np
from astropy.coordinates import SkyCoord
import subprocess

def foreground_ne(cat,dist):
    '''
    Args:
        cat ([SkyCoord]): SkyCoord coordinate list
        dist ([float]): Distance to DM
        
    Return:
        ret ([array]): foreground ne for the corresponding position
    '''
    ret=np.zeros((len(cat.galactic)))
    for i in enumerate(cat.galactic):
        check = subprocess.run(['/home/jackl/YMW_electron/ymw16_ne',
                               '-d','/home/jackl/YMW_electron/', f'{i[1].l.deg}', f'{i[1].b.deg}' ,f'{dist}' ,'2'],
                      stdout=subprocess.PIPE)
        ret[i[0]]=str(check.stdout).split('ne:')[-1].split('\\n')[0]
    return ret

def find_DM(cat,dist):
    '''
    Args:
        cat ([SkyCoord]): SkyCoord coordinate list
        dist ([float]): Distance to DM
        
    Return:
        ret ([array]): [foreground DM for the corresponding position, the modelled DM for the cloud]
    '''
    ret=np.zeros((len(cat.galactic),2))
    for i in enumerate(cat.galactic):
        check = subprocess.run(['/home/jackl/YMW_electron/ymw16',
                               '-d','/home/jackl/YMW_electron/','MC', f'{i[1].l.deg}', f'{i[1].b.deg}' ,f'{dist}' ,'2'],
                      stdout=subprocess.PIPE)
        
        DM_MW=str(check.stdout).split()[8]
        DM_MC=str(check.stdout).split()[10]
        
        ret[i[0]]=(float(DM_MW),float(DM_MC))
    return ret
