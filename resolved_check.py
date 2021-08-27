import numpy as np
import path
from astropy.io import fits
from spectral_cube import SpectralCube
from radio_beam import Beam
from scipy import interpolate
import matplotlib.pyplot as plt

def beam_axes(location):
    '''
    Determine the beam axes of a source
    
    Args:
        location ([str]): path to .dat files of Stokes Spectrum
        
    Returns:
        ret ([list]): RA and DEC size of the beam
    '''
    name = path.Path(location).parent.parent+'/'+path.Path(location).parent.parent.name+'.i.single.fits'
    header = fits.getheader(name)
    beam = Beam.from_fits_header(header)
    maj=float(str(beam).split('=')[1].split(' ')[0])
    mino=float(str(beam).split('=')[2].split(' ')[0])
    pa=np.deg2rad(float(str(beam).split('=')[3].split(' ')[0]))
    ra=maj*abs(np.sin(pa))+mino*abs(np.cos(pa))
    dec=maj*abs(np.cos(pa))+mino*abs(np.sin(pa))
    return ra,dec

def det_res(location,plot=False):
    '''
    Checks if a source is resolved by inspecting Stokes I 
    
    Args:
        location ([str]): path to .dat files of Stokes Spectrum
        plot ([bool]): enable plotting of stokes I peak (Default = False)
    Returns:
        ret ([array]): 2D array that has a 1 if the source dimension is resolved and 0 if unresolved
    
    '''
    ret = np.zeros(2)
    pix=(int(path.Path(location).name.split('-')[-2]),int(path.Path(location).name.split('-')[-1].split('.')[0]))
    file=str(path.Path(location).parent.parent+'/'+path.Path(location).parent.parent.name+'.i.smooth.fits')
    convert=fits.getheader(file.split('.smooth')[0]+'.single.fits')['CDELT2']*3600
    beam=np.array(beam_axes(location))*convert
    cube=SpectralCube.read(file)
    subcube_y=cube.subcube(ylo=pix[0]-30,yhi=pix[0]+30)
    subcube_x=cube.subcube(xlo=pix[1]-30,xhi=pix[1]+30)
    x=np.arange(0,60)
    y=np.arange(0,60)
    subcube_y=cube.subcube(ylo=pix[0]-30,yhi=pix[0]+30)[0][:,pix[1]]
    subcube_x=cube.subcube(xlo=pix[1]-30,xhi=pix[1]+30)[0][pix[0],:]

    fwhm_x=np.max(subcube_x.value)/2
    fwhm_y=np.max(subcube_y.value)/2


    xfunc=interpolate.interp1d(x, subcube_x.value)
    yfunc=interpolate.interp1d(y, subcube_y.value)
    bem=np.mean(beam)
    ret[0]=(((xfunc(30-bem)>1.*fwhm_x) or (xfunc(30+bem)>1.*fwhm_x)))
    ret[1]=(((yfunc(30-bem)>1.*fwhm_y) or (yfunc(30+bem)>1.*fwhm_y)))
    if plot == True:
        fig = plt.figure(figsize=(5,3),dpi=200)
        ax = fig.subplots(nrows=1,ncols=2,sharey=True)
        ax[0].plot(x,subcube_x)
        ax[0].axhline(np.max(subcube_x.value)/2,ls='--',label='fwhm',color='red')
        ax[1].axhline(np.max(subcube_y.value)/2,ls='--',label='fwhm',color='red')
        ax[1].plot(y,subcube_y)
        ax[0].set_xlabel('RA')
        ax[1].set_xlabel('DEC')
        ax[0].axvline(30+np.mean(beam),ls='--',color='green',label='beam')
        ax[0].axvline(30-np.mean(beam),ls='--',color='green')
        ax[1].axvline(30+np.mean(beam),ls='--',color='green',label='beam')
        ax[1].axvline(30-np.mean(beam),ls='--',color='green')
        ax[0].set_ylabel('Intensity [Jy/beam]')
        ax[1].legend(fontsize=7)
    return ret
