import numpy as np
import path
import astropy.constants as c
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import json
from astropy.coordinates import SkyCoord 
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

modeltypes={2:'1 Screen with EFD',
            3:'1 Screen with EFD and IFD',
            4:'2 Screens with different EFD',
            5:'2 Screens with EFD and IFD',
            6:'3 Screens'}

# Collects all the models I want which in this case are m2 and m4
qumodels = {}
import importlib.util
for m in [2,3,4,5,6]:
    spec = importlib.util.spec_from_file_location(f"qumodels.m{m}",
                                                  f"/home/jackl/RM/RMtools_1D/models_ns/m{m}.py")
    m2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m2)
    qumodels[m] = m2

def snr(location):
    '''
    Gets the dictionary of parameters for each model for each object
    
    Input:
    location = path where the model.json file is (single path)
    
    Output:
    results = the model parameters
    '''
    with open(location, 'r') as f:
        return json.load(f)['snrPIfit']
    
def fwhm(location):
    '''
    Gets the dictionary of parameters for each model for each object
    
    Input:
    location = path where the model.json file is (single path)
    
    Output:
    results = the model parameters
    '''
    with open(location, 'r') as f:
        return json.load(f)['fwhmRMSF']
    
def phi(location=None):
    '''
    Gets the dictionary of parameters for each model for each object
    
    Input:
    location = path where the model.json file is (single path)
    
    Output:
    results = the model parameters
    '''
    with open(location, 'r') as f:
        return json.load(f)['phiPeakPIfit_rm2']
    
    
def stokes_plot(location):
    '''
    Plots suite of stokes parameters
    
    Input:
    location = path where the stokes information file is (single path with .dat)
    
    Ouput:
    Plot 1 = Stokes I vs Frequency
    Plot 2 = Stokes q,u,p vs lamdba^2
    Plot 3 = Stokes q vs Stokes u
    Plot 4 = chi vs lamdbda^2 (RM slope)
    '''
    s_path=path.Path(location)
    freq=(np.loadtxt(s_path)[:,0])
    stokes_i=np.loadtxt(s_path)[:,1]
    stokes_q=np.loadtxt(s_path)[:,2]
    stokes_u=np.loadtxt(s_path)[:,3]
    stokes_v=np.loadtxt(s_path)[:,5]
    
    lamsq=(c.c.value/freq)**2
    
    dpsi=np.rad2deg((abs(stokes_v/stokes_q+1j*stokes_v/stokes_u
                 )*(stokes_u/stokes_q))*((stokes_u/stokes_q)**2+1)**-1)
    
    fig = plt.figure(figsize=(15,15))
    ax2 = plt.subplot(412)
    ax2.set_title(f'{s_path.name}')
    ax2.errorbar(freq/1e9,stokes_i,yerr=stokes_v*np.sqrt(2),fmt='.',color='k')
    ax2.set_xlabel('Frequency [GHz]')
    ax2.set_ylabel('Stokes I [Jy/beam]')

    ax3 = plt.subplot(425)
    ax3.errorbar(lamsq,stokes_q/stokes_i,yerr=stokes_v/stokes_i,fmt='.',color='blue',label='Stokes q')
    ax3.errorbar(lamsq,stokes_u/stokes_i,yerr=stokes_v/stokes_i,fmt='.',color='red',label='Stokes u')
    ax3.errorbar(lamsq,abs(stokes_q/stokes_i+1j*stokes_u/stokes_i),yerr=np.sqrt(2)*stokes_v/stokes_i,
                 fmt='.',color='black',label='Stokes p')
    ax3.set_ylabel('Stokes q,u,p [pol. fraction]')
    ax3.set_xlabel('$\lambda^{2}\,[\mathrm{m}^{2}]$')
    
    ax4 = plt.subplot(414)
    ax4.errorbar(lamsq,0.5*np.rad2deg(np.arctan(stokes_u/stokes_q)),yerr=dpsi,fmt='.',color='black',label='$\chi$');
    ax4.set_xlabel('$\lambda^{2} [\mathrm{m}^{2}]$')
    ax4.set_ylabel('$\chi$ [deg]')
    ax4.set_ylim(-90,90)
    
    ax5 = plt.subplot(426)
    ax5.scatter(stokes_q/stokes_i,stokes_u/stokes_i,marker='.',c=lamsq,cmap='gist_rainbow',zorder=3)
    ax5.set_xlabel('Stokes q')
    ax5.set_ylabel('Stokes u')
    fig.tight_layout()
    
def FDF_plot(location):
    '''
    Plots FDF of a source
    
    Input:
    location = path where the stokes information file is (single path with .dat)
    
    Ouput:
    Plot 1 = FDF of source
    '''
    s_path=path.Path(location)
    base_nam=s_path.split(".")[0]
    cl_data=np.loadtxt(f'{base_nam}_FDFclean.dat')
    di_data=np.loadtxt(f'{base_nam}_FDFdirty.dat')
    rmsf_data=np.loadtxt(f'{base_nam}_RMSF.dat')
    snr_r=snr(f'{base_nam}_RMclean.json')
    phi_r=phi(f'{base_nam}_RMclean.json')
    fwhm_r=fwhm(f'{base_nam}_RMsynth.json')
    FDF_name=f'{base_nam}_FDFclean.dat'
    
    if 
      
        h=fits.getheader(f'{path.Path(FDF_name).parent.parent}/{path.Path(FDF_name).parent.parent.name}.i.smooth.fits')
        w=WCS(h)
        if len(path.Path(FDF_name).name.split('_')) == 5:
                pixx=int(path.Path(FDF_name).name.split('_')[2])
                pixy=int(path.Path(FDF_name).name.split('_')[3])
                posies=[w.wcs_pix2world(pixx,pixy,0,0)[0],w.wcs_pix2world(pixx,pixy,0,0)[1]]
                #print(np.asarray(posies))
        else:
                pixx=int(path.Path(FDF_name).name.split('F')[0].split('-')[-2])
                pixy=int(path.Path(FDF_name).name.split('F')[0].split('-')[-1][:-1])
                posies=[w.wcs_pix2world(pixx,pixy,0,0)[0],w.wcs_pix2world(pixx,pixy,0,0)[1]]
                #print(np.asarray(posies))
        s_pos=np.asarray(posies)
        cat=SkyCoord(s_pos[0],s_pos[1],unit='deg')
        s_ra=cat.ra.to_string(unit=u.hourangle, sep='', precision =2, pad=True)
        s_dec=cat.dec.to_string(sep='', precision=2, alwayssign=True, pad=True)
        title_name='J{0}{1}'.format(s_ra,s_dec)
    
    #loadin data
    cl_phis=cl_data[:,0]
    cl_real=cl_data[:,1]
    cl_imag=cl_data[:,2]
    
    rmsf_phis=rmsf_data[:,0]
    rmsf_real=rmsf_data[:,1]
    rmsf_imag=rmsf_data[:,2]
    
    di_phis=di_data[:,0]
    di_real=di_data[:,1]
    di_imag=di_data[:,2]
    
    #bit of maths
    cl_total=abs(cl_real+1j*cl_imag)
    rmsf_total=abs(rmsf_real+1j*rmsf_imag)
    di_total=abs(di_real+1j*di_imag)
    
    cutoff=7*(np.max(cl_total)/snr_r)
    
    #peak fitting
    peaks, _ = find_peaks(cl_total, height=cutoff)
    
    #plotting
    fig = plt.figure(figsize=(15,15))
    ax1 = plt.subplot(411)
    ax1.set_title(f'{title_name} [RM = {phi_r:.0f} $\pm$ {fwhm_r/snr_r:.0f} '+'$\mathrm{rad}\,\mathrm{m}^{-2}$]')
    ax1.plot(cl_phis,cl_total,label='F($\phi$) clean',zorder=3)
    ax1.plot(di_phis,di_total,label='F($\phi$) dirty')
    ax1.plot(rmsf_phis+phi_r,rmsf_total*np.max(cl_total),label='RMSF')
    ax1.scatter(cl_phis[peaks],cl_total[peaks],marker='x',color='k')
    ax1.axhline(cutoff,ls='--',color='red',label='7 $\sigma$ cutoff')
    ax1.set_xlabel('$\phi\,[\,\mathrm{rad}\,\mathrm{m}^{-2}]$')
    ax1.set_ylabel('|F($\phi$)|')
    ax1.set_xlim(min(cl_phis),max(cl_phis))
    ax1.legend()
    fig.tight_layout()
    
def full_plot(location):
    '''
    Plots suite of stokes parameters and FDF of a source
    
    Input:
    location = path where the stokes information file is (single path with .dat)
    
    Ouput:
    Plot 1 = FDF of source
    Plot 2 = Stokes I vs Frequency
    Plot 3 = Stokes q,u,p vs lamdba^2
    Plot 4 = Stokes q vs Stokes u
    Plot 5 = chi vs lamdbda^2 (RM slope)
    '''
    s_path=path.Path(location)
    freq=(np.loadtxt(s_path)[:,0])
    stokes_i=np.loadtxt(s_path)[:,1]
    stokes_q=np.loadtxt(s_path)[:,2]
    stokes_u=np.loadtxt(s_path)[:,3]
    stokes_v=np.loadtxt(s_path)[:,5]
    
    lamsq=(c.c.value/freq)**2
    
    dpsi=np.rad2deg((abs(stokes_v/stokes_q+1j*stokes_v/stokes_u
                 )*(stokes_u/stokes_q))*((stokes_u/stokes_q)**2+1)**-1)
    base_nam=s_path.split(".")[0]
    cl_data=np.loadtxt(f'{base_nam}_FDFclean.dat')
    di_data=np.loadtxt(f'{base_nam}_FDFdirty.dat')
    rmsf_data=np.loadtxt(f'{base_nam}_RMSF.dat')
    snr_r=snr(f'{base_nam}_RMclean.json')
    phi_r=phi(f'{base_nam}_RMclean.json')
    fwhm_r=fwhm(f'{base_nam}_RMsynth.json')
    FDF_name=f'{base_nam}_FDFclean.dat'
  
    h=fits.getheader(f'{path.Path(FDF_name).parent.parent}/{path.Path(FDF_name).parent.parent.name}.i.smooth.fits')
    w=WCS(h)
    if len(path.Path(FDF_name).name.split('_')) == 5:
            pixx=int(path.Path(FDF_name).name.split('_')[2])
            pixy=int(path.Path(FDF_name).name.split('_')[3])
            posies=[w.wcs_pix2world(pixx,pixy,0,0)[0],w.wcs_pix2world(pixx,pixy,0,0)[1]]
            #print(np.asarray(posies))
    else:
            pixx=int(path.Path(FDF_name).name.split('F')[0].split('-')[-2])
            pixy=int(path.Path(FDF_name).name.split('F')[0].split('-')[-1][:-1])
            posies=[w.wcs_pix2world(pixx,pixy,0,0)[0],w.wcs_pix2world(pixx,pixy,0,0)[1]]
            #print(np.asarray(posies))
    s_pos=np.asarray(posies)
    cat=SkyCoord(s_pos[0],s_pos[1],unit='deg')
    s_ra=cat.ra.to_string(unit=u.hourangle, sep='', precision =2, pad=True)
    s_dec=cat.dec.to_string(sep='', precision=2, alwayssign=True, pad=True)
    title_name='J{0}{1}'.format(s_ra,s_dec)
    
    #loadin data
    cl_phis=cl_data[:,0]
    cl_real=cl_data[:,1]
    cl_imag=cl_data[:,2]
    
    rmsf_phis=rmsf_data[:,0]
    rmsf_real=rmsf_data[:,1]
    rmsf_imag=rmsf_data[:,2]
    
    di_phis=di_data[:,0]
    di_real=di_data[:,1]
    di_imag=di_data[:,2]
    
    #bit of maths
    cl_total=abs(cl_real+1j*cl_imag)
    rmsf_total=abs(rmsf_real+1j*rmsf_imag)
    di_total=abs(di_real+1j*di_imag)
    
    cutoff=7*(np.max(cl_total)/snr_r)
    
    #peak fitting
    peaks, _ = find_peaks(cl_total, height=cutoff)
    
    #plotting
    fig = plt.figure(figsize=(15,15))
    ax1 = plt.subplot(411)
    ax1.set_title(f'{title_name} [RM = {phi_r:.0f} $\pm$ {fwhm_r/snr_r:.0f} '+'$\mathrm{rad}\,\mathrm{m}^{-2}$]')
    ax1.plot(cl_phis,cl_total,label='F($\phi$) clean',zorder=3,color='blue')
    ax1.plot(di_phis,di_total,label='F($\phi$) dirty',color='orange')
    ax1.plot(rmsf_phis+phi_r,rmsf_total*np.max(cl_total),label='RMSF',color='green')
    ax1.scatter(cl_phis[peaks],cl_total[peaks],marker='x',color='k')
    ax1.axhline(cutoff,ls='--',color='red',label='7 $\sigma$ cutoff')
    ax1.set_xlabel('$\phi\,[\,\mathrm{rad}\,\mathrm{m}^{-2}]$')
    ax1.set_ylabel('|F($\phi$)|')
    ax1.set_xlim(min(cl_phis),max(cl_phis))
    ax1.legend()
    
    ax2 = plt.subplot(412)
    ax2.errorbar(freq/1e9,stokes_i,yerr=stokes_v*np.sqrt(2),fmt='.',color='k')
    ax2.set_xlabel('Frequency [GHz]')
    ax2.set_ylabel('Stokes I [Jy/beam]')

    ax3 = plt.subplot(425)
    ax3.errorbar(lamsq,stokes_q/stokes_i,yerr=stokes_v/stokes_i,fmt='.',color='blue',label='Stokes q')
    ax3.errorbar(lamsq,stokes_u/stokes_i,yerr=stokes_v/stokes_i,fmt='.',color='red',label='Stokes u')
    ax3.errorbar(lamsq,abs(stokes_q/stokes_i+1j*stokes_u/stokes_i),yerr=np.sqrt(2)*stokes_v/stokes_i,
                 fmt='.',color='black',label='Stokes p')
    ax3.set_ylabel('Stokes q,u,p [pol. fraction]')
    ax3.set_xlabel('$\lambda^{2}\,[\mathrm{m}^{2}]$')
    
    ax4 = plt.subplot(414)
    ax4.errorbar(lamsq,0.5*np.rad2deg(np.arctan(stokes_u/stokes_q)),yerr=dpsi,fmt='.',color='black',label='$\chi$');
    ax4.set_xlabel('$\lambda^{2} [\mathrm{m}^{2}]$')
    ax4.set_ylabel('$\chi$ [deg]')
    ax4.set_ylim(-90,90)
    
    ax5 = plt.subplot(426)
    ax5.scatter(stokes_q/stokes_i,stokes_u/stokes_i,marker='.',c=lamsq,cmap='gist_rainbow',zorder=3)
    ax5.set_xlabel('Stokes q')
    ax5.set_ylabel('Stokes u')
    fig.tight_layout()
    
def QU_plot(location,QU_dict):
    '''
    Plots suite of stokes parameters and FDF of a source with QU fitted parameters
    
    Input:
    location = path where the stokes information file is (single path with .dat)
    QU_dict = path where the sources's best fitting model is located (limited to m2 and m4 models)
    
    Ouput:
    Plot 1 = FDF of source
    Plot 2 = Stokes I vs Frequency
    Plot 3 = Stokes q,u,p vs lamdba^2
    Plot 4 = Stokes q vs Stokes u
    Plot 5 = chi vs lamdbda^2 (RM slope)
    '''
    s_path=path.Path(location)
    freq=(np.loadtxt(s_path)[:,0])
    stokes_i=np.loadtxt(s_path)[:,1]
    stokes_q=np.loadtxt(s_path)[:,2]
    stokes_u=np.loadtxt(s_path)[:,3]
    stokes_v=np.loadtxt(s_path)[:,5]
    
    lamsq=(c.c.value/freq)**2
    
    dpsi=np.rad2deg((abs(stokes_v/stokes_q+1j*stokes_v/stokes_u
                 )*(stokes_u/stokes_q))*((stokes_u/stokes_q)**2+1)**-1)
    base_nam=s_path.split(".")[0]
    cl_data=np.loadtxt(f'{base_nam}_FDFclean.dat')
    di_data=np.loadtxt(f'{base_nam}_FDFdirty.dat')
    rmsf_data=np.loadtxt(f'{base_nam}_RMSF.dat')
    snr_r=snr(f'{base_nam}_RMclean.json')
    phi_r=phi(f'{base_nam}_RMclean.json')
    
    #print(total_dictionaries[i])
    read_dictionary = np.load(QU_dict,allow_pickle='TRUE').item()
    mod_number=int(list(read_dictionary.keys())[0].split('m')[-1])
    modq=qumodels[mod_number].model(list(read_dictionary.values())[0],lamsq).real
    modu=qumodels[mod_number].model(list(read_dictionary.values())[0],lamsq).imag
    
    if mod_number == 2:
        RM_radm2=list(read_dictionary.values())[0]['RM_radm2']
    if mod_number == 4:
        RM_radm2=[list(read_dictionary.values())[0]['RM1_radm2'],list(read_dictionary.values())[0]['RM2_radm2']]
    if mod_number == 6:
        RM_radm2=[list(read_dictionary.values())[0]['RM1_radm2'],list(read_dictionary.values())[0]['RM2_radm2']]
    
    #loadin data
    cl_phis=cl_data[:,0]
    cl_real=cl_data[:,1]
    cl_imag=cl_data[:,2]
    
    rmsf_phis=rmsf_data[:,0]
    rmsf_real=rmsf_data[:,1]
    rmsf_imag=rmsf_data[:,2]
    
    di_phis=di_data[:,0]
    di_real=di_data[:,1]
    di_imag=di_data[:,2]
    
    #bit of maths
    cl_total=abs(cl_real+1j*cl_imag)
    rmsf_total=abs(rmsf_real+1j*rmsf_imag)
    di_total=abs(di_real+1j*di_imag)
    
    cutoff=7*(np.max(cl_total)/snr_r)
    
    #peak fitting
    peaks, _ = find_peaks(cl_total, height=cutoff)
    
    #plotting
    fig = plt.figure(figsize=(15,15))
    ax1 = plt.subplot(411)
    ax1.set_title(f'{s_path.name}')
    ax1.plot(cl_phis,cl_total,label='FDFclean',zorder=3,color='blue')
    ax1.plot(di_phis,di_total,label='FDFdirty',color='orange')
    ax1.plot(rmsf_phis+phi_r,rmsf_total*np.max(cl_total),label='RMSF',color='green')
    ax1.scatter(cl_phis[peaks],cl_total[peaks],marker='x',color='k')
    ax1.axhline(cutoff,ls='--',color='red',label='7 $\sigma$ cutoff')
    ax1.set_xlabel('$\phi\,[\,\mathrm{rad}\,\mathrm{m}^{-2}]$')
    ax1.set_ylabel('|F($\phi$)|')
    ax1.set_xlim(min(cl_phis),max(cl_phis))
    
    if mod_number == 2:
        ax1.axvline(RM_radm2,color='k',ls='-.',
                         label='$\mathrm{RM}_{1} =$'+f'{RM_radm2:.0f}'+'$\,\mathrm{rad}\,\mathrm{m}^{-2}$')
    if mod_number == 4:
        ax1.axvline(RM_radm2[0],color='k',ls='-.',
                         label='$\mathrm{RM}_{1} =$'+f'{RM_radm2[0]:.0f}'+'$\,\mathrm{rad}\,\mathrm{m}^{-2}$')
        ax1.axvline(RM_radm2[1],color='k',ls='-.',
                         label='$\mathrm{RM}_{2} =$'+f'{RM_radm2[1]:.0f}'+'$\,\mathrm{rad}\,\mathrm{m}^{-2}$')
        
    if mod_number == 6:
        ax1.axvline(RM_radm2[0],color='k',ls='-.',
                         label='$\mathrm{RM}_{1} =$'+f'{RM_radm2[0]:.0f}'+'$\,\mathrm{rad}\,\mathrm{m}^{-2}$')
        ax1.axvline(RM_radm2[1],color='k',ls='-.',
                         label='$\mathrm{RM}_{2} =$'+f'{RM_radm2[1]:.0f}'+'$\,\mathrm{rad}\,\mathrm{m}^{-2}$')
        ax1.axvline(-346.41549683,color='k',ls='-.',
                         label='$\mathrm{RM}_{2} =$'+f'{-346.41549683:.0f}'+'$\,\mathrm{rad}\,\mathrm{m}^{-2}$')
    
    ax2 = plt.subplot(412)
    ax2.errorbar(freq/1e9,stokes_i,yerr=stokes_v*np.sqrt(2),fmt='.',color='k')
    ax2.set_xlabel('Frequency [GHz]')
    ax2.set_ylabel('Stokes I [Jy/beam]')

    ax3 = plt.subplot(425)
    ax3.errorbar(lamsq,stokes_q/stokes_i,yerr=stokes_v/stokes_i,fmt='.',color='dodgerblue',label='Stokes q')
    ax3.errorbar(lamsq,stokes_u/stokes_i,yerr=stokes_v/stokes_i,fmt='.',color='tomato',label='Stokes u')
    ax3.errorbar(lamsq,abs(stokes_q/stokes_i+1j*stokes_u/stokes_i),yerr=np.sqrt(2)*stokes_v/stokes_i,
                 fmt='.',color='black',label='Stokes p')
    ax3.set_ylabel('Stokes q,u,p [pol. fraction]')
    ax3.set_xlabel('$\lambda^{2}\,[\mathrm{m}^{2}]$')
    
    ax3.plot(lamsq,abs(modq+1j*modu),zorder=4,linewidth=5.0,color='lime',ls='--')
    ax3.plot(lamsq,modu,zorder=3,linewidth=5.0,color='red',ls='--')
    ax3.plot(lamsq,modq,zorder=3,linewidth=5.0,color='blue',ls='--')
    
    ax4 = plt.subplot(414)
    ax4.errorbar(lamsq,0.5*np.rad2deg(np.arctan(stokes_u/stokes_q)),yerr=dpsi,fmt='.',color='black',label='$\chi$');
    ax4.set_xlabel('$\lambda^{2} [\mathrm{m}^{2}]$')
    ax4.set_ylabel('$\chi$ [deg]')
    ax4.set_ylim(-90,90)
    ax4.plot(lamsq,0.5*np.rad2deg(np.arctan(modu/modq))
                 ,ls='--',color='lime',linewidth=5,zorder=3)
    
    ax5 = plt.subplot(426)
    ax5.scatter(stokes_q/stokes_i,stokes_u/stokes_i,marker='.',c=lamsq,cmap='gist_rainbow',zorder=3)
    ax5.set_xlabel('Stokes q')
    ax5.set_ylabel('Stokes u')
    
    ax5.plot(modq,modu,zorder=3,ls='--',color='k',linewidth=5)
    
    fig.tight_layout()
