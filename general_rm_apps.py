import numpy as np
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats
import rm_misc
import rm_plotting
import milky_way_foreground
from run_ymw16_ne import find_DM,foreground_ne
import resolved_check
import path
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from regions import DS9Parser,write_ds9
from matplotlib import colors
import astropy.units as u
import astropy.constants as const

def stats_board(data_set,confidence=2/3,num_iterations=10000):
    '''
    Collects the stats of an array into a dictionary 
    
    Args:
        data_set ([array]): 1D array of input data
        confidence ([float]): confidence of result as fraction (Default = 2/3 or 66%)
        num_iterations ([float]): number of iterations for bootstrapping (Default = 10000)
    Returns:
        stats ([dictionary]): dictionary of stats in the form of 
        [length of array, max of array, min of array, 
        (mean, lower bound, upper bound), (median, lower bound, upper bound), (std, lower bound, upper bound)]
    
    '''
    alpha=1-confidence
    samples=np.ma.masked_invalid(data_set).compressed()
    
    mean_res=bs.bootstrap(samples,stat_func=bs_stats.mean,alpha=alpha,num_iterations=num_iterations)
    mean=mean_res.value #mean using bootstrapping
    mean_ler=mean-mean_res.lower_bound #lower bound on mean
    mean_uer=mean_res.upper_bound-mean #upper bound on mean
    
    median_res=bs.bootstrap(samples,stat_func=bs_stats.median,alpha=alpha,num_iterations=num_iterations)
    median=median_res.value #median using bootstrapping
    median_ler=median-median_res.lower_bound #lower bound on median
    median_uer=median_res.upper_bound-median #upper bound on median
    
    std_res=bs.bootstrap(samples,stat_func=bs_stats.std,alpha=alpha,num_iterations=num_iterations)
    std=std_res.value #std using bootstrapping
    std_ler=std-std_res.lower_bound #lower bound on std
    std_uer=std_res.upper_bound-std #upper bound on std
    
    stats=dict()
    stats['n']=len(data_set)
    stats['max']=round(np.nanmax(data_set),2)
    stats['min']=round(np.nanmin(data_set),2)
    stats['mean']=(round(mean,2),round(mean_ler,2),round(mean_uer,2))
    stats['median']=(round(median,2),round(median_ler,2),round(median_uer,2))
    stats['std']=(round(std,2),round(std_ler,2),round(std_uer,2))
    
    return stats

import matplotlib.pyplot as plt
import rm_misc

def second_moment(location,cutoff=7,plot=False):
    '''
    Calculate the second moment of a given FDF spectrum (for use with RM-tools)
    
    Args:
        location ([path]): path where the stokes information file is (single path with .dat)
        cutoff ([float]): cutoff of number of times the noise to use to calculate the second moment (Default = 7)
        plot ([bool]): enable plotting of FDF and second moment (Default = False)
    Returns:
        ret ([array]): the second moment and the error of the measurement in units of rad m^-2
    '''
    s_path=path.Path(location)
    base_nam=s_path.split(".")[0]
    FDF_name=f'{base_nam}_FDFclean.dat'
    s_nam=path.Path(FDF_name).name
    data=np.loadtxt(FDF_name)
    base_FDF_nam=path.Path(FDF_name).split("FDF")[0]
    snr_r=rm_misc.snr(f'{base_FDF_nam}RMclean.json')
    fwhm_r=rm_misc.fwhm(f'{base_nam}_RMsynth.json')
    
    phis=data[:,0]
    amp=abs(data[:,2]+1j*data[:,1])
    noise_r=(np.max(amp)/snr_r)
    m1_mask=np.where(amp>noise_r*cutoff)[0]
    m1 = np.sum(phis[m1_mask]*amp[m1_mask])/np.sum(amp[m1_mask])
    m1_err = np.std(phis[m1_mask]*amp[m1_mask])/np.sum(amp[m1_mask])
    peaks, proper = find_peaks(amp, height=noise_r*cutoff)
    h=noise_r*cutoff
    
    if plot==True:
        plt.figure(figsize=(15,10))
        plt.plot(phis,amp)
        plt.plot(phis[peaks], amp[peaks], "x")
        plt.axhline(h,color='C1')
        plt.axhline(0.0,color='gray',ls='--')
        plt.title(f'{s_nam}')
        
    if np.isnan(np.std(phis[peaks])) == False:
        ret = (np.sqrt(np.sum((phis[peaks]-m1)**2*amp[peaks])/np.sum(amp[peaks])),(np.sqrt(m1_err**2+(fwhm_r/snr_r)**2)))
    else:
        ret = (0.0,fwhm_r/snr_r)
    return ret

def sort_xy(pos):
    '''
    Orders coordinates in clockwise order
    Args:
        pos ([array,Nx2]): The positions of vertices in the region in y,x format
        
    Returns:
        array,Nx2: Positions ordered clockwise
    '''
    x,y=pos.T
    x0 = np.mean(x)
    y0 = np.mean(y)

    r = np.sqrt((x-x0)**2 + (y-y0)**2)

    angles = np.where((y-y0) > 0, np.arccos((x-x0)/r), 2*np.pi-np.arccos((x-x0)/r))

    mask = np.argsort(angles)

    x_sorted = x[mask]
    y_sorted = y[mask]

    return np.array([x_sorted, y_sorted]).T

def name_hex(c_name):
    col_rgb = np.array(colors.hex2color(colors.cnames[f'{c_name}']))
    r,g,b=col_rgb*255
    return f"#{int(r):02x}{int(g):02x}{int(b):02x}"

def ds9_poly(pos,color,text,frame='galactic',out_path):
    '''
    Creates ds9 region file that can be used with APLpy
    Args:
        pos ([array,Nx2]): The positions of vertices in the region in y,x format
        color ([str]): The color of the chosen region
        text ([str]): The text associated with the region
        frame ([str]): The frame of the positions (Default = galactic)
        out_path ([str]): Path to save region file (Default = /home/jackl/MC)
        
    Output:
        ds9 region: Creates ds9 region file saved as the txt name given
        
    '''
    pos_list=sort_xy(pos)
    k=''
    for i in pos_list:
        k=k+f'{i[1]},{i[0]},'
    coord_str=k
    centr_pos=pos.sum(axis=0)/len(pos)
    centr_str=f'{centr_pos[1]},{centr_pos[0]}'
    color_name=name_hex(color)
    #print(centr_str)
    poly_str=f'{frame} \npolygon({coord_str}) # color={color_name} width=2 dash=1 dashlist=8 3,'
    form_str=f'{frame} \ntext({centr_str}) # color={color_name} text={text} font="times 14 bold italic" width=2'
    reg_string=poly_str+form_str
    parser = DS9Parser(reg_string)
    regions = parser.shapes.to_regions()
    filename = f'{out_path}_{text}.reg'
    write_ds9(regions, filename)
    
def faraday_cap(freq):
    '''
    Calculates the RM-synthesis capabilities based on the input frequency range
    Args:
        freq ([list]) = input array of frequencies as a list between the highest and lowest frequencies in units of Hertz
    
    Output:
        info ([dict]) = contains the max and min frequencies, frequency step, max and min lambda^2, lamda^2 step, resolution
        in phi space, maximum measurable scale in phi space, and maximum measurable range of phi
    '''
    fmin=np.min(freq/1e9)
    fmax=np.max(freq/1e9)
    fstep=np.mean(abs((np.diff(freq))))/1e6
    fchannels=(((fmax * u.GHz)-(fmin * u.GHz)).to(u.MHz)/(fstep*u.MHz))
    lamsqmax=((const.c/(fmin * u.GHz))**2).to(u.m*u.m)
    lamsqmin=((const.c/(fmax * u.GHz))**2).to(u.m*u.m)
    deltalamsq=(lamsqmax-lamsqmin)
       
    info={'f_max':(fmax*u.GHz).to(u.MHz),'f_min':(fmin*u.GHz).to(u.MHz),'f_step':(fstep*u.MHz).to(u.kHz),
                    'lsq_max':lamsqmax,'lsq_min':lamsqmin,
                   'delta_lsq':deltalamsq,'phi_resolution':((3.79)/deltalamsq)*u.rad,
                   'max_phi_scale':(np.pi*(lamsqmin**-1))*u.rad,
                    'max_phi_range':((np.sqrt(3)*(fchannels*u.rad)/((deltalamsq))))}
    return info
