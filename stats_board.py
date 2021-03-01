import numpy as np
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats

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
    
    mean_str=str(
        bs.bootstrap(samples,stat_func=bs_stats.mean,alpha=alpha,num_iterations=num_iterations)).split(' ')
    mean=float(mean_str[0]) #mean using bootstrapping
    mean_ler=mean-float(mean_str[-2][1:-1]) #lower bound on mean
    mean_uer=float(mean_str[-1][:-1])-mean #upper bound on mean
    
    median_str=str(
        bs.bootstrap(samples,stat_func=bs_stats.median,alpha=alpha,num_iterations=num_iterations)).split(' ')
    median=float(median_str[0]) #median using bootstrapping
    median_ler=median-float(median_str[-2][1:-1]) #lower bound on median
    median_uer=float(median_str[-1][:-1])-median #upper bound on median
    
    std_str=str(
        bs.bootstrap(samples,stat_func=bs_stats.std,alpha=alpha,num_iterations=num_iterations)).split(' ')
    std=float(std_str[0]) #std using bootstrapping
    std_ler=std-float(std_str[-2][1:-1]) #lower bound on std
    std_uer=float(std_str[-1][:-1])-std #upper bound on std
    
    stats=dict()
    stats['n']=len(data_set)
    stats['max']=round(np.nanmax(data_set),2)
    stats['min']=round(np.nanmin(data_set),2)
    stats['mean']=(round(mean,2),round(mean_ler,2),round(mean_uer,2))
    stats['median']=(round(median,2),round(median_ler,2),round(median_uer,2))
    stats['std']=(round(std,2),round(std_ler,2),round(std_uer,2))
    
    return stats
