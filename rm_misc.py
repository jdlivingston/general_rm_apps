import json

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
