import numpy as np
from dotted_dict import DottedDict

# This function creates a dictionay where the coordinated of the bunch will be saved every turn. Careful
# The values of the dictionary are different in every turn.
# If later you want to analyze the turn by turn data you need to save them during the tracking.

def create_bunch(particles):
    bunch = DottedDict()
    bunch.x = np.zeros(particles)
    bunch.px = np.zeros(particles)
    bunch.y = np.zeros(particles)
    bunch.py = np.zeros(particles)
    return bunch


# #--- <x^2> ---
def mean2(numb):
    return np.nanmean( (numb - np.nanmean(numb))**2 )

# #--- <xx'> ---
def mean3(numbx , numbpx):
    return np.nanmean( (numbx - np.nanmean(numbx)) * (numbpx - np.nanmean(numbpx)) )

# #--- sqrt(<x^2> * <px^2> - <xx'>^2) --- compute statistical emittance
def cmp_emit(position, angle):
    return  np.sqrt(mean2(position) * mean2(angle) - mean3(position,angle)**2) # geometrical emittance

# # Get the parameters used by the RF map and the drift element, for the inclusion of the longitudinal motion.
# # Dispersion included. 
def prepere_longitudinal_motion(m0, c, C0, E_0, gamma_tr, sigma_delta_madx, rf): 
    
    # 1. Estimate the parameters of the relativistic particles
    E_rest = m0 # [eV]
    P0_C = np.sqrt(E_0**2-E_rest**2)  # Acutally it is P0*C --> reference momentum times the speed of ligth --> [eV]
    gamma_0 =  E_0/E_rest # gamma realtivistic of the reference particle  # crosscheckd with mad
    beta_0 = np.sqrt(1-1/gamma_0**2) # beta realtivistic of the reference particle
    f_rev_0 = (beta_0*c)/C0 # revolution frequency [Hz]
    T_rev_0 = 1/f_rev_0 # revolution period [s]
    
    
    # estimate the harmonic number
    h_float  = rf.f_RF/ f_rev_0
    h = int(rf.f_RF/ f_rev_0) 
    
    # 2. Estimate synchrotron parameters
    alpha_c = (1/(gamma_tr**2)) # compaction factor
    eta = alpha_c - 1/(gamma_0**2)  # slip factor
    
    if eta < 0 : 
        phi_s = np.pi
    else:
        phi_s = 0
    
    # 3. Estimate the harmonic number
    h_float  = rf.f_RF/ f_rev_0
    h = int(rf.f_RF/ f_rev_0) 
   
    
    # 4. Calculate theoretically the synchrotron tune
    Q_s = np.sqrt(rf.V_RF*eta*h/(2*np.pi*E_0*(beta_0**2)) * np.cos(phi_s))
    
    # 5. Matching
    sigma_delta = sigma_delta_madx/(beta_0**2) # rms bunch length in [m]
    beta_z = C0*eta/(2*np.pi*Q_s) #[m]
    sigma_z = beta_z * sigma_delta
    
    
    my_dict = DottedDict()
    my_dict.eta = eta
    my_dict.P0_C = P0_C
    my_dict.beta_0 = beta_0
    my_dict.sigma_z = sigma_z
    my_dict.sigma_delta = sigma_delta
    my_dict.h = h
    
    return my_dict

def get_rf_bucket(c, beta_0, eta, h, E_0, rf): # stationary rf bucket, no acceleration
    bucket_length = c/rf.f_RF # [m]
    max_bucket_height = beta_0*np.sqrt((2*rf.V_RF)/(np.pi*h*eta*E_0)) # no units
    z_left = - bucket_length/2. # [m]
    z_right =  bucket_length/2. # [m]
    return z_left, z_right