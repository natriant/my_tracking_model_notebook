from math import *
import random
import numpy as np 

def rotation_no_twiss(Qx_rot, Qy_rot, twiss , x, px, y, py):
    x1 = cos(Qx_rot)*x+sin(Qx_rot)*px
    px1 = -sin(Qx_rot)*x+cos(Qx_rot)*px
    y1 = cos(Qy_rot)*y+sin(Qy_rot)*py
    py1 = -sin(Qy_rot)*y+cos(Qy_rot)*py
    return x1, px1, y1, py1

def rotation(Qx_rot, Qy_rot, twiss, x, px, y, py):
    x1 = (cos(Qx_rot) + twiss.alpha_x*sin(Qx_rot))*x + twiss.beta_x*sin(Qx_rot)*px
    px1 = -twiss.gamma_x*sin(Qx_rot)*x+ (cos(Qx_rot)-twiss.alpha_x*sin(Qx_rot))*px
    y1 = (cos(Qy_rot) + twiss.alpha_y*sin(Qy_rot))*y + twiss.beta_y*sin(Qy_rot)*py
    py1 = -twiss.gamma_y*sin(Qy_rot)*y+ (cos(Qy_rot)- twiss.alpha_y*sin(Qy_rot))*py
    return x1, px1, y1, py1

def rotation_with_detuners(particles_tunes_x, particles_tunes_y, Qx_rot, Qy_rot, twiss, x, px, y, py):
    Qx_new = Qx_rot+ 2*np.pi*particles_tunes_x
    Qy_new = Qy_rot + 2*np.pi*particles_tunes_y
    x1 = (np.cos(Qx_new) + twiss.alpha_x*np.sin(Qx_new))*x + twiss.beta_x*np.sin(Qx_new)*px
    px1 = -twiss.gamma_x*np.sin(Qx_new)*x+ (np.cos(Qx_new)-twiss.alpha_x*np.sin(Qx_new))*px
    y1 = (np.cos(Qy_new) + twiss.alpha_y*np.sin(Qy_new))*y + twiss.beta_y*np.sin(Qy_new)*py
    py1 = -twiss.gamma_y*np.sin(Qy_new)*y+ (np.cos(Qy_new)- twiss.alpha_y*np.sin(Qy_new))*py
    
    return x1, px1, y1, py1

def octupole_map(k3, x, px, y, py):
    k3_kick = k3/6.
    x1 = x
    px1 = px - k3_kick*(x**3-3*x*y**2)
    y1 = y
    py1 = py - k3_kick*(y**3-3*y*x**2)
    return x1, px1, y1, py1
    
def noise_map(Delta, sigmapx, sigmapy, x, px, y, py) :
    z_px = random.gauss(0,1) 
    x1 = x
    px1 = px + Delta*sigmapx*z_px
    y1 = y
    py1 = py + Delta*sigmapy*z_px
    return x1, px1, y1, py1

def BB_4D_map(xi, sigmax, sigmapx, sigmay, sigmapy, x, px, y, py):    
    # Note: The formula of Lebedev for the BB kick is for normalised coordinates. 
    # As here physical coordinates are used, the BB parameters and the kick itself need 
    # to be normalised with the sigmax and sigmapx. 
    
    r_sq = x**2/(sigmax**2) + y**2/(sigmay**2) # r_sq is r^2 uncomment when V motion is included

    x1 = x
    factor_x = ((8.*np.pi*xi*x1/sigmax)/r_sq)*sigmapx
    px1 = px + factor_x*(1.-np.exp(-r_sq/2.))
   
    y1 = y
    factor_y = ((8.*np.pi*xi*y1/sigmay)/r_sq)*sigmapy # uncomment when V motion is included
    py1 = py + factor_y*(1.-np.exp(-r_sq/2.))
   
    return x1, px1, y1, py1

def feedback_system_map(gain, sigmapx, sigmapy, x, px, y, py):
    px_average, py_average = np.nanmean(px/sigmapx), np.nanmean(py/sigmapy)
    
    x1 = x
    px1 = px - gain*sigmapx*px_average
    y1 = y
    py1 = py - gain*sigmapy*py_average
    return x1, px1, y1, py1

def aperture_limits_x_px_y_py(max_aperture_value, x, px, y, py):
    # Function to manage particle losses. An aperture is defined as a condition on the phase space coordinates.
    # If one of x, px, y, py for each particle doesn't meet the condition is considered lost and its coordinates are set as NaN.
    x1, px1, y1, py1 = x, px, y, py # type: nd.array()
    for i in range(len(x1)): # all the arrays have the same length, x1 was chosen arbitrary
        if x1[i]>max_aperture_value or px1[i]>max_aperture_value or y1[i]>max_aperture_value or py1[i]>max_aperture_value:
            x1[i] = px1[i]= y1[i] = py1[i] = np.nan 
            
    return x1, px1, y1, py1

def amplitude_detuning(k3_equivalent_x, k3_equivalent_y, twiss, x, px, y, py):
    # Detuning with amplitude , see documentation apo pyheadtail
    x1, px1, y1, py1 = x, px, y, py # type: nd.array()
    
    x1_norm, px1_norm = x1/sqrt(twiss.beta_x), px1*sqrt(twiss.beta_x) # later on include alpha and gamma to be correct in every case
    Jx = (x1_norm**2 + px1_norm**2)/2. # actions for each particle 
    
    y1_norm, py1_norm = y1/sqrt(twiss.beta_y), py1*sqrt(twiss.beta_y) # later on include alpha and gamma to be correct in every case
    Jy = (y1_norm**2 + py1_norm**2)/2. # actions for each particle 
        
    
    a_xx = (1/(8*np.pi))*k3_equivalent_x*(twiss.beta_x**2) # the detuning coefficient
    a_yy = (1/(8*np.pi))*k3_equivalent_y*(twiss.beta_y**2) # the detuning coefficient
    
    particles_tunes_x = a_xx*Jx/2. # DQ(Jx)
    particles_tunes_y = a_yy*Jy/2. # DQ(Jy)
    
    return particles_tunes_x, particles_tunes_y