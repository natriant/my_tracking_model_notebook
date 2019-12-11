from dotted_dict import DottedDict
from my_elements import *


def my_machine(Qx_init, Qy_init, segments, k3, turn, bunch, twiss, rf, long_motion_parameters, flag_oct, flag_noise, flag_BB, flag_feedback, max_aperture_value, Delta = 0., ksi = 0., gain = 0., sigmax = 0., sigmapx = 0., sigmay=0., sigmapy = 0.):
    if turn == 1:
        print('flag_oct {}, flag_noise {}, flag_BB {}, flag_feedback {}'.format(flag_oct, flag_noise, flag_BB, flag_feedback))
        print('aperture limit is {} [m]'.format(max_aperture_value))
    for segment in range(segments):
        bunch.x, bunch.px, bunch.y, bunch.py = rotation(Qx_rot = Qx_init/segments, Qy_rot = Qy_init/segments, twiss = twiss, x = bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
        if flag_oct:
            bunch.x, bunch.px, bunch.y, bunch.py = octupole_map(k3/segments, x = bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
            if turn == 1 and (segment == segments-1):
                print('octupole kick', k3)
        if flag_noise and (segment == segments-1):
            if turn ==1 :
                print('Delta',Delta)
            bunch.x, bunch.px, bunch.y, bunch.py = noise_map(Delta, sigmapx, sigmapy, x = bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
        if flag_BB and (segment == segments-1):
            if turn ==1:
                print('ksi{}, sigmax{}, sigmapx{}, sigmay{}, sigmapy{}'.format(ksi, sigmax, sigmapx, sigmay, sigmapy))
            bunch.x, bunch.px, bunch.y, bunch.py = BB_4D_map(ksi, sigmax, sigmapx, sigmay, sigmapy, x= bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
        
        # Aperture limitations
        bunch.x, bunch.px, bunch.y, bunch.py = aperture_limits_x_px_y_py(max_aperture_value, x= bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
    
    
        if flag_feedback and (segment == segments-1):
            if turn ==1:
                print('g',gain)
            bunch.x, bunch.px, bunch.y, bunch.py = feedback_system_map(gain, sigmapx, sigmapy, x= bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
            
    # Longitudinal motion, as long as there is no dispersion, it can be seperetely.
    bunch.z, bunch.delta = drift(long_motion_parameters.eta, long_motion_parameters.C0, bunch.z, bunch.delta)
    bunch.z, bunch.delta = RF_map(long_motion_parameters.P0_C, long_motion_parameters.m0, long_motion_parameters.beta_0, rf.f_RF, long_motion_parameters.c, rf.lag, rf.V_RF, long_motion_parameters.E_0, bunch.z, bunch.delta)
            
    return bunch

def my_machine_with_detuners(Qx_init, Qy_init, k3_equivalent_x, k3_equivalent_y, turn, bunch, twiss, rf, long_motion_parameters, flag_detuner, flag_noise, flag_BB, flag_feedback, max_aperture_value, Delta = 0., ksi = 0., gain = 0., sigmax = 0., sigmapx = 0., sigmay=0., sigmapy = 0.):
    if turn == 1:
        print('flag_detuner {}, flag_noise {}, flag_BB {}, flag_feedback {}'.format(flag_detuner, flag_noise, flag_BB, flag_feedback))
        print('aperture limit is {} [m]'.format(max_aperture_value))
   
    # Rotation
    #if turn == 1: # first turn, use the original working point
    #    bunch.x, bunch.px, bunch.y, bunch.py = rotation(Qx_rot = Qx_init, Qy_rot = Qy_init, twiss = twiss, x = bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
    #else:    
    if flag_detuner:
        if turn == 1 :
                print('k3 equivalent x {}, k3 equivalent y {}'.format(k3_equivalent_x, k3_equivalent_y))
        particles_tunes_x, particles_tunes_y = amplitude_detuning(k3_equivalent_x, k3_equivalent_y, twiss,x = bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
        bunch.x, bunch.px, bunch.y, bunch.py = rotation_with_detuners(particles_tunes_x, particles_tunes_y, Qx_rot = Qx_init, Qy_rot = Qy_init, twiss = twiss, x = bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
        
    if flag_noise:
        if turn ==1 :
            print('Delta',Delta)
        bunch.x, bunch.px, bunch.y, bunch.py = noise_map(Delta, sigmapx, sigmapy, x = bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
   
    if flag_BB :
        if turn ==1:
            print('ksi{}, sigmax{}, sigmapx{}, sigmay{}, sigmapy{}'.format(ksi, sigmax, sigmapx, sigmay, sigmapy))
        bunch.x, bunch.px, bunch.y, bunch.py = BB_4D_map(ksi, sigmax, sigmapx, sigmay, sigmapy, x= bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
        
    # Aperture limitations
    bunch.x, bunch.px, bunch.y, bunch.py = aperture_limits_x_px_y_py(max_aperture_value, x= bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)
    
    
    if flag_feedback:
        if turn ==1:
            print('g',gain)
        bunch.x, bunch.px, bunch.y, bunch.py = feedback_system_map(gain, sigmapx, sigmapy, x= bunch.x, px = bunch.px, y = bunch.y, py=bunch.py)        
        
        
        
    # Longitudinal motion, as long as there is no dispersion, it can be seperetely.
    bunch.z, bunch.delta = drift(long_motion_parameters.eta, long_motion_parameters.C0, bunch.z, bunch.delta)
    bunch.z, bunch.delta = RF_map(long_motion_parameters.P0_C, long_motion_parameters.m0, long_motion_parameters.beta_0, rf.f_RF, long_motion_parameters.c, rf.lag, rf.V_RF, long_motion_parameters.E_0, bunch.z, bunch.delta)

    return bunch