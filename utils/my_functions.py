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
