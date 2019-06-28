# -*- coding: utf-8 -*-

#Library with all the functions to detect signal's minima and interpolate around them.

# Main differences from the previous version:
# Now "min_analysis" returns a dictionary whose elements are
# 1) 'min_time': a list of lists, eacn one containing the
# un-quantized minima time for a particoular pixel.
# 2) 'params': a list of arrays containing the quadratic
# fit parameters of every minimum in each pixel.
# eg. min_analysis['min_time'][998][3] would be the time at wich the
# third minimum of the 998th (19th row, 48th column) pixel occurs
# while min_analysis['params'][998][3] would be its quadratic fit parameters


#Sets output values of divisions as floating points, as in Python 3

from __future__ import division

#----------------------------------------------------

def find_minimum(clean_signal, t_min, t_max):

    min = []
    #iter on time
    for t in range(t_min, t_max):
        if ((clean_signal[t-1] > clean_signal[t])
        and (clean_signal[t] < clean_signal[t+1])):
            #min[] contains indexes and clean_signals values where it has a minimum
            min.append(t)

    return min

#-----------------------------------------------------

def parabola(x, a, b, c):
        return a*x**2 + b*x + c

#-----------------------------------------------------

def min_interpol(clean_signals, min_collection, points, threshold):
# Threshold defines the valor of excitability above which a minumun is identified as noise

    import numpy as np

    around_min = np.zeros(points)
    t = np.zeros(points)
    interpol_min_collection = []
    quadratic_params = []

    # Iterates
    for i in range(0, len(min_collection)):
        k = int(i/np.size(clean_signals,1))
        l = i%(np.size(clean_signals,1))
        pixel_minima = []
        pixel_params = []

        for j in min_collection[i][2]:

            #Sets t as an array of lenght 'points' around the local minimum time and ar_min as it's relative intensities
            t0 = j-int(points/2)

            for n in range(0, points):
                t[n] = t0 + n

                around_min[n] = clean_signals[k, l, int(t0 + n)]

                #Fits the local minimum with a parabola

            params = np.polyfit(t, around_min, 2)
            t_min = -params[1]/(2*params[0])

            # Adds parabola's minimum and parameters to the proper pixel's list

            if params[2] > threshold:
                pixel_params.append(params)
                pixel_minima.append(t_min)

        #Adds each pixel's list to a list of lists
        quadratic_params.append(pixel_params)
        interpol_min_collection.append(pixel_minima)

    #Returns a dictionary of two elements

    return {'min_time' : interpol_min_collection, 'params' : quadratic_params}


#------------------------------------------

def min_analysis(clean_signals, points, t_min, t_max, threshold = 0):
#Suggested values for a quadratic fit on a 1000 frames image set are
# point = 5, t_min = 3, t_max = 997

    import numpy as np

    min_collection = []

    #iter on macro-pixels

    for x in range(0,np.size(clean_signals,0)):
        for y in range(0,np.size(clean_signals,1)):

            clean_signals_xy = clean_signals[x][y]
            time = find_minimum(clean_signals_xy, t_min, t_max)

            min_collection.append([x, y, time])

    #Be careful the function returns a DICTIONARY in which the element
    # 'min_time' is a collection of all the interpolated minima
    # 'params' are the parameters of each parabola
    return min_interpol(clean_signals, min_collection, points, threshold)

#-----------------------------------------------

#Creates a .mat variable containing min_collection on chosen path

def min_to_matlab(min_collection, number, time, name):

    import scipy.io

    #REMEMBER to change the path manually here!
    #path = 'C:/Users/Marco/Desktop/WaveScales/SlowStudents/Analysis/Topo_' + number + '/t' + time + '/'

    scipy.io.savemat(name , mdict={name: min_collection})

#-------------------------------------------------------

#Produces a plot of the minima's wave at a particular time t
#only considering the pixels for which a minimum occurs in a given
#range 'dist' around t
#The value of each pixel, if a minimum is near enough, is
#the value that the interpolated parabola takes at time t

# Typically params will be something as " min_analysis['params'] "
# min_time something as " min_analysis['min_time'] "
def plot_wave(min_time, params, t, dist, rows = 40, columns = 50, save = False, ret = False, time_depth = False):
# If plot = True the function saves the image under the name "wave_'t'_'dist'.png"
# In the same path where this file is

# If time_depth = True the only information used in calculating the depth
# of minima is the temporal distance for the observational time t
    import matplotlib.pyplot as plt
    import numpy as np

    min_signal = np.zeros((rows, columns))

    for x in range(0, rows):
        for y in range(0, columns):

            index = x*columns + y

            if params[index] == []:
                min_signal[x, y] = np.nan

            else:
                time_dist = [abs(z - t) for z in min_time[index]]
                i = np.argmin(time_dist)

                # Mi assicuto che il minimo ci sia stato prima del tempp t di osservazione
                # seguo solo neuroni che si stanno accendendo, NON spegenendo
                # Inoltre l'accensione deve essere occorsa in una finestra temporale 'dist'
                # precedente al tempo t
                if (min_time[index][i] < t) & (time_dist[i] < dist):
                    if time_depth:
                    # Divido le profondità in 5 fasce temporali precedenti a t
                    # se voglio N fasce ho 20 -> 100/N e 0.5 -> 0.N
                        rel_dist = (time_dist[i]/dist)
                        min_signal[x, y] = - 1 + rel_dist
                        # min_signal[x, y] = (int(int(rel_dist*100)/20)) / 10 - 0.5
                    else:
                    #Voglio allineare i minimi delle parabole alla stessa altezza
                    #Calcolo il valore che assume al minimo
                        y_min = parabola(min_time[index][i], *params[index][i])
                    # Traslo il valore che la parabola assume in t di modo che
                    # il suo minimo sia in -0.005 (vlore arbitrario)
                        min_signal[x, y] = parabola(t, *params[index][i]) - (0.005 + y_min)
                else:
                    min_signal[x, y] = 0

    if ret:
        return min_signal
    elif save:
        plt.imshow(min_signal, interpolation = 'nearest')
        plt.text(0.05, 2, 'Observation time = ' + "{:.2f}".format(t * 0.04) + 's')
        plt.text(0.05, 39, 'Time window = ' + str(dist * 0.04) + 's')
        plt.colorbar()
        plt.savefig('wave_' + str(t) + '_' + str(dist) + '.png')
    else:
        plt.imshow(min_signal, interpolation = 'nearest')
        plt.text(0.05, 2, 'Observation time = ' + "{:.2f}".format(t * 0.04) + 's')
        plt.text(0.05, 39, 'Time window = ' + str(dist * 0.04) + 's')
        plt.colorbar()
