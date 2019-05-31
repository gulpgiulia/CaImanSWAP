# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from pylab import *
from collections import OrderedDict
import matplotlib.pyplot as plt
# Needed to handle exits properly
import sys, traceback
import math
import os
from joblib import Parallel, delayed

import SetData
#from memory_profiler import LogFile
#------------------------------PARAMETER DEFINITION-----------------------------
DIM_X = int(SetData.DIM_X/SetData.MACRO_PIXEL_DIM)
DIM_Y = int(SetData.DIM_Y/SetData.MACRO_PIXEL_DIM)

# Check if a directory for 'Results' exists; if note, create it
#if os.path.exists(SetData.ANALYSIS_DIR  + 'RESULTS/')== False:
#    os.makedirs(SetData.ANALYSIS_DIR +'RESULTS/')
#if os.path.exists(SetData.ANALYSIS_DIR  + SetData.IMG_TYPE +'data_analysis/')== False:
#    os.makedirs(SetData.ANALYSIS_DIR  + SetData.IMG_TYPE +'data_analysis/')
#if os.path.exists(SetData.ANALYSIS_DIR +'RESULTS/'+SetData.IMG_TYPE +'/') == False:
#    os.makedirs(SetData.ANALYSIS_DIR +'RESULTS/'+SetData.IMG_TYPE +'/')

ResultsDir = SetData.ANALYSIS_DIR  + SetData.IMG_TYPE +'data_analysis/'
if os.path.exists(ResultsDir)== False:
    os.makedirs(ResultsDir)
print 'ResultsDir: ', ResultsDir

#-------------------------------------------------------------------------------
def excitability(params, exc, height, width, im_grid):
    #excitability = []
    # Here we define the image grid

    for i in range(0, len(params)):
        col = i % width
        row = i // width
        #col = i % height
        #row = i // height
        for j in range (0, len(params[i])):
            exc.append(params[i][j][0])
            #im_grid[row, col] = params[i][j][0]/len(params[i])
            im_grid[row, col] = params[i][j][0]


    return exc, im_grid
#------------------------------------
def readParams(filename):
    pix = []

    with open(filename, "r") as file:
        for line in file:
            if len(line) > 1:
                #print [[(i) for i in x.split(" ")] for x in line[:-2].split(";")]
                #print ciao
                pix.append([ np.array([float(i) for i in x.split(" ")]) for x in line[:-2].split(";")])
            else:
                pix.append([])

    file.close()
    return pix
#------------------------------------
# Stack Overflow algorithms
def bisection(array,value):

    n = len(array)
    if (value < array[0]):
        return -1
    elif (value > array[n-1]):
        return n
    jl = 0 # Initialize lower
    ju = n-1 # and upper limits.
    while (ju-jl > 1): # If we are not yet done,
        jm=(ju+jl) >> 1 # compute a midpoint with a bitshift
        if (value >= array[jm]):
            jl=jm # and replace either the lower limit
        else:
            ju=jm # or the upper limit, as appropriate.
        # Repeat until the test condition is satisfied.
    if (value == array[0]):# edge cases at bottom
        return 0
    elif (value == array[n-1]): # and top
        return n-1
    else:
        return jl
#------------------------------------
def min_to_txt(min_time, path):
    file = open(path + 'min_time.txt', 'w')
    for i in range(0, len(min_time)):
        file.write(str(i) + ' ' + str(min_time[i]) + '\n')
    file.close()
#------------------------------------
def lowhigh_loading():
    print 'Lowcut and highcut loading...'
    lowcut_tot = []
    highcut_tot = []
    zone_tot = []
    for set in SetData.num_measures:
        SET_NUMBER = str(set)
        print'set=' + SET_NUMBER
        #np_data = np.zeros(10)
        filename = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/' + 'low_high.txt'
        if os.path.isfile(filename)== False:
            print('ERROR: low_high file does not exist. EXIT')
            sys.exit()
            #sys.exit('ERROR: low_high does not exists' )
        np_data = np.genfromtxt(filename, delimiter=';',comments= '#', dtype = str)
        lowcut_set = np_data[0][0].split(',')
        highcut_set = np_data[1][0].split(',')
        zone_set = np_data[2][0].split(',')
        del np_data
        #print zone_set
        for i in range(0,len(lowcut_set)):
            try:
                lowcut_set[i] = float(lowcut_set[i])
            except:
                del lowcut_set[i]
        for i in range(0,len(highcut_set)):
            try:
                highcut_set[i] = float(highcut_set[i])
            except:
                del highcut_set[i]
        #print highcut_set
        for i in range(0, len(zone_set)):
            #print zone_set[i]
            if zone_set[i] == '':
                del zone_set[i]
        #print zone_set
        lowcut_tot.append(lowcut_set)
        highcut_tot.append(highcut_set)
        zone_tot.append(zone_set)
    del lowcut_set, highcut_set, zone_set
    for i in range(0, len(lowcut_tot)):
        if len(lowcut_tot[i]) == 0:
            del lowcut_tot[i]
    for i in range(0, len(highcut_tot)):
        if len(highcut_tot[i]) == 0:
            del highcut_tot[i]
    for i in range(0, len(zone_tot)):
        if len(zone_tot[i]) == 0:
            del zone_tot[i]
    #print zone_tot
    print 'Lowcut and highcut loaded!'
    return lowcut_tot, highcut_tot, zone_tot
#===============================================================================

def Data_analysis(z, lowcut_tot, highcut_tot, zone_tot, DIM_X, DIM_Y):

    print 'Computing ' + z + '...'
    # Create the directory where results are going to be saved
    #SavePath = SetData.ANALYSIS_DIR +'RESULTS/' + SetData.IMG_TYPE + '/' + z + '/'
    SavePath = ResultsDir
    if not os.path.exists(SavePath):
        os.makedirs(SavePath)

    #variables setting
    params = []
    exc = []
    im_grid_exc = np.zeros((DIM_X, DIM_Y))

    im_grid = np.zeros((DIM_X, DIM_Y))
    #mean_velocity = []
    #X_mouses_vel = np.zeros(len(SetData.num_measures)*200)
    #Y_mouses_vel = np.zeros(len(SetData.num_measures)*200)
    modulo_vel = []

    for index, set in enumerate(SetData.num_measures):

        for elem in range(0,len(zone_tot[index])):

            if zone_tot[index][elem] == z:
                #print 'entrato'
                # Computation of the excitability loading the file params where fit parameters are saved
                SET_NUMBER = str(set)
                lowcut = lowcut_tot[index][elem]
                highcut = highcut_tot[index][elem]
                #============================ EXCITABILIY ============================
                #print '         computing  exictability ...'
                #params loading...
                #filename = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + '/t' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'params_'+ SetData.IMG_TYPE + '_' + str(set) + '[' + str(lowcut) + '_' + str(highcut) + ']'+'.txt'
		filename = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'params.txt'
                print 'File with parameters: ', filename

                if os.path.exists(filename)== False:
		    print('ERROR: params file does not exist. EXIT')
                    sys.exit()

                params = readParams(filename)

                # Computing excitability
                exc_temp, im_grid_exc_temp = excitability(params, exc, DIM_X, DIM_Y, im_grid_exc)
                exc += exc_temp
                im_grid_exc += im_grid_exc_temp
                del params
                #im_grid_exc = im_grid_exc/np.nanmax(im_grid_exc[:,:])
                #elem += 1

                #=============================== ORIGIN of WAVES ==========================
                #print '         waves surce computing... '
                #variables initialization
                #image grid is re-initialized

                #min_points loading
                filename = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'min_points.txt'
                print 'File with min_points: ', filename
                if os.path.isfile(filename)== False:
                    print('ERROR: min_points file does not exist. EXIT')
                    sys.exit()

                # Here we load the file
                with open(filename) as f:
                    # We built is as a dict with key = pixel-ID ; value = activation-time
                    trans_dict = OrderedDict()
                    last_up_time = -1
                    first_up_time = 1E9

                    last_id = -1
                    first_id = -1

                    for line in f:
                        pos = line.find('[')
                        key = int(line[:pos])
                        try:
                            value = [float(x) for x in line[pos + 1:-3].split(",")]

                            # We keep track of the first and last up transition
                            if last_up_time < max(value):
                                last_up_time = max(value)
                                last_id = np.argmax(value)
                                key_last = key

                            if first_up_time > min(value):
                                first_up_time = min(value)
                                first_id = np.argmin(value)
                                key_first = key


                        except ValueError as err:
                            value = []

                        trans_dict[key] = value

                print "        Dictionary is loaded"

                if trans_dict.keys()[-1] + 1 != DIM_X * DIM_Y:
                    print "\nError on initiFromFile: pixels-IDs do not math simulation dimentions"
                    print "Max pixels-ID is %d while Width x Height is %d\n" % (trans_dict.keys()[-1] + 1, DIM_X * DIM_Y)
                    sys.exit(-1)

                #print "        First up transition detected at time T = %lf for ID = %d" % (first_up_time, key_first)
                #print "        Last up transition detected at time T = %lf for ID = %d" % (last_up_time, key_last)

                # Here we build two arrays: a time array in which all activation times for all
                # the pixels are stored (and ordered) and an idx-array in which the corresponding
                # sorted pixels idx values are saved.

                up_times = []
                up_idx = []
                for key in trans_dict.keys():
                    up_times.extend(trans_dict[key])
                    up_idx.extend([key for i in range(len(trans_dict[key]))])
                del trans_dict

                # Now we order our time and idx arrays
                ord_up_times = [x for x,_  in sorted(zip(up_times, up_idx))]
                ord_idx = [y for _, y in sorted(zip(up_times, up_idx))]
                del up_times
                del up_idx

                # Here we manually define a set of time for wave segmentation
                # Loading BeginTime.txt
                #timefile = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + '/t' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'BeginTime_' + SetData.IMG_TYPE + '_t' + str(set)  + '.txt'
                timefile = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'BeginTime.txt'
                if os.path.isfile(filename)== False:
                    print('ERROR: BeginTime file does not exist. EXIT')
                    sys.exit() 
                    #sys.exit('ERROR: BeginTime does not exists' )
                #loading
                f = open(timefile, 'r')
                wave_times = []
                for line in f:
                    wave_times.append(float(line))
                f.close()

                # Here we find the closest elements of ord_time to wave_times[i]
                for trigger in wave_times:
                    j_match = bisection(ord_up_times, trigger)

                    # Here we grab the corresponding idx of pixels
                    source_idx = np.array(ord_idx[j_match : j_match + SetData.k_cluster])

                    # Here we deduce the rows and colums idxs
                    colums = source_idx % DIM_Y
                    rows = source_idx // DIM_Y
                    #rows = source_idx % DIM_X
                    #colums = source_idx // DIM_X
                    #for row, col in zip(rows, colums):
                    #    im_grid[row][col] += 1
                    for row in rows:
                        for col in colums:
                            im_grid[row][col] += 1
                    del colums
                    del rows

                #===========================WAVE VELOCITY===========================
                #print '         waves velocity computing...'
                #load end_times of each wave
                # Loading of EndTime
                #timefile = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + '/t' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + '/EndTime_' + SetData.IMG_TYPE + '_t' + str(set) + '.txt'
                timefile = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + '/EndTime.txt'                
                if os.path.isfile(filename)== False:
                    print('ERROR: EndTime file does not exist. EXIT')
                    sys.exit()

                end_times = []
                f = open(timefile, 'r')
                for line in f:
                    end_times.append(float(line))

                X_velocity = np.zeros(len(wave_times))
                Y_velocity = np.zeros(len(wave_times))

                #ciclo su tutte le onde di una data raccolta dati
                for w in range(0,len(wave_times)):

                    ord_up_times =[round (elem,6) for elem in ord_up_times]
                    #print ord_up_times[0:20]
                    begin = ord_up_times.index(wave_times[w])
                    end = ord_up_times.index(end_times[w])


                    T = ord_up_times[begin:end]
                    #print len(T)
                    indici = ord_idx[begin:end]

                    grid = np.zeros((DIM_Y,DIM_X))

                    #creo la griglia
                    for k in range(0,len(T)):
                        x = indici[k]%DIM_Y
                        y = indici[k]//DIM_Y
                        grid [x, y] = T[k]

                    contatore = 0
                    #velocity = 0
                    velocity = []
                    Tx_temp = 0
                    Ty_temp = 0
                    #begin_time = wave_times[w]
                    #begin_position = ord_idx[begin]
                    #Tx = [] #derivata in x
                    #Ty = [] #derivata in y
                    #flag = 0

                    for x in range(1,DIM_Y-1):
                        for y in range(1,DIM_X-1):


                            if ( grid[x+1,y] != 0  and grid[x-1,y] != 0 and grid[x,y+1] != 0  and grid[x,y-1] != 0) and ( grid[x+1,y] != nan  and grid[x-1,y] != nan and grid[x,y+1] != nan  and grid[x,y-1] != nan):
                                #print 'entrato!'
                                Tx_temp=((grid[x+1,y] - grid[x-1,y])/2)
                                Ty_temp=((grid[x,y+1] - grid[x,y-1])/2)
                                velocity.append(1/sqrt(Tx_temp**2+Ty_temp**2))

                                contatore +=1
                                #print 'contatore aumentato!'


                    if contatore != 0:
                        mean_velocity = mean(velocity)
                        #qui aggiustiamo le unità di misura
                        modulo_vel.append(mean_velocity*SetData.MACRO_PIXEL_DIM*SetData.PIXEL_SIZE/(SetData.SAMPLING_TIME))
                del velocity, Tx_temp, Ty_temp, mean_velocity, grid, wave_times, end_times


    print '     Exitability plotting...'
    #excitability map and histogram plotting
    #print im_grid_exc
    for i in range(DIM_X):
        for j in range(DIM_Y):
            if im_grid_exc[i,j] == 0:
                im_grid_exc[i, j] = np.nan
                im_grid[i,j] = np.nan


    #im_grid plotting (color map for excitability)
    im_grid_exc = im_grid_exc/np.nanmax(im_grid_exc[:,:])
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    plt.Normalize()
    plt.imshow(im_grid_exc, cmap='plasma')
    #name_fig = SetData.IMG_TYPE + '_' + z
    name_fig = z
    title_fig = 'ExcitabilityMap'
    plt.colorbar()
    #plt.title(title_fig)
    #plt.suptitle(name_fig, fontsize=16)
    #plt.ylabel('uhm')
    #plt.xlabel('s^{-2}')
    name_saved = SavePath + title_fig + '_' + name_fig + '.png'
    plt.savefig(name_saved)
    #plt.show()
    plt.close()
    del im_grid_exc
    #excitability histogram plotting

    a = min(exc)
    b = max(exc)
    bins = arange(a, b, (b-a)/50)
    hist, bin, bho = plt.hist(exc, bins = bins,  color = SetData.color, density =False)
    #plt.hist(exc, bins = bins,  color = c, density =False)
    histogram = np.zeros(len(hist)+1)
    for i in range(0,len(hist)):
        histogram[i] = hist[i]
    plt.close()
    plt.bar(bins, histogram/np.max(hist), width=(b-a)/50,  color = SetData.color)
    #name_fig = SetData.IMG_TYPE + '_' + z
    name_fig = z
    title_fig = 'Excitability'
    #plt.title(title_fig)
    #plt.suptitle(name_fig, fontsize=16)
    plt.ylabel('Number of minima')
    plt.xlabel(r'$s^{-2}$', usetex = True)
    name_saved = SavePath + title_fig + '_' + name_fig + '.png'
    plt.savefig(name_saved)
    #plt.show()
    plt.close()
    del histogram, hist, bin, bho, bins, a, b
    print '     Excitability plotted!'

    print '     Waves origin plotting...'

    #waves origin plotting
    #print np.nanmax(im_grid[:,:])
    im_grid = im_grid/np.nanmax(im_grid[:,:])
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    plt.imshow(im_grid, cmap='plasma')
    #plt.Normalize(vmin = -1, vmax = 1)
    #name_fig = SetData.IMG_TYPE + z
    name_fig = z
    title_fig = 'Waves_origin'
    plt.colorbar()
    #plt.title(title_fig)
    #plt.suptitle(name_fig, fontsize=16)
    #plt.ylabel('uhm')
    #plt.xlabel('uhm')
    name_saved = SavePath + title_fig + '_' + name_fig + '.pdf'
    plt.savefig(name_saved)
    #plt.show()
    plt.close()
    del im_grid

    print '     Wave origin plotted!'

    print len(modulo_vel)
    print '     Waves velocity plotting...'
    a = min(modulo_vel)
    b = max(modulo_vel)
    bins = arange(a, b, (b-a)/50)
    #plt.xlim(0,50)
    hist, bin, bho = plt.hist(modulo_vel, bins = bins,  color = SetData.color, density = False)
    histogram = np.zeros(len(hist)+1)
    for i in range(0,len(hist)):
        histogram[i] = hist[i]
    plt.close()
    plt.bar(bins, histogram/np.max(hist), width=(b-a)/50,  color = SetData.color)
    plt.ylabel('Number of waves')
    plt.xlabel('v(mm/s)')
    #name_fig = SetData.IMG_TYPE + z
    name_fig = z    
    title_fig = 'Waves_Velocity'
    name_saved = SavePath + title_fig + '_' + name_fig + '.pdf'
    plt.savefig(name_saved)
    #plt.show()
    plt.close()
    del histogram, hist, bin, bho, bins, a, b

    print '     Waves velocity plotted!'



    #Save Results
    filename = SavePath + 'Results' + z + '.txt'
    with open(filename, "w") as file:
        file.write('Mean excitability = %s' % mean(exc))
        file.write('\n')
        file.write('Excitability standard deviation = %s' % std(exc))
        file.write('\n')
        file.write('Mean velocity = %s' % mean(modulo_vel))
        file.write('\n')
        file.write('Velocity standard deviation = %s' % std(modulo_vel))

    print('Mean excitability = %s' % mean(exc))
    print('Excitability standard deviation = %s' % std(exc))
    print('Mean velocity = %s' % mean(modulo_vel))
    print('velocity standard deviation = %s' % std(modulo_vel))

    print 'Results plotted!'
    del modulo_vel, exc
    return 0

# -------------------------------------------------------------------
# -------------------------------------------------------------------

# -------- Load passband parameters (lowcut, highcut, zone) ---------
lowcut_tot, highcut_tot, zone_tot = lowhigh_loading()
# Create zone_elem
zone_elem = []

# ... for each spectrum zone to be analyzed:
for set in range(0,len(zone_tot)):
    for elem in zone_tot[set]:
        if not elem in zone_elem:
            zone_elem.append(elem)

# ---------------------------- DATA ANALYSIS ------------------------
#Data_analysis('Zone_2', lowcut_tot, highcut_tot, zone_tot, DIM_X, DIM_Y)
Parallel(n_jobs=SetData.nprocs)(delayed(Data_analysis)(z, lowcut_tot, highcut_tot, zone_tot, DIM_X, DIM_Y) for z in zone_elem)
#sys.stdout.write = LogFile('memory_profile_log')

print 'Data Analysis completed!'

