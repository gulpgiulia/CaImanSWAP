#===============================================================================
#                                       MAIN
#===============================================================================

import os
import numpy as np
from pylab import *
from scipy.signal import butter, lfilter, filtfilt, freqz
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy.ma as ma

import json
import os.path
from joblib import Parallel, delayed

#import psutil

import skimage
from skimage import data, io, filters, measure
from skimage import img_as_float, img_as_uint

from frequency_analysis import bandpass_filter, multiplot
from new_min_analysis import min_analysis, min_to_matlab, plot_wave

import SetData

#------------------------------FUNCTION DEFINITION-----------------------------
def Spectrum_plotting(img_spectrum_freq, spectrum_mean, set, path, lowcut = None, highcut = None, show = True, save = False):
    print('     Plotting the spectrum...')

    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    plt.plot(img_spectrum_freq[1:SetData.spectrum_time]*1000, (spectrum_mean[1:SetData.spectrum_time, 0, 0].real)**2, linewidth=0.45, c = SetData.color)
    #name_fig = SetData.IMG_TYPE + '_' + str(set)
    name_fig = 'plot'
    title_fig = 'Image_spectrum'
    #plt.suptitle(name_fig, fontsize=16)

    if lowcut is not None and highcut is None:
        plt.xlim(lowcut, )
        plt.title('Filter:' + str(lowcut) + '- ?' +'Hz')

    if lowcut is None and highcut is not None:
        plt.xlim (0., highcut)
        plt.title('Filter:' + '? ' + '-' + str(highcut) +'Hz')


    if lowcut is not None and highcut is not None:
        plt.xlim (lowcut, highcut)
        plt.title('Filter:' + str(lowcut) + '-' + str(highcut) +'Hz')

    plt.yscale(value='log')

    plt.ylabel(r'$\mathcal{A}$')
    plt.xlabel(r'$\nu (Hz)$')
    name_saved = path + title_fig + 'Log_' + name_fig + '.png'
    print('Spectrum_plotting: name_saved', name_saved)
    if show:
        plt.show()
    if save:
        plt.savefig(name_saved)
    plt.close()


    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    plt.plot(img_spectrum_freq[1:SetData.spectrum_time]*1000, (spectrum_mean[1:SetData.spectrum_time, 0, 0].real)**2, linewidth=0.45, c = SetData.color)
    #name_fig = SetData.IMG_TYPE + '_' + str(set)
    name_fig = 'plot'
    title_fig = 'Image_spectrum'

    if lowcut is not None and highcut is None:
        plt.xlim(lowcut, )
        plt.title('Filter : ' + str(lowcut) + '- ?' +'Hz')
        y_lim = find_nearest(img_spectrum_freq[1:SetData.spectrum_time], lowcut / 1000.)
        plt.ylim(0., np.max((spectrum_mean[y_lim : SetData.spectrum_time, 0, 0].real)**2))

    if lowcut is None and highcut is not None:
        plt.xlim (0, highcut)
        plt.title('Filter :'  + '? ' + '-' + str(highcut) +'Hz')
        y_lim = find_nearest(img_spectrum_freq[1:SetData.spectrum_time], highcut / 1000.)
        plt.ylim(0., np.max((spectrum_mean[1:y_lim, 0, 0].real)**2))

    if lowcut is not None and highcut is not None:
        plt.xlim (lowcut, highcut)
        plt.title('Filter : ' + str(lowcut) + '-' + str(highcut) +'Hz')
        y_lim_low = find_nearest(img_spectrum_freq[1:SetData.spectrum_time], lowcut / 1000.)
        y_lim_high = find_nearest(img_spectrum_freq[1:SetData.spectrum_time], highcut / 1000.)
        plt.ylim(0., np.max((spectrum_mean[y_lim_low:y_lim_high, 0, 0].real)**2))


    #plt.ylim(0., np.max((spectrum_mean[1:SetData.spectrum_time, 0, 0].real)**2))
    plt.ylabel(r'$\mathcal{A}$')
    plt.xlabel(r'$\nu (Hz)$')

    if show:
        plt.show()
    if save:
        name_saved = path + title_fig + '_' + name_fig + '.png'
        plt.savefig(name_saved)

    plt.close()
    print('     Spectrum plotted!')
#------------------------------------------------------------
def RawSignal_plotting(img_coll_norm, set, lowcut, highcut, path):
    print('     Plotting the raw signal...')
    to_plot = np.zeros(SetData.TIME_SAMPLES)
    for i in range(0, SetData.TIME_SAMPLES):
        to_plot[i] = img_coll_norm[i, SetData.PIXEL_SAMPLE_x, SetData.PIXEL_SAMPLE_y]
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    x = np.linspace(0, SetData.TIME_SAMPLES*SetData.SAMPLING_TIME, SetData.TIME_SAMPLES)
    plt.plot(x, to_plot, linewidth = 0.45, c = SetData.color)
    #name_fig = SetData.IMG_TYPE + '_' + str(set)
    name_fig = 'plot'
    title_fig = 'Raw_Signal'
    plt.title('Filter:' + str(lowcut) + '-' + str(highcut) +'Hz')
    plt.ylabel(r'$\frac{\mathcal{F}}{\mathcal{F}_{max}}$')
    plt.xlabel('t(s)')
    name_saved = path + title_fig +'_'+ name_fig + '.png'
    plt.savefig(name_saved)
    plt.close()

    print('     Raw signal plotted!')
#------------------------------------------------------------
def CleanSignal_plotting(clean_signals_norm, set, lowcut, highcut, path):
    print('     Plotting the cleaned signal...')

    to_plot = np.zeros((SetData.TIME_SAMPLES))
    for i in range(0, SetData.TIME_SAMPLES):
        to_plot[i] = clean_signals_norm[SetData.PIXEL_SAMPLE_x, SetData.PIXEL_SAMPLE_y, i]

    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')

    x = np.linspace(0, SetData.TIME_SAMPLES*SetData.SAMPLING_TIME, SetData.TIME_SAMPLES)
    plt.plot(x, to_plot, linewidth = 0.45, c = SetData.color)
    #name_fig = SetData.IMG_TYPE + '_' + str(set)
    name_fig = 'plot'
    title_fig = 'Cleaned_signal'
    plt.title('Filter:' + str(lowcut) + '-' + str(highcut) +'Hz')
    #plt.suptitle(name_fig, fontsize=16)
    #plt.title(title_fig)
    #plt.suptitle(name_fig, fontsize=16)
    #plt.ylim(ymin=-10**16, ymax = 10**16)
    #plt.yscale(value='log')
    plt.xlim(xmin=(SetData.TIME_SAMPLE_min)*SetData.SAMPLING_TIME, xmax=(SetData.TIME_SAMPLE_max)*SetData.SAMPLING_TIME)
    plt.ylabel(r'$\frac{\mathcal{F}}{\mathcal{F}_{max}}$')
    plt.xlabel('t(s)')
    name_saved = path + title_fig + '_' + name_fig + '.png'
    print('CleanSignal_plotting: name_saved', name_saved)

    plt.savefig(name_saved)
    #plt.show()
    plt.close()

    print('     Cleaned Signal plotted!')
#------------------------------------------------------------
def saveParams(filename, params):
    with open(filename, "w") as file:
        for i, pix in enumerate(params):
            #file.write(str(i)+':')
            for min in pix:
                msg = str(min[0]) + " " + str(min[1]) + " " + str(min[2]) + ";"
                file.write(msg)

            file.write("\n")

    file.close()
#------------------------------------------------------------
def saveMin(filename, min_time):
    with open(filename, "w") as file:
        for i, pix in enumerate(min_time):
            file.write(str(i)+' ')
            file.write(str(pix))
            #for min in pix:
                #msg = str(min[0]) + " " + str(min[1]) + " " + str(min[2]) + ";"
            #    file.write(min)

            file.write("\n")

    file.close()
#------------------------------------------------------------
def find_nearest(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()
#------------------------------------------------------------
def getZone():
    zone_ans = True;
    print("Which ZONE does this band belong to?")

    while zone_ans:
        zone = safeRawInput('Insert Zone name: ')

        isZone = False

        for SD_zone in SetData.zone:
            if SD_zone is not None:
                if zone in SD_zone:
                    isZone = True
        '''
        if isZone:
            zone_ans = False

        else:
            print(zone + " is not currently a known Zone. Do you want to add it? (Y/N)")
            add = input ()
            if add == 'Y':
                zone_ans = False
        '''
        zone_ans=False
    return zone
#------------------------------------------------------------
def safeInput(msg):
    ans = True
    while ans:
        try:
            temp = input(msg)
            ans = False

        except SyntaxError:
            print("Invalid input. Please re-insert.")

        except NameError:
            print("Invalid input. Please re-insert.")

    return temp
#------------------------------------------------------------
def safeRawInput(msg, valid_ans = None):
    ans = True
    while ans:
        try:
            temp = input(msg)

            if valid_ans is not None:
                if not (temp in valid_ans):
                    raise SyntaxError

            ans = False

        except SyntaxError:
            print("Invalid input. Please re-insert.")

        except NameError:
            print("Invalid input. Please re-insert.")

    return temp
#------------------------------------------------------------
print('---------------------------Image_analysis-------------------------------')

DIM_X_REDUCED = int(SetData.DIM_X / SetData.MACRO_PIXEL_DIM)
DIM_Y_REDUCED = int(SetData.DIM_Y / SetData.MACRO_PIXEL_DIM)

lowcut_tot = []
highcut_tot = []
zone_tot = []

#===============================================================================

#-------------------------------------------------------------------------------

def Spectrum_evaluation(set):

    SET_NUMBER = str(set)
    print('Set number = ', SET_NUMBER)

    #-----------------------------IMAGE SET LOADING---------------------------------
    # In this section, images are loaded from a txt file
    filename = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/' + 'initialized_images_t' + SET_NUMBER + '.txt'
    #filename = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + '/t' + SET_NUMBER + '/' + 'initialized_images_' + SetData.IMG_TYPE + '_t' + SET_NUMBER + '.txt'
    print('input filename: ', filename)

    # check the file existence
    if os.path.isfile(filename)== False:
        print('ERROR: initialized_images does not exist' )
        sys.exit(-1)

    #print('     Loading data...')
    Images2D = np.loadtxt(filename)
    img_collection_reduced = np.reshape(Images2D,
                                        (SetData.TIME_SAMPLES+1,
                                         DIM_X_REDUCED,
                                         DIM_Y_REDUCED),
                                        order='C')
    del Images2D

    # ad hoc settings for specific datasets (isoflurane)
    if SetData.IMG_TYPE == 'ISO_180322':
        img_collection_reduced[:,29,10] = np.nan
        img_collection_reduced[:,42,22] = np.nan
        img_collection_reduced[:,26,15] = np.nan

    #plt.matshow(img_collection_reduced[0])
    #plt.show()

    #filename = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + '/t' + SET_NUMBER + '/' + 'background_' + SetData.IMG_TYPE + '_t' + SET_NUMBER + '.txt'
    #Background2D = np.loadtxt(filename)
    #img_background = np.reshape(Background2D, ((DIM_X_REDUCED, DIM_Y_REDUCED)),order='C')
    #img_collection_reduced = img_collection_reduced - img_background
    #print('     Background loaded!')

    #------------------------------SIGNAL SPECTRUM------------------------------

    path = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/'
    print('path = ', path)

    #print('     Evaluating the spectrum...')

    img_coll_norm = np.zeros((SetData.TIME_SAMPLES, DIM_X_REDUCED, DIM_Y_REDUCED))

    # *** HERE - NORMALIZATION of the dataset with the maximum value ***
    img_coll_norm = img_collection_reduced/np.nanmax(img_collection_reduced)
    # The fluctuations of the Ca signal are normalized to the max fluorescence of each pixel,
    # in contrast to the standard normalization to the baseline or median of the signal.
    # [almost by definition the max fluorescence is more affected by noise than any reasonable estimate of the signal baseline]
    # --> REPLACE the normalization factor with the MEDIAN
    #img_coll_norm = img_collection_reduced/np.nanmedian(img_collection_reduced)

    del img_collection_reduced

    # Next we compute the Fourier transform of the input along the temporal axis
    img_spectrum = np.fft.rfftn(img_coll_norm, axes = [0])
    img_spectrum_freq = np.fft.rfftfreq(img_coll_norm.shape[0], d = 1. / SetData.SAMPLING_TIME)

    # We take the mean of the spectrum
    spectrum_mean = measure.block_reduce(img_spectrum, (1, np.size(img_coll_norm,1), np.size(img_coll_norm,2)), np.nanmean)
    #print('spectrum_mean size: ', size(spectrum_mean))
    #print('img_spectrum size: ', size(img_spectrum))
    #print(spectrum_mean)
    del img_spectrum

    #print('     Spectrum evaluated!')

    #------------ Save the total spectrum
    #print('     Plotting the total spectrum...')
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    plt.plot(img_spectrum_freq[1:SetData.spectrum_time]*1000, (spectrum_mean[1:SetData.spectrum_time, 0, 0].real)**2, linewidth=0.45, c = SetData.color)
    #name_fig = SetData.IMG_TYPE + 't' + str(set)
    name_fig = 'mean-across-pixels'
    title_fig = 'ImageSpectrum'
    plt.title(title_fig+'LogScale')
    plt.yscale(value='log')
    plt.ylabel(r'$\mathcal{A}$')
    plt.xlabel(r'$\nu (Hz)$')
    name_saved = path + title_fig + 'Log_' + name_fig + '.png'
    #plt.show()
    plt.savefig(name_saved)
    plt.close()

    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    plt.plot(img_spectrum_freq[1:SetData.spectrum_time]*1000, (spectrum_mean[1:SetData.spectrum_time, 0, 0].real)**2, linewidth=0.45, c = SetData.color)
    #name_fig = SetData.IMG_TYPE + '_' + str(set)
    name_fig = 'mean-across-pixels'
    title_fig = 'ImageSpectrum'
    plt.title(title_fig)
    plt.ylabel(r'$\mathcal{A}$')
    plt.xlabel(r'$\nu (Hz)$')
    name_saved = path + title_fig + '_' + name_fig + '.png'
    plt.savefig(name_saved)
    plt.close()
    #print(img_spectrum_freq)
    returning = [spectrum_mean, img_spectrum_freq, img_coll_norm]

    del spectrum_mean
    del img_spectrum_freq
    del img_coll_norm

    return returning
#------------------------------------------------------------

print('Loading pre-processed data and evaluating the frequency spectrum...')
spectrum_mean_measures = []
img_spectrum_freq_measures = []
img_coll_norm_measures = []
returning = []
returning.extend(Parallel(n_jobs=int(SetData.nprocs))(delayed(Spectrum_evaluation)(set) for set in SetData.num_measures))
#returning.extend((Spectrum_evaluation)(set) for set in SetData.num_measures) # NO PARALLE EXEC
#print('Ret:' + str(len(returning)))

for i in range(0, len(SetData.num_measures)):
    spectrum_mean_measures.append(returning[i][0])
    img_spectrum_freq_measures.append(returning[i][1])
    img_coll_norm_measures.append(returning[i][2])
    #del returning[i]

print('Spectrum evaluated!')

#===============================================================================

#-------------------------------------------------------------------------------

for index, set in enumerate(SetData.num_measures):
    #lowcut_set[index] = np.zeros(max(len(SetData.lowcut), len(SetData.highcut)))
    #highcut_set[index] = np.zeros(max(len(SetData.lowcut), len(SetData.highcut)))

    SET_NUMBER = str(set)
    print('Set number = ', SET_NUMBER)
    img_coll_norm = img_coll_norm_measures[index]
    spectrum_mean = spectrum_mean_measures[index]
    img_spectrum_freq = img_spectrum_freq_measures[index]

    #-----------------------------------------------Notch filter---------------

    answer = safeRawInput('Do you want to apply a notch filter? (Y/N):', valid_ans = ['Y', 'N'])
    if answer == 'Y':
        f0 = safeInput('Frequency to be eliminated:')

    #-----------------------------------------------Lowcut and highcut selection

    # cycle over the zones of interest of the spectrum
    if SetData.lowcut[index] is None:
        answer = 'Y'
        cont = 0
        lowcut_set = []
        highcut_set = []
        zone_set = []

        while answer == 'Y':
            # plot of the total spectrum
            path = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + '/t' + SET_NUMBER + '/'
            Spectrum_plotting(img_spectrum_freq, spectrum_mean, set, path)

            #input
            ans = True;
            while ans:
                ans = False
                lowcut = safeInput('Lowcut: ')
                highcut = safeInput('Highcut: ')

                if (lowcut >= highcut):
                    ans = True
                    print("WARNING: Lowcut = " + str(lowcut) + " Highcut = " + str(highcut))
                    print("Lowcut is larger than highcut. Please re-insert both values.")

            zone = getZone()

            lowcut_set.append(lowcut)
            highcut_set.append(highcut)
            zone_set.append(zone)
            #lowcut_set[index][cont] = lowcut

            # path were things are saved
            path = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut) + ']Hz/'
            if os.path.exists(path)== False:
                os.makedirs(path)
            # plotting...
            Spectrum_plotting(img_spectrum_freq, spectrum_mean, set, path, lowcut, highcut, show = False, save = True)
            #del spectrum_mean
            #del spectrum_mean_measures[index]

            RawSignal_plotting(img_coll_norm, set, lowcut, highcut, path)

            #------------------------------- BANDPASS FILTER ---------------------------
            # Clean the signal through a bandpass-pass filter
            print('     Filtering the signal...')

            #clean_signals = bandpass_filter(img_collection_reduced, SetData.order, SetData.fs, lowcut, highcut)
            clean_signals_norm = bandpass_filter(img_coll_norm, SetData.order, SetData.fs, lowcut, highcut)

            CleanSignal_plotting(clean_signals_norm, set, lowcut, highcut, path)
            #del img_coll_norm

            #--------------------------- INTERPOLATION of the MINIMA -------------------
            # In this section, the minima are found and interpolated with a quadratic function
            # Then, minima are saved in a matlab and/or txt file

            # Identify min and max values
            min_analysis1 = min_analysis(clean_signals_norm, SetData.points, SetData.t_min, SetData.t_max)

            #Save file in MATLAB
            #print('Do you want to save the minima collection in a matlab file? Y/N')
            #name = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + '/t' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut) + ']Hz/' + 'min_points_matlab'+ SetData.IMG_TYPE + '_' + str(set) + '[' + str(lowcut) + '_' + str(highcut) + ']'
            #min_to_matlab(min_analysis1['min_time'], number, time, name)

            #SAVE FILE TXT
            #print('Do you want to save the minima collection in a txt file? Y/N')
            name = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'min_points_'+ SetData.IMG_TYPE + '_' + str(set) + '[' + str(lowcut) + '_' + str(highcut) + ']' + '.txt'
            saveMin(name, min_analysis1['min_time'])
            name = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'params_'+ SetData.IMG_TYPE + '_' + str(set) + '[' + str(lowcut) + '_' + str(highcut) + ']' + '.txt'
            #np.savetxt(name, min_analysis1['params'], fmt = '%s', delimiter = ' ')
            saveParams(name, min_analysis1['params'])

            answer = safeRawInput('Are you interested in a different spectrum range? (Y/N):', valid_ans = ['Y', 'N'])
            cont += 1

        lowcut_tot.append(lowcut_set)
        highcut_tot.append(highcut_set)
        zone_tot.append(zone_set)
        del img_coll_norm
        del spectrum_mean
        del img_spectrum_freq

    #----- else... passband is specified in the SetData.py
    else: #if the user already know what filter parameters should be used
        lowcut_set = []
        highcut_set = []
        zone_set = []

        for t in range(0, len(SetData.lowcut[index])):
            path = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/'
            print('path HERE: ', path)
            #Highcut and lowcut setting
            if not SetData.lowcut[index][t] is None:
                lowcut = SetData.lowcut[index][t]
                lowcut_set.append(lowcut)

            else:
                highcut = SetData.highcut[index][t]
                #plot spectrum
                Spectrum_plotting(img_spectrum_freq, spectrum_mean, set, path, highcut = highcut)

                #take lowcut as input
                ans = True;
                while ans:
                    ans = False
                    lowcut = safeInput('Lowcut:')

                    if (lowcut >= highcut):
                        ans = True
                        print("WARNING: Lowcut = " + str(lowcut) + " Highcut = " + str(highcut))
                        print("Lowcut is greater than highcut. Please re-insert lowcut value.")

                SetData.lowcut[index][t] = lowcut
                lowcut_set.append(lowcut)
                #print(lowcut)
                #lowcut_set[index][t] = lowcut
                print(SetData.lowcut[index][t])

            if not SetData.highcut[index][t] is None:
                highcut = SetData.highcut[index][t]
                highcut_set.append(highcut)

            else: #if highcut is NONE

                #plot spectrum
                Spectrum_plotting(img_spectrum_freq, spectrum_mean, set, path, lowcut = lowcut)

                #take lowcut as input
                ans = True;
                while ans:
                    ans = False
                    highcut = safeInput('Highcut:')

                    if (lowcut >= highcut):
                        ans = True
                        print("Lowcut is greater than highcut. Please re-insert highcut value.")


                highcut_set.append(highcut)

            if SetData.zone[index][t] is None:
                zone = getZone()

                zone_set.append(zone)

            else:
                zone_set.append(SetData.zone[index][t])

            #lowcut_set = np.asarray(lowcut_set)
            #highcut_set = np.asarray(highcut_set)
            #path were things are saved
            path = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut) + ']Hz/'
            if os.path.exists(path)== False:
                os.makedirs(path)

            Spectrum_plotting(img_spectrum_freq, spectrum_mean, set, path, lowcut, highcut, show = False, save = True)
            #del spectrum_mean
            #del spectrum_mean_measures[index]

            RawSignal_plotting(img_coll_norm, set, lowcut, highcut, path)

            #------------------------------- BANDPASS FILTER ---------------------------
            # Clean the signal through a bandpass-pass filter
            print('     Filtering the signal...')

            #clean_signals = bandpass_filter(img_collection_reduced, SetData.order, SetData.fs, lowcut, highcut)
            clean_signals_norm = bandpass_filter(img_coll_norm, SetData.order, SetData.fs, lowcut, highcut)
            #del img_coll_norm

            CleanSignal_plotting(clean_signals_norm, set, lowcut, highcut, path)

            #--------------------------- INTERPOLATION of the MINIMA -------------------
            # In this section, minima are found and interpolated with a quadratic function
            # Then, minimuma are saved in a matlab and/or txt file

            # Identify min and max values
            min_analysis1 = min_analysis(clean_signals_norm, SetData.points, SetData.t_min, SetData.t_max)

            #Save file in MATLAB
            #print('Do you want to save the minima collection in a matlab file? Y/N')
            #print('Lowcut:' +str(lowcut))
            #print('Highcut:'+str(highcut))

            #name = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + '/t' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut) + ']Hz/' + 'min_points_matlab'+ SetData.IMG_TYPE + '_' + str(set) + '[' + str(lowcut) + '_' + str(highcut) + ']'
            #min_to_matlab(min_analysis1['min_time'], number, time, name)

            #SAVE FILE TXT
            #print('Do you want to save the minima collection in a txt file? Y/N')

            #name = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'min_points_'+ SetData.IMG_TYPE + '_' + str(set) + '[' + str(lowcut) + '_' + str(highcut) + ']' + '.txt'
            name = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'min_points.txt'
            saveMin(name, min_analysis1['min_time'])

            #name = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'params_'+ SetData.IMG_TYPE + '_' + str(set) + '[' + str(lowcut) + '_' + str(highcut) + ']' + '.txt'
            name = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/[' + str(lowcut) + '_' + str(highcut)  + ']Hz/' + 'params.txt'
            #np.savetxt(name, min_analysis1['params'], fmt = '%s', delimiter = ' ')
            saveParams(name, min_analysis1['params'])

        lowcut_tot.append(lowcut_set)
        highcut_tot.append(highcut_set)
        zone_tot.append(zone_set)
        del img_coll_norm
        del spectrum_mean
        del img_spectrum_freq


    # Here we save lowcut and highcut im a txt file
    filename = SetData.ANALYSIS_DIR + SetData.IMG_TYPE + 't' + SET_NUMBER + '/' + 'low_high.txt'
    print('output filename for passband: ', filename)

    for lc in lowcut_tot:
        low = ""
        for elem in lc:
            low += str(elem) + ','
        low += ';'

    for lc in highcut_tot:
        high = ""
        for elem in lc:
            high += str(elem) + ','
        high += ';'

    for lc in zone_tot:
        z = ""
        for elem in lc:
            z += str(elem) + ","
        z += ';'

    with open(filename, "w") as file:
        file.write('#In this file lowcut and highcut parameters are saved. PLEASE DO NOT MODIFY;')
        file.write('\n')
        file.write(low)
        file.write('\n')
        file.write(high)
        file.write('\n')
        file.write(z)

    img_coll_norm_measures[index]=0
    spectrum_mean_measures[index]=0
    img_spectrum_freq_measures[index]=0
#================================================================================#
#================================================================================#
#print(psutil.cpu_count())
#psutil.virtual_memory()

print('Image analysis completed!')
