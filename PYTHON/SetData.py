#Here you need to set few parameters that describe your dataset.
#Please, write the value under each definition.

import sys
import os

#------------------------------ Image information ------------------------------
#DIM_X: X dimension of your image (in pixels)
DIM_X = 100
#DIM_Y: Y dimension of your image (in pixels)
DIM_Y = 100
#MACRO_PIXEL_DIM: dimension of new 'macro-pixel' after the spacial smoothing
MACRO_PIXEL_DIM = 2
### --> Evaluate to suppress the macro-pixel smoothing
#SAMPLING_TIME: time rate of the acquisition process (in seconds)
SAMPLING_TIME = 0.04
#TIME_MIN: first time frame number (string)
TIME_MIN = 1
#TIME_MAX: last time frame number (string)
TIME_MAX = 250
TIME_SAMPLES = TIME_MAX-TIME_MIN
PIXEL_SIZE = 0.05 #mm
#number of processor to use during the analysis
nprocs = 2

#----------------------------- Dataset information -----------------------------
#images should be stored as
#<IMG_PATH>/<IMG_TYPE>/t<SET_NUMBER>/<IMG_NAME><SET_NUMBER>_<TIME_FRAME>.<EXTENSION>

#IMG_PATH (string): path to the dataset
DATA_DIR = '/home/rgutzen/Sciebo/own/Data/WaveScalES/LENS/'
#IMG_TYPE (string)
IMG_TYPE = '170110_mouse2_deep/'
#IMG_NAME (string)
IMG_NAME = 'provevideo'
#EXTENSION
EXTENSION = '.tif'
#num_measures --> number of datasets to analyse
num_measures = [1] #numbers associated to<SET_NUMBER>

###DATA_DIR = '/Users/giuliadebonis/GoogleDriveSapienza/LENS-INFN/DATA/Ketamine/ketamine_25Hz/170110/mouse2/deep/t1/'

#----------------------------- Analysis information ----------------------------
#Contour_Limit is the input parameter of the 'find_contours' function in the measure.find_contours module of the scikit-image Python package, https://scikit-image.org.
Contour_Limit = [0.197]
#ANALYSIS PATH (string): path of the directory where results have to be saved
ANALYSIS_DIR = '/home/rgutzen/ProjectsData/wavescales/CaImaging/'

#*** FILTER ***
#FILTER_TYPE = 'Butterworth'
#Butterworth order
order = 6
### --> Evaluate other types of filters (Bessel, FIR, ...)

# Passband
#lowcut: lowcut of the Butterworth filter. If not sure of this number type None
lowcut = [[0.5]]
#highcut: highcut of the Butterworth filter. If not sure of this number type None
highcut = [[3.0]]
#zone: tag given to the portion of the frequency spectrum selected by the filter
zone = [['Zone_1']]
### len(lowcut), len(highcut), len(zone) == len(num_measures) (see Consistency Checks)
### The filter settings can be differentiate for each dataset
### If len(lowcut[i]) is >1 for dataset i, more than 1 band-pass filter is applied (whose settings have to be specified), and thus more than 1 Zone can be identified in the spectrum.

#*** FIT ***
#POLYFIT_ORDER = 2
### polynomial fit, order 2 (parabola)

#points: how many points should be used?
ZOOM = 4 #ZOOM ???
points = 5

#*** Visualization Options ***

#The represented spectrum is computed over a time interval of length 'spectrum_time' (in samples)
spectrum_time = 1000
#fs = sampling frequency
fs = 25.00 #[Hz]

# Here we define the number of pixels closest in time to the wave born time, that are consdered as part of the wave birth itself
k_cluster = 5
t_min = 5
t_max = 995

### SamplePixel (far from bounderies and blood vessels) for illustrating results:
#PIXEL_SAMPLE_x: x-coordinate of the SamplePixel (already reduced!)
PIXEL_SAMPLE_x = 30
#PIXEL_SAMPLE_y: y-cohordinate of the SamplePixel (already reduced!)
PIXEL_SAMPLE_y = 30
TIME_SAMPLE_min =0
TIME_SAMPLE_max =TIME_SAMPLES
#graph color
color = 'blue'

#----------------------------- Consinstency Checks -------------------------------
#-----------------------------PLEASE DO NOT MODIFY--------------------------------
#check of Contour_Limit and num_measures
if len(Contour_Limit) != len(num_measures):
    print('ERROR: Contour_Limit and num_measures dimensions do not match')
    sys.exit(-1)

#check of lowcut and highcut with num_measures
if len(lowcut) != len(highcut):
    print('ERROR: Wrong lowcut and highcut dimensions')
    sys.exit(-1)

if len(zone) != len(lowcut):
    print ('ERROR: Wron zone and passband dimensions')
    sys.exit(-1)

if len(lowcut) != len(num_measures):
    print('ERROR: Wrong passband and num_measures dimensions')
    sys.exit(-1)

#check of lowcut and highcut components

for i in range(0,len(highcut)):
    try:
        if len(lowcut[i]) != len(highcut[i]):
            print('ERROR: Non-matching dimension for the %d element of lowcut and highcut' % i)
            sys.exit(-1)

        if len (zone[i]) != len (lowcut[i]):
            print ('ERROR: Non-matching dimension for %d element of zone and lowcut' % i)
            print ('Len zone is %d while Len lowcut is %d' % (len(zone[i]), len(lowcut[i])))
            sys.exit(-1)

        for l_el, h_el in zip (lowcut[i], highcut[i]):
            if l_el is None or h_el is None:
                continue

            if l_el > h_el:
                print ('ERROR: Lowcut[' + str(i) + '] element ' + str(l_el) + \
                       ' is greater than highcut[' + str(i) +'] element ' + str(h_el))
                sys.exit(-1)

    except TypeError:
        if lowcut[i] is None and highcut[i] is None and zone[i] is None:
            continue

        else:
            print('ERROR: Wrong syntax element %d' % i)
            sys.exit(-1)

#check directory
if os.path.exists(DATA_DIR) == False:
    print('ERROR:The selected DATA_DIR does not exist')
    sys.exit(-1)
