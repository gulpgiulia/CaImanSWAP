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
#SAMPLING_TIME: time rate of the acquisition process (in seconds)
SAMPLING_TIME = 0.04
#TIME_MIN: first time frame number (string)
TIME_MIN = 1
#TIME_MAX: last time frame number (string)
TIME_MAX = 1000
TIME_SAMPLES = TIME_MAX-TIME_MIN
PIXEL_SIZE = 0.05 #mm
#number of processor to use during the analysis
nprocs = 2
#----------------------------- Dataset information -----------------------------
#images should be stored as
#<IMG_PATH>/<IMG_TYPE>/t<SET_NUMBER>/<IMG_NAME><SET_NUMBER>_<TIME_FRAME>.<EXTENSION>

#IMG_PATH (string): path to the dataset
DATA_DIR = '../../DATA/'
#IMG_TYPE (string)
IMG_TYPE = '170110'
#IMG_NAME (string)
IMG_NAME = 'provevideo'
#EXTENSION
EXTENSION = '.tif'
#num_measures
num_measures = [1] #numbers associated to<SET_NUMBER>


#----------------------------- Analysis information ----------------------------
Contour_Limit = [0.197, 0.258, 0.2715, 0.25, 0.195, 0.185]
#ANALYSIS PATH (string): path of the directory where results have to be saved
ANALYSIS_DIR = '../../ANALYSIS/'
#lowcut: lowcut of the butterworth filter. If not sure of this number type None
lowcut = [[0.5],[0.5],[0.5],[0.5],[0.5],[0.5]]
zone = [['Zone_1'], ['Zone_1'], ['Zone_1'], ['Zone_1'], ['Zone_1'], ['Zone_1']]
#highcut: highcut of the butterworth filter. If not sure of this number type None
highcut = [[2.5],[3.5],[3.5],[3.5],[3.5],[3.5]]
#highcut = [[3.5, 3.3], [2.5], [3.5], [3.5, 5.0], [3.5, 5.0], [3.5, 5.0]]
spectrum_time = 1000
#Butterworth order
order = 6
fs = 25.00 #Hz
#fit parameter: how many points shoul be used?
ZOOM = 4
points = 5
# Here we defined the number of elements (of pixels idx's) past-to the closest value
# of wave born time to be consider as part of the wave birth itself
k_cluster = 5
t_min = 5
t_max = 995
#PIXEL_SAMPLE_x: x cohordinate of the Pixel studied (already reduced!)
PIXEL_SAMPLE_x = 30
#PIXEL_SAMPLE_y: y cohordinate of the Pixel studied (already reduced!)
PIXEL_SAMPLE_y = 30
TIME_SAMPLE_min =0
TIME_SAMPLE_max =TIME_SAMPLES
#graph colo
color = 'blue'

#----------------------------- Consinstecy check -------------------------------
#----------------------------PLEASE DO NOT MODIFY-------------------------------
#check of Contour_Limit and num_measures
if len(Contour_Limit) != len(num_measures):
    print('ERROR: Contour_Limit and num_measures dimensions not match')
    sys.exit(-1)

#check of lowcut and highcut with num_measures
if len(lowcut) != len(highcut):
    print('ERROR: Wrong lowcut and highcut dimensions')
    sys.exit(-1)

if len (zone) != len (lowcut):
    print ('ERROR: Wron zone and lowcut dimensions')
    sys.exit(-1)

if len(lowcut) != len(num_measures):
    print('ERROR: Wrong lowcut/highcut and num_measures dimensions')
    sys.exit(-1)

#check of lowcut and highcut components

for i in range(0,len(highcut)):
    try:
        if len(lowcut[i]) != len(highcut[i]):
            print('ERROR: Non matching dimension for the %d element of lowcut and highcut' % i)
            sys.exit(-1)

        if len (zone[i]) != len (lowcut[i]):
            print ('ERROR: Non matching dimension for %d element of zone and lowcut' % i)
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
