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
SAMPLING_TIME = 0.01
#TIME_MIN: first time frame number (string)
TIME_MIN = 3
#TIME_MAX: last time frame number (string)
TIME_MAX = 10000
TIME_SAMPLES = TIME_MAX-TIME_MIN
PIXEL_SIZE = 0.01
nprocs = 4
#----------------------------- Dataset information -----------------------------
#images should be stored as
#<IMG_PATH>/<IMG_TYPE>/t<SET_NUMBER>/<IMG_NAME><SET_NUMBER>_<TIME_FRAME>.<EXTENSION>

#IMG_PATH (string): path to the dataset
DATA_DIR = '../../DATA/'
#IMG_TYPE (string)
IMG_TYPE = 'ISO_180322'
#IMG_NAME (string)
IMG_NAME = 'image'
#EXTENSION
EXTENSION = '.tif'
#num_measures = [24]
#num_measures = [2,3,4] #numbers associated to<SET_NUMBER>
num_measures = [2, 3, 4, 5, 11, 12, 13, 16, 17, 18, 22, 23, 24]


#----------------------------- Analysis information ----------------------------
#Contour_Limit = [None]
Contour_Limit = [0.14, 0.15, 0.14, 0.13,
                0.13, 0.11, 0.11, 0.105,
                0.105, 0.11, 0.10, 0.095,
                0.101]

#ANALYSIS PATH (string): path of the directory where results have to be saved
ANALYSIS_DIR = '../../ANALYSIS/'
#lowcut: lowcut of the butterworth filter. If not sure of this number type None
#lowcut = [None]
#lowcut = [[0.4, 0.8, 1.7],[0.25, 0.8, 1.7],[0.25, 0.8, 1.7]]
lowcut = [[0.4, 1.1, 2.2, 3.2, 3.9, 4.3],[0.4, 1.1, 2.1, 3.0, 3.9, 4.3],[0.4, 1.1, 2.1, 3.2],[0.4, 1.1, 2.3, 3.4, 3.9, 4.3],[0.4, 1.0, 2.0, 3.2, 3.9, 4.2],[0.4, 1.0, 2.0, 3.0, 3.9],[0.4, 0.9, 1.9, 2.8, 3.8],[0.25, 0.9, 1.7, 2.6, 3.9],[0.25, 0.8, 1.8, 2.6, 3.9],[0.3, 0.9, 1.7, 2.7, 3.8],[0.25, 0.9, 2.03, 3.0, 3.9],[0.25, 1.0, 2.1, 3.2, 3.9, 4.3],[0.25, 1.0, 2.2, 3.3, 3.9, 4.5]]
#zone = [None]
#zone = [['Zone_1', 'Zone_2', 'Zone_3'],['Zone_1', 'Zone_2', 'Zone_3'], ['Zone_1', 'Zone_2', 'Zone_3']]
zone = [['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5', 'Zone_6'],['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5', 'Zone_6'],['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4'],
        ['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5', 'Zone_6'],['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5', 'Zone_6'],['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5'],
        ['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5'], ['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5'],['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5'],
        ['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5'],['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5'],['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5', 'Zone_6'],
        ['Zone_1', 'Zone_2', 'Zone_3', 'Zone_4', 'Zone_5', 'Zone_6']]
#highcut: highcut of the butterworth filter. If not sure of this number type None
#highcut = [None]
#highcut = [[0.8, 1.7, 2.5],[0.8, 1.7, 2.5],[0.8, 1.7, 2.5]]
highcut = [[0.6, 1.3, 2.4, 3.5, 4.2, 4.5],[0.6, 1.3, 2.4, 3.5, 4.2, 4.5],[0.6, 1.3, 2.5, 3.5],[0.6, 1.3, 2.5, 3.8, 4.2, 4.5],[0.6, 1.2, 2.2, 3.5, 4.2, 4.4],[0.6, 1.3, 2.4, 3.4, 4.2],[0.6, 1.1, 2.2, 3.2, 4.2], [0.6, 1.1, 1.9, 2.9, 4.2],[0.5, 1.0, 2.0, 3.0, 4.2],[0.6, 1.0, 2.0, 3.0, 4.2],[0.3, 1.2, 2.25, 3.3, 4.1],[0.6, 1.2, 2.3, 3.4, 4.1, 4.5],[0.5, 1.2, 2.4, 3.5, 4.1, 4.6]]

spectrum_time = 10000
#Butterworth order
order = 4
fs = 100.00 #Hz
#fit parameter: how many points shoul be used?
ZOOM = 4
points = 5
# Here we defined the number of elements (of pixels idx's) past-to the closest value
# of wave born time to be consider as part of the wave birth itself
k_cluster = 5
t_min = 5
t_max = 9995
#PIXEL_SAMPLE_x: x cohordinate of the Pixel studied (already reduced!)
PIXEL_SAMPLE_x = 20
#PIXEL_SAMPLE_y: y cohordinate of the Pixel studied (already reduced!)
PIXEL_SAMPLE_y = 20
TIME_SAMPLE_min =0
TIME_SAMPLE_max =TIME_SAMPLES
#graph colo
color = 'orange'

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
