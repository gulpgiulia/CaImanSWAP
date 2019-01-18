
This program has been developed by Marco Celotto, Chiara De Luca and Paolo Muratore within the HBP_WAVESCALES2 project under the supervision of Pier Stanislao Paolucci and Giulia De Bonis.

The code implements an analysis pipeline that has been developed to analyze '.tif' cortical images obtained from anestethized mice through GCaMP6f fluorescence. The goal of this code is to identify and reconstruct cortical slow delta waves (0.5-4Hz band) from images and to take measures of this phenomena. Further details about this work can be found in the following article: arXiv:1811.11687 


--------------------------STRUCTURE OF THE PROGRAM-------------------------------

Except for where differently specified, the code is written in Python2.7

#ATTENTION: To correctly execute the program, the user doesn't have to modify the directories structure that the code automatically generates to manages the data.

The program is devided in 4 phases:
  - Image_initialization.py: In this phase the raw images that has to be analyzed is acquired as an input. To only select the effective area of images from where the signal of interest comes, an optimal mask is constructed and applied on the raw data. The background of the images set is then removed and a spatial smoothing is made to lower the noise. At the end of this phase a series of '.txt' files containing numpy arrays used to represent the processed images is generated. Each data set can be initalized just once and then used as a starting point for different analysis trials (e.g. for various regions of the spectrum that the user wants to analyze, as we will furtherly explain).  
 
  - Image_analysis.py: In this phase the mean spectrum of a set of images is computed and it is presented to the user. The user can select one or more regions of interests that have to be analyzed separately. In such a way the waves analysis can be performed on different subregions of the delta frequency band. Once these regions have been declared, they are selected via a band-pass Butterworth filter. Minima of signals (passages of the waves down-to-up states wavefront) are identified for each spectrum region and a parabolic interpolation is computed for each minimum. In such a way minima times are obtained with greater precision than the one imposed by the images sampling frequency. Parabulas parameters are stored for each minimum in a proper data structure as well as minima times for each pixel (minima times collection). 

#ATTENTION: This next pahse has been implemented in MATLAB starting from a code written by Giulia De Bonis and Maurizio Mattia. 

  - DatasetLounch.mat: Starting from minima times collection waves begin and end times are identified. A 'Begin_Time.txt' and 'End_time.txt' files couple containing these information is produced for each requested dataset spectrum region. 

  - Data_analysis.py: In this last phase the following measures are taken: neuronal population excitability (both in the form of an excitability distribution and of a cortical map of excitabilities), speed of waves and origin of waves. All the results are saved in the RESULTS directory.



--------------------------PROGRAM USAGE--------------------------------

The preliminary step to use the program is to compile the SetData.py script: it contains all the parameters relevant to the analysis.

  -SetData.py rules: variables names in this script are the actual names of the variables used in the entire Python2.7 code. Spectrum regions of interests for each dataset can here be declared a priori (e.g. known as a result of previous executions of the code) or left blank whether the user wants to decide frequency ranges, or even the number of regions itself, during the analysis. These decision can be taken by the user and written in SetData.py through the syntax explained in the following example:

      # We have 4 dataset to analyze
      num_measures = [1,2,3,4]

      # First of all we set the parameter that controls the mask cut level. There have to be as many level parameters as the numer of sets we are analyzing. For the sets where this parameter is unknown, insert 'None' in the correspondat space.       

      Contour_Limit = [0.18, None, None, 0.10]

      # Let's make some assumptions about our previous knowledge on the spectrum region of interest of the dataset. We will then report a slice of code where these informations are reported in the SetData.py syntax
      # In set 1 we know that just one spectrum region has to be analyzed: we know, a priori, that 1.9Hz is the superior limit of this region but we don't know the inferior limit and we want to dinamically select it during the execution of the program.  
      # In set 2 we still know one spectrum region has to be analyzed but here is the inferior limit to be known (0.2Hz) and the superior is unknown. 
      # In set 3 no information is known about the number of spectrum regions we want to analyze and, consequently, about their limits. 
      # In set 4 we know that we want to separately analize two regions of the spectrum, yet the inferior limit of the second one is unknown.

      lowcut = [[None], [0.2], None, [0.2, None]]
      highcut = [[1.9], [None], None, [1.0, 1.3]]
      zone = [['Zone_1'], ['Zone_1'], None, ['Zone_2', 'Zone_3']]

An example of execution in LINUX/IOS:
cd PYTHON
python2 Auto_contours.py
python2 Image_analysis.py
cp SetData.py ../MATLAB/SetData.txt
cd ..
cd MATLAB
 /bin/Matlab_path -nodesktop -nosplash -r 'DatasetLounch.mat'
cd ..
cd PYTHON
python2 DataAnalisys
