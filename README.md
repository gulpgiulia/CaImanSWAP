# CaImanSWAP
CaImanSWAP: Analysis Pipeline for Optical Data (wide-field fluorescence microscopy + calcium imaging)

------------------------------ INTRODUCTION ------------------------------

CaImanSWAP (Calcium-Imaging Slow Waves Analysis Pipeline) is the analysis pipeline for the identification of the slow-wave activity (delta band, [0.5-4.0] Hz) in Ca-Imaging datasets, as the images obtained with wide-field fluorescence microscopy on transgenic mouse model expressing a calcium indicator.

The program was originally developed by Marco Celotto, Chiara De Luca and Paolo Muratore within the WaveScalEs project (HBP-SP3) under the supervision of Pier Stanislao Paolucci and Giulia De Bonis (Istituto Nazionale di Fisica Nucleeare, INFN Rome, Italy), and applied to datasets produced by Francesco Resta and Anna Letizia Allegra Mascaro (HBP-SP1) using the LENS facilities (European Laboratory for Non-linear Spectroscopy, University of Florence, Italy) on GCaMP6f mouse model.
Further details in arXiv:1811.11687 [q-bio.NC].

The Human Brain Project HBP has received funding from the European Union's Horizon 2020 Framework Programme for Research and Innovation under the Specific Grant Agreements No. 785907 (Human Brain Project SGA2). HBP-SP1 is the sub-projects devoted to Mouse Brain Organization; HBP-SP3 is the sub-project devoted to Systems and Cognitive Neuroscience.  

---------------------- STRUCTURE OF THE PROGRAM --------------------------

Except for where differently specified, the code is written in Python2.7

#ATTENTION: To correctly execute the program, the user should not modify the directory structure that the code automatically generates to process the data.

The program is made up of 4 phases:
  1) Image_initialization.py: In this phase, the raw images to be analyzed are acquired as an input and translated in a matrix of floats. In order to select only the effective area where the signal of interest comes from, an optimal mask is constructed and applied to the raw data. Pixel-by-pixel background subtraction follows and a spatial smoothing can be made to further reduce the noise. At the end of this phase, '.txt' files with numpy arrays representing the processed images are generated. Usually, the image_initalization phase is unique for each set of raw data, resulting in a set of processed data that can be used as a starting point for different analysis trials (e.g. inspecting different regions of the spectrum, or adopting different filter families, see next).  
 
  2) Image_analysis.py: In this phase, the frequency spectrum of an image set is computed. The user can identify one or more regions of interests, that are analyzed separately to investigate the features of the delta frequency band. Once these regions are declared, they are selected via a band-pass Butterworth filter. The minima of the filtered signal (for each selected frequency band) are identified, and a parabolic interpolation is computed for each minimum. The minima are interpreted as the onset of the down-to-up transition, induced by the passage of a slow wave across the cortex. The parabolic interpolation is adopted to estimate the time of the minimum (i.e. the time of transition) with greater precision than the one imposed by the sampling frequency. The parabolic fit parameters are stored for each minimum in a proper data structure (collection of minima for each pixel).

  3) DatasetLaunch.mat (implemented in MATLAB starting from a code written by Giulia De Bonis and Maurizio Mattia): With the collection of minima as an input, the starting and ending times of each wave are identified; the pair of files 'Begin_Time.txt' and 'End_time.txt' is produced for each dataset (for each spectral region of interest), resulting in a collection of waves.

  4) Data_analysis.py: In this phase, the following measures are taken on the identified collection of waves: excitability of the neuronal population (represented both as a distribution and as a cortical map), speed of waves and origin of waves; results are saved in the RESULTS directory.

-------------------------- PROGRAM USAGE --------------------------------

**SETTINGS (SetData.py)**

The preliminary step to use the program is to compile the script 'SetData.py' that contains all the parameters relevant to the analysis.
NOTE: the names of variables in this script are the actual names used in the code, please do not modify them in the script. 

  -SetData.py:  Spectral regions of interest for each dataset can be declared a priori (e.g. if known as a result of previous executions of the code, or if some assumptions can be made) or left blank('None') if the user wants to decide frequency ranges, or even the number of regions itself, during the analysis. The following example illustrates how to specify the different options in the SetData.py.

      # If we have 4 datasets to analyze:
      num_measures = [1,2,3,4] 

      # First of all, we set the parameter that controls the mask cut level. We need to specify as many level parameters as the number of datasets we are analyzing. If the level value is unknown, insert 'None' in the corresponding list entry.       
      Contour_Limit = [0.18, None, None, 0.10]
      # ('Contour_Limit' is the input parameter of the 'find_contours' function in the measure.find_contours module of the scikit-image Python package, https://scikit-image.org).	

      # Concerning the selection of the spectral regions of interest:
      # - in dataset 1, we know that just one spectrum region has to be analyzed: we assume, a priori, that 1.9 Hz is the upper limit of this band but we don't know the lower limit and we want to dinamically determine it during the execution of the program; 
      # - in dataset 2, we still know that only one spectral region has to be selected, but here the lower limit is assumed (0.2Hz) and the upper is unknown; 
      # - in dataset 3, no information is available about the number of spectral regions we want to analyze and, consequently, about their limits; 
      # - in dataset 4, we want to separately analyze two regions of the spectrum, yet the lower limit of the second one is unknown.

      lowcut = [[None], [0.2], None, [0.2, None]]
      highcut = [[1.9], [None], None, [1.0, 1.3]]
      zone = [['Zone_1'], ['Zone_1'], None, ['Zone_2', 'Zone_3']]

**EXAMPLE EXECUTION (Linux/IOS)**

cd PYTHON
<edit SetData.py>
python2 Image_initialization.py
python2 Image_analysis.py
cp SetData.py ../MATLAB/SetData.txt
cd ..
cd MATLAB
 /your/path/to/bin/matlab -nodesktop -nosplash -r "run ('DatasetLaunch.m')" 
cd ..
cd PYTHON
python2 Data_analysis.py
