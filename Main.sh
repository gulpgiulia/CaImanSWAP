cd PYTHON
python2 Auto_contours.py
python2 Image_analysis.py
cp SetData.py ../MATLAB/SetData.txt
cd ..
cd MATLAB
 /Applications/MATLAB_R2015b.app/bin/matlab -nodesktop -nosplash -r 'DatasetLounch.mat'
cd ..
cd PYTHON
python2 DataAnalisys
