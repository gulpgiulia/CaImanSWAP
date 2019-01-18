% setParamsAndOptions.m
%
%   Load analysis params and information about recording dataset.
%
%   Copyright 2015 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Feb. 18, 2015
%   Version: 2.0 - Giu. 29, 2018

%% Optical experiments parameters

%width = 50;
%height = 50; %40 for 700110 50 for 700111
%sampling_time = 0.001; % [s]
%pixel_size = 0.0001; %[m]
path = DataDir; 

%% 

debug=0; % debug level
test=0; % for testing, to not overwrite results

%% Load analysis params and information about recording dataset.
Options.PeriodToAnalyze = [0 500];

Options.SaveMUA = 0;

Options.LogMUA.FreqBand = [200 1500];
Options.LogMUA.SmoothingWindow = 0.010; % It depends on the sampling rate.
Options.UD.ThresholdModulation = 0.6;   % Between 0 and 1. 0.5 is the mean of the Up and Down levels.

Options.LFP.SmoothingWindow = Options.LogMUA.SmoothingWindow; % s...
Options.LFP.SubsamplingRate = Options.LogMUA.FreqBand(1);   % Hz...

Options.LFPMUAplot.MUArange = [-0.75 3.25];
Options.LFPMUAplot.LFPrange = [-2000 2000];

% N.B. struct UpTimeLags depends on the electrode array %CAMBIO!
Options.UpTimeLags.MaxAbsTimeLag = 0.800;       % Maximum reasonable time lag between electrodes...
% Options.UpTimeLags.ReferenceChannel = 11;        % Channel whose triggers are used as reference. If not defined use channel close to the middle is automatically selected.
Options.UpTimeLags.SortingChannels = [05 06 07]; % A guess to sorting timelags array. Usually the channels belonging to the same cluster.
Options.UpTimeLags.PCNumOfSTD = 3;              % Num of st.dev. in the PC space to avoid outliers in the PCA.
% % Options.UpTimeLags.InvertElectrodePos = 1;      % If exist and is 1 invert the order of the electrodes position.
Options.UpTimeLags.NumOfClusters = 10;          % Number of clusters of wavefronts.
Options.UpTimeLags.WavefrontTimeGap = 20;       % Time gap between wavefront in ms.

% if ~exist('loop','var')
%     disp('loop NOT EXIST')
%     User='local'; % 'Maurizio', 'Giulia' (on ikkio), 'local' (Giulia on local nemo22) ...
%     setPath % setting of the workspace directories
% else
%     if ~loop
%         disp('loop IS FALSE')
%         User='local'; % 'Maurizio', 'Giulia' (on ikkio), 'local' (Giulia on local nemo22) ...
%         setPath % setting of the workspace directories
%     end
% end


% -------------------------------------------------------------------------
% ELECTRODE ARRAY - 32 channels
% -------------------------------------------------------------------------
% --> distinguish between RIGHT (RH) and LEFT (LH) Hemisphere
% --> define a reference system
% --> enable the opportunity to select a specific Cortical Area:
%   V = Visual Cortex
%   R = Retrosplenial Cortex 
%   P = Parietal Association Area (PtA)
%   S = Somatosensory Cortex
%   M = Motor Cortex

%disp(['*** FileName: ' FileName]);
fprintf('\n*** FileName: %s\n', FileName);
hemisphere = FileName(end-1:end); % 'RH/LH'

SelectedArea = 'A'; % A = all; or 'V', or 'R', or 'P' ... (see above)
evoked=0;

%RecordingSet = [1:nCh];
%if strcmpi(hemisphere,'RH')

%     
%     n = 1;
%     RecordingSet(n).label ='01.Ch27';
%     RecordingSet(n).port = '27';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;    %High-Pass Filter cut-off frequency
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 0*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1;
%     RecordingSet(n).label ='02.Ch18';
%     RecordingSet(n).port = '18';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 0*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1;
%     RecordingSet(n).label ='03.Ch28';
%     RecordingSet(n).port = '28';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 1*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1;
%     RecordingSet(n).label ='04.Ch19';
%     RecordingSet(n).port = '19';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;    
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 1*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1; 
%     RecordingSet(n).label ='05.Ch10';
%     RecordingSet(n).port = '10';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 1*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1; 
%     RecordingSet(n).label ='06.Ch29';
%     RecordingSet(n).port = '29';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 2*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1; 
%     RecordingSet(n).label ='07.Ch20';
%     RecordingSet(n).port = '20';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 2*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1; 
%     RecordingSet(n).label ='08.Ch11';
%     RecordingSet(n).port = '11';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 2*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1;
%     RecordingSet(n).label ='09.Ch30';
%     RecordingSet(n).port = '30';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 3*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%  
%     n = n + 1;
%     RecordingSet(n).label ='10.Ch21';
%     RecordingSet(n).port = '21';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 3*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1;
%     RecordingSet(n).label ='11.Ch12';
%     RecordingSet(n).port = '12';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 3*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='12.Ch4';
%     RecordingSet(n).port = '4';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 3*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='13.Ch31';
%     RecordingSet(n).port = '31';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 4*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='14.Ch22';
%     RecordingSet(n).port = '22';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration = 0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 4*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='15.Ch13';
%     RecordingSet(n).port = '13';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 4*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='16.Ch5';
%     RecordingSet(n).port = '5';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 4*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='17.Ch32';
%     RecordingSet(n).port = '32';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 5*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='18.Ch23';
%     RecordingSet(n).port = '23';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration = 0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 5*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='19.Ch14';
%     RecordingSet(n).port = '14';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 5*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='20.Ch6';
%     RecordingSet(n).port = '6';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 5*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='21.Ch33';
%     RecordingSet(n).port = '33';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 6*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='22.Ch24';
%     RecordingSet(n).port = '24';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 6*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='23.Ch15';
%     RecordingSet(n).port = '15';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration = 0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 6*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='24.Ch7';
%     RecordingSet(n).port = '7';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 6*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='25.Ch3';
%     RecordingSet(n).port = '3';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 4*0.550;
%     RecordingSet(n).YPos = 6*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='26.Ch34';
%     RecordingSet(n).port = '34';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 7*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='27.Ch25';
%     RecordingSet(n).port = '25';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 7*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='28.Ch16';
%     RecordingSet(n).port = '16';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 7*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='29.Ch8';
%     RecordingSet(n).port = '8';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration = 0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 7*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='30.Ch26';
%     RecordingSet(n).port = '26';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 8*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='31.Ch17';
%     RecordingSet(n).port = '17';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 8*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='32.Ch9';
%     RecordingSet(n).port = '9';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 8*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%             
% elseif strcmpi(hemisphere,'LH')
%     
%     n = 1;
%     RecordingSet(n).label ='01.Ch10';
%     RecordingSet(n).port = '10';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 0*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1;
%     RecordingSet(n).label ='02.Ch19';
%     RecordingSet(n).port = '19';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 0*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1;
%     RecordingSet(n).label ='03.Ch9';
%     RecordingSet(n).port = '9';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 1*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1;
%     RecordingSet(n).label ='04.Ch18';
%     RecordingSet(n).port = '18';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;    
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 1*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1; 
%     RecordingSet(n).label ='05.Ch27';
%     RecordingSet(n).port = '27';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 1*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1; 
%     RecordingSet(n).label ='06.Ch8';
%     RecordingSet(n).port = '8';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 2*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1; 
%     RecordingSet(n).label ='07.Ch17';
%     RecordingSet(n).port = '17';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 2*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1; 
%     RecordingSet(n).label ='08.Ch26';
%     RecordingSet(n).port = '26';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 2*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1;
%     RecordingSet(n).label ='09.Ch7';
%     RecordingSet(n).port = '7';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 3*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%  
%     n = n + 1;
%     RecordingSet(n).label ='10.Ch16';
%     RecordingSet(n).port = '16';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 3*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
% 
%     n = n + 1;
%     RecordingSet(n).label ='11.Ch25';
%     RecordingSet(n).port = '25';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 3*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='12.Ch33';
%     RecordingSet(n).port = '33';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 3*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='13.Ch6';
%     RecordingSet(n).port = '6';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 4*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='14.Ch15';
%     RecordingSet(n).port = '15';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration = 0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 4*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='15.Ch24';
%     RecordingSet(n).port = '24';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 4*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='16.Ch32';
%     RecordingSet(n).port = '32';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 4*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='17.Ch5';
%     RecordingSet(n).port = '5';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 5*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='18.Ch14';
%     RecordingSet(n).port = '14';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration = 0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 5*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='19.Ch23';
%     RecordingSet(n).port = '23';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 5*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='20.Ch31';
%     RecordingSet(n).port = '31';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 5*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='21.Ch4';
%     RecordingSet(n).port = '4';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 6*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='22.Ch13';
%     RecordingSet(n).port = '13';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 6*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='23.Ch22';
%     RecordingSet(n).port = '22';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration = 0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 6*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='24.Ch30';
%     RecordingSet(n).port = '30';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 6*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='25.Ch34';
%     RecordingSet(n).port = '34';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 4*0.550;
%     RecordingSet(n).YPos = 6*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='26.Ch3';
%     RecordingSet(n).port = '3';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 0*0.550;
%     RecordingSet(n).YPos = 7*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='27.Ch12';
%     RecordingSet(n).port = '12';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 7*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='28.Ch21';
%     RecordingSet(n).port = '21';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 7*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='29.Ch29';
%     RecordingSet(n).port = '29';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration = 0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 7*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='30.Ch11';
%     RecordingSet(n).port = '11';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 1*0.550;
%     RecordingSet(n).YPos = 8*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='31.Ch20';
%     RecordingSet(n).port = '20';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 2*0.550;
%     RecordingSet(n).YPos = 8*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%     
%     n = n + 1;
%     RecordingSet(n).label ='32.Ch28';
%     RecordingSet(n).port = '28';
%     RecordingSet(n).AbsoluteThreshold = 0.40;
%     RecordingSet(n).MinStateDuration =  0.080;
%     RecordingSet(n).HPFCutOffFreq = 0.1;
%     RecordingSet(n).XPos = 3*0.550;
%     RecordingSet(n).YPos = 8*0.550;
%     % RecordingSet(n).PCAforUpDownDetection = 1;
%        
% else 
%     disp('Please, specify if Left Hemisphere or Right Hemisphere')
% end
%     
% if evoked
%     EventArray.port = 'Stim';
% end
% 
% % --- SelectedArea ---
% if SelectedArea == 'V'     % V = Visual Cortex
%     ChannelSet = [18 19 20 22 23 24 25 27 28 29 30 31 32];
% elseif SelectedArea == 'R' % R = Retrosplenial Cortex 
%     ChannelSet = [13 17 21 26];
% elseif SelectedArea == 'P' % P = Parietal Association Area (PtA)
%     ChannelSet = [14 15 16];
% elseif SelectedArea == 'S' % S = Somatosensory Cortex
%     ChannelSet = [4 5 7 8 10 11 12];
% elseif SelectedArea == 'M' % M = Motor Cortex
%     ChannelSet = [1 2 3 6 9];
% elseif SelectedArea == 'A' % A = All the electrodes
%     ChannelSet = [1:32];
% else
%     disp('Please, specify the Selected Area of the Cortex')
% end
% 

%% LOCALITY criteria
%

% %%Geometry only
arrayMask = [];

h = height;
w = width;

for (i=1:1:(h*w)) 
    
    %Non ho condizioni periodiche al bordo, sto individuando i primi
    %vicini!    
    vicini = [];
    
    %controllo a dx
    if (mod((i+1),w) ~= 1) 
        vicini = [vicini i+1];    
    end    
    
    if (mod((i+w+1),w) ~= 1 & ((i+w+1)<(w*h+1))) 
        vicini = [vicini i+w+1];
    end
    
    if (mod((i-w+1),w) ~= 1 & ((i-w+1)>1))
        vicini = [vicini i-w+1];
    end 
    
    %controllo a sx
    if (mod((i-1),w) ~= 0) 
        vicini = [vicini i-1];
    end 
    
    if ((mod((i-w-1),w) ~= 0) & ((i-w-1)>0)) 
        vicini = [vicini i-w-1];
    end
    
    if ((mod((i+w-1),w) ~= 0) & ((i+w-1)<(w*h)))
        vicini = [vicini i+w-1];
    end
    
    %spigoli
    if ((i+w) < (w*h+1)) 
        vicini = [vicini i+w];
    end
    
    if ((i-w) > 0)
        vicini = [vicini i-w];
    end     
        
    arrayMask = [arrayMask {[i, {vicini}]}];
    
end
