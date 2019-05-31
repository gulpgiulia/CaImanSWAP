
%% --- load data settings
clear
%BaseDir = '../..';
%ScriptDir = [BaseDir '/SRC/MATLAB/'];
ScriptDir = pwd;
%cd(ScriptDir)
%import data
file = readtable('SetData.txt', 'delimiter', '\n');

%% --- search for interesting lines
for i=1:height(file)
    string = file{i,:};
    if startsWith(string, 'num_measures')
        line_measures = string{1};
    end
    if startsWith(string, 'IMG_TYPE')
        line_type = string{1};
    end
    if startsWith(string, 'DIM_X')
        dimX = string{1};
    end
    if startsWith(string, 'DIM_Y')
        dimY = string{1};
    end
    if startsWith(string, 'MACRO_PIXEL_DIM')
        macro_dim = string{1};
    end
    if startsWith(string, 'SAMPLING_TIME')
        sampling_time = string{1};
    end
    if startsWith(string, 'PIXEL_SIZE')
        pixel_size = string{1};
    end
    if startsWith(string, 'ANALYSIS_DIR')
        analysis_dir = string{1};
    end
end

%% --- load parameters
%% PIXEL_SIZE
for idx = 1:length(pixel_size)
    if pixel_size(idx) == '='
        pixel_size(idx);
        inizio = idx+1;
    end
    if pixel_size(idx) == '#'
        fine = idx-1;
    end
end
pixel_size = pixel_size(inizio:fine);
pixel_size = strrep(pixel_size, ' ', '');
pixel_size = str2num(pixel_size);

%% MACROPIXEL_DIM
for idx = 1:length(macro_dim)
    if macro_dim(idx) == '='
        inizio = idx+1;
    end
end
fine = idx;
macro_dim = macro_dim(inizio:fine);
macro_dim = strrep(macro_dim, ' ', '');
macro_dim = str2num(macro_dim);

%% SAMPLING_TIME
for idx = 1:length(sampling_time)
    if sampling_time(idx) == '='
        inizio = idx+1;
    end
end
fine = idx;
sampling_time = sampling_time(inizio:fine);
sampling_time = strrep(sampling_time, ' ', '');
sampling_time = str2num(sampling_time);

%% DIM_X and DIM_Y
for idx = 1:length(dimX)
    if dimX(idx) == '='
        inizio = idx+1;
    end
end
fine = idx;
dimX = dimX(inizio:fine);
dimX = strrep(dimX, ' ', '');
dimX = str2num(dimX);
dimX = dimX / macro_dim;

for idx = 1:length(dimY)
    if dimY(idx) == '='
        inizio = idx+1;
    end
end
fine = idx;
dimY = dimY(inizio:fine);
dimY = strrep(dimY, ' ', '');
dimY = str2num(dimY);
dimY = dimY / macro_dim;

%% DIM_X and DIM_Y
flag = 0;
for idx=1:length(line_type)
    if line_type(idx) == "'"
        if flag == 0
            flag = 1;
            inizio = idx+1;
        else
            fine = idx-1;
        end
    end
end
IMG_TYPE = line_type(inizio:fine);

%% NUM_MEASURES
fine = 0;
for idx=1:length(line_measures)
    if line_measures(idx) == ']'
        fine = idx -1;
    end
    if line_measures(idx) == '['
        inizio = idx+1;
    end
end
if fine == 0
    fine = idx;
end
num_measures = line_measures(inizio:fine);
num_measures = strrep(num_measures, ',', ' ');
num_measures= str2num(num_measures);

%% ANALYSIS_DIR
for idx = 1:length(analysis_dir)
    if analysis_dir(idx) == '='
        inizio = idx+3;
    end
end
fine = idx-1;
analysis_dir = analysis_dir(inizio:fine);
analysis_dir = strrep(analysis_dir, ' ', '');

%%
AnalysisDir = [analysis_dir IMG_TYPE];

FileName = 'Spontaneous_RH';

%cd(ScriptDir);

SetNum = length(num_measures);
topo = IMG_TYPE;
TopoDir = [AnalysisDir];

%cycle on t-sets. 
for num = num_measures
    TempDir = [TopoDir strcat('t', num2str(num), '/')];
    %here we need to read lowcut and highcut
    file = importdata(strcat(TempDir, 'low_high.txt'));
    lowcut = file{2};
    lowcut = strrep(lowcut, ',', ' ');
    lowcut = strrep(lowcut, ';', '');
    low_precision = split(lowcut, ' ');
    for i=1:length(low_precision)
        index = strfind(low_precision{i}, '.');
        low_precision{i} = length(low_precision{i})-index;
    end
    lowcut = str2num(lowcut);
    highcut = file{3};
    highcut = strrep(highcut, ',', ' ');
    highcut = strrep(highcut, ';', '');
    high_precision = split(highcut, ' ');
    for i=1:length(high_precision)
        index = strfind(high_precision{i}, '.');
        high_precision{i} = length(high_precision{i})-index;
    end
    highcut = str2num(highcut);
    
    %ok, now we are ready to begin the analysis
    %cycle on the lowcut/highcut directory of the analysed set
    for idx = 1:length(lowcut)
        idx;
        lowcut;
        highcut;
        if idx==1
            DataDir = [TempDir strcat('[', num2str(lowcut(idx), strcat('%.',num2str(low_precision{idx}),'f')), '_', num2str(highcut(idx), strcat('%.', num2str(high_precision{idx}), 'f')), ']Hz/')]
        else
            DataDir = [strcat('[', num2str(lowcut(idx), strcat('%.',num2str(low_precision{idx}),'f')), '_', num2str(highcut(idx), strcat('%.', num2str(high_precision{idx}), 'f')), ']Hz/')]
        end
        cd(DataDir);
        mintimeFile = dir('min_points*.txt');
        min_time = mintimeFile.name;
        
        min_time = importdata(mintimeFile.name);
        min_time{586};
        for i = 1:length(min_time)
            min_time{i} = split(min_time{i}, ' ');
            for j = 1:length(min_time{i})
                min_time{i}{j} = strrep(min_time{i}{j}, ',', '');
                min_time{i}{j} = strrep(min_time{i}{j}, '[', '');
                min_time{i}{j} = strrep(min_time{i}{j}, ']', '');
                min_time{i}{j} = str2num(min_time{i}{j});
            end
            min_time{i} = min_time{i}(2:end);
        end
        for i=1:length(min_time)
            min_time{i} = cell2mat(min_time{i});
            min_time{i} = reshape(min_time{i}, [1,length(min_time{i})]);
        end
        min_time = reshape(min_time, [1,length(min_time)]);
        DataDir = [TempDir strcat('[', num2str(lowcut(idx), strcat('%.',num2str(low_precision{idx}),'f')), '_', num2str(highcut(idx), strcat('%.', num2str(high_precision{idx}), 'f')), ']Hz/')]
        %DataDir = [TempDir strcat('[', num2str(lowcut(idx)), '_', num2str(highcut(idx)), ']Hz/')];
        %cd(strcat('../../', ScriptDir));
        cd(ScriptDir)
        WaveHuntTimes; % -- from Trigger Times to Wave Collection
        cd(TempDir)
    end
    cd('../')
end
