
clear
BaseDir = '../..';

ScriptDir = [BaseDir '/SRC/MATLAB/'];
cd(ScriptDir)
%import data
file = readtable('SetData.txt', 'delimiter', '\n');

%search for interesting lines
for i=1:height(file)
    string = file{i,:};
    if startsWith(string, 'num_measures')
        line_measures = string{1};
    end
    if startsWith(string, 'IMG_TYPE')
        line_type = string{1};
    end
    if startsWith(string, 'DIM_X')
        width = string{1};
    end
    if startsWith(string, 'DIM_Y')
        height = string{1};
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
end

%selecting parameters
for idx = 1:length(pixel_size)
    if pixel_size(idx) == '='
        pixel_size(idx)
        inizio = idx+1;
    end
end
fine = idx;
pixel_size = pixel_size(inizio:fine);
pixel_size = strrep(pixel_size, ' ', '');
pixel_size = str2num(pixel_size);

for idx = 1:length(macro_dim)
    if macro_dim(idx) == '='
        inizio = idx+1;
    end
end
fine = idx;
macro_dim = macro_dim(inizio:fine);
macro_dim = strrep(macro_dim, ' ', '');
macro_dim = str2num(macro_dim);
macro_dim;

for idx = 1:length(sampling_time)
    if sampling_time(idx) == '='
        inizio = idx+1;
    end
end
fine = idx;
sampling_time = sampling_time(inizio:fine);
sampling_time = strrep(sampling_time, ' ', '');
sampling_time = str2num(sampling_time);
sampling_time;


for idx = 1:length(width)
    if width(idx) == '='
        inizio = idx+1;
    end
end
fine = idx;
width = width(inizio:fine);
width = strrep(width, ' ', '');
width = str2num(width);
width = width / macro_dim;

for idx = 1:length(height)
    if height(idx) == '='
        inizio = idx+1;
    end
end
fine = idx;
height = height(inizio:fine);
height = strrep(height, ' ', '');
height = str2num(height);
height = height / macro_dim;

%selecting IMG_TYPE
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

%selecting num_measures
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


AnalysisDir = [BaseDir '/ANALYSIS/' IMG_TYPE '/'];

FileName = 'Spontaneous_RH';

cd(ScriptDir);

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
        cd(strcat('../../', ScriptDir));
        WaveHuntTimes;
        cd(TempDir)
    end
    cd('../')
end
