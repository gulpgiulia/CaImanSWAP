% readPixelData.m
%
% (v0) first attempt to read the outcome of Python (50x40 grid of pixels)

%% LOAD input data

nTrans=[];
for i =1:1:numel(min_time)
    min_time{i};
   nTrans(i)=numel(min_time{i});
end
%nTrans
pixelID=find(nTrans ~=0);
nPixel=numel(find(nTrans ~=0));
nCh = nPixel;
RecordingSet = [1:nCh];
%% CREATE UpTrans
%
UpTrans = [];
PixelLabel = [];

for idx = pixelID
   UpTrans = [UpTrans min_time{idx}];
   PixelLabel = [PixelLabel repmat(idx, 1, numel(min_time{idx}))];
end

[UpTrans,ndx] = sort(UpTrans);
PixelLabel = PixelLabel(ndx);

%figure
%plot(UpTrans,PixelLabel,'k.');
%plot(UpTrans*sampling_time,PixelLabel,'k.');
%xlabel('Time (s)');
%ylabel('Pixel');
%print('-deps2c', [OutputDir 'SequenceOfUpTrans.eps']);  % path for the results
