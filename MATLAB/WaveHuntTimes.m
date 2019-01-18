%% Load analysis params and information about recording dataset.
%
version = 'v2';
OutputDir = [DataDir 'WAVES_' version '/'];

setParamsAndOptions;
readPixelData;

%readPixelData;
ChLabel = PixelLabel;
ChannelSet = [1:nCh];

%% [GDB] Set the Path for the results
%
if exist(OutputDir,'dir') == 0
    mkdir(OutputDir);
end

FilesToDelete = [OutputDir 'ClusteredWave_*'];
if ~isempty(dir([OutputDir 'ClusteredWave_*']))
     delete(FilesToDelete);
    % (to be sure that the new execution replace the previous results)
end

fid = fopen([OutputDir FileName '.out'], 'wt');
fprintf(fid,'\n*** FileName: %s\n', FileName);

%% [GDB] EXCLUDE the OUTLIER CHANNELS
%
nCh= nPixel;

%% Setting additional constants used throughout the script (colormap)
%
MAX_ABS_TIMELAG = Options.UpTimeLags.MaxAbsTimeLag;

TL_RANGE = [-300 300];

RED_BLUE = 1;
CM = jet();
if RED_BLUE
   CM = gray(32);
   CM(:,3) = 1;
   Red = flipud(gray(32));
   Red(2:end,1) = 1;
   CM = [CM; Red(2:end,:)];
end
CM(1,:) = [1 1 1]*0.25;



%% ExpectedTrans i.e. estimate of the Number of Waves
%ExpectedTrans used to estimate/optimize IWI
TransPerCh = nTrans;
%ExpectedTrans = median(TransPerCh);
%"TransPerCh(find(TransPerCh))" to account for the case of a subset of electrodes

eExpectedTrans(1) = median(TransPerCh(find(TransPerCh)));

[Active,nA] = sort(TransPerCh(find(TransPerCh)),'descend');
eExpectedTrans(2) = median(Active(1:8)); % consider THE 8 MOST ACTIVE CHANNELS only

eExpectedTrans(3) = mean(TransPerCh(find(TransPerCh)))+std(TransPerCh(find(TransPerCh)));
% mix median and std? maybe mean...?

sel1=1;
ExpectedTrans=eExpectedTrans(sel1);

fprintf('\n--- Estimate of Expected Transitions ---');
for i =1:1:numel(eExpectedTrans)
    if i==sel1; comment='******'; else comment=''; end;
    fprintf('\nExpectedTrans(%d): %d\t%s', i,round(eExpectedTrans(i)),comment);
end
fprintf('\n\n')


%% Compute the time lags matrix looking for optimal MAX_ABS_TIMELAG...
%  (depends on the distance between electrodes in the array)
%
DeltaTL = diff(UpTrans);
OneMoreLoop = 1;
while OneMoreLoop
   WnW = DeltaTL<=MAX_ABS_TIMELAG; % parameter initialized in SetParametersAndOptions
   ndxBegin = find(WnW==1,1,'first');
   nw = 0;
   clear('Wave','WaveUnique','WaveSize','WaveTime');
   while ndxBegin < numel(DeltaTL)
      ndxEnd = find(WnW(ndxBegin:end)==0,1,'first')+ndxBegin-1; % isolated transitions are EXCLUDED
      if isempty(ndxEnd)
         ndxEnd = numel(UpTrans);
      end
      nw = nw + 1;
      Wave(nw).ndx = ndxBegin:ndxEnd;
      WaveUnique(nw) = numel(Wave(nw).ndx) == numel(unique(ChLabel(Wave(nw).ndx)));
      WaveSize(nw) = numel(Wave(nw).ndx);
      WaveTime(nw) = mean(UpTrans(Wave(nw).ndx));
      if ndxEnd == numel(UpTrans)
         ndxBegin = numel(DeltaTL);
      else
         ndxBegin = find(WnW(ndxEnd:end)==1,1,'first')+ndxEnd-1;
      end
   end
   OneMoreLoop = 0;
   if min(WaveUnique) == 0 % the first scan requests that all the waves are UNIQUE
      %fprintf('MaxAbsTimelag too large: %f -> %f\n',MAX_ABS_TIMELAG,MAX_ABS_TIMELAG*0.75);
      %fprintf(fid,'MaxAbsTimelag too large: %f -> %f\n',MAX_ABS_TIMELAG,MAX_ABS_TIMELAG*0.75);
      MAX_ABS_TIMELAG = MAX_ABS_TIMELAG*0.75;
      OneMoreLoop = 1;
   else
%       if max(WaveSize) < numel(RecordingSet)
%          fprintf('MaxAbsTimelag too small: %f -> %f\n',MAX_ABS_TIMELAG,MAX_ABS_TIMELAG*1.25);
%          MAX_ABS_TIMELAG = MAX_ABS_TIMELAG*1.125;
%          OneMoreLoop = 1;
%       end
   end
end % while OneMoreLoop


%% Wave Hunt -- step 1: Compute the Time Lags, Find Unique Waves --> WaveCollection1
% [GDB]
%
WaveColl(1).nWaves=length(WaveSize);

X=[1:1:nCh];
dX = X(2) - X(1);
for i=X
    WaveColl(1).NumberOf(i)=length(find(WaveSize==i));
end

%% IWI distribution and histogram
%
IWI = diff(WaveTime);          % IWI = Inter-Wave-Interval...

%%% histogram(X,edges) sorts X into bins with the bin edges specified by 
%%% the vector, edges. Each bin includes the left edge, but does not include 
%%% the right edge, except for the last bin which includes both edges.
%figure
%hIWI=histogram(IWI,[0:0.250:round(max(IWI))]); 
% 250ms is the MINIMUM IWI to be compatible with a frequency within 4Hz
%eWC2=sum(hIWI.Values(2:end)); % estimates the number of slow waves (frequency <4Hz)
% This number should be comparable with the number of elements in WaveCollection(2)


%% Recollect small-waves in full waves (involving a wider area).

MAX_IWI = max(IWI)*0.5
%MAX_IWI = median(IWI); % ... another estimate (NOT TAKING INTO ACCOUNT the 4Hz LIMIT)
medianIWI = median(IWI(find(IWI>=0.250))); % estimate TAKING INTO ACCOUNT the 4Hz LIMIT

%%ExpectedTrans=eWC2;
ACCEPTABLE_REJECTION_RATE = 0.1;

OneMoreLoop = 1;
while OneMoreLoop
   WnW = IWI<=MAX_IWI;
   ndxBegin = find(WnW==1,1,'first');
   ndxBegin;
   numel(IWI);
   nw = 0;
   clear('FullWave','FullWaveUnique','FullWaveSize','FullWaveTime','FullWaveUniqueSize');
   while ndxBegin < numel(IWI)
      ndxEnd = find(WnW(ndxBegin:end)==0,1,'first')+ndxBegin-1;
      if isempty(ndxEnd)
         ndxEnd = numel(WaveTime);
      end
      nw = nw + 1;
      FullWave(nw).ndx = [Wave(ndxBegin:ndxEnd).ndx];
      FullWaveUnique(nw) = numel(FullWave(nw).ndx) == numel(unique(ChLabel(FullWave(nw).ndx)));
      FullWaveSize(nw) = numel(FullWave(nw).ndx);
      FullWaveUniqueSize(nw) = numel(unique(ChLabel(FullWave(nw).ndx)));
      FullWaveTime(nw) = mean(UpTrans(FullWave(nw).ndx));
      FullWaveUnique;
      if ndxEnd == numel(WaveTime)
         ndxBegin = numel(IWI);
      else
         ndxBegin = find(WnW(ndxEnd:end)==1,1,'first')+ndxEnd-1;
         for j = ndxEnd+1:ndxBegin-1
            nw = nw + 1;
            FullWave(nw).ndx = Wave(j).ndx;
            FullWaveUnique(nw) = numel(FullWave(nw).ndx) == numel(unique(ChLabel(FullWave(nw).ndx)));
            FullWaveSize(nw) = numel(FullWave(nw).ndx);
            FullWaveUniqueSize(nw) = numel(unique(ChLabel(FullWave(nw).ndx)));
            FullWaveTime(nw) = mean(UpTrans(FullWave(nw).ndx));
         end
      end
   end
   BadWavesNum = numel(find(FullWaveUnique==0));
   fprintf('Max wave size: %d, Bad waves: %d (%.3g%%), Good waves: %d\n',...
       max(FullWaveSize),BadWavesNum,BadWavesNum/numel(FullWaveUnique)*100,numel(FullWaveUnique)-BadWavesNum);
   fprintf('Num. waves: %d, max. non-unique ch.: %d\n',...
       numel(FullWaveUnique),max(FullWaveSize-FullWaveUniqueSize));
   fprintf(fid,'Max wave size: %d, Bad waves: %d (%.3g%%), Good waves: %d\n',...
       max(FullWaveSize),BadWavesNum,BadWavesNum/numel(FullWaveUnique)*100,numel(FullWaveUnique)-BadWavesNum);
   fprintf(fid,'Num. waves: %d, max. non-unique ch.: %d\n',...
       numel(FullWaveUnique),max(FullWaveSize-FullWaveUniqueSize));
   
   OneMoreLoop = 0;
   if numel(FullWaveUnique) <= ExpectedTrans % If not we have an artifactual amplification of small waves...
      if BadWavesNum/numel(FullWaveUnique) > ACCEPTABLE_REJECTION_RATE;
         if min(FullWaveUnique) == 0 % at lest a Wave non-unique
            fprintf('MaxIWI too large: %f -> %f\n',MAX_IWI,MAX_IWI*0.75);
            fprintf(fid,'MaxIWI too large: %f -> %f\n',MAX_IWI,MAX_IWI*0.75);   
            MAX_IWI = MAX_IWI*0.75;
            OneMoreLoop = 1;
         else % only unique waves
            %if max(WaveSize) < numel(RecordingSet) % at least one wave MUST BE GLOBAL (i.e. involving the full set of electrodes)
            if max(WaveSize) < numel(ChannelSet) % at least one wave MUST BE GLOBAL (i.e. involving the full set of electrodes)
               fprintf('MaxIWI too small: %f -> %f\n',MAX_IWI,MAX_IWI*1.25);
               fprintf(fid,'MaxIWI too small: %f -> %f\n',MAX_IWI,MAX_IWI*1.25); 
               MAX_IWI = MAX_IWI*1.25;
               OneMoreLoop = 1;
            end
         end
      end %else...no more loop
   end %else...no more loop
end % while OneMoreLoop

fprintf('\nExpected waves     : %d\n', round(ExpectedTrans));
fprintf('Reconstructed waves: %d\n', numel(FullWaveUnique)); 
fprintf(fid,'\nExpected waves     : %d\n', round(ExpectedTrans));
fprintf(fid,'Reconstructed waves: %d\n', numel(FullWaveUnique)); 
%... but still several non-unique waves, i.e. waves with repeated channels

totFullWaveTrans=sum(FullWaveSize); 


%% Remove from waves nonunique channels...
%

FullWaveBackup = FullWave;
TOTNonUnique=numel(find(FullWaveUnique==0));

fprintf('\tWaves to clean     : %d\n', numel(find(FullWaveUnique==0)));
fprintf('\t(Max. non-unique ch.: %d)\n', max(FullWaveSize-FullWaveUniqueSize));
fprintf(fid,'\tWaves to clean     : %d\n', numel(find(FullWaveUnique==0)));
fprintf(fid,'\t(Max. non-unique ch.: %d)\n', max(FullWaveSize-FullWaveUniqueSize));

listOfnw=[]; % List of "critical" waves: non-unique but possibly hiding multiple unique waves
nAdded=0; % total number of new waves added to the collection 
          % (useful also to point back to the original wave indexing in FullWaveBackup)
nNew=[];  
nSegments=[]; % (how many segments non-unique waves are made of)
nSize=[]; % (size of the segments)

alert(1:6)=0;

%THR = 3*wmDelta;

% *** TWO THRESHOLDS ***
% THR1 = maxRisingTime, estimated as 1/2 of mean(UpStateDuration)
% THR2 = 250ms = min IWI (frequency limit) 

meanUP=0;
%for i = ChannelSet
%    meanUP=meanUP+results(i).UpState.mean_duration;
%end
%meanUP=meanUP/numel(ChannelSet);
meanUP = 0.200
THR1 = 0.5*meanUP; %RisingTime
THR2 = 0.250; %4Hz LIMIT

% Save plot of "critical" non-unique waves, before & after the treatment
% Create the folder for the storage of plots if non-existing, and empty if existing
FilesToDelete = [OutputDir 'nW/nW*'];
if exist([OutputDir 'nW/'],'dir') == 0
    mkdir([OutputDir 'nW/']);
elseif ~isempty(dir([OutputDir 'nW/']))
     delete(FilesToDelete);    
end

nw=1;
NonUnique=0; % COUNTER for NonUnique Waves
while nw<=numel(FullWaveUnique)
%for nw = 1:numel(FullWaveUnique) 
    % numel(FullWaveUnique) = number of waves in the collection
    % FullWaveUnique==0 if the waves contains repeated channels
    
    NewFullWave=[];NewFullWaveUnique=[];NewFullWaveSize=[];NewFullWaveTime=[];seg=[]; % empty arrays
    nCHS=[];
    
    if FullWaveUnique(nw) == 0
        NonUnique=NonUnique+1;
              
        chs = ChLabel(FullWave(nw).ndx);
        nCHS(1)=numel(chs);
        wndx = 1:1:numel(chs);
% FullWave(nw).ndx = index of the transitions that make up the wave(nw)
% ChLabel() = the corresponding channel for each transition 
% --> in this list there are repeated channels
        nch = hist(chs,1:numel(RecordingSet)); % WARNING! not numel(ChannelSet) even in the case of a subset
 
% --- (1) First CHANNEL CLEANING: check the TIME DISTANCE BETWEEN REPETITIONS
        rep=find(nch>1);
        
        k=1;
        while k<=numel(rep)
        %for i=rep
            i=rep(k);
        
            timeStamp=UpTrans(FullWave(nw).ndx(chs==i));
            idx=FullWave(nw).ndx(chs==i);
            if diff(timeStamp)<0.125
                % open a time-window around each occurrence of i and check
                % how many members of the clan are in the time-window 
                fprintf('Channel %d has close repetitions\n',i);
                nClan=[];
                for j=1:1:nch(i)
                    window=[timeStamp(j)-0.125 timeStamp(j)+0.125];
                    %window=[timeStamp(j)-0.125 timeStamp(j)+0.125];
                    inWindow=FullWave(nw).ndx(UpTrans(FullWave(nw).ndx)>window(1) & UpTrans(FullWave(nw).ndx)<window(2)); 
                    inWindow=setdiff(inWindow,idx(j));
                    nClan(j)=numel(intersect(ChLabel(inWindow),arrayMask{i}{2}));
                end
                
                if ~isempty(find(nClan==0)) % to be more generic, consider the case min(nClan)
                    [FullWave(nw).ndx,a]=setdiff(FullWave(nw).ndx,idx(find(nClan==0)));
                    b=setdiff(wndx,a);
                    chs(b)=[];
                    fprintf('%d repetition(s) deleted\n',numel(b));
                    nch=hist(chs,1:numel(RecordingSet));
                    if(nch(i)>1) 
                        fprintf('...Now, re-do channel %d\n',i);
                    else
                        k=k+1;
                    end
                else
                    fprintf('CloseRepProblem for channel i = %d NOT SOLVED\n',i);
                    k=k+1;
                end
            else
                k=k+1;
            end
            
        end
        nCHS(2)=numel(chs); % update the number of transitions (to keep count of removed transitions)

        
        
        
        
%[GDB] Check the time difference between transitions 
        delta=diff(UpTrans(FullWave(nw).ndx));
        delta = cast(delta, 'double');
        mD=mean(delta);
        stdD=std(delta);
       
% % --- PLOT non-unique wave, in order to CHECK how it is treated ---
%         f=figure('visible','off');
%         
%         % --- deltaHistogram
%         subplot(1,4,1) 
%         nBin=round(numel(delta)/3);
%         hDelta=histogram(delta,nBin);
%         xlabel('delta (s)');
%         hold on
%         line([mD mD],[0 max(hDelta.Values)],'Color','r')
%         line([mD+stdD mD+stdD],[0 max(hDelta.Values)],'Color',[0.75 0 0],'LineStyle','--')
%         line([mD+2*stdD mD+2*stdD],[0 max(hDelta.Values)],'Color',[0.5 0 0],'LineStyle','--')
%         line([mD+3*stdD mD+3*stdD],[0 max(hDelta.Values)],'Color',[0.25 0 0],'LineStyle','--')
%         hold off
%         set(gca,'Position',[0.05 0.11 0.2 0.8150]); %[left bottom width height]
%                
%         % --- UpTrans
%         subplot(1,4,2:4)        
%         plot(UpTrans(FullWave(nw).ndx),ChLabel(FullWave(nw).ndx),'k.')
%         set(gca, 'YLim',[0 nCh]);
%         hold on
%         xlabel('Time (s)');
%         ylabel('Channel');
%         set(gca,'YLim',[0 nCh]); 
%         set(gca,'Position', [0.35 0.11 0.6 0.8150]); %[left bottom width height]
%         
%         dim0 = [0.735 0.9 0.1 0.1];
%         str0 = {['nW = ' num2str(nw) '   [ ' num2str(NonUnique) '/' num2str(TOTNonUnique) ' ]']};
%         annotation('textbox',dim0,'String',str0,'FitBoxToText','on','BackgroundColor','none');

%             XRange = [min(UpTrans(FullWaveBackup(i).ndx)) min(UpTrans(FullWaveBackup(i).ndx))+2*MAX_IWI];
%             XRange = XRange + [-1 +1]*diff(XRange)*0.1;
%             set(gca,'XLim',XRange)            
%             print('-dpdf', [OutputDir 'nW/nW' num2str(i)])
% -----------------------------------------------------------------               


% Look for CANDIDATE SEGMENTS  
        THR=mD+3*stdD;
        B=find(delta>THR); % larger sepration between transitions 
                           % --> possible BREAK of the non-unique wave
        
        if ~isempty(B) % --- IDENTIFY SEGMENTS
            % SEGMENTS i.e. candidate NEW waves
            segments=numel(B)+1;
            istart=1;
            i=1;
            while i<segments
                seg(i).ndx=(istart:B(i));
                seg(i).chs=chs(istart:B(i));
                istart=B(i)+1;
                i=i+1;
            end
            seg(i).ndx=(istart:numel(chs));
            seg(i).chs=chs(istart:end);

            % --- CLEAN SEGMENTS---
            delSeg=[];
            for i=1:1:segments % CHECK if SEGMENTS have to be cleaned
                
                % 1) non-unique segments
                if(numel(unique(seg(i).chs))~=numel(seg(i).chs)) % the subset of transition is NOT unique 
                    % --> 'CLEAN', i.e. scan the transition sequence and
                    % find repetition of channels
                    fprintf('ALERT(1) nw=%d, seg=%d\n', nw,i);
                    fprintf(fid,'ALERT(1) nw=%d, seg=%d\n', nw,i);
                    alert(1)=alert(1)+1;
                    
                    %j=1;
                    %while j<=numel(seg(i).chs) % while-loop beacause the numel(chs) is updated in the loop
                    delCh=[];
                    for j=1:numel(seg(i).chs)
                        listNdx=(find(seg(i).chs==seg(i).chs(j)));
                        if(numel(listNdx)~=1)
%                         % Keep the first occurrance, remove the others (interpreted as noise)
%                             seg(i).chs(listNdx(2:end))=[];
%                             seg(i).ndx(listNdx(2:end))=[];

%                         % Keep the occurrance which is the closest to the other occurances in the "clan"
                            t0Clan=0; nClan=0;
                            for k = arrayMask{seg(i).chs(j)}{2}
                  
                                %disp(k)
                                %disp(UpTrans(FullWave(nw).ndx(seg(i).ndx(seg(i).chs==k))))
                                tClan=UpTrans(FullWave(nw).ndx(seg(i).ndx(seg(i).chs==k)));
                                if ~isempty(tClan)
                                    nClan=nClan+1;
                                    t0Clan=t0Clan+mean(tClan);  
                                end                               
                            end
                            t0Clan=t0Clan/nClan; 
                            %PROVIAMO AD AGGIUNGERLA NOI... ATTENZIONE!
                            %t0Clan = mean(t0Clan);
                            tCh=UpTrans(FullWave(nw).ndx(seg(i).ndx(seg(i).chs==seg(i).chs(j))));
                            [t0Ch,index]=min(abs(tCh-t0Clan));
                            delCh=[delCh setdiff(listNdx,listNdx(index))];                            
                        end 
                        %j=j+1;
                    end 
                    seg(i).chs(unique(delCh))=[];
                    seg(i).ndx(unique(delCh))=[];
                end 
                
                % 2) channels non-LocallyConnected (see arrayMask)
                if numel(seg(i).chs)<=5 % 5 is min(numel(arrayMask{:}{2})
                    %for j=seg(i).chs
                    delList=[];
                    for j=1:1:numel(seg(i).chs)
                        k=seg(i).chs(j);
                        %disp(['k = ' num2str(k)])
                        %disp(['intersect: ' num2str(intersect(arrayMask{k}{2},setdiff(seg(i).chs,k)))])
                        if isempty(intersect(arrayMask{k}{2},setdiff(seg(i).chs,k)))
                          delList=[delList j];  
                        end 
                    end
                    if ~isempty(delList)
                        fprintf('ALERT(2) nw=%d, seg=%d\n', nw,i);
                        alert(2)=alert(2)+1;
                        seg(i).chs(delList)=[];
                        seg(i).ndx(delList)=[];
                    end
                end
                
                % PREPARE TO REMOVE EMPTY SEGMENTS
                if isempty(seg(i).ndx)
                   fprintf('ALERT(3) nw=%d, seg=%d\n', nw,i);
                   fprintf(fid,'ALERT(3) nw=%d, seg=%d\n', nw,i);
                   alert(3)=alert(3)+1;
                   delSeg=[delSeg i];
                end
            
            end
                
            % REMOVE EMPTY SEGMENTS
            if ~isempty(delSeg)
                seg(delSeg)=[];
            end
            segments=numel(seg); % update the value in 'segments' = number of segments

            % coalescence of segments if no repetitions with adjacent(s) one(s)
            % N.B. a "small" intersection is admitted
            %for i=1:1:(segments-1)   
            i=1;
            while i<=(segments-1)
                %if isempty(intersect(seg(i).chs,seg(i+1).chs))
                if numel(intersect(seg(i).chs,seg(i+1).chs))<=floor(1/4*min(numel(seg(i).chs),numel(seg(i+1).chs)))
                    % CANDIDATE SEGMENTS for COALESCENCE
                    % check also if distance between segments'border is smaller than 250ms = 1/4Hz
                    distance=UpTrans(FullWave(nw).ndx(seg(i+1).ndx(1))) - UpTrans(FullWave(nw).ndx(seg(i).ndx(end)));
                    
                    if distance>=0.250
                        % FREQUENCY ALERT: distance compatible with SWA, the two segments should be kept separated
                        fprintf('>>> ALERT(4) (Frequency Alert)\n nw=%d, seg=%d AND seg=%d\n', nw,i,i+1 );
                        fprintf(fid,'>>> ALERT(4) (Frequency Alert)\n nw=%d, seg=%d AND seg=%d\n', nw,i,i+1 );
                        alert(4)=alert(4)+1;
                        disp(['--- ' num2str(seg(i).chs)]);
                        disp(['--- ' num2str(seg(i+1).chs)]);
                        fprintf('TOT Transitions: %d\n',numel(seg(i).chs)+numel(seg(i+1).chs));
                        disp(['intersect: ' num2str(intersect(seg(i).chs,seg(i+1).chs))]);
                        fprintf('numel(intesect): %d\n',numel(intersect(seg(i).chs,seg(i+1).chs)));
                        fprintf('threshold: %d\n',floor(1/4*min(numel(seg(i).chs),numel(seg(i+1).chs))));
                        i=i+1; %increment the pointer only if no coalescence is made
                    else
                        % COALESCENCE
                        % The two segments are close enough that can be merged into a single wave
                        % (consider them separated waves would mean the SWA frequency is larger than 4Hz)
                        
                        disp('*** Empty/small intersection for consecutive segments ==> coalescence of segments');
                        fprintf('ALERT(5) nw=%d, seg=%d AND seg=%d\n', nw,i,i+1);
                        fprintf(fid,'*** Empty/small intersection for consecutive segments ==> coalescence of segments\n');
                        fprintf(fid,'ALERT(5) nw=%d, seg=%d AND seg=%d\n', nw,i,i+1);
                        alert(5)=alert(5)+1;
                        disp(['--- ' num2str(seg(i).chs)]);
                        disp(['--- ' num2str(seg(i+1).chs)]);
                        fprintf('TOT Transitions: %d\n',numel(seg(i).chs)+numel(seg(i+1).chs));
                        disp(['intersect: ' num2str(intersect(seg(i).chs,seg(i+1).chs))]);
                        fprintf('numel(intesect): %d\n',numel(intersect(seg(i).chs,seg(i+1).chs)));
                        fprintf('threshold: %d\n',floor(1/4*min(numel(seg(i).chs),numel(seg(i+1).chs))));

                        % COALESCENCE of consecutive SEGMENTS
                        mergedCHS=[seg(i).chs seg(i+1).chs];
                        mergedNDX=[seg(i).ndx seg(i+1).ndx];
                        
                        % CHECK for REPETITIONS (and treat them as usual...
                        % looking at the meanTime in the Clan)
                        delCh=[];
                        for j=1:numel(mergedCHS)
                            listNdx=(find(mergedCHS==mergedCHS(j)));
                            if(numel(listNdx)~=1)
    %                         % Keep the first occurrance, remove the others (interpreted as noise)
    %                             seg(i).chs(listNdx(2:end))=[];
    %                             seg(i).ndx(listNdx(2:end))=[];

    %                         % Keep the occurrance which is the closest to the other occurances in the "clan"
                                t0Clan=0; nClan=0;
                                for k = arrayMask{mergedCHS(j)}{2}
                                    %disp(k)
                                    %disp(UpTrans(FullWave(nw).ndx(seg(i).ndx(seg(i).chs==k))))
                                    tClan=UpTrans(FullWave(nw).ndx(mergedNDX(mergedCHS==k)));
                                    if ~isempty(tClan)
                                        nClan=nClan+1;
                                        t0Clan=t0Clan+tClan;
                                    end
                                end
                                t0Clan=t0Clan/nClan;
                                tCh=UpTrans(FullWave(nw).ndx(mergedNDX(mergedCHS==mergedCHS(j))));
                                [t0Ch,index]=min(abs(tCh-t0Clan));
                                delCh=[delCh setdiff(listNdx,listNdx(index))];                            
                            end 
                        end 
                        
                        mergedCHS(unique(delCh))=[];
                        mergedNDX(unique(delCh))=[];
                        % N.B. unique(delCh) is empty if delCh is empty -->
                        % no cancellation is made in mergedCHS and mergedNDX
                        
                        seg(i).chs=mergedCHS; % COALESCENCE
                        seg(i).ndx=mergedNDX; % COALESCENCE  
                        seg(i+1)=[]; % coalesced segments are at index i, segment at index i+1 is REMOVED 
                        segments=segments-1;
                        
                        % if segments are merged, do not increment the pointer but
                        % check if the new segment can be merged with the following one
                                                                     
                    end % END IF (FREQUENCY ALERT)
                    
                else % consecutive segments intersect too much...
                    i=i+1; %increment the pointer only if no coalescence is made
                end
            end
                    
            if(segments~=numel(seg)) disp('ERROR - Number of Segments');end; % (for doubt's sake)
                
            % $$$$$ N.B. the number of segments has to be updated    
            for i=1:1:segments      
                NewFullWave(i).ndx = FullWave(nw).ndx(seg(i).ndx);
                NewFullWaveUnique(i) = 1; % update the logical value (...we are "cleaning" the waves)
                NewFullWaveSize(i) = numel(NewFullWave(i).ndx);
                NewFullWaveTime(i) = mean(UpTrans(NewFullWave(i).ndx));
                nSize=[nSize NewFullWaveSize(i)];
            end
            nSegments=[nSegments segments];
            
        else % NO SEGMENTS identified -->
             % repeated chiannels are due to noise (and not to the presence of more than one wave)
             % CLEAN the wave, i.e. keep only the first channel occurrance
            fprintf('ALERT(6) nw=%d\n', nw); 
            fprintf(fid,'ALERT(6) nw=%d\n', nw);
            alert(6)=alert(6)+1;
            
            % ..no more so sure that the while loop is correct when elemets are deleted (and not inserted) in the array
%             j=1;
%             while j<=numel(chs) % while-loop beacause the numel(chs) is updated in the loop 
%                 listNdx=(find(chs==chs(j)));
%                 if(numel(listNdx)~=1)
%                 % Keep the first occurrance, remove the others (interpreted as noise)
%                     chs(listNdx(2:end))=[];
%                     FullWave(nw).ndx(listNdx(2:end))=[];
%                 end
%                 j=j+1;
%             end 
            
            delCh=[];
            for j=1:numel(chs)
                listNdx=(find(chs==chs(j)));
                if(numel(listNdx)~=1)
% Keep the occurrance which is the closest to the other occurances in the "clan"
                    t0Clan=0; nClan=0;
                    for k = arrayMask{chs(j)}{2}
                        %disp(k)
                        tClan=UpTrans(FullWave(nw).ndx(chs==k));
                        if ~isempty(tClan)
                            %nClan=nClan+1;
                            %t0Clan=t0Clan+tClan;
                            nClan=nClan+numel(tClan); % take into account the case the CLAN has repetions
                            t0Clan=t0Clan+sum(tClan); % take into account the case the CLAN has repetions  
                        end                               
                    end
                    t0Clan=t0Clan/nClan; 
                    tCh=UpTrans(FullWave(nw).ndx(chs==chs(j)));
                    [t0Ch,index]=min(abs(tCh-t0Clan));
                    delCh=[delCh setdiff(listNdx,listNdx(index))];                            
                end 
            end 
            chs(unique(delCh))=[];
            FullWave(nw).ndx(unique(delCh))=[];
       
            % wave is "cleaned"; store and plot the updated wave
            NewFullWave.ndx = FullWave(nw).ndx;
            NewFullWaveUnique = 1; % update the logical value (...we are "cleaning" the waves)
            NewFullWaveSize = numel(NewFullWave.ndx);
            NewFullWaveTime = mean(UpTrans(NewFullWave.ndx));
            
            nSize=[nSize NewFullWaveSize];
            nSegments=[nSegments 1];
        end

        
        % --- REPLACE CurrentWave with NewWave(s) 
        % [its segments or its 'cleaned' version]
        if nw~=1 
            FullWave=[FullWave(1:nw-1) NewFullWave FullWave(nw+1:end)];
            FullWaveUnique = [FullWaveUnique(1:nw-1) NewFullWaveUnique FullWaveUnique(nw+1:end)];
            FullWaveSize = [FullWaveSize(1:nw-1) NewFullWaveSize FullWaveSize(nw+1:end)];
            FullWaveTime = [FullWaveTime(1:nw-1) NewFullWaveTime FullWaveTime(nw+1:end)];
        else
            FullWave=[NewFullWave FullWave(nw+1:end)];
            FullWaveUnique = [NewFullWaveUnique FullWaveUnique(nw+1:end)];
            FullWaveSize = [NewFullWaveSize FullWaveSize(nw+1:end)];
            FullWaveTime = [NewFullWaveTime FullWaveTime(nw+1:end)];
        end
        
% --- PLOT the wave after the treatment --- 
        %mask={[1 0 0] [0 1 0] [0 0 1]}; % red, green, blue [alternance of 3 colors]
        % C={'r','b','g','m','c','y','k'};
        %for j=1:1:numel(NewFullWave)
            %c=mask{mod(j,3)+1}; % selection from a set of 3 colors)
            %plot(UpTrans(FullWave(nw-1+j).ndx),ChLabel(FullWave(nw-1+j).ndx),...
                %'color',c,'marker','o','linestyle','none');
        %end
        %print('-dpdf', [OutputDir 'nW/nW' num2str(nw)])
        %close(f);
% -----------------------------------------        

        % --- INCREMENT the pointer
        if numel(NewFullWave)>1 % SEGMENTS ARE NEW WAVES          
            nAdded=nAdded+numel(NewFullWave)-1;
            nw=nw+numel(NewFullWave); % increment (point at the next wave)
        else % no segments identified, the current wave is a New Wave, because it has been cleaned
            nw=nw+1; % increment (point at the next wave)
        end
                  
        
    else nw=nw+1; % increment (point at the next wave) [current wave is already unique]      
    end
    
    if ~isempty(NewFullWave)
        nNew=[nNew numel(NewFullWave)-1];
    end
    
end % END WHILE LOOP  
fprintf('alert = [%d,%d,%d,%d,%d,%d]\n', alert); 
fprintf(fid,'alert = [%d,%d,%d,%d,%d,%d]\n', alert); 


%% Remove from waves nonunique channels...
%
if exist('OLDv','var')
    FullWaveBackup = FullWave;

    fprintf('\tWaves to clean     : %d\n', numel(find(FullWaveUnique==0)));
    fprintf('\t(Max. non-unique ch.: %d)\n', max(FullWaveSize-FullWaveUniqueSize));
    fprintf(fid,'\tWaves to clean     : %d\n', numel(find(FullWaveUnique==0)));
    fprintf(fid,'\t(Max. non-unique ch.: %d)\n', max(FullWaveSize-FullWaveUniqueSize));

    listOfnw=[]; % List of "critical" waves: non-unique but possibly hiding multiple unique waves
    for nw = 1:numel(FullWaveUnique) 
        % numel(FullWaveUnique) = number of waves in the collection
        % FullWaveUnique==0 if the waves contains repeated channels
       if FullWaveUnique(nw) == 0
          chs = ChLabel(FullWave(nw).ndx);
    % FullWave(nw).ndx = index of the transitions that make up the wave(nw)
    % ChLabel() = the corresponding channel for each transition 
    % --> in this list there are repeated channels

          nch = hist(chs,1:numel(RecordingSet)); % WARNING! not numel(ChannelSet) even in the case of a subset

          ndx = find(nch>1); % non-unique channels

    %[GDB] CHECK the amount of channels and transitions suppressed from the non-unique wave      

          C=length(ndx); % amount of suppressed channels from non-unique wave
          L1=length(chs);
          [chs,ndx] = setdiff(chs,ndx); % unique channels
    % "ndx" now contains the index of transition in the wave corresponding to the ordered list of unique channels
          L2=length(chs);
          T=L1-L2; % amount of suppressed transitions from non-unique wave

          if C>ceil(nCh/2)
              fprintf('\n*** nW = %d\tsuppressed channels = %d \tsuppressed transitions = %d',nw,C,T);
              fprintf(fid,'\n*** nW = %d\tsuppressed channels = %d \tsuppressed transitions = %d',nw,C,T);    
              listOfnw=[listOfnw, nw];
          end

          if T>ceil(nCh/2)
              fprintf('\n### nW = %d\tsuppressed channels = %d \tsuppressed transitions = %d',nw,C,T);
              fprintf(fid,'\n### nW = %d\tsuppressed channels = %d \tsuppressed transitions = %d',nw,C,T);
              listOfnw=[listOfnw, nw];
          end

          FullWave(nw).ndx = FullWave(nw).ndx(ndx);
    % FullWave(nw) is now a "reduced" wave with only transitions occurred in unique (ie. non repeated) channels
    % N.B. REPEATED CHANNELS ARE SUPPRESSED ! ! !  (and their transitions removed)

          FullWaveUnique(nw) = 1; % update the logical value (...we are "cleaning" the waves)
          FullWaveSize(nw) = numel(FullWave(nw).ndx);
          FullWaveTime(nw) = mean(UpTrans(FullWave(nw).ndx));
       end
    end

    % Save plot of "critical" non-unique waves = hidden unique waves
    % Create the folder for the storage of plots if non-existing, and empty if existing
    FilesToDelete = [OutputDir 'nW/nW*'];
    if exist([OutputDir 'nW/'],'dir') == 0
        mkdir([OutputDir 'nW/']);
    elseif ~isempty(dir([OutputDir 'nW/']))
         delete(FilesToDelete);    
    end

    listOfnw=unique(listOfnw); % remove duplicate indexes
    fprintf('\n\nCritical waves: %d\n',numel(listOfnw));
    fprintf(fid,'\n\nCritical waves: %d\n',numel(listOfnw));
    if(~isempty(listOfnw))
        for i = listOfnw
%           figure
%           plot(UpTrans(FullWaveBackup(i).ndx),ChLabel(FullWaveBackup(i).ndx),'k.')
%           hold on
%           plot(UpTrans(FullWave(i).ndx),ChLabel(FullWave(i).ndx),'ro')
%           xlabel('Time (s)');
%           ylabel('Channel');
%           dim0 = [0.2175 0.825 0.1 0.1];
%           str0 = {['nW = ' num2str(i)]};
%           annotation('textbox',dim0,'String',str0,'FitBoxToText','on','BackgroundColor','none');
    %       XRange = [min(UpTrans(FullWaveBackup(i).ndx)) min(UpTrans(FullWaveBackup(i).ndx))+2*MAX_IWI];
    %       XRange = XRange + [-1 +1]*diff(XRange)*0.1;
    %       set(gca,'XLim',XRange)
          print('-dpdf', [OutputDir 'nW/nW' num2str(i)])
        end
    else
        for i = [1:1:numel(FullWave)]
            if mod(i,20)==0
%                 figure
%                 plot(UpTrans(FullWaveBackup(i).ndx),ChLabel(FullWaveBackup(i).ndx),'k.')
%                 hold on
%                 plot(UpTrans(FullWave(i).ndx),ChLabel(FullWave(i).ndx),'ro')
%                 xlabel('Time (s)');
%                 ylabel('Channel');
%                 dim0 = [0.2175 0.825 0.1 0.1];
%                 str0 = {['nW = ' num2str(i)]};
%                 annotation('textbox',dim0,'String',str0,'FitBoxToText','on','BackgroundColor','none');
    %             XRange = [min(UpTrans(FullWaveBackup(i).ndx)) min(UpTrans(FullWaveBackup(i).ndx))+2*MAX_IWI];
    %             XRange = XRange + [-1 +1]*diff(XRange)*0.1;
    %             set(gca,'XLim',XRange)            
                print('-dpdf', [OutputDir 'nW/nW' num2str(i)])
            end
        end
    end
    % HOW TO SOLVE THE PROBLEM of CRITICAL NON-UNIQUE WAVES? Possible solutions:
    % - reduce IWI (for critical waves only, of for the full waves collection?)
    % - keep the firts occurrance of each channel (to preserve at least one wave)
    % - check the time distance between transitions at the same channel. 
    %   If larger than a threshold, try to "cut" the sequence in more than one wave

end

%% Wave Hunt -- step 2: Coalescence of 'Short' Waves, Rejection of 'Small' Waves --> WaveCollection2
% [GDB]
%
WaveColl(2).nWaves=length(FullWaveSize);

X=[1:1:nCh];
dX = X(2) - X(1);
for i=X
    WaveColl(2).NumberOf(i)=length(find(FullWaveSize==i));
end

value=sum(WaveColl(2).NumberOf(ceil(nCh/2):nCh))/WaveColl(2).nWaves;


%% Find Maxima and Minima in the histogram of FullWaveSize

[hh,xx]=hist(FullWaveSize,(1:1:nCh));

for i=1:1:numel(hh)-1
    if abs(hh(i)-hh(i+1))==1
        hh(i+1)=hh(i);
    end   
end
RI=diff(hh);

bMax=[];
bMin=[];
for i=1:1:numel(RI)-1
    %disp([RI(i), RI(i+1)])
    if sign(RI(i))*sign(RI(i+1))~=1
        if sign(RI(i))~=0
            if sign(RI(i))>sign(RI(i+1))
                %fprintf('BIN %d is a MAXIMUM\n',i+1);
                bMax=[bMax i+1];
            else
                %fprintf('BIN %d is a MINIMUM\n',i+1);
                bMin=[bMin i+1];
            end
        else
            %fprintf('BIN %d has the same value than %d\n',i+1,i);
            if ~isempty(intersect(i,bMin))
                bMin=[bMin i+1];
            end
            if ~isempty(intersect(i,bMax))
                bMin=[bMin i+1];
            end               
        end  
    end
end
if hh(end)>hh(end-1)
    fprintf('BIN %d is a MAXIMUM\n',numel(hh));
    bMax=[bMax i+1];
end


%% Remove small waves and those rejected...
%
%MIN_CH_NUM = 4; 
% If a sub-area is selected, this number could be too large...

if SelectedArea == 'R'
    MIN_CH_NUM = 3;
elseif SelectedArea == 'P'
    MIN_CH_NUM = 2;
% maybe a parameter ad hoc for each cortical area? ...to be OPTIMISED
else
    %MIN_CH_NUM = 4; % too few for 32 electrodes?
    %MIN_CH_NUM = 28; 
%--- [GDB] introduce a requirement of GLOBALITY
%-- at least half of the array:
    if value >= 0.5
        MIN_CH_NUM = ceil(nCh/2);
    % add a requirement concerning an ACCEPTABLE_REJECTION_RATE = 0.5 
    % i.e. at least 50% of the waves in the collection must involve the 50%
    % or more of the electrode array
    else
        R=[ceil(nCh/2):-1:1];
        for i=R
            value=sum(WaveColl(2).NumberOf(i:nCh))/WaveColl(2).nWaves;
            if value >= 0.5
                MIN_CH_NUM = i;
                break
            end
        end
    end    
end
    
if version=='v2' %CAMBIO
    MIN_CH_NUM = 300;
end
    
ndx = find(FullWaveUnique==1 & FullWaveSize >= MIN_CH_NUM);
RejectedWaves = numel(FullWaveUnique) - numel(ndx); % rejected beacuse too small 
% (i.e. involving too few channels)
Wave = FullWave(ndx);
WaveUnique = FullWaveUnique(ndx);
WaveSize = FullWaveSize(ndx);
WaveTime = FullWaveTime(ndx);

fprintf('\nMinimum number of channels: %d\n',MIN_CH_NUM);
fprintf('Rejected waves: %d (too small, i.e. involve less than minimum number of channels)\n',RejectedWaves);
fprintf('Accepted waves: %d\n',numel(Wave));
fprintf(fid,'\nMinimum number of channels: %d\n',MIN_CH_NUM);
fprintf(fid,'Rejected waves: %d (too small, i.e. involve less than minimum number of channels)\n',RejectedWaves);
fprintf(fid,'Accepted waves: %d\n',numel(Wave));
WaveColl(3).nWaves=numel(Wave);

TimeLagMatrix = NaN(numel(Wave),numel(RecordingSet)); % TLM initialized with NaN
UpTransNdxMatrix = NaN(numel(Wave),numel(RecordingSet)); % UpTransMatrix initialized with NaN
for k = 1:numel(Wave)
   PXs = [];
   TLMs = UpTrans(Wave(k).ndx);
   CLs = ChLabel(Wave(k).ndx);
   for i=1:numel(CLs)
    PXs(i) = find(pixelID == CLs(i)); 
   end
   %TimeLagMatrix(k,CLs) = TLMs - mean(TLMs); % each wave is centered at the mean time
   TimeLagMatrix(k, PXs) = TLMs - mean(TLMs); % each wave is centered at the mean time
   UpTransNdxMatrix(k,CLs) = Wave(k).ndx;
end
TimeLagRange = [min(TimeLagMatrix(isnan(TimeLagMatrix)==0)) ...
                max(TimeLagMatrix(isnan(TimeLagMatrix)==0))]; % max duration of the waves
TLMtoPlot = TimeLagMatrix;
TLMtoPlot(isnan(TLMtoPlot)==1) = -1; % '-1' replace 'NaN' in the plot

%% Duration of Waves
%

WaveDuration=[];
for i = 1:1:numel(FullWave)
    WaveDuration = [WaveDuration UpTrans(FullWave(i).ndx(end))-UpTrans(FullWave(i).ndx(1))];
end

%% Create an array with the beginning time of each wave
% Saves the array in a 'BeginTime.txt' file in the path folder
BeginTime = [];
EndTime = [];
BeginTimeFile = [DataDir strcat('BeginTime_', topo, '_t', num2str(num),'.txt')]
EndTimeFile = [DataDir strcat('EndTime_', topo, '_t', num2str(num),'.txt')];
F = fopen(BeginTimeFile, 'wt');
G = fopen(EndTimeFile, 'wt');
% Check Begin and End Times files are properly opened
if (F > 0 && G > 0)
       fprintf('Begin and End Times files have been successfully opened\n');
end
for i = 1:1:numel(Wave)
       BeginTime = [BeginTime UpTrans(Wave(i).ndx(1))];
       EndTime = [EndTime UpTrans(Wave(i).ndx(end))];
       fprintf(F, '%f \n', BeginTime(i));
       fprintf(G, '%f \n', EndTime(i));
end
fclose('all');