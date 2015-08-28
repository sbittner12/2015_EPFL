function [triggeredAmplitudes,frequencies,time,indices,sessions]=trigeredAmplitude(data,norma,triggers,offsetBefore,offsetAfter,winSize,step,srate,windowFunc,chosenSessions,useTrilwiseNorm)
% function trigeredAmplitude - calculate amplitudes for sessioned and
% trigered multichanneled data and concatinates them into a single
% multidimensional array
%
% syntax:
% triggeredAmplitudes = baselineAmplitude(data,triggers,offsetBefore,offsetAfter,winSize,step,srate,windowFunc,chosenSessions)
% [triggeredAmplitudes,frequencies] = baselineAmplitude(data,triggers,offsetBefore,offsetAfter,winSize,step,srate,windowFunc,chosenSessions)
% [triggeredAmplitudes,frequencies,time]=trigeredAmplitude(data,triggers,offsetBefore,offsetAfter,winSize,step,srate,windowFunc,chosenSessions)
%
% data      - matrix or cell of matrices of data (channel-wise recordings)
% norma     - matrix or cell of normalization values
% trigers   - vector or cell of vectors containing triger possitions
% offsetBefore- distance of taken data to the next triger
% offsetAfter- distance of taken data to the previous triger
% winSize   - size of the TRFFT window
% step      - step to move the window
% srate     - sampling rate
% windowFunc- kernel function used for the TRFFT
% chosenSessions - choosen sessions to concatinate data from
%
% triggeredAmplitudes - output amplitudes; first dimension: frequencies, second
% dimension: channels, third dimension: trials
% frequencies - vector containing frequencies
% time        - vector containing time point in respect to the trigger
%

% Tomislav Milekovic, 07/16/2008

if (nargin<10 || isempty(chosenSessions))
    chosenSessions=1:length(data);
end

if (nargin<11 || isempty(useTrilwiseNorm))
    useTrilwiseNorm = false;
end

% Calculate number of channels, number of frequencies and number of time
% bins
noOfChannels=size(data{1},2);
noOfFrequncies=ceil(winSize/2+1);
noOfTimeBins=ceil((offsetBefore+offsetAfter-winSize+2)/step);
time=[];
frequencies=[];

if (~iscell(norma))
    tmpNorma=norma;
    norma=cell(size(data));
    for ii=1:length(data)
        norma{ii}=tmpNorma;
    end
end

% Normalize the window function
windowFunc=windowFunc./norm(windowFunc);

triggers_r=cell(size(triggers));

% Counting the number of trials
noOfTrials=0;
for ii=chosenSessions
    if (~isempty(triggers{ii}))
        goodSt=find(triggers{ii}-offsetBefore>0,1,'first');
        goodEnd=find(triggers{ii}+offsetAfter<size(data{ii},1),1,'last');
        triggers_r{ii}=triggers{ii}(goodSt:goodEnd);
        noOfTrials=noOfTrials+length(triggers_r{ii});
    else
        triggers_r{ii}=[];
    end
end

% Allocating
triggeredAmplitudes=zeros([noOfFrequncies,noOfTimeBins,noOfChannels,noOfTrials]);
indices=zeros([noOfTimeBins noOfTrials]);
sessions=zeros([1 noOfTrials]);

% Preforming the TRFFT amplitude calculation
currentTrial=0;
startTrial = 0;
endTrial = 0;
for sess=chosenSessions
    if (isempty(triggers_r{sess}))
        continue;
    end

    startTrial=currentTrial+1;

    for jj=1:length(triggers_r{sess})
        trialInd = triggers_r{sess}(jj)-offsetBefore:triggers_r{sess}(jj)+offsetAfter;
        if (trialInd(1)<1 || trialInd(end)>size(data{sess},1))
            continue;
        else
            currentTrial=currentTrial+1;
            trailData=data{sess}(trialInd,:);
            [trailFFT,frequencies,time]=trFFT(trailData,winSize,step,srate,windowFunc);

            triggeredAmplitudes(:,:,:,currentTrial)=abs(trailFFT);
            indices(:,currentTrial)=triggers_r{sess}(jj)+time-offsetBefore-1;
            sessions(currentTrial)=sess;
        end
    end

    endTrial=currentTrial;
    
    if (~isempty(norma{sess}))
        noOfTrials = endTrial-startTrial+1;
        if (useTrilwiseNorm)
            for tr = 1:noOfTrials
                normMat=repmat(permute(norma{sess},[1 3 2 4]),[1 noOfTimeBins 1 1]);
                triggeredAmplitudes(:,:,:,startTrial+tr-1)=triggeredAmplitudes(:,:,:,startTrial+tr-1)./sqrt(srate)./normMat;
            end                
        else
            normMat=repmat(permute(norma{sess},[1 3 2 4]),[1 noOfTimeBins 1 noOfTrials]);

            triggeredAmplitudes(:,:,:,startTrial:endTrial)=triggeredAmplitudes(:,:,:,startTrial:endTrial)./sqrt(srate)./normMat;
        end
    else
        triggeredAmplitudes(:,:,:,startTrial:endTrial)=triggeredAmplitudes(:,:,:,startTrial:endTrial)./sqrt(srate);
    end
end

triggeredAmplitudes(:,:,:,endTrial+1:end)=[];
    

if (~isempty(time))
    time=(time-offsetBefore)/srate*1000.0;
end














