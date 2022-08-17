function spS=vIRt_SpikingSpectrum(spikeTimes,dataMask) %

if nargin>1
    %% Data masking
    wEpochs.behav=bwconncomp(dataMask.behav);
    % mask epochs with short whisking bouts
    durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=2000;
    dataMask.behav(vertcat(wEpochs.behav.PixelIdxList{~durationThd}))=false;
    wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(durationThd);
    cumulDur=cumsum(cellfun(@numel, wEpochs.behav.PixelIdxList)/1000);
    timeLimitIdx=find(cumulDur>=30,1); %keep only 30s or so of wisking
    if isempty(timeLimitIdx) %then keep all
        timeLimitIdx=numel(wEpochs.behav.PixelIdxList);
    end
   
    % apply mask 
%     wEpochs.behav.PixelIdxList=vertcat(wEpochs.behav.PixelIdxList{1:timeLimitIdx});
    wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(1:timeLimitIdx);
%     if ~isstruct(spikeTimes)
%         spikeTimes=struct('times',spikeTimes);
%     end
%     for epochNum=numel(wEpochs.behav.PixelIdxList):-1:1
%         epochIdx=spikeTimes(1).times>=wEpochs.behav.PixelIdxList{epochNum}(1)/1000 & ...
%             spikeTimes(1).times<=wEpochs.behav.PixelIdxList{epochNum}(end)/1000;
%         spikeTimes(epochNum).times=spikeTimes(1).times(epochIdx);
%     end
    for epochNum=1:numel(wEpochs.behav.PixelIdxList)-1
        interEpochIdx=spikeTimes>wEpochs.behav.PixelIdxList{epochNum}(end)/1000 & ...
            spikeTimes<wEpochs.behav.PixelIdxList{epochNum+1}(1)/1000;
        spikeTimes(interEpochIdx)=NaN;
    end
    spikeTimes(spikeTimes>wEpochs.behav.PixelIdxList{end}(end)/1000)=NaN;
    spikeTimes=spikeTimes(~isnan(spikeTimes));
else % keep first 30s
%     spikeTimes=spikeTimes(spikeTimes<=30);
end

%% Power spectrum
% set parameters
params.Fs=1000; % sampling frequency
params.fpass=[3 30]; % % band of frequencies to be kept
params.NW=min([100 floor(max(spikeTimes))*2.5]); %3
params.tapers=[params.NW params.NW*2-1]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[2 0.01];
params.trialave=0;
movingwin=[0.5 0.05];

% tic
[spS.spectrumVals,spS.freqVals,spS.R,spS.Serr]=mtspectrumpt(spikeTimes,params);
% toc
% [spS.spectrumVals,spS.times,spS.freqVals,spS.R,spS.Serr]=mtspecgrampt_optimized(spikeTimes,movingwin,params);

% convert to dB (power spectral density)
spS.spectrumValsPSD= 10*log10(spS.spectrumVals);
spS.SerrPSD(1,:)=10*log10(spS.Serr(1,:));
spS.SerrPSD(2,:)=10*log10(spS.Serr(2,:));
spS.RPSD=10*log10(spS.R);
 
spS.StatSigIdx=spS.SerrPSD(1,:)<=spS.RPSD & spS.SerrPSD(2,:)>=spS.RPSD;

%% figures
if false
    figure;
    plot(spS.freqVals,spS.spectrumVals)
    
    figure
    plot(spS.freqVals,10*log10(spS.spectrumVals),spS.freqVals,10*log10(spS.Serr(1,:)),...
        spS.freqVals,10*log10(spS.Serr(2,:)));
    line(get(gca,'xlim'),[10*log10(spS.R) 10*log10(spS.R)]);
    hold on
    plot(spS.freqVals(spS.StatSigIdx),spS.spectrumValsPSD(spS.StatSigIdx),'k')
    
end

