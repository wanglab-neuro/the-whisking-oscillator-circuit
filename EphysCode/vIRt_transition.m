function [ISIs,spS]=vIRt_transition(wAngle,ephysData,dataMask,uIdx,savePlots)
% burstiness plots

if nargin < 5
    savePlots =false;
end
fileName=ephysData.recInfo.sessionName;

%% Data masking: look at each whisking epoch
wEpochs.behav=bwconncomp(dataMask.behav);
% mask epochs with short whisking bouts
durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=3000;

%% keep epochs with significant phase tuning (or just duration threshold if Slow)
epochID=find(durationThd);
durationThd(epochID(~dataMask.epochIdx))=false;

% % mask epochs with too much preceding whisking
% keepEpoch=false(numel(epochID),1);
% for epochNum=1:numel(epochID)
%     epochWindow=wEpochs.behav.PixelIdxList{epochID(epochNum)};
%     if epochWindow(1)<501; continue; end
%     if mean(abs(wAngle(epochWindow(1)-501:epochWindow(1)-1)))<1
%         keepEpoch(epochNum)=true;
%     end
% %     mean(abs(wAngle(epochWindow(1):epochWindow(1)+499)))
% end
% durationThd(epochID(~keepEpoch))=false;

dataMask.behav(vertcat(wEpochs.behav.PixelIdxList{~durationThd}))=false;
wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(durationThd);
wEpochs.behav.NumObjects=sum(durationThd);
% wEpochs.behav.PixelIdxList={vertcat(wEpochs.behav.PixelIdxList{:})};
% wEpochs.behav.NumObjects=1;

% % do the same for ephys data
wEpochs.ephys=bwconncomp(dataMask.ephys);
dataMask.ephys(vertcat(wEpochs.ephys.PixelIdxList{~durationThd}))=false;
wEpochs.ephys.PixelIdxList=wEpochs.ephys.PixelIdxList(durationThd);
wEpochs.ephys.NumObjects=sum(durationThd);
% wEpochs.ephys.PixelIdxList={vertcat(wEpochs.ephys.PixelIdxList{:})};
% wEpochs.ephys.NumObjects=1;

spikeRasters = ephysData.rasters(ephysData.selectedUnits,:);
numEpochs=wEpochs.behav.NumObjects;

for unitNum=1:size(spikeRasters,1)
    %     [mVals,wb_mVals]=deal(cell(numEpochs,1));
    
    unitSpikeEvent=spikeRasters(unitNum,:) ;%wEpochs.ephys.PixelIdxList{wEpochNum});
    %     ISI=diff(find(unitSpikeEvent))/1000;
    %       figure;  histogram(ISI,30)
    %         wb_ISI=ISI(ISI<=0.015); %within bursts
    
    spikeTimes=struct('interW',[],'intraW',[],'transition',[]);
    for wEpochNum=1:numEpochs
        epochIdx=wEpochs.ephys.PixelIdxList{wEpochNum};
        if epochIdx(1)<1000 || wEpochs.behav.PixelIdxList{wEpochNum}(1)<3000; continue; end
        
        %compare histogram inter and intra-epoch
        spikeTimes(wEpochNum).interW=find(unitSpikeEvent(epochIdx(1)-3000:epochIdx(1)-1))';
        spikeTimes(wEpochNum).intraW=find(unitSpikeEvent(epochIdx(1:3000)))';
        spikeTimes(wEpochNum).transition=find(unitSpikeEvent(epochIdx(1)-2000:epochIdx(3000)))'/1000;
        %         interISI=diff(spikeTimes(wEpochNum).interW);
        %         intraISI=diff(spikeTimes(wEpochNum).intraW);
        %             if numel(interISI)<20 || numel(intraISI)<20; continue; end
        %             figure; hold on
        %             histogram(interISI,2:15:150)
        %             histogram(intraISI,2:15:150)
    end
    wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(~cellfun(@isempty, {spikeTimes.interW}));
    wEpochs.ephys.PixelIdxList=wEpochs.ephys.PixelIdxList(~cellfun(@isempty, {spikeTimes.interW}));
    spikeTimes=spikeTimes(~cellfun(@isempty, {spikeTimes.interW}));
    [numEpochs,wEpochs.behav.NumObjects,wEpochs.ephys.NumObjects]=deal(size(spikeTimes,2));
    
    % isi with Chronux
    try
        [interISI_hist,interBins]=isi(rmfield(spikeTimes,'intraW'));
        [intraISI_hist,intraBins]=isi(rmfield(spikeTimes,'interW'));
        
        if savePlots
            figure('Color','white','name',...
                [fileName ' Unit' num2str(ephysData.selectedUnits(unitNum))]); hold on
            subplot(1,4,4); hold on
            bar(intraBins,intraISI_hist,2)
            bar(interBins,interISI_hist,2)
            legend('intra-whisking ISI','inter-whisking ISI')
        end
        
        ISIs.inter.vals=interISI_hist;
        ISIs.inter.bins=interBins;
        ISIs.intra.vals=intraISI_hist;
        ISIs.intra.bins=intraBins;
    catch
        ISIs=[];
    end
    
    for wEpochNum=1:numEpochs
        spikeTimes(wEpochNum).interW=spikeTimes(wEpochNum).interW/1000;
        spikeTimes(wEpochNum).intraW=spikeTimes(wEpochNum).intraW/1000;
    end
    
    %% compute continuous spectrum
    params.Fs=1000; % sampling frequency
    params.fpass=[1 25]; % % band of frequencies to be kept
    params.NW=2; %3
    params.tapers=[params.NW params.NW*2-1]; % taper parameters
    params.pad=2; % pad factor for fft
    params.err=[2 0.01];
    params.trialave=1;
    movingwin=[0.5 0.05];
    
    [spS.spectrumVals,spS.time,spS.freqVals]=...
        mtspecgrampt(rmfield(spikeTimes,{'intraW','interW'}),movingwin,params);
    
    if savePlots
        subplot(1,4,1:3); hold on
        imagesc(spS.time,spS.freqVals,10*log10(spS.spectrumVals)');
        %     imagesc(spS.time,spS.freqVals,spS.spectrumVals');
        xline(2,'LineWidth',2,'Color',[0 0 0 0.5])
        axis('tight');box off;
        xlabel('Time (ms)')
        ylabel('Frequency (Hz)');
        set(gca,'xlim',[0 5],'ylim',params.fpass,'xticklabels',-2:0.5:3,...
            'FontSize',10,'FontName','Helvetica','TickDir','out','Color','white');%'Calibri'
        cbH=colorbar;
        set(cbH,'ydir','normal','TickDir','out','box','off')
        ylabel(cbH, 'Power (dB)','fontsize',10);% if "raw" PSD (\muV^2/Hz)

        figDir='D:\Vincent\Figures\vIRt';
        savefig(gcf,fullfile(figDir, 'transition', ['Spectrogram - Cell' num2str(uIdx) fileName ' - Unit' num2str(unitNum) ' - all whisking epochs average']));
        print(gcf,fullfile(figDir, 'transition', ['Spectrogram - Cell' num2str(uIdx) fileName ' - Unit' num2str(unitNum) ' - all whisking epochs average']),'-dpng');
        close(gcf)
    end
end
end

