function [latency,jitter]=vIRt_PhotoTagPlots(ephysData,pulses,uIdx,savePlots)

if nargin < 4
    savePlots =false;
end

%% variables
fileName=ephysData.recInfo.sessionName;
TTLs.start=pulses.TTLTimes;
pulseDur=pulses.duration;
IPI=mode(diff(TTLs.start));
delay=0.005;
preAlignWindow=0.050;
postAlignWindow=0.20;
SRR=ephysData.recInfo.SRratio;
traceExcerpt.excerptSize=SRR;

if islogical(ephysData.selectedUnits) %logical array
    ephysData.selectedUnits=find(ephysData.selectedUnits);
end

%% compute rasters
if isfield(ephysData,'rasters')
    spikeRasters=ephysData.rasters;
else
    spikeRasters=EphysFun.MakeRasters(ephysData.spikes.times,ephysData.spikes.unitID,...
        1,int32(size(ephysData.traces,2)/ephysData.spikes.samplingRate*1000)); %ephysData.spikes.samplingRate
end
if size(spikeRasters,1)>1
    spikeRasters=spikeRasters(ephysData.selectedUnits,:);
end
alignedRasters=EphysFun.AlignRasters(spikeRasters,TTLs.start,preAlignWindow,postAlignWindow,1000);

%% compute spike density functions
% spikeRate=EphysFun.MakeSDF(spikeRasters);

%% Figures
%if need to load ephys data:
spikeSortingDir=[ephysData.recInfo.dirName filesep 'SpikeSorting' filesep ephysData.recInfo.sessionName];
LoadSpikeData(fullfile(spikeSortingDir, [ephysData.recInfo.sessionName '_export_res.mat'])) ;

for cellNum=1:size(ephysData.selectedUnits,1)

    figure('Position',[214   108   747   754],'Color','w','name',...
        [fileName ' Unit' num2str(ephysData.selectedUnits(cellNum))] ); %Ch' num2str(spikeData.selectedUnits(cellNum))

    %% raw trace
    subplot(3,3,7:9); hold on
    if ~isfield(ephysData.recInfo,'SRratio')
        SRR=double(ephysData.spikes.samplingRate/1000);
    end
    traceExcerpt.location=TTLs.start(1)*SRR;
    if exist('traceData','var') && isa(traceData,'memmapfile')
        winIdxStart=(traceExcerpt.location-traceExcerpt.excerptSize)*traceData.traceInfo.numChan+1;
        winSize=2; %default 1 pulse
        winIdxEnd=winIdxStart+(winSize*2*traceExcerpt.excerptSize*traceData.traceInfo.numChan);
    else
        winIdxStart=(traceExcerpt.location-traceExcerpt.excerptSize); %*traceData.traceInfo.numChan+1;
        winIdxEnd=traceExcerpt.location+traceExcerpt.excerptSize;
    end
    excerptWindow=int32(winIdxStart:winIdxEnd-1);%-SRR;
    if exist('traceData','var') && isa(traceData,'memmapfile')
        traceExcerpt.data=traceData.allTraces.Data(excerptWindow);
        traceExcerpt.data=reshape(traceExcerpt.data,[traceData.traceInfo.numChan traceExcerpt.excerptSize*2*winSize]);
        preprocOption={'CAR','all'};
        traceExcerpt.data=PreProcData(traceExcerpt.data,30000,preprocOption);
        traceExcerpt.data=traceExcerpt.data(channelNum,:);%     figure; plot(dataExcerpt(11,:))
    else
        prefElec=double(ephysData.spikes.preferredElectrode(ismember(...
            ephysData.spikes.unitID,ephysData.selectedUnits(cellNum))));
        try
            keepTrace=mode(prefElec(ephysData.spikes.times>(traceExcerpt.location-...
                traceExcerpt.excerptSize)/SRR));
        catch
            keepTrace=mode(prefElec);
        end
        if excerptWindow(end) > size(ephysData.traces,2)
            excerptWindow = excerptWindow(1):size(ephysData.traces,2);
        end
        traceExcerpt.data=ephysData.traces(keepTrace,excerptWindow);
    end

    excerptTTLtimes=double(TTLs.start(TTLs.start>(traceExcerpt.location-...
        traceExcerpt.excerptSize)/SRR &...
        TTLs.start<(traceExcerpt.location+traceExcerpt.excerptSize)/SRR)-...
        (traceExcerpt.location-traceExcerpt.excerptSize)/...
        SRR)*SRR;

    try
        excerptSpikeTimes={double(ephysData.spikes.times(ephysData.spikes.times>(traceExcerpt.location-...
            traceExcerpt.excerptSize)/SRR &...
            ephysData.spikes.times<(traceExcerpt.location+traceExcerpt.excerptSize)/SRR)-...
            (traceExcerpt.location-traceExcerpt.excerptSize)/...
            SRR)*SRR};
    catch
        excerptSpikeTimes={NaN};
    end
    %     figure; plot(traceExcerpt.data)
    OptoRawTrace(traceExcerpt,excerptSpikeTimes,...
        SRR,excerptTTLtimes,pulseDur,'',gca)

    %% waveforms
    if isfield(ephysData.spikes,'wF')
        waveForms=ephysData.spikes.wF(ephysData.selectedUnits(cellNum)).spikesFilt;
        keepTrace=range(mean(waveForms,3))==max(range(mean(waveForms,3)));
        waveForms=squeeze(waveForms(:,keepTrace,:))'*ephysData.recInfo.bitResolution; %ephysData.recInfo.channelMap==keepTrace
        ephysData.spikes.waveforms=ephysData.spikes.waveforms(:,1:size(waveForms,2));
        ephysData.spikes.waveforms(ephysData.spikes.unitID==ephysData.selectedUnits(cellNum),:)=waveForms;
    elseif isfield(ephysData.spikes,'waveform') && ~isempty(ephysData.spikes.waveform)
        % all good
    else
        spikesTimes=ephysData.spikes.times(ephysData.spikes.unitID==ephysData.selectedUnits(cellNum));
        waveForms=NaN(size(spikesTimes,1),50);
        %         electrodesId=unique(spikes.preferredElectrode);
        waveForms=ExtractChunks(ephysData.traces(keepTrace,:),... %foo = PreProcData(foo,30000,{'bandpass',[300 3000]});
            spikesTimes*ephysData.recInfo.samplingRate,50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
        % scale to resolution
        waveForms=waveForms.*ephysData.recInfo.bitResolution;
        ephysData.spikes.waveforms(ephysData.spikes.unitID==ephysData.selectedUnits(cellNum),:)=waveForms;
    end
    subplot(3,3,5); hold on
    onSpikes=OptoWaveforms(ephysData.spikes,TTLs.start,ephysData.selectedUnits(cellNum),delay,gca);

    %% Jitter
    subplot(3,3,2); hold on
    latency=OptoJitter(ephysData.spikes,TTLs.start,ephysData.selectedUnits(cellNum),pulseDur,gca);

    % %% ISI
    subplot(3,3,3); hold on
    OptoISI(ephysData.spikes,TTLs.start,ephysData.selectedUnits(cellNum),pulseDur,gca)
    %
    % %% ACG
    subplot(3,3,6); hold on
    OptoACG(ephysData.spikes,TTLs.start,ephysData.selectedUnits(cellNum),pulseDur,gca)

    %% rasters
    subplot(3,3,4);
    if ~iscell(alignedRasters); alignedRasters={alignedRasters}; end
    OptoRasters(alignedRasters(cellNum),preAlignWindow*1000,pulseDur,IPI,gca);
    % title(['Channel ' num2str(channelNum) ', Neuron ' num2str(spikeData.selectedUnits(cellNum))],'FontName','Cambria');

    %% SDF
    subplot(3,3,1)
    OptoSDF(alignedRasters(cellNum),preAlignWindow*1000,pulseDur*1000,IPI*1000,gca)

    if savePlots
        %         if ~exist(fullfile(cd, 'Figures'),'dir')
        %             mkdir('Figures')
        %         end
        figDir='D:\Vincent\Figures\vIRt';
        savefig(gcf,fullfile(figDir, 'PT', ['Cell' num2str(uIdx) fileName '_Unit' num2str(cellNum) '_PT.fig']));
        print(gcf,fullfile(figDir, 'PT', ['Cell' num2str(uIdx) fileName '_Unit' num2str(cellNum) '_PT']),'-dpng');
        close(gcf)
    end
end
end

