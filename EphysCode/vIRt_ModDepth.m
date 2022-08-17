function modDepth=vIRt_ModDepth(wPhase,ephysData,dataMask)

% Modulation depth of the spike rates of vIRt neurons. 
% The depth is calculated from the spike rates as <(minimum - maximum)/average>,
% where the averaging is on a per cycle basis. 
% Data plotted as a function of the phase in the whisk cycle at which the spike rate is maximal.

%% Data masking: look at each whisking epoch
wEpochs.behav=bwconncomp(dataMask.behav);
% mask epochs with short whisking bouts
durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=3000;

%% keep epochs with significant phase tuning
epochID=find(durationThd);
durationThd(epochID(~dataMask.epochIdx))=false;

dataMask.behav(vertcat(wEpochs.behav.PixelIdxList{~durationThd}))=false;
wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(durationThd);
wEpochs.behav.NumObjects=sum(durationThd);
% wEpochs.behav.PixelIdxList={vertcat(wEpochs.behav.PixelIdxList{:})};
% wEpochs.behav.NumObjects=1;

% do the same for ephys data
wEpochs.ephys=bwconncomp(dataMask.ephys);
dataMask.ephys(vertcat(wEpochs.ephys.PixelIdxList{~durationThd}))=false;
wEpochs.ephys.PixelIdxList=wEpochs.ephys.PixelIdxList(durationThd);
wEpochs.ephys.NumObjects=sum(durationThd);
% wEpochs.ephys.PixelIdxList={vertcat(wEpochs.ephys.PixelIdxList{:})};
% wEpochs.ephys.NumObjects=1;

spikeRasters = ephysData.rasters(ephysData.selectedUnits,:);
numEpochs=wEpochs.behav.NumObjects;
modDepth=cell(size(spikeRasters,1),numEpochs);

for unitNum=1:size(spikeRasters,1)    
    for wEpochNum=1:numEpochs   
        unitSpikeEvent=spikeRasters(unitNum,wEpochs.ephys.PixelIdxList{wEpochNum});
        
        wEpochPhase=wPhase(wEpochs.behav.PixelIdxList{wEpochNum});
        [~,wPhaseCycles] = findpeaks(wEpochPhase);
        spikeRate=EphysFun.MakeSDF(unitSpikeEvent,20);
        cycleModDepth=nan(numel(wPhaseCycles)-1,1);
        for wCycle=1:numel(wPhaseCycles)-1
            cycleSpikeRate=spikeRate(wPhaseCycles(wCycle):wPhaseCycles(wCycle+1));
            cycleModDepth(wCycle)=abs((min(cycleSpikeRate)-max(cycleSpikeRate))/nanmean(cycleSpikeRate));
        end
        modDepth{unitNum,wEpochNum}=nanmean(cycleModDepth);
    end
end

modDepth=mean([modDepth{:}]);
