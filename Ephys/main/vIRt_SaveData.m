function ephys=vIRt_SaveData(sessBaseName,baseDir,...
                whiskers,bWhisk,whiskingEpochs,whiskingEpochsList,wEpochMask,...
                breathing,ephys,sessData)
            
[ephys.selectedUnits,keep]=deal(ephys.keepList);
rmfield(ephys,'keepList')

%% organize selectedUnits by depth
ephys=vIRt_SortByDepth(ephys);
% traces=ephys.traces; % If traces are needed, get them from original location listed in recInfo.sessFiles
recInfo=ephys.recInfo;
recInfo.sessFiles=sessData.sessFiles;
ephys=rmfield(ephys,{'traces';'recInfo'});
pulses=sessData.pulses;

if ~exist(fullfile(baseDir,'Analysis','Data',sessBaseName),'dir')
    mkdir(fullfile(baseDir,'Analysis','Data',sessBaseName))
end

%% save session data: ephys, traces, whisking, breathing
% save(fullfile(baseDir,'Analysis','Data',sessBaseName,[sessBaseName '_recInfo.mat']),...
%     'traces');

save(fullfile(baseDir,'Analysis','Data',sessBaseName,[sessBaseName '_recInfo.mat']),...
    'recInfo');
save(fullfile(baseDir,'Analysis','Data',sessBaseName,[sessBaseName '_pulses.mat']),...
    'pulses');

save(fullfile(baseDir,'Analysis','Data',sessBaseName,[sessBaseName '_behavior.mat']),...
    'whiskers','bWhisk','whiskingEpochs','whiskingEpochsList','wEpochMask',...
    'breathing');

for unitNum=1:numel(keep)
    unitIdx=ephys.spikes.unitID==keep(unitNum);
    
    unitId=keep(unitNum);
    spikeTimes=ephys.spikes.times(unitIdx);
    waveForms=ephys.spikes.waveforms(unitIdx,:);
    preferredEl=ephys.spikes.preferredElectrode(unitIdx,:);
    bitResolution=ephys.spikes.bitResolution;
    samplingRate=ephys.spikes.samplingRate;
    rasters=ephys.rasters(keep(unitNum));
    spikeRate=ephys.spikeRate(keep(unitNum));
    coordinates=ephys.unitCoordinates(unitNum,:);
  
    save(fullfile(baseDir,'Analysis','Data',sessBaseName,[sessBaseName '_Unit'...
        num2str(unitId) '.mat']),...
        'unitId','spikeTimes','waveForms','preferredEl','rasters','spikeRate',...
        'bitResolution','samplingRate','coordinates')
end

ephys=rmfield(ephys,{'spikes'});
save(fullfile(baseDir,'Analysis','Data',sessBaseName,[sessBaseName '_ephys.mat']),...
    'ephys');

