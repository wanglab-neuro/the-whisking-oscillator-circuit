function CV2=vIRt_CV2(wAngle,ephysData,dataMask)
% local coefficient of variation

%% Data masking: look at each whisking epoch
wEpochs.behav=bwconncomp(dataMask.behav);
% mask epochs with short whisking bouts
durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=3000;

%% keep epochs with significant phase tuning
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

% do the same for ephys data
wEpochs.ephys=bwconncomp(dataMask.ephys);
dataMask.ephys(vertcat(wEpochs.ephys.PixelIdxList{~durationThd}))=false;
wEpochs.ephys.PixelIdxList=wEpochs.ephys.PixelIdxList(durationThd);
wEpochs.ephys.NumObjects=sum(durationThd);
% wEpochs.ephys.PixelIdxList={vertcat(wEpochs.ephys.PixelIdxList{:})};
% wEpochs.ephys.NumObjects=1;

spikeRasters = ephysData.rasters(ephysData.selectedUnits,:);
numEpochs=wEpochs.behav.NumObjects;
CV2=struct('mVals',[],'wb_mVals',[],'CV2_all_ISI',[],'CV2_withinBurst_ISI',[]);

for unitNum=1:size(spikeRasters,1)
    
    [mVals,wb_mVals]=deal(cell(numEpochs,1));
    
    for wEpochNum=1:numEpochs
        unitSpikeEvent=spikeRasters(unitNum,wEpochs.ephys.PixelIdxList{wEpochNum});
        ISI=diff(find(unitSpikeEvent))/1000;
%       figure;  histogram(ISI,30)
        wb_ISI=ISI(ISI<=0.015); %within bursts
        
        % Formula: m=2*(|ISI−ISI'|/ISI+ISI')
        mVals{wEpochNum}=2*[0 abs(diff(ISI))]./movsum(ISI,2);
        wb_mVals{wEpochNum}=2*[0 abs(diff(wb_ISI))]./movsum(wb_ISI,2); %within bursts

        % check that there are enough m-values. Reliability of the estimates
        % requires a minimum of 20 m-values per window across all trials
        if numel(mVals{wEpochNum})<20
            [mVals{wEpochNum},wb_mVals{wEpochNum}]=deal(NaN);
        else
            mVals{wEpochNum}=mVals{wEpochNum}(2:end);    
            wb_mVals{wEpochNum}=wb_mVals{wEpochNum}(2:end);    
        end
              
    end
    CV2(unitNum).mVals=[mVals{:}];
    CV2(unitNum).wb_mVals=[wb_mVals{:}];
    CV2(unitNum).CV2_all_ISI=mean(CV2(unitNum).mVals);
    CV2(unitNum).CV2_withinBurst_ISI=mean(CV2(unitNum).wb_mVals);
end


