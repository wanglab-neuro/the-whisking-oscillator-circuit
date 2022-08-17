function [spikeLatency,ISI]=OptoJitter(spikeData,TTLtimes,keepCell,pulseDur,axisHandle)
%% plots spike latency standard deviation
% Input: alignedRasters - cell array with stim aligned rasters
%         alignTime - time used to align responses
%         pulseDur - pulse duration

for cellNum=1:length(keepCell)
    %find(mean(meanChan)==max(mean(meanChan)));
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761]); hold on
    end
    
    if size(spikeData.rasters,1)==1
        spikeTimes=(find(spikeData.rasters)'-0.5)/1000;
    else
        spikeTimes=spikeData.times(spikeData.unitID==keepCell(cellNum),:);
    end
    
    %get wich spike times occur during TTL
    spikeLatency=nan(numel(TTLtimes),1);
    pulseIdx=false(size(spikeTimes,1),size(TTLtimes,1));
    
    for TTLNum=1:length(TTLtimes)
        pulseIdx(:,TTLNum)=spikeTimes>TTLtimes(TTLNum) & spikeTimes<TTLtimes(TTLNum)+pulseDur;
        try
            spikeLatency(TTLNum)=spikeTimes(find(pulseIdx(:,TTLNum),1))-TTLtimes(TTLNum);
        catch
            spikeLatency(TTLNum)=NaN;
        end
    end
    spikeLatency=spikeLatency(~isnan(spikeLatency))*1000;
    onSpikes=any(pulseIdx,2);
    ISI=diff(spikeTimes(~onSpikes))*1000;
    
    if isgraphics(axisHandle) || ~isnan(axisHandle)
        plot(mean(spikeLatency),std(spikeLatency),'marker','o','MarkerFaceColor',...
            [0.3 0.75 0.93], 'MarkerEdgeColor','k','MarkerSize',8);
        plot(mean(ISI),std(ISI),'marker','o','MarkerFaceColor',...
            [1.0000    0.6784    0.0980], 'MarkerEdgeColor','k','MarkerSize',8); %[0.5 0.5 0.5]
        
        box off; grid('on');
        set(gca,'xlim',[1 ceil(max(get(gca,'xlim'))/10)*10],...
            'ylim',[0.1 ceil(max(get(gca,'ylim'))/10)*10]);
        set(gca,'xscale','log','yscale','log','GridAlpha',0.25,'MinorGridAlpha',1)
        set(gca,'XTick',[1,10, 100],'XTickLabel',[1,10, 100],...
            'YTick',[0.1 1 10],'YTickLabel',[0.1 1 10]);
        set(gca,'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
        %     legend('Latency','ISI','location','southeast','box','off')
        
        xlabel('Latency vs ISI (ms)')
        ylabel('Latency jitter, SD (ms)');
        hold off
    end
       
end
