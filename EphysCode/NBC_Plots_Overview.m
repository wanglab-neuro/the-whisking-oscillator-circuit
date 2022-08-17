function NBC_Plots_Overview(whisker,whiskingEpochs,breathing, ephys,TTLtimes,zoomin,saveFig)
%% define figure colormap
cmap=lines;cmap=[cmap(1:7,:);(lines+flipud(copper))/2;autumn];

% find base setpoint and redress values if needed
baseSP=mode(round(whisker.setPoint/10)*10);
if baseSP<0
    whisker.setPoint=whisker.setPoint+180;
    whisker.angle=whisker.angle+180;
end

%% allocate
rasters=ephys.rasters(ephys.selectedUnits,:);
scaleF=ceil(max(whisker.angle));
if ~exist('TTLtimes','var') | isempty(TTLtimes)
    TTLtimes=[]; 
else
    pulseDur=0.01; 
%     pulseDur=mode(diff(TTLtimes)); 
end

%% open figure
figure('color','white','position',[62,167,1307,727],...
    'name',['Overview_' ephys.recInfo.baseName],'NumberTitle', 'off');
angleAxH = axes('Position',[0.1 0.5 0.8 0.4]); hold on;
spikesAxH = axes('Position',[0.1 0.1 0.8 0.4]); hold on;

%% plot whisker angle and setpoint
axes(angleAxH)
pH{1}=plot(whisker.timestamp,whisker.angle,'color',cmap(1,:),'linewidth',1.2);
if ~isempty(breathing)
    baseSP=mode(round(whisker.setPoint/10)*10);
    plot(whisker.timestamp,whisker.setPoint,'color',[cmap(2,:) 0.5],'linewidth',1.2);
    pH{2}=plot(breathing.ts,breathing.data+baseSP,'color',[cmap(3,:) 0.5],'linewidth',1.2);
else
    pH{2}=plot(whisker.timestamp,whisker.setPoint,'color',[cmap(2,:) 0.5],'linewidth',1.2);
end
for pulseNum=1:size(TTLtimes,1)
    patch([TTLtimes(pulseNum,1), TTLtimes(pulseNum,1),...
        TTLtimes(pulseNum,1)+pulseDur, TTLtimes(pulseNum,1)+pulseDur], ...
        [get(gca,'ylim') fliplr(get(gca,'ylim'))], ...
        [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.7);
end

% plot Phase
patchColor=[222, 207, 103]./255; %180, 225, 228
protractPeriods=bwconncomp(whisker.phase<=0 & whiskingEpochs);
proPIdx=cellfun(@(proP) [whisker.timestamp(proP(1)) whisker.timestamp(proP(1))...
    whisker.timestamp(proP(end)) whisker.timestamp(proP(end))],...
    protractPeriods.PixelIdxList,'un',0); proPIdx=[proPIdx{:}];
pH{3}=patch('Faces',reshape(1:protractPeriods.NumObjects*4,[4,protractPeriods.NumObjects])',...
    'Vertices',[proPIdx',...
    repmat([scaleF-10;scaleF;scaleF;scaleF-10],protractPeriods.NumObjects,1)],...
    'FaceVertexCData',repmat([0;0;6;6],protractPeriods.NumObjects,1),...
    'FaceColor',patchColor,'EdgeColor','none','FaceAlpha',0.5);
set(gca,'xtick',[],'xcolor','none','TickDir','out','box','off');
ylabel('whisker angle')
if ~isempty(breathing)
    legend([pH{:}],{'Whisker angle','Breathing air flow','Protraction phase'},'location','southeast')
else
    legend([pH{:}],{'Whisker angle','Whisker setpoint','Protraction phase'},'location','southeast')
end
legend('boxoff')

%% plot rasters
axes(spikesAxH)
EphysFun.PlotRaster(rasters,ephys.timestamps','lines',[],'k')

% plot Phase
scaleF=size(rasters,1);
patch('Faces',reshape(1:protractPeriods.NumObjects*4,[4,protractPeriods.NumObjects])',...
    'Vertices',[proPIdx',...
    repmat([0;scaleF;scaleF;0],protractPeriods.NumObjects,1)],...
    'FaceVertexCData',repmat([0;0;6;6],protractPeriods.NumObjects,1),...
    'FaceColor',patchColor,'EdgeColor','none','FaceAlpha',0.5);
for pulseNum=1:size(TTLtimes,2)
    patch(spikesAxH,[TTLtimes(1,pulseNum), TTLtimes(1,pulseNum),...
        TTLtimes(1,pulseNum)+pulseDur, TTLtimes(1,pulseNum)+pulseDur], ...
        [get(gca,'ylim') fliplr(get(gca,'ylim'))], ...
        [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.7);
end
set(gca,'xtick',0:size(rasters,2),... %:1000
    'xticklabel',(0:size(rasters,2)),...
    'ytick',0.5:size(rasters,1)-0.5,...
    'yticklabel',1:size(rasters,1),...
    'TickDir','out','box','off');
if isfield(ephys,'selectedUnits')
    set(gca,'yticklabel',ephys.selectedUnits);
end
xlabel('Time (s)')
ylabel('unit number')

% axes(spikesAxH)
% % yyaxis right
% axes('YAxisLocation','right')
% ylabel('X/Y Probe Geometry')
% set(gca,'ytick',0.5:size(rasters,1)-0.5,...
%         'yticklabel',ephys.unitCoordinates(:,2),...
%         'TickDir','out','box','off');

%link axes
linkaxes([angleAxH,spikesAxH],'x')

% zoom in whisking epochs
if zoomin
    wEpochs=bwconncomp(whiskingEpochs);
    viewStartPoint=round(wEpochs.PixelIdxList{1, 9}(1)/1000)*1000;
    set(gca,'xlim',[viewStartPoint viewStartPoint+2000])
end
%% save figure
if saveFig
    savefig(gcf,['Overview_Whole_' ephys.recInfo.sessionName '.fig'])
    %     saveas(gcf,['Overview_Whole_' ephys.recName '.png'])
end