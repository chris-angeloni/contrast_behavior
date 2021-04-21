function f = plot_single_cell_psych(cellID,res,spikeData,sessionData,ops)

% index this cell id in the res struct
rI = contains({res.single_cell.cellID},cellID);

if isempty(rI)
    error(sprintf(...
        'ERROR in plot_single_cell_psych.m: %s does not exist in res...\n',...
        cellID));
    f = [];
else

    clf;
    
    % extract res
    r = res.single_cell(rI);
    
    % index all units and sessions
    uI = contains(spikeData.cellinfo(:,7),cellID);
    sI = spikeData.cellinfo{uI,2} == [sessionData.session.sessionID];

    % make trial psth
    [~,spikes,trials,PSTHs] = makePSTH(spikeData.spikes{uI},...
                                       sessionData.events(sI).targOn,...
                                       ops.edges,5);
    tt = sessionData.behavior(sI).trialType;
    [sortt,sorti] = sortrows(tt);
    [~,spkSort] = ismember(trials,sorti);

    pp = sessionData.session(sI).plot;
    pp.colors1 = gen_colors(sessionData.session(sI).cond,.5,.8,0);

    % plot target raster
    nrows = 6;
    subplot(nrows,1,[1 2]); hold on;
    cv = pp.colors1(sortt(spkSort)+1,:);
    scatter(spikes,spkSort,15,cv,'.');
    ylim([1 length(tt)]); xlim([ops.edges(1) ops.edges(end)]);
    plot([0 0],ylim,'k','LineWidth',1.5);
    plotPrefs;
    ylabel('Trials');


    % plot mean psth
    subplot(nrows,1,3); hold on;
    uv = unique(tt(:,1));
    for i = 1:length(uv)
        plot(ops.time,mean(PSTHs(tt(:,1)==uv(i),:),1),...
             'color',pp.colors1(i,:),'LineWidth',1);
    end
    xlim([ops.edges(1) ops.edges(end)]);
    plot([0 0],ylim,'k','LineWidth',1.5);
    plotPrefs;
    ylabel('spks/s');

    % plot rt histogram
    subplot(nrows,1,4); hold on;
    [n,x] = hist(sessionData.behavior(sI).RT,ops.edges);
    stairs(x,n,'color',pp.contrastColor(end,:),'LineWidth',1);
    xlim([ops.edges(1) ops.edges(end)]);
    plot([0 0],ylim,'k','LineWidth',1.5);
    plotPrefs;
    ylabel('Lick Count');
    
    % plot auc
    subplot(nrows,1,[5 6]); hold on;
    plot([r.vols(2) r.vols(end)]+[-2 2],[.5 .5],'k');
    errorBars(r.vols(2:end),r.auc,pp.colors(end,:),[],[],r.auc_pct,[],...
                          'LineWidth',1,'Marker','.','MarkerSize',20);
    axis tight; plotPrefs;
    xlabel('Target Volume (dB SNR)');
    ylabel('AUC');
    
end


