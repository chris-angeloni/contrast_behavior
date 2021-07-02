function [res,r] = run_choice(spikeData,sessionData,ops)

% function res = run_choice(spks,cellinfo,sessionData,ops)

sig_cells = ops.sig_neurons;
resFile = ops.resFile;

% get neurons in the task that have good waveforms
included_cells = ops.include & contains(spikeData.cellinfo(:,end),ops.task);
spikes = spikeData.spikes(included_cells);
cellInfo = spikeData.cellinfo(included_cells,:);

% reaction time analysis
sI = contains({sessionData.session.task},'psychometric');
s = sessionData.session(sI);
b = sessionData.behavior(sI);

for i = 1:numel(b)
    
    % for each session, get reaction times for each volume
    rtI = b(i).RT > .05 & b(i).RT < 1;
    rt = b(i).RT(rtI);
    tt = b(i).lickTT(rtI);
    snrI = [0 1 2 3 4 5 6];
    currentSNRS = ismember(snrI,tt);
    RT(i,:) = nan(size(currentSNRS));
    snr(i,:) = [-inf s(i).stimInfo.targetDBShift];
    RT(i,currentSNRS) = grpstats(rt,tt);
    snri(i,:) = snrI;
    contrast(i) = contains(s(i).cond,'lohi');
    
end
cc = repmat(contrast',1,7);
grpRT_mean = grpstats(RT(:),{snri(:),cc(:)},'mean');
grpRT_sem = grpstats(RT(:),{snri(:),cc(:)},'sem');
grps = grpstats([snri(:),cc(:)],{snri(:),cc(:)});
grpSNR = grpstats(snr(:),{snri(:),cc(:)});

color = {'b','r'};
hold on
for i = 1:2
    x = grpSNR(grps(:,2)==(i-1));
    ym = grpRT_mean(grps(:,2)==(i-1));
    ye = grpRT_sem(grps(:,2)==(i-1));
    errorbar(x(2) - mean(diff(x(2:end))),ym(1),ye(1),...
             color{i},'linewidth',1,'marker','o');
    errorbar(x(2:end),ym(2:end),ye(2:end),...
             color{i},'linewidth',1,'marker','o');
    
end
xlim([-7 27]); 
ylabel('RT (s)'); xlabel('Target Volume (mean dB SNR)');
plotPrefs;

keyboard

if ~exist(fullfile(ops.resDir,resFile),'file')
    
    % set up figure for internal use
    fh = figure(1234); clf;
    set(fh,'Visible','on');
    

    % first run analysis of single neurons
    t0 = tic;
    fprintf('CHOICE ANALYSIS: single neurons\n');
    for i = 1:length(spikes)
    end
end

