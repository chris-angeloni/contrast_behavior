clear all; close all;
%addpath(genpath('~/chris-lab/code_general/'));
addpath(genpath('./_functions/'));

%% setup
dataDir = './_data'; % location of data files

% load all spike and session data
fprintf('Loading spike data... '); tic;
spikeData = load(fullfile(dataDir,'spikeData.mat')); toc;
fprintf('Loading session data... '); tic;
sessionData = load(fullfile(dataDir,'sessionData.mat')); toc;

% set the seed
ops.seed = 1989;
rng(ops.seed);


%% waveform analysis
res_wave = run_wave(spikeData);
%f1 = figure(1); clf; plot_waveforms(res_wave);


%% cell inclusion criteria
%  1) Spike rate > 1Hz
%  2) good waveform shapes (not inverted, good pvalue, not super wide)
wave_include = ones(size(spikeData.cellinfo,1),1);
wave_include(~isnan(vertcat(spikeData.waveform.FWHM))) = res_wave.include;
include = [spikeData.cellinfo{:,8}]'>1 & wave_include;
ops.include = include;

%% Figure 1: run and plot the normative model simulations
% (this takes 20-30 minutes the first time, several minutes after
%  saving the initial simulatuion results)
run_normative_model;


%% Figure 3: run and plot behavior
% (requires behavior folder: ~/gits/gain-gonogo/')
stats_beh = run_behavior;


%% Figure 4: run and plot muscimol results
% (requires behavior folder: ~/gits/gain-gonogo/')
% (this will take roughly 20 minutes first run, then ~5 min rerunning)
stats_muscimol = run_muscimol;


%% Figure 5: psych
% (this will take quite a while the first time, several hours)

% options
ops.resDir = './_data';
ops.bin = .005;
ops.baseline = [-.1 0];
ops.target = [0 .1];
ops.response = [.1 1];
ops.edges = ops.target(1)+ops.baseline(1) : ...
    ops.bin : ...
    ops.target(2)+ops.response(2)-ops.baseline(1);
ops.time = ops.edges(1:end-1) + mean(diff(ops.edges));
ops.targetLevel = [5 6];
ops.noiseLevel = [0];
ops.smooth = 2;
ops.timeInd = 1;
ops.sig_neurons = true;
ops.resFile = '_res_psych.mat';

% run
[res_psych,r_psych] = run_psych(spikeData,sessionData,ops);

% plot summary
stats_psych = plot_psych_summaries(r_psych,ops);

% plot example neuron
f41 = figure(41); clf;
%res_psych.single_cell(find(contains({res_psych.single_cell.cellID},'CA118_2007070944'))).cellID
cellID = 'CA118_2007070944_09_040_mu';
plot_single_cell_psych(cellID,res_psych,spikeData,sessionData,ops);
saveFigPDF(f41,[300 475],'./_plots/_psych_ex_neuron_loc.pdf');

% plot example sessions
sz = [500 520];
f42 = figure(42); clf; set(f42,'Position',[0 0 sz]);
sessID = 2007070944; % ca118 - low contrast
plot_session_summary_psych_fig(sessID,res_psych,spikeData,sessionData,ops,...
                               r_psych,'critp_adj');
saveFigPDF(f42,sz,'./_plots/_psych_ex_session_loc.pdf');

f43 = figure(43); clf; set(f43,'Position',[0 0 sz]);
sessID = 2006110945; % ca118 - high contrast
plot_session_summary_psych_fig(sessID,res_psych,spikeData,sessionData,ops,...
                               r_psych,'critp_adj');
saveFigPDF(f43,sz,'./_plots/_psych_ex_session_hic.pdf');




%% Figure 5: offset

% options
ops.resDir = './_data';
ops.bin = .005;
ops.baseline = [-.1 0];
ops.target = [0 .1];
ops.response = [.1 1];
ops.edges = ops.target(1)+ops.baseline(1) : ...
    ops.bin : ...
    ops.target(2)+ops.response(2)-ops.baseline(1);
ops.time = ops.edges(1:end-1) + mean(diff(ops.edges));
ops.targetLevel = [2];
ops.noiseLevel = [0];
ops.smooth = 2;
ops.timeInd = 1;
ops.include = include;
ops.sig_neurons = false;
ops.resFile = '_res_offset.mat';

% run
[res_off,r_off] = run_offset(spikeData,sessionData,ops);

% plot summary
stats_off = plot_offset_summaries(r_off);

% example neuron
plot_single_cell_offset(...
    spikeData.cellinfo{26,7},res_off,spikeData,sessionData,ops);

keyboard



%% figure 6: LN model
% (fitting each cell takes some time... ~2-3 days for the initial run)

% model options
ops.bin = .025;
ops.fs = 1/ops.bin;
ops.w = 0.3;
ops.pad = ops.w;
ops.f = []; % frequencies depend on the data
ops.t = (-ops.w:1/ops.fs:0);
ops.nbins = 50;
ops.weight = true;
ops.sigma = .01;
ops.includeZero = false;
ops.modelType = 'exponential';
ops.kfolds = 10;
ops.models = {'sta','NRC','glm'};
ops.modelStrings = {'Spike-Triggered Average',...
                    'Normalized Reverse Correlation',...
                    'Regularized GLM'};
ops.targetExclude = .05;
ops.trialCutoff = 4.5 * ops.fs;
ops.resDir = './_data';
ops.resFile = '_res_lnmodel.mat';
ops.fig_visible = 'off';

% run
[res_ln] = run_lnmodel(spikeData,sessionData,ops);

% plots and stats
stats_ln = plot_lnmodel_summaries(res_ln,res_psych,r_psych,ops)


%% figure ?: STRF analysis

% options
ops.bins = 8;
ops.strf_its = 1000;
ops.resFile = '_res_strf_glm.mat';


% run
res_strf = run_strf(spikeData,sessionData,ops)

% find and remove empty cells
empty = arrayfun(  @(a)isempty(a.sessionIndex),  res_strf  );
res_strf(empty) = [];

% compute frequency and time sparseness
tI = ops.t > -0.1;
for i = 1:length(res_strf)
    
    % compute overall SNR for each STRF method
    SNR_glm(i) = std(res_strf(i).strf(:,tI),1,[1 2]) ./ ...
            std(res_strf(i).strf(:,~tI),1,[1 2]);
    SNR_sta(i) = std(res_strf(i).sta(:,tI),1,[1 2]) ./ ...
            std(res_strf(i).sta(:,~tI),1,[1 2]);
    
    clear strf fk tk;
    for j = 1:2
        % frequencies
        strf = res_strf(i).strf_c(:,tI,j);
        fk = mean(strf,2) .* std(strf,1,2);

        % time
        strf = res_strf(i).strf_c(:,:,j);
        tk = mean(strf,1) .* std(strf,1,1);

        % sparseness over frequencies and time (doesn't like
        % negative numbers)
        Sf(i,j) = sparseness(fk+min(fk));
        St(i,j) = sparseness(tk+min(tk));
        Sa(i,j) = sparseness(strf(:) + min(strf(:)));
        
        % compute SNR (signal being the std in the first 100ms
        % vs the std in the remaining field)
        SNR(i,j) = std(strf(:,tI),1,[1 2]) ./ ...
            std(strf(:,~tI),1,[1 2]);
        
% $$$         % compute SNR of shuffled STRF
% $$$         for ii = 1:1000
% $$$             S = reshape(strf(randperm(numel(strf))),size(strf));
% $$$             SNR_shuff(i,j,ii) = std(S(:,tI),1,[1 2]) ./ ...
% $$$                 std(S(:,~tI),1,[1 2]);
% $$$             
% $$$         end
    end
end

f3 = figure(33); clf;
% plot STRFs with different SNRs
% SNR2use = [1.501 max(SNR(:))];
SNR2use = [min(SNR(:)) max(SNR(:))];
nstrf = 25;
targetpoints = linspace(min(SNR2use),max(SNR2use),nstrf);

% find strfs close to each target point
for i = 1:numel(targetpoints)
    [~,strfI(i)] = min(mean(abs(targetpoints(i) - SNR),2));
end
strfI = unique(strfI,'stable');

% plot the high contrast strfs
nrows = ceil(sqrt(nstrf)); ncols = nrows;
for i = 1:numel(strfI)
    subplot(nrows,ncols,i);
    RF = res_strf(strfI(i)).strf_c(:,:,1);
    RF = (RF - mean(RF(:))) ./ std(RF(:));
    imagesc(RF,[-5 5])
    title(sprintf('%3.2f',SNR(strfI(i),2)));
    if i < numel(targetpoints)
        set(gca,'xticklabels',[],'yticklabels',[]);
    end
    axis square;
end



% look at strfs with high SNRs
I = find(all(SNR > 1.75,2));
for i = 1:length(I)
    
    fprintf('%g ',i);
    
    f = res_strf(1).ops.f / 1000;
    t = 1:13;
    
    clf; subplot(1,2,1);
    imagesc(t,f,res_strf(I(i)).strf_c(:,:,1)); colorbar;
    title(sprintf('%s\nSNR: %3.2f, BF: %3.2f',...
                  res_strf(I(i)).cellID,...
                  SNR(I(i),1),...
                  res_strf(I(i)).bf(1)/1000),...
          'interpreter','none');
    
    subplot(1,2,2);
    imagesc(t,f,res_strf(I(i)).strf_c(:,:,2)); colorbar;
    title(sprintf('SNR: %3.2f, BF: %3.2f',...
                  SNR(I(i),2),...
                  res_strf(I(i)).bf(2)/1000));
    
    pause;
    
end

f32 = figure(32); clf;
cells2use = [5 12 18 20 57];
for i = 1:length(cells2use)
    f = res_strf(1).ops.f / 1000;
   t = 0:-.025:-.3;
    
    subplot(length(cells2use),2,(i-1)*2 + 1);
    imagesc(t,f,fliplr(res_strf(I(cells2use(i))).strf_c(:,:,1))); ...
        colorbar;
    title(sprintf('%s',res_strf(I(cells2use(i))).cellID),...
          'interpreter','none')
    text(-.29,10,sprintf('SNR: %3.2f\nBF: %3.2f',...
                  SNR(I(cells2use(i)),1),...
                  res_strf(I(cells2use(i))).bf(1)/1000),...
         'color','w');
    axis square;
    plotPrefs;
    
    subplot(length(cells2use),2,(i-1)*2 + 2);
    imagesc(t,f,fliplr(res_strf(I(cells2use(i))).strf_c(:,:,2))); colorbar;
    text(-.29,10,sprintf('SNR: %3.2f\nBF: %3.2f',...
                  SNR(I(cells2use(i)),2),...
                  res_strf(I(cells2use(i))).bf(2)/1000),...
         'color','w');
    axis square;
    plotPrefs;
    if i == length(cells2use)
        xlabel('Time (s)'); ylabel('Frequency (kHz)');
    end
    
end
saveFigPDF(f32,[700 800],'./_plots/_strf_ex.pdf',.2)
    

f123 = figure(123); clf;

% select the data
SNRcut = mean(SNR(:)) + 2*std(SNR(:))
I = find(all(SNR > SNRcut,2));
BF = cat(1,res_strf(I).bf) / 1000;
LAG = cat(1,res_strf(I).lag);
RANGE = cat(1,res_strf(I).range);
sz = mean(SNR(I),2)

% plot params
sz = 10;
alpha = .01;
nrows = 4;
ncols = 1;

% selection from SNR dist
subplot(nrows,ncols,1); hold on;
histogram(SNR,50,'facecolor','k')
plot([SNRcut SNRcut],ylim,'r');
xlabel('STRF SNR'); ylabel('n cells');
plotPrefs; axis tight square;

% best frequency
subplot(nrows,ncols,2); hold on;
plot([min(BF(:)) max(BF(:))],[min(BF(:)) max(BF(:))],...
     'color',[.5 .5 .5]);
w = 1;
scatter(BF(:,1) + unifrnd(-w,w,length(LAG),1),...
        BF(:,2) + unifrnd(-w,w,length(LAG),1),...
        sz,'ok','markerfacecolor','k','markerfacealpha',alpha);
xlabel('Low Contrast BF (kHz)');
ylabel('High Contrast BF (kHz)');
plotPrefs; axis tight square;

[p,h,stats] = signrank(BF(:,1),BF(:,2))
title(sprintf('signrank: p=%03.2f, z = %3.2f',...
              p,stats.zval));

% lag
subplot(nrows,ncols,3); hold on;
w = 0.025 / 5;
plot([0 .1],[0 .1],'color',[.5 .5 .5]);
scatter(-LAG(:,1) + unifrnd(-w,w,length(LAG),1),...
        -LAG(:,2)+ unifrnd(-w,w,length(LAG),1),...
        sz,'ok','markerfacecolor','k','markerfacealpha',alpha);
xlabel('Low Contrast Lag (s)');
ylabel('High Contrast Lag (s)');
plotPrefs; axis tight square;

[p,h,stats] = signrank(-LAG(:,1),-LAG(:,2))
title(sprintf('signrank: p=%03.2f, z = %3.2f',...
              p,stats.zval));

% range
subplot(nrows,ncols,4); hold on;
plot([min(RANGE(:)) max(RANGE(:))],[min(RANGE(:)) max(RANGE(:))],'color',[.5 .5 .5]);
scatter(RANGE(:,1),RANGE(:,2),sz,'ok','markerfacecolor','k','markerfacealpha',alpha);
xlabel('Low Contrast Range (au)');
ylabel('High Contrast Range (au)');
plotPrefs; axis tight square;

[p,h,stats] = signrank(RANGE(:,1),RANGE(:,2))
title(sprintf('signrank: p=%03.2f, z = %3.2f',...
              p,stats.zval));
saveFigPDF(f123,[400 800],'./_plots/_strf_res.pdf',.2)




%% figure ?: normative model with "readout" noise
run_readout_noise(0:15);


%% figure ?: fit lohse model
ops = [];
run_lohse_model(ops);






    










































if 2+2 == 5

    %% Figure 7: choice

    % options
    ops.resDir = './_data';
    ops.bin = .005;
    ops.baseline = [-.1 0];
    ops.target = [0 .1];
    ops.response = [.1 1];
    ops.edges = ops.target(1)+ops.baseline(1) : ...
        ops.bin : ...
        ops.target(2)+ops.response(2)-ops.baseline(1);
    ops.time = ops.edges(1:end-1) + mean(diff(ops.edges));
    ops.timeStep = .05;
    ops.timeWindow = ops.timeStep*2;
    ops.timeCent = (-.1 + ops.timeWindow/2):ops.timeStep:...
        (max(1) - ops.timeWindow/2);
    ops.smooth = 2;
    ops.timeInd = 1;
    ops.include = include;
    ops.sig_neurons = false;
    ops.task = 'psychometric';

    % results filename
    ops.resFile = 'res_choice.mat';

    % run
    [res_choice,r_choice] = run_choice(spikeData,sessionData,ops);

    % plot summary
    stats_choice = plot_choice_summary(res_choice,res_ln,sessionData, ...
                                       ops);

end


%% noise adaptation
run_noise_adapt = false
if run_noise_adapt
    ops.bin = .005;
    ops.baseline = [-.1 0]; 
    ops.edges = ops.baseline(1) : ops.bin : 5-ops.baseline(1);
    ops.time = ops.edges(1:end-1) + mean(diff(ops.edges));
    ops.tsne_bins = 50;
    ops.tsne_smooth = 2.5;
    ops.kpcs = 7;
    ops.time_index = ops.time > 2 & ops.time < 4;
    ops.include = include;

    % run
    res_adapt = run_adapt(spikeData,sessionData,ops);

    % plot summary
    f2 = figure(2); clf; plot_noise_clusters(res_adapt,ops);
    saveFigPDF(f2,[800 1000],'./_plots/noise_FR_clustering.pdf')
end








