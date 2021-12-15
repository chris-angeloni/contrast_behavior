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










%% simulate inactivation
% run regularly
sigma_n = 1;
sigma_c = [1 2];
gc = [];
f = .25;
run_normative_model_inactivation(sigma_n,sigma_c,gc,f)

% run it without gain control
sigma_n = 1;
sigma_c = [1 2];
gc = 'off';
f = .25;
run_normative_model_inactivation(sigma_n,sigma_c,gc,f)

% run it with low gain in contrast
sigma_n = 1;
sigma_c = [1 2];
gc = 'low';
f = .25;
run_normative_model_inactivation(sigma_n,sigma_c,gc,f)

% run it with med noise in contrast
sigma_n = 2;
sigma_c = [1 2];
gc = [];
f = .25;
run_normative_model_inactivation(sigma_n,sigma_c,gc,f)

% run it with med noise in contrast
sigma_n = 3;
sigma_c = [1 2];
gc = [];
f = .25;
run_normative_model_inactivation(sigma_n,sigma_c,gc,f)

% run it with high noise in contrast
sigma_n = 5;
sigma_c = [1 2];
gc = [];
f = .25;
run_normative_model_inactivation(sigma_n,sigma_c,gc,f)


% simulate silent conditions by scaling the background
silent_scale = 0.1;

% run it with regular noise in silence
sigma_n = 1;
sigma_c = [1 1];   % background is only a product of neural noise
gc = [];
f = 2 * .25;       % high contrast targets (scaled as before)
run_normative_model_inactivation(sigma_n,sigma_c,gc,f,silent_scale)

% run it with med noise in silence
sigma_n = 2;
sigma_c = [1 1];  
gc = [];
f = 2 * .25;      
run_normative_model_inactivation(sigma_n,sigma_c,gc,f,silent_scale)

% run it with med noise in silence
sigma_n = 3;
sigma_c = [1 1];   % background is only a product of neural noise
gc = [];
f = 2 * .25;      
run_normative_model_inactivation(sigma_n,sigma_c,gc,f,silent_scale)

% run it with high noise in silence
sigma_n = 5;
sigma_c = [1 1];   % background is only a product of neural noise
gc = [];
f = 2 * .25;       
run_normative_model_inactivation(sigma_n,sigma_c,gc,f,silent_scale)


%  % run it with regular noise in silence
%  sigma_n = 1;       % neural noise = 1
%  sigma_c = [1 1];   % background is only a product of neural noise
%  gc = [];
%  f = 2 * .25;       % high contrast targets (scaled as before)
%  run_normative_model_inactivation(sigma_n,sigma_c,gc,f)
%  
%  % run it with med noise in silence
%  sigma_n = 2;       % neural noise = 5 to simulate muscimol
%  sigma_c = [1 1];   % background is only a product of neural noise
%  gc = [];
%  f = 2 * .25;       % high contrast targets (scaled as before)
%  run_normative_model_inactivation(sigma_n,sigma_c,gc,f)
%  
%  % run it with med noise in silence
%  sigma_n = 3;       % neural noise = 5 to simulate muscimol
%  sigma_c = [1 1];   % background is only a product of neural noise
%  gc = [];
%  f = 2 * .25;       % high contrast targets (scaled as before)
%  run_normative_model_inactivation(sigma_n,sigma_c,gc,f)
%  
%  % run it with high noise in silence
%  sigma_n = 5;       % neural noise = 5 to simulate muscimol
%  sigma_c = [1 1];   % background is only a product of neural noise
%  gc = [];
%  f = 2 * .25;       % high contrast targets (scaled as before)
%  run_normative_model_inactivation(sigma_n,sigma_c,gc,f)










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








