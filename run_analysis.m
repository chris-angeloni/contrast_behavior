clear all; close all;
addpath(genpath('~/chris-lab/code_general/'));
addpath(genpath('./_functions/'));

%% setup
% load all spike and session data
fprintf('Loading spike data... '); tic;
spikeData = load('~/data/gain_behavior/spikeData.mat'); toc;
fprintf('Loading session data... '); tic;
sessionData = load('~/data/gain_behavior/sessionData.mat'); toc;


%% waveform analysis
res_wave = run_wave(spikeData);
f1 = figure(1); clf; plot_waveforms(res_wave);



%% inclusion criteria
%  1) Spike rate > 1Hz
%  2) good waveform shapes (not inverted, good pvalue, not super wide)
wave_include = ones(size(spikeData.cellinfo,1),1);
wave_include(~isnan(vertcat(spikeData.waveform.FWHM))) = res_wave.include;
include = [spikeData.cellinfo{:,8}]'>1 & wave_include; 



%% noise adaptation
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




%% psych

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
ops.include = include;
ops.sig_neurons = true;

% results filename
if ops.sig_neurons
    ops.resFile = 'res_psych_sigcells.mat';
else
    ops.resFile = 'res_psych_allcells.mat';
end

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

keyboard

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



%% offset

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

% results filename
if ops.sig_neurons
    ops.resFile = 'res_offset_sigcells.mat';
else
    ops.resFile = 'res_offset_allcells.mat';
end

% run
[res_off,r_off] = run_offset(spikeData,sessionData,ops);

% plot summary
stats_off = plot_offset_summaries(r_off);

% example neuron
plot_single_cell_offset(...
    spikeData.cellinfo{26,7},res_off,spikeData,sessionData,ops);






