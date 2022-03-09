clear all; close all;
%addpath(genpath('~/chris-lab/code_general/'));
addpath(genpath('~/chris-lab/projects/contrast_glm/_functions/'));

%% setup
% load options
load('./_data/sim_tau001.mat');

% generate neural data
[stim,contrast] = stimgen(ops);
obj.tau = 10;
[y,l,h,g,x0] = simulate(obj,stim,contrast,ops);

% adjust silent periods so that they are at the default stimulus value
stim(stim==0) = obj.x0;

% make the stim predictors and fit the STRF
X = lagDesignMatrix(stim-obj.operating_point,length(ops.t))';

%[coeffs, dev, stats] = glmfit(X,y,'poisson');
%nt = length(ops.t); nf = length(ops.f);
%subplot(1,2,1);
%imagesc(reshape(res.strf(2:end),nf,[]));
%subplot(1,2,2);
%imagesc(reshape(coeffs(2:end),nf,[]));

% set up index for fitting the strf
ns = (ops.blockLength*2+ops.pad) * ops.fs;
t = ones(1,ns);
t(1:ops.gain_lags) = 0;
t(ops.blockLength*ops.fs:ops.blockLength*ops.fs+ops.gain_lags) = 0;
t(ops.blockLength*2*ops.fs:end) = 0;
ops.strf_include = logical(repmat(t,1,ops.ntrials*ops.exemplars));
ops.strf_include = logical(ones(size(y)));

% gain fitting index
t = ones(1,ns);
t(1:(ops.blockLength-.1)*ops.fs) = 0;
t((ops.blockLength*2-ops.pad)*ops.fs+1:end) = 0;
ops.gain_include = logical(repmat(t,1,ops.ntrials*ops.exemplars));

% transition index
t = ones(1,ns);
t(1:ops.blockLength*ops.fs) = 0;
t((ops.blockLength*2-ops.pad)*ops.fs+1:end) = 0;
ops.transition = logical(repmat(t,1,ops.ntrials*ops.exemplars));


%% fitting
% fit STRF in inclusion index
fprintf('fitting strf...'); tic;
[coeffs, dev, stats] = glmfit(X,y,'poisson'); toc;

fprintf('simulating taus...'); tic;
nt = ops.exemplars*ops.ntrials;
taus = [10 1 .5 .1 .05];
clear err y_pred yp
for i = 1:length(taus)
    [err(i),y_pred(i,:)] = gain_fun(taus(i),...
                                    y,...
                                    X,...
                                    contrast,...
                                    coeffs,...
                                    obj,...
                                    ops);
    yp(:,:,i) = reshape(y_pred(i,ops.gain_include),[],nt)';
end
toc;

% spike raster
yt = reshape(y(ops.gain_include),[],nt)';

% time
ind = logical(mean(reshape(ops.gain_include,[],nt)'));
time = (1:ns) ./ ops.fs - 1/ops.fs;

f1 = figure(1); clf;
subplot(2,3,1);
imagesc(obj.beta); title('STRF True');
axis square;
subplot(2,3,4);
imagesc(reshape(coeffs(2:end),length(ops.f),[]));
axis square;
title('STRF Estimate');
subplot(2,3,2);
imagesc(time(ind),1:nt,yt)
xlabel('Time (steps)'); ylabel('Trial');
title('Trial Activity');
subplot(2,3,5); hold on;
ph = plot(time(ind),squeeze(mean(yp,1)));
plot(time(ind),mean(yt),'k','LineWidth',1);
plot([2 2],ylim,'r--');
xlabel('Time (steps)');
ylabel('spks/s');
legend(ph,num2str(taus'))
title(sprintf('FR\ntrue tau = %g (black line)',obj.tau));
axis tight;
subplot(2,3,6); hold on
plot(taus,err,'k');
plot([obj.tau obj.tau],ylim,'r--');
plot(taus(err==min(err)),min(err),'ro')
legend('error','true value','best fit')
set(gca,'xscale','log');
xlabel('Tau Estimate'); ylabel('error');
title('Error'); 

addpath(genpath('~/chris-lab/code_general/plotting/'));
saveFigPDF(f1,[700 500],sprintf('./_plots/tau-estimator_true%04d.pdf',obj.tau*100))

