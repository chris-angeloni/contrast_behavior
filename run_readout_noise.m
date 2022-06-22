function run_readout_noise(nlevels)

if nargin < 1
    nlevels = 0:15;
end

% load simulation data
load('./_data/normative_model.mat')

% parameters
muTlow = 0:.25:3;
sigmaLow = 1;
sigmaHigh = 2;
f = .25;
nLevels = 16;
times = 1:50;
iT = 1:length(times);
bins = -.5:(nLevels-.5);
bb = [bins,nLevels+.5,nLevels+1.5]-1;
mu_lo = 1.5;
mu_hi = 2.25;
cInd = [100 50];

% histogram bins
hC = 0:nLevels-1;

fprintf('Computing target distributions with readout noise... '); tic;

% for each noise level, add noise to the buffers
for kk = 1:length(nlevels)
    
    % for each target volume
    for j = 1:length(rr)
        
        % each noise level randomly adds or subtracts a number
        % of spikes based on a uniform sample from +-k
        readout_noise_t = randi([-nlevels(kk) nlevels(kk)], ...
                                size(rr(j).buffer.yT));
        readout_noise_b = randi([-nlevels(kk) nlevels(kk)],...
                                size(rr(j).buffer.yT));

        
        % targets and background with readout noise
        yT = rr(j).buffer.yT + readout_noise_t;
        yB = rr(j).buffer.yB + readout_noise_b;
        
        % fix negative spike rates and rates above max
        yT(yT < 0) = 0; yT(yT > nLevels-1) = nLevels-1;
        yB(yB < 0) = 0; yB(yB > nLevels-1) = nLevels-1;

        % for each time index
        for i = 1:numel(iT)
            
            % for each contrast
            for cc = 1:2
                
                % extract responses to targets and background
                T = yT(times(i)+cInd(cc),:);
                B = yB(times(i)+cInd(cc),:);
                
                % generate a histogram
                cnT = histc(T,hC) + 1;
                cnB = histc(B,hC) + 1;
                cnT = cnT ./ sum(cnT);
                cnB = cnB ./ sum(cnB);
                
                % save histograms
                histT(j,i,cc,kk,:) = cnT;
                histB(j,i,cc,kk,:) = cnB;
                
                % KL divergence
                indsP = cnB > 0 & cnT > 0;
                klD(j,i,cc,kk) = 0.5 * (sum(cnB(indsP) .* log2(cnB(indsP) ./ cnT(indsP))) + sum(cnT(indsP) .* log2(cnT(indsP) ./ cnB(indsP))));
                
                % distribution overlap
                koverlap(j,i,cc,kk) = 1 - sqrt(cnB * cnT');
                
            end
        end
    end
end

% at the desired time index, fit psychometric functions for each
% noise level
T = 26; % time step near full adaptation
for i = 1:numel(nlevels)
    for j = 1:2
        x = cat(1,rr.muT_low);
        y = squeeze(koverlap(:,T,j,i));
        [prms(j,i,:),mdl] = fitLogGrid(x,y);
        mxslope(j,i) = max(diff(y)./diff(x));
    end
end


% set up color gradient for different noise levels
cgrad{1} = [linspace(0,.8,numel(nlevels))' ...
            linspace(0,.8,numel(nlevels))' ...
            ones(numel(nlevels),1)];
cgrad{2} = [ones(numel(nlevels),1) ...
            linspace(0,.8,numel(nlevels))' ...
            linspace(0,.8,numel(nlevels))'];

% plot the effect of noise on the response histograms
f42 = figure(42); clf;
vol = 10;
for i = 1:numel(nlevels)
    subplot(2,2,1); hold on;
    bar(hC,squeeze(histB(vol,25,1,i,:)),...
        'facecolor',cgrad{1}(i,:),...
        'linestyle','-',...
        'barwidth',1);
    legend(num2str(nlevels'),'location','northwestoutside')
    xlabel('k (spikes)'); ylabel('p(k)'); plotPrefs;
    ylim([0 .5]);

    subplot(2,2,2); hold on;
    bar(hC,squeeze(histT(vol,25,1,i,:)),...
        'facecolor',cgrad{1}(i,:),...
        'linestyle','-',...
        'barwidth',1);
    legend(num2str(nlevels'),'location','northwestoutside')
    xlabel('k (spikes)'); ylabel('p(k)'); plotPrefs;
    ylim([0 .5]);

    subplot(2,2,3); hold on;
    bar(hC,squeeze(histB(vol,25,2,i,:)),...
        'facecolor',cgrad{2}(i,:),...
        'linestyle','-',...
        'barwidth',1);
    legend(num2str(nlevels'),'location','northwestoutside')
    xlabel('k (spikes)'); ylabel('p(k)'); plotPrefs;
    ylim([0 .5]);

    subplot(2,2,4); hold on;
    bar(hC,squeeze(histT(vol,25,2,i,:)),...
        'facecolor',cgrad{2}(i,:),...
        'linestyle','-',...
        'barwidth',1);
    legend(num2str(nlevels'),'location','northwestoutside')
    xlabel('k (spikes)'); ylabel('p(k)'); plotPrefs;
    ylim([0 .5]);

end

saveFigPDF(f42,[800 500],'./_plots/_readoutNoiseHistograms.pdf',.2)


% for each level of readout noise, plot psychometric performance
f32 = figure(32); clf;
T = 26; % time step near full adaptation
vols = cat(1,rr.muT_low);
xf = linspace(min(vols),max(vols),100);
for i = 1:numel(nlevels)
    subplot(1,4,1); hold on;
    scatter(vols,squeeze(koverlap(:,T,1,i)),20,cgrad{1}(i,:),...
            'filled','markerfacealpha',.5);
    plot(xf,mdl(squeeze(prms(1,i,:)),xf),...
         'color',cgrad{1}(i,:),'linewidth',1);
    xlabel('Target Level'); ylabel('Discriminability');
    ylim([.55 .9]); plotPrefs;
    
    subplot(1,4,2); hold on;
    scatter(vols,squeeze(koverlap(:,T,2,i)),20,cgrad{2}(i,:),...
            'filled','markerfacealpha',.5);
    plot(xf,mdl(squeeze(prms(2,i,:)),xf),...
         'color',cgrad{2}(i,:),'linewidth',1);
    xlabel('Target Level'); ylabel('Discriminability');
    ylim([.55 .9]); plotPrefs;

end

subplot(1,4,3); hold on;
plot(nlevels,prms(1,:,1)./prms(1,:,2),'b');
plot(nlevels,prms(2,:,1)./prms(2,:,2),'r');
xlabel('Readout Noise'); ylabel('Threshold');
plotPrefs;

subplot(1,4,4); hold on;
plot(nlevels,mxslope(1,:),'b');
plot(nlevels,mxslope(2,:),'r');
xlabel('Readout Noise'); ylabel('Slope');
plotPrefs;

saveFigPDF(f32,[800 200],'./_plots/_readoutNoisePsych.pdf',.2)

% alternative plot of curves for each readout noise level
f31 = figure(31); clf
for i = 1:numel(nlevels)
    subplot(1,numel(nlevels),i); hold on;
    scatter(vols,squeeze(koverlap(:,T,1,i)),20,cgrad{1}(i,:),...
            'filled','markerfacealpha',.5);
    plot(xf,mdl(squeeze(prms(1,i,:)),xf),...
         'color',cgrad{1}(i,:),'linewidth',1);
    scatter(vols,squeeze(koverlap(:,T,2,i)),20,cgrad{2}(i,:),...
            'filled','markerfacealpha',.5);
    plot(xf,mdl(squeeze(prms(2,i,:)),xf),...
         'color',cgrad{2}(i,:),'linewidth',1);
    xlabel('Target Level'); ylabel('Discriminability');
    ylim([.55 .9]); plotPrefs;
    title(sprintf('readout noise = %d',nlevels(i)));
    
end

saveFigPDF(f31,[2400 200],'./_plots/_readoutNoisePsychPlots.pdf',.2)
