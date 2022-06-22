function [res] = run_strf(spikeData,sessionData,ops)

%% function [res] = run_strf(spikeData,sessionData,ops)
resFile = ops.resFile;
resFN = fullfile(ops.resDir,resFile);

if ~exist(resFN,'file')

    fprintf('Checking stim files... '); tic;
    % make a big cell array with the stimulus for each cell
    incl = find(ops.include);
    clear stim;
    for i = 1:length(incl)
        % find the session for this cell
        sI = find(vertcat(sessionData.session.sessionID) == ...
                  spikeData.cellinfo{incl(i),2});
        fp = strsplit(sessionData.session(sI).stimInfo.stim,'\');
        
        % check if the raw stimulus file exists
        stimFile = fp{end};
        stim{i} = fullfile('./_stimuli',stimFile);
        stimExist(i) = exist(stim{i});
        
        % check if the spectrogram file exists
        specFile{i} = sprintf('./_data/_spectrograms/%s_%d-spec.mat',...
                              sessionData.session(sI).mouse,...
                              sessionData.session(sI).sessionID);
        specExist(i) = exist(specFile{i});
        
    end
    stimcells = false(size(ops.include));
    stimcells(incl) = stimExist>0 | specExist > 0;
    toc;
    
    % load all the specfiles to check they are not corrupt
    fprintf('Checking spectrogram files...'); tic;
    specFILES = unique(specFile(specExist>0));
    for i = 1:length(specFILES)
        try
            a = load(specFILES{i});
        catch
            keyboard
        end
    end
    toc;

    % of included neurons, find ones with the right stimulus
    included_cells = stimcells;
    spikes = spikeData.spikes(included_cells);
    cellInfo = spikeData.cellinfo(included_cells,:);

    % session for the first neuron
    last_sI = 0;

    % models to run
    models = ops.models;

    % figure handling
    sz = [1500 400];
    f1 = figure(1); set(f1,'visible',ops.fig_visible,'Position',[0 0 sz]);
    cols = {'b','r'};

    %% single neurons
    clear sI;
    t1 = 0;
    START = 1 % 2074; %% CHANGE ME BACK TO 1
    for c = START:length(spikes)

        % c = 122;
        % c = 8

        % cell info
        fprintf('CELL %d/%d... \n',c,length(spikes));
        fn = sprintf('./_data/_lnmodel/%s.mat',cellInfo{c,7});
        sI(c) = find(vertcat(sessionData.session.sessionID) == ...
                     cellInfo{c,2});

        t0 = tic;

        if ~exist(fn,'file')


            % format spikes and stimuli (only if cell is in a new session)
            stimInfo = sessionData.session(sI(c)).stimInfo;
            behavior = sessionData.behavior(sI(c));
            ops.f = stimInfo.freqs;
            stimInfo.offsets = sessionData.session(sI(c)).offsets;
            events = sessionData.events(sI(c)).trialOn;
            [stim,spec,y0,ops] = formatStimAndSpikes(...
                spikes{c},events,stimInfo,behavior,ops,last_sI==sI(c));
            
            % high contrast index
            [~,hiC] = max(stimInfo.sd); 

            % index for matched stimuli
            I = stim.index(7,:) == 1;
            if sum(y0(I)) < 100
                fprintf('\tLess than 100 spikes, skipping... \n');
                continue;
            end

            % design matrix
            X0 = normalize(lagDesignMatrix(stim.db(:,I),length(ops.t)));
            
            % glmfit
            % b = glmfit(X0',y0(I)','poisson')
            % strf = fliplr(reshape(b(2:end),length(ops.f),[]));
            
            % cross-validated fit with glmnet
            options.alpha = 0;
            nfolds = 10;
            family = 'gaussian';
            cvfit = cvglmnet(X0',y0(I)',family,options,[],nfolds);
            coeffs = cvglmnetCoef(cvfit,'lambda_1se');
            STRF = fliplr(reshape(coeffs(2:end),length(ops.f),[]));
            
            % spike triggered average
            X = stim.db - mean(stim.db,2);
            STA = genSTA(find(y0>0 & I),X,ops.w,ops.fs);
            STA = STA - mean(STA,2); % subtract mean over time
            
            % compute strf for each contrast
            clear sta sta_c strf_c cfi cf lag strf_range cstrf cstrf_all
            for j = 1:2
                
                if j == 1      % low contrast first
                    [~,cI] = min(stimInfo.sd);
                elseif j == 2  % high contrast second
                    [~,cI] = max(stimInfo.sd);
                end
                
                I = stim.index(2,:) == cI & stim.index(7,:) == 1;
                included_spikes = find(y0>0 & I);
                
                
                fprintf('\tfitting strf... '); tic;
                % sta (compute cause its quick)
                sta = genSTA(included_spikes,X,ops.w,ops.fs);
                sta_c(:,:,j) = sta - mean(sta,2); % subtract mean over time
                
                % compute cross-validated STRF
                X0 = normalize(lagDesignMatrix(stim.db(:,I),length(ops.t)));
                cvfit = cvglmnet(X0',y0(I)',family,options,[],nfolds);
                coeffs = cvglmnetCoef(cvfit,'lambda_1se');
                strf_c(:,:,j) = fliplr(reshape(coeffs(2:end),length(ops.f),[]));
                toc;
                
                if (2+2) == 5
                    for i = 1:ops.strf_its
                        % calculate STA from each half of data and correlate
                        I1 = randsample(included_spikes,round(length(included_spikes)/2));
                        I2 = included_spikes(~ismember(included_spikes,I1));
                        sta1 = genSTA([I1 y0(I1)],stim.db,ops.w,ops.fs);
                        sta2 = genSTA([I2 y0(I2)],stim.db,ops.w,ops.fs);
                        r(i,j) = corr(sta1(:),sta2(:));
                    end
                end
                
                % centered STRF
                sf = strf_c(:,:,j);
                nnf = length(ops.f);
                nt = length(ops.t);
                sf(sf<0) = nan;
                
                % compute time and frequency means, weighted by the standard
                % deviation in each bin (this will place more weight onto
                % "active" channels)
                pkt = mean(sf,1,'omitnan') .* std(sf,1,'omitnan');
                pkf = mean(sf,2,'omitnan') .* std(sf,[],2,'omitnan');
                if all(pkt == 0)
                    pkt(2) = 1;
                end
                if all(pkf == 0)
                    pkf(ops.bins+1) = 1;
                end
                [~,cfi(j)] = max(pkf); 
                cf(j) = ops.f(cfi(j));
                lag(j) = ops.t(pkt==max(pkt));
                strf_range(j) = range(sf(:));
                
            end
            for j = 1:2
                
                % index around center frequency in high contrast
                csfi = [cfi(2)-ops.bins : cfi(2)+ops.bins];
                csi = [1:ops.bins*2+1];
                csi = csi(csfi <= nnf & csfi > 0);
                csfi = csfi(csfi <= nnf & csfi > 0);
                
                % extract positive strf values around the cf
                cstrf(:,:,j) = nan(ops.bins*2+1,nt);
                cstrf(csi,:,j) = sf(csfi,:);
                
                % extract all strf values around cf
                stf = strf_c(:,:,j);
                cstrf_all(:,:,j) = nan(ops.bins*2+1,nt);
                cstrf_all(csi,:,j) = stf(csfi,:);
                
            end
            
            res(c).bf = cf;
            res(c).lag = lag;
            res(c).cstrf = cstrf;
            res(c).cstrf_all = cstrf_all;
            res(c).range = strf_range;
            res(c).sta_c = sta_c;
            res(c).sta = STA;
            res(c).strf_c = strf_c;
            res(c).strf = STRF;
            
            %% noise power
            % make a raster
            trialEdges = -.1:.005:5.1;
            trialTime = trialEdges(1:end-1) + .005/2;
            [noiseSortT,noiseSortI] = ...
                sortrows([behavior.trialType(:,3) behavior.trialType(:,1)]);
            [trialPSTH,trialRaster,trialTrials,trialPSTHsm] = ...
                makePSTH(spikes{c},events,trialEdges,5);
            [~,noiseSortTrials] = ismember(trialTrials,noiseSortI);
            
            % find max trial length to use for average PSTH
            nMax = max(spec.trial_actual_length + 2*spec.pad_samps);

            % compute noise power, excluding transition periods and
            % excluding target trials
            timeI = -ops.pad : 1/ops.fs : (nMax/ops.fs-ops.pad);
            timeI = timeI(1:end-1) + 1/ops.fs/2;
            contrastI = timeI > 3;
            switchI = timeI < 1.5 | ...
                      timeI > 3 & timeI < ops.pad + 3;
            I = behavior.trialType(:,1) == 0;
            figure(123); clf;
            for i = 1:2
                
                if i == 1 % low contrast
                    [~,cI] = min(stimInfo.sd);
                elseif i == 2 % high contrast
                    [~,cI] = max(stimInfo.sd);
                end
                tI = contrastI+1 == cI & ~switchI;
                NR(:,i) = responsePower(trialPSTH(I,tI),behavior.trialType(I,3));
                
                fI = find(I);
                [~,sorti] = sort(behavior.trialType(I,3));
                subplot(4,2,[1 2]); hold on
                imagesc(timeI(tI),1:length(fI),trialPSTH(fI(sorti),tI),...
                        [0 1]);
                plot(timeI,tI*(length(fI)+2),cols{i},'linewidth',3);
                ylabel('Time');
                
            end
            axis tight;
            res(c).noiseRatio = NR;
            NR(NR < 0 | isinf(NR)) = nan;
            res(c).NR = median(NR,1,'omitnan');

            % add some other stuff
            res(c).ops = ops;
            res(c).maxTrialSamps = nMax;
            res(c).cellInfo = cellInfo(c,:);
            res(c).cellID = cellInfo{c,7};
            res(c).sessionID = sessionData.session(sI(c)).sessionID;
            res(c).cond = sessionData.session(sI(c)).cond;
            res(c).sessionIndex = sI(c);
            res(c).spec_file = spec.spec_file;
            res(c).timeIndex = timeI;
            res(c).contrastIndex = contrastI;
            
        end
        
        t1(c) = toc(t0);
        fprintf('\tRuntime: %g\n',t1(c)); 
        
        if c == START
            sta_1 = zeros(size(res(c).cstrf_all(:,:,1)));
            sta_2 = zeros(size(res(c).cstrf_all(:,:,1)));
        end
        sta_1 = sum(cat(3,sta_1,res(c).cstrf_all(:,:,1)),3,'omitnan');
        sta_2 = sum(cat(3,sta_2,res(c).cstrf_all(:,:,2)),3,'omitnan');
        
        subplot(4,2,3);
        imagesc(res(c).strf_c(:,:,1)); colorbar; axis square;
        title(res(c).NR(1));
        subplot(4,2,4);
        imagesc(res(c).strf_c(:,:,2)); colorbar; axis square;
        title(res(c).NR(2));
        subplot(4,2,5);
        imagesc(sta_1./c); colorbar; axis square;
        subplot(4,2,6);
        imagesc(sta_2./c); colorbar; axis square;
        subplot(4,2,[7 8]);
        plot(1:c,t1); xlim([0 length(spikes)]);
        xlabel('cell'); ylabel('runtime (s)');
        title(sprintf('%d: %s (%s)',c,res(c).cellID,res(c).cond),'interpreter','none');
        drawnow;
        

    end

    save(resFN,'res');

else
    
    load(resFN);
    
end
