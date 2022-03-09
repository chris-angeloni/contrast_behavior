function [res] = run_strf(spikeData,sessionData,ops)

%% function [res] = run_strf(spikeData,sessionData,ops)
resFile = ops.resFile;
resFN = fullfile(ops.resDir,resFile);

fprintf('Checking stim files... '); tic;
% make a big cell array with the stimulus for each cell
incl = find(ops.include);
for i = 1:length(incl)
    % find the session for this cell
    sI = find(vertcat(sessionData.session.sessionID) == ...
              spikeData.cellinfo{incl(i),2});
    fp = strsplit(sessionData.session(sI).stimInfo.stim,'\');
    stimFile = fp{end};
    stim{i} = fullfile('./_stimuli',stimFile);
    stimExist(i) = exist(stim{i});
end
stimcells = false(size(ops.include));
stimcells(incl) = stimExist>0;
toc;

% of included neurons, find ones with the right stimulus and with
% higher fr
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


%% single neurons
clear sI;
for c = 1:length(spikes)

    %c = 122;

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
        
        % high contrast indec
        [~,hiC] = max(stimInfo.sd); 

        % index for matched stimuli
        I = stim.index(7,:) == 1;
        if sum(y0(I)) < 100
            fprintf('\tLess than 100 spikes, skipping... \n');
            continue;
        end

        % design matrix (not currently used)
        % S = stim.db - mean(stim.db,2);
        % S = S - mean(S(:));
        % X0 = lagDesignMatrix(S,length(ops.t));
        
        % spike triggered average for everything
        X = stim.db - mean(stim.db,2);
        STA = genSTA(find(y0>0 & I),X,ops.w,ops.fs);
        STA = STA - mean(STA,2); % subtract mean over time
        
        % compute STA for each contrast
        clear sta sta_c cfi cf lag strf_range cstrf cstrf_all
        for j = 1:2
            I = stim.index(2,:) == j & stim.index(7,:) == 1;
            included_spikes = find(y0>0 & I);
            
            % sta
            sta = genSTA(included_spikes,X,ops.w,ops.fs);
            sta_c(:,:,j) = sta - mean(sta,2); % subtract mean over time
            
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
            sf = sta_c(:,:,j);
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
            csfi = [cfi(hiC)-ops.bins : cfi(hiC)+ops.bins];
            csi = [1:ops.bins*2+1];
            csi = csi(csfi <= nnf & csfi > 0);
            csfi = csfi(csfi <= nnf & csfi > 0);
            
            % extract positive strf values around the cf
            cstrf(:,:,j) = nan(ops.bins*2+1,nt);
            cstrf(csi,:,j) = sf(csfi,:);
            
            % extract all strf values around cf
            stf = sta_c(:,:,j);
            cstrf_all(:,:,j) = nan(ops.bins*2+1,nt);
            cstrf_all(csi,:,j) = stf(csfi,:);
            
        end
        
        res(c).bf = cf;
        res(c).lag = lag;
        res(c).cstrf = cstrf;
        res(c).cstrf_all = cstrf_all;
        res(c).range = strf_range;
        res(c).strf_c = sta_c;
        res(c).strf = STA;
        
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
        I = behavior.trialType(:,1) == 0;
        for i = 1:2
            NR(:,i) = responsePower(trialPSTH(I,contrastI+1 == i), ...
                                    behavior.trialType(I,3));
        end
        res(c).noiseRatio = NR;

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
    
    fprintf('\tRuntime: '); toc(t0);
    
end

keyboard
