function run_lohse_model(ops)

addpath(genpath('./_functions/'));
addpath(genpath('../contrast_glm/'));
addpath(genpath('~/chris-lab/code_general/'));

resPath = '~/data/gain_opto/glm_res';
dc = [0 1 2];
models = {'ln-static','ln-gc','glm','glm-symm','glm-asymm(1)', ...
          'glm-asymm(2)'};
ops.mdls = models;
%tag = '';
tag = 'spline-3-3_';
%tag = 'spline-5-3_';
%tag = 'spline-7-3_';
%tag = 'spline-10-3_';
%tag = 'spline-15-3_';
resfile = sprintf('./acute_res_%s2tau_lohse.mat',tag);

% setup figure
plot_on = true;
plot_visible = false;
sz = [900 900];
f1 = figure(1); set(f1,'visible',plot_visible,'Position',[0 0 sz]);

if ~exist(resfile,'file')
    % get file list for one model (all neurons)
    fileList = dir(fullfile(resPath,sprintf('acute_%sdc%d*.mat',tag,dc(1))));
    
    for i = 1:length(fileList)
        
        fprintf('%d/%d... ',i,length(fileList)); tic;
        
        for m = 1:length(dc)
            
            % get file list for this model
            fileList = dir(fullfile(resPath,sprintf('acute_%sdc%d*.mat',tag,dc(m))));

            load(fullfile(fileList(i).folder,fileList(i).name));
            
            % common variabiles
            if m == 1
                
                % centered strf with lag/bf
                ops.bins = 7;
                [cstrf,bf,lag] = centeredSTRF(res.strf,ops);
                
                % activity/strf
                r(i).cellID = res.cellID;
                r(i).sessionID = res.sessionID;
                r(i).y = res.y;
                r(i).c = res.c;
                r(i).scene = res.scene;
                r(i).strf = res.strf;
                r(i).cstrf = cstrf;
                r(i).bf = bf;
                r(i).lag = lag;
                r(i).nr = median(res.nr,1,'omitnan');
                
                % ln results
                r(i).static_ln_gain = res.static_ln.ahat(3);
                r(i).static_ln_mxslope = res.static_ln.maxslope;
                r(i).static_ln_meanslope = res.static_ln.maxslope;
                r(i).static_ln_beta = res.static_ln.beta;
                r(i).pred(1,:) = res.static_ln.pred;
                r(i).gc_ln_gain = res.gc_ln.ahat(:,3);
                r(i).gc_ln_mxslope = res.gc_ln.maxslope;
                r(i).gc_ln_meanslope = res.gc_ln.maxslope;
                r(i).gc_ln_beta = res.gc_ln.beta;
                r(i).pred(2,:) = res.gc_ln.pred;
                
                % static strf result (exponent)
                r(i).pred(3,:) = exp(res.strf_fit.pred + ...
                                     res.strf_fit.coeffs(1));
                
                % gain parameters for gc_ln
                r(i).gc_ln.params = res.gc_ln.ahat;
                
                % fit a time-varying exponential
                OPS.fs = 1/.025;
                OPS.nbins = 50;
                OPS.weighting = true;
                OPS.modelType = 'exponential';
                y_tau = ln_tau(r(i).gc_ln.params,tau,OPS);
                keyboard
                                
            end
            
            % model results
            r(i).beta0{m} = res.beta0;
            r(i).beta1{m} = res.beta1;
            r(i).beta2{m} = res.beta2;
            r(i).beta3{m} = res.beta3;
            r(i).w(m,:) = res.w;
            r(i).pred(3+m,:) = res.gain_fit.pred;
            
            % reshape w (by 4 trials to get all transitions with padding)
            r(i).ws(m,:) = mean(reshape(res.w,ops.blockLength*4*ops.fs,[])',1);
            
            if isfield(res,'basis')
                nb = size(res.basis,2);
                if m == 1
                    r(i).b2_kernel{m} = res.basis*res.beta2;
                    r(i).b3_kernel{m} = res.basis*res.beta3;
                elseif m == 2
                    r(i).b2_kernel{m} = [res.basis*res.beta2(1:nb),...
                                        res.basis*res.beta2(nb+1:end)];
                    r(i).b3_kernel{m} = res.basis*res.beta3;
                elseif m == 3
                    r(i).b2_kernel{m} = [res.basis*res.beta2(1:nb),...
                                        res.basis*res.beta2(nb+1:end)];
                    r(i).b3_kernel{m} = [res.basis*res.beta3(1:nb),...
                                        res.basis*res.beta3(nb+1:end)];
                end
            end
            
        end
        
        % reshape y
        r(i).ys = reshape(r(i).y,ops.blockLength*ops.fs,[])';
        
        % compute corr and r2 for each model for each contrast
        for ii = 1:size(r(i).pred,1)
            r(i).ps(ii,:,:) = reshape(r(i).pred(ii,:),ops.blockLength*ops.fs,[])';
            for jj = 1:2
                I = ops.order_r(:,1) == ops.contrast(jj);
                r(i).corr(ii,jj) = ...
                    corr(mean(squeeze(r(i).ys(I,:)))',...
                         mean(squeeze(r(i).ps(ii,I,:)))');
                r(i).r2(ii,jj) = r(i).corr(ii,jj).^2;
            end
        end
        
        % compute corr and r2 across contrasts
        yt = reshape(r(i).y,ops.blockLength*ops.fs*2,[])';
        for ii = 1:size(r(i).pred,1)
            pt = reshape(r(i).pred(ii,:),ops.blockLength*ops.fs*2,[])';
            r(i).corr_all(ii) = corr(mean(yt)',mean(pt)');
            r(i).r2_all(ii) = r(i).corr(ii).^2;
        end
        
        % compute excluding transition window
        I = repmat([false(1,ops.fs) true(1,ops.fs*(ops.blockLength-1))],1,2);
        for ii = 1:size(r(i).pred,1)
            pt = reshape(r(i).pred(ii,:),ops.blockLength*ops.fs*2,[])';
            r(i).corr_not(ii) = corr(mean(yt(:,I))',mean(pt(:,I))');
            r(i).r2_not(ii) = r(i).corr_not(ii).^2;
        end
        
        % cell stats
        r(i).mfr = u.mfr;
        r(i).trough_peak = u.trough_peak;
        r(i).peak_inflect = u.peak_inflect;
        r(i).FWHM = u.FWHM;
        
        % plot individual results
        if plot_on
            blks = 4; % number of blocks to repeat
            twin = blks * ops.stimInfo.baseNoiseD;
            ops.mdls = models;
            plot_acute_neuron_mdls(res,r(i),u,s,ops,twin);
            fn = sprintf('./_plots/%s.pdf',[tag r(i).cellID]);
            saveFigPDF(f1,sz,fn);
            clf(f1);
        end
        
        toc;
        
    end
    
    % model labels
    ops.mdls = models;
    ops.dc = dc;

    save(resfile,'r','ops');

else
    
    load(resfile);

end