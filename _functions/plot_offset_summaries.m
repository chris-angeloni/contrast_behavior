function res = plot_offset_summaries(r)


%% analysis options
% col indices to average
clear ind;
ind{1} = [1 2]; % [1 2];
ind{2} = [3 4];
ind{3} = [5 6];
ind{4} = [7 8];
ind{5} = [9 10];

% inclusions
sig_auc_thresh = sum(squeeze(r.auc_sig(:,1,:)),2,'omitnan');
sig_auc_high = sum(squeeze(r.auc_sig(:,2,:)),2,'omitnan');
bad_mouse1 = (contains(r.mouse,'CA121') & r.contrastI==1);

% exclude CA122 lo->hi (couldn't learn), include only population
% data withmore than 2 significant volumes, and sessions with fewer
% than 25% false alarm rate
% ~bad_mouse1 & sig_auc_thresh > 1 & r.fa < .25; <- good svm tau results/badfdr
% ~bad_mouse1 & sig_auc_thresh > 2 & r.fa < .5; <- good critp tau/fdr
% ~bad_mouse1 & sig_auc_thresh > 2 & r.fa < .25; <- all tau ns/good fdr
% ~bad_mouse1 & sig_auc_thresh > 3 & r.fa < .2; <- okay for sig decoders
include = ~bad_mouse1 & sig_auc_thresh > 2 & r.fa < .5;

% format sessions as strings
r.session = cellstr(num2str(r.sessionID));

% options
x = combindstats(r.offs,ind);
grpvar = 'mouse'; %'session';
lags = [1:length(x)-1];
ugrp = unique(r.(grpvar));
bootstrap = false;



%% setup
% plot stuff
nrows = 8;
ncols = 3;
yl = [.5 1];


% set up plot colors
colors{1}{1} = [0.5 0.5 0.5];
colors{1}{2} = [0.7 0.7 1.0];
colors{1}{3} = [0.0 0.0 1.0];
colors{2}{1} = [0.5 0.5 0.5];
colors{2}{2} = [1.0 0.7 0.7];
colors{2}{3} = [1.0 0.0 0.0];

% colors rgb
cols = [0 0 1; 1 0 0];

% pre compute PC for desired fields
fields = {'beh_rate','critp','svm_rate','mean_critp'};
for i = 1:length(fields)
    fn = sprintf('%s_PC',fields{i});
    r.(fn) = squeeze(PDtoPC(r.(fields{i})(:,[2 3],:),r.(fields{i})(:,1,:)));
end

% compute stats for all fields for each contrast
fields = {'beh_rate_PC','critp_PC','svm_rate_PC','mean_critp_PC',...
         'auc','mean_auc','prop_sig_neurons'};
labels = {'PC_{behavior}','PC_{population criterion}','PC_{svm}',...
          'PC_{sig. neurons criterion}','PC_{population auc}',...
          'PC_{sig. neurons auc}',['Proportion Significant ' ...
                    'Neurons']};

clear res;
for i = 1:length(fields)
    % volume
    for j = 1:2
        
        prms = nan(2,length(ugrp),3);
        taus = nan(2,length(ugrp));
        n2 = nan(2,length(ugrp));
        mean_perf = nan(2,length(ugrp),length(x));
        
        % contrast
        for k = 1:2
            
            fprintf('%s: volume %d, contrast %d... ',fields{i},j,k); tic;
            
            % column mean and error
            d = r.(fields{i})(r.contrastI==(k-1) & include,j,:);
            [ym,ye,dat] = combindstats(d,ind);
            res(j).(fields{i}).mean(k,:) = ym;
            res(j).(fields{i}).sem(k,:) = ye;
            
            % sign rank test of first time point against all others
            tmp = squeeze(mean(cat(3,dat{:}),2,'omitnan'));
            for ii = 1:4
                [p(k,ii),~,stat(k,ii)] = signrank(...
                    tmp(:,1),tmp(:,ii+1));
            end
            n(k) = length(tmp);
            
            % fit taus for each group
            grps = r.(grpvar)(r.contrastI==(k-1) & include);
            mtmp = grpstats(tmp(:,lags),grps);
            this_grp = unique(grps);
            n2(k,contains(ugrp,this_grp)) = mode(grpstats(tmp(:,lags),grps,'numel'),2);
            for ii = 1:size(this_grp,1)
                mI = contains(ugrp,this_grp(ii));
                [prm,mdl,tau] = fitExpGrid(x(lags),mtmp(ii,:));
                prms(k,mI,:) = prm;
                taus(k,mI) = tau;
            end
            mean_perf(k,contains(ugrp,this_grp),:) = grpstats(tmp,grps);
            
            if bootstrap

                % tau on mean
                [prm,mdl,tau] = fitExpGrid(x(lags),mean(tmp(:,lags),1,'omitnan'));
                res(j).(fields{i}).exp_mean_prms(k,:) = prm;
                res(j).(fields{i}).exp_mean_tau(k) = tau;
                
                % bootstrapped
                fprintf('bootstrapping... ');
                parfor its = 1:500
                    samp = randsample(size(tmp,1),size(tmp,1),true);
                    [~,~,tau_boot(its)] = fitExpGrid(x(lags),mean(tmp(samp,lags),1,'omitnan'));
                end
                res(j).(fields{i}).exp_boot_tau(k,:) = tau_boot;

            end
            toc;
            
        end
        % fdr correction for all tests and save timing stats
        res(j).(fields{i}).h_fdr = fdr_bh(p);
        res(j).(fields{i}).p = p;
        res(j).(fields{i}).stats = stat;
        res(j).(fields{i}).test = 'signrank';
        res(j).(fields{i}).n = n;
        res(j).(fields{i}).label = labels{i};
        
        % stats on exponential fits
        if contains(grpvar,'mouse')
            [p,h,stats] = signrank(taus(1,:),taus(2,:));
            test = 'signrank';
        else
            [p,h,stats] = ranksum(taus(1,~isnan(taus(1,:))),...
                                  taus(2,~isnan(taus(2,:))));
            test = 'ranksum';
        end
        
        res(j).(fields{i}).groups = ugrp;
        res(j).(fields{i}).mean_perf = mean_perf;
        res(j).(fields{i}).exp_prm = prms;
        res(j).(fields{i}).exp_tau = taus;
        res(j).(fields{i}).exp_n = sum(~isnan(taus),2);
        res(j).(fields{i}).exp_p = p;
        res(j).(fields{i}).exp_h = h;
        res(j).(fields{i}).exp_stats = stats;
        res(j).(fields{i}).exp_test = test;
        res(j).(fields{i}).exp_sess_count = n2;
        disp(res(j).(fields{i}).exp_p);
        
    end
end


%% plot individual mice
f12 = figure(12); clf;
xf = linspace(min(x),max(x),100);
cnt = 1;
for j = 1:size(ugrp,1)
    for i = 1:length(fields)
        subplot(size(ugrp,1),length(fields),cnt); hold on;
        if j == 1; title(fields{i},'interpreter','none'); end
        if i == 1; ylabel(ugrp{j}); end
        for k = 1:2
            plot(x,squeeze(res(1).(fields{i}).mean_perf(k,j,:)),...
                 'o','color',colors{k}{end});
            plot(xf,mdl(squeeze(res(1).(fields{i}).exp_prm(k,j,:)),xf),...
                 'color',colors{k}{end})
            %ylim([.5 1]); 
            plotPrefs;
            
        end
        
        cnt = cnt + 1;
    end
end
saveFigPDF(f12,[1200 1400],'./_plots/_offset_summary_allmice.pdf');


    


%% plot each field (volume == new fig)


% volume
for j = 1:2
    fh(j) = figure(j*10); clf;
    
    % field
    for i = 1:length(fields)
        
        subplot(nrows,ncols,[i i+1]+2*(i-1)); hold on
        
        % contrast
        clear eh;
        for k = 1:2
            
            % error plot
            eh(k) = errorBars(x,res(j).(fields{i}).mean(k,:),...
                              colors{k}{end},[],[],...
                              res(j).(fields{i}).sem(k,:),[],...
                              'Marker','.','MarkerSize',20, ...
                              'LineWidth',1);
        end
        YL = ylim;
        
        % plot the test results
        cnt = 0;
        for k = 1:2
            for ii = 1:4
                if res(j).(fields{i}).h_fdr(k,ii)
                    cnt = cnt + 1;
                    plot([x(1) x(ii+1)],...
                         [YL(2)+.05*cnt*diff(YL) YL(2)+.05*cnt*diff(YL)],...
                         'color',colors{k}{end},'linewidth',.75)
                end
            end
        end
        axis tight;
        eh(end+1) = plot([0 0],ylim,'k'); plotPrefs;
        if i == 1
            legend(eh,'H->L','L->H','switch','location','se');
        end
        set(gca,'xtick',x); set(gca,'xticklabels',num2str(round(x'*1000)));
        ylabel(res(j).(fields{i}).label);  xlim([-.05 1.05]);
        
        subplot(nrows,ncols,[i+2]+2*(i-1)); hold on
        plot(res(j).(fields{i}).exp_tau*1000,'k');
        xt = [ones(1,size(res(j).(fields{i}).exp_tau,2)); ...
              ones(1,size(res(j).(fields{i}).exp_tau,2))*2];
        scatter(xt(:),res(j).(fields{i}).exp_tau(:)*1000,100,...
                cols(xt(:),:),'.');
        YL = ylim;
        if res(j).(fields{i}).exp_h
            plot([1 2],[YL(2) YL(2)]+.05*range(YL),'k','linewidth',1);
        end
        title(sprintf('p=%3.5f',res(j).(fields{i}).exp_p));
        set(gca,'xtick',[1 2]);
        set(gca,'xticklabels',{'Low','High'});
        xlim([.75 2.25]); ylabel('\tau (ms)');
        
    end
    
    % firing rate plot is different
    % mean fr
    lw = [.5 1];
    subplot(nrows,ncols,[22 23]); hold on;
    for k = 1:2
        mkfill = {'none',colors{k}{end}};
        [ym,ye] = combindstats(...
            r.mean_fr(r.contrastI==(k-1) & include,1,:),ind);
        eh = errorBars(x,ym,colors{k}{end},[],[],ye,[],...
                       'Marker','o','MarkerFaceColor',mkfill{1}, ...
                       'LineWidth',lw(1));
        [ym,ye] = combindstats(...
            r.mean_fr(r.contrastI==(k-1) & include,j+1,:),ind);
        eh = errorBars(x,ym,colors{k}{end},[],[],ye,[],...
                       'Marker','o','MarkerFaceColor',mkfill{2},'LineWidth',lw(2));
    end
    plot([0 0],ylim,'k'); plotPrefs; xlim([-.05 1.05]);
    set(gca,'xtick',x); set(gca,'xticklabels',num2str(round(x'*1000)));
    ylabel('FR_{sig neurons}'); xlabel('Time (ms)');
    
    
    subplot(nrows,ncols,[1 2]);
    title(sprintf('Volume %d',j));
    saveFigPDF(fh(j),[500 1000],sprintf('./_plots/_offset_summary_vol%d.pdf',j));
end

   


for i = 1:length(fields)
    
    [i i+1]+2*(i-1)
    
end