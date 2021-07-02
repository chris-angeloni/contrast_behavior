function [mouseSNR,mousePerf,mouseContrast] = plotMouseCurves(t,mouseList,mdl,colors,ms)

mouseSNR = []; mousePerf = []; mouseContrast = [];
for i = 1:numel(mouseList)
    
    subplot(5,5,i); hold on;
    plot([-10 30],[.5 .5],'k');
    for j = 1:2
        
        % extract data for this mouse and contrast
        mI = contains(t.mouse,mouseList{i});
        if size(mI,1) == 1
            mI = mI';
        end
        I = t.contrast == j-1 & mI;
        
        if sum(I) > 0
            x = round(t(I,:).vols,1);
            xf = linspace(min(x(:)),max(x(:)),100);
            y = t(I,:).pc;
            
            % save out mean performance per SNR for later
            usnrs = grpstats(x(:),x(:),'mean');
            mouseSNR = [mouseSNR; usnrs];
            mousePerf = [mousePerf; grpstats(y(:),x(:),'mean')];
            mouseContrast = [mouseContrast; repmat(j-1,length(usnrs),1)];
            
            % plot it
            plot(x(:),y(:),'.','color',colors(j,:)+(colors(j,:)==0)*.7,...
                 'MarkerSize',ms)
            mp = median(t.prms(I,:)); mt = mp(1)/mp(2);
            plot(xf,mdl(mp,xf),'color',colors(j,:),'LineWidth',1);
            plot([mt mt],[.5 mdl(mp,mt)],'color',colors(j,:));
            plotPrefs;
        end
        
    end
    title(mouseList{i});
    axis tight; ylim([0 1]);
    
    if i == numel(mouseList)
        xlabel('Target Volume (dB SNR)');
        ylabel('Percent Correct');
    end
end