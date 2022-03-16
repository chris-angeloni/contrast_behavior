function h = plotPval(p,x,y,units)

YL = ylim;
y = y + range(YL)*.1;

if ~exist('units','var') | isempty(units)
    units = 'normalized';
end

hold on;
[psym,pval] = pvalStr(p);
h = text(x,y,sprintf('p=%s\n%s',pval,psym),...
         'units',units,...
         'horizontalAlignment','center');
hold off;