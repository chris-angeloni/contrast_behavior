function r = run_pseudopop(res_psych,ops)

keyboard

% preprocess the data
dat = res_psych.single_cell;
FR = cat(1,dat.fr);
FRz = (FR - mean(FR(:,1),2));
[sortfr,sorti] = sort(mean(FRz(:,2:end),2));

% visualize all of the data
clim = [-20 20]; colormap(divergentColors(100,'brewer'));
imagesc(FRz(sorti,:),clim);
caxis(clim); ch = colorbar;
set(get(ch,'Title'),'String','\DeltaFR')
set(gca,'ydir','normal');




