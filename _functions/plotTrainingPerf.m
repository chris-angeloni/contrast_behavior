function plotTrainingPerf(t,mouseList,colors)

% for each mouse, find the date of the first high contrast and low
% contrast session, then make a vector of performance per session
for i = 1:length(mouseList)
    for j = 1:2
        
        % first date and date vector
        mI = contains(t.mouse,mouseList{i});
        if size(mI,1) == 1
            mI = mI';
        end
        I = t.contrast == j-1 & mI;
        firstDate{i,j} = t.date(find(I,1,'first'));
        expDays{i,j} = round(days(t.date(I) - firstDate{i,j}));
        perfDays{i,j} = t.mean_pc(I);
        n(i,j) = sum(I);
        
    end
end

% make this padded with nans to have better plots
dayLabel = [0 max(n(:))];
perfMat = nan(max(n(:))+1,length(mouseList),2);
for i = 1:length(mouseList)
    for j = 1:2
        perfMat([expDays{i,j}+1]',i,j) = perfDays{i,j};
    end
end
perfMat(perfMat==0) = nan;
firstHighContrast = days(vertcat(firstDate{:,1}) - vertcat(firstDate{:,2})) ...
    > 0;

% plot
hold on;
for i = 1:2
    I = firstHighContrast==(i-1);
    plot(perfMat(:,I,i),'.','Color',colors(i,:)+(colors(i,:)==0)*.8);
    ph(i) = plot(movmean(nanmedian(perfMat(:,I,i),2),7,'omitnan'),...
                 'Color',colors(i,:),'LineWidth',1);
    text(.8,.1+(.1*(i-1)),sprintf('n=%d',sum(I)),...
         'units','normalized','color',colors(i,:));
end
ch = plot(xlim,[.5 .5],'k'); axis tight;
uistack(ph,'top');
xlabel('Time (days rel. task exposure)');
ylabel('Percent Correct'); plotPrefs;
hold off;