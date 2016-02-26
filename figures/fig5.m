function fig5()


    load('../CuratedCells.mat');
    markersize = 10;
    
    ddIndex = [1,2]; cellIndex = {[1,33],[1,33]};
    fig = figure();
    dds = [5,24];
    for dd = 1:length(ddIndex);
        currentCells = Base(ddIndex(dd)).Cells;
        subplot(1,2,dd); hold on;
        
        % Add bisection line
        plot([0,0],[-2,2],'k --','linewidth',2);
        plot([-2,2],[0,0],'k --','linewidth',2);
        
        acReg = cat(1,currentCells.regAc);
        acSlope = acReg(:,2);
        hmReg= cat(1,currentCells.regHm);
        hmSlope = hmReg(:,2);
        
        plot(acSlope,hmSlope,'o','markerfacecolor',[0.1,0.6,0.15],'markeredgecolor','k',...
            'linewidth',1,'markersize',markersize);
        
        
        
        xlabel('Anticorrelated slope');
        
        set(gca,'ytick',-1:0.1:1,'xtick',-1:0.5:1);
        ylim([-0.1,0.4]); xlim([-1.25,1]);
        
        if dd == 2;
            set(gca,'yticklabel',[]);
        else
            ylabel('Half-matched slope');
        end
        
        [r,p] = corr(hmSlope,acSlope);
        %text(-0.15,0.75,sprintf('r = %.2f, P = %.3f',r,p),'fontsize',15);
        fprintf('r = %.2f, p = %.3f\n',r,p);
        specialCells = cellIndex{dd};
        plot(acSlope(specialCells(1)),hmSlope(specialCells(1)),'s','markerfacecolor',[0.1,0.1,0.8],'markeredgecolor','k',...
            'linewidth',1,'markersize',markersize*1.5);
        plot(acSlope(specialCells(2)),hmSlope(specialCells(2)),'s','markerfacecolor',[0.8,0.1,0.8],'markeredgecolor','k',...
            'linewidth',1,'markersize',markersize*1.5);
        
        
        title(sprintf('%i%% dot density',dds(dd)));
    end
    
    set_plot_params(fig);
    set(gcf,'position',[200,200,1300,500]);
    savefig('fig5.fig');
end