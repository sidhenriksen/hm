function fig3()


    load('../CuratedCells.mat');    
    
    ddIndex = [1,2];
    fig = figure();
    cellIndex = [1,1];
    dds = [5,24];
    for dd = 1:length(ddIndex);
        currentCells = Base(ddIndex(dd)).Cells;
        subplot(1,2,dd); hold on;
        
        % Add bisection line
        plot([-1,2],[0,0],'k --','linewidth',2);
        
        DDI = [currentCells.DDI];
        hmReg= cat(1,currentCells.regHm);
        hmSlope = hmReg(:,2);
        Phm = [currentCells.Phm];
        
        for k = 1:length(DDI);
            if ~any(cellIndex==k);
                if Phm(k) > 0.01
                    col = [0.1,0.6,0.15];
                    s = 'o';
                else
                    col = [0.8,0.1,0.1];
                    s = '^';
                end
                plot(DDI(k),hmSlope(k),s,'markerfacecolor',col,'markeredgecolor','k',...
                    'linewidth',1,'markersize',12);
            end
        end
        
        
        
        xlabel('Disparity Discrimination Index');
        
        set(gca,'ytick',[-1:0.1:1],'xtick',0:0.25:1);
        ylim([-0.1,0.4]); xlim([0,1]);
        
        if dd == 2;
            set(gca,'yticklabel',[]);
        else
            ylabel('Half-matched slope');
        end
        
        plot(DDI(cellIndex(dd)),hmSlope(cellIndex(dd)),'s','markerfacecolor',[0.1,0.1,0.8],'markeredgecolor','k',...
            'linewidth',1,'markersize',15);
        [r,p] = corr(DDI',hmSlope);
        fprintf('Density %i: N=%i\n',dd,sum(Phm<0.01));
        
        title(sprintf('%i%% dot density',dds(dd)));
        
    end

    
    set_plot_params(fig);
    set(gcf,'position',[200,200,1300,500]);
    savefig('fig3.fig');    

end