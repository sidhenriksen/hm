function fig2()
    % Example tuning curve for both densities
    % with scatter plots...
    
    % Location of half-matched data    
    load('../CuratedCells.mat')
    
    markersize = 8;
    lw = 3;
    
    densityIndex = [1,2];
    cellIndex = [1,1];
    dds = [5,24];
    fig=figure();
    for dd = 1:2;
        currentCell = Base(densityIndex(dd)).Cells(cellIndex(dd));
        dx = currentCell.Dx;
        
        % Multiply by 2.5 since trial duration is 400 ms (to get
        % spikes/sec)
        c = currentCell.correlatedResponse*2.5;
        hm = currentCell.halfmatchedResponse*2.5;
        ac = currentCell.anticorrelatedResponse*2.5;
        
        subplot(2,2,dd); hold on;
        E1 = errorbar(dx,c,currentCell.correlatedSEM*2.5);
        E2 = errorbar(dx,hm,currentCell.halfmatchedSEM*2.5);
        E3 = errorbar(dx,ac,currentCell.anticorrelatedSEM*2.5);
        
        set(E1,'color','r','linewidth',lw,'markersize',markersize,'linestyle','-',...
            'markerfacecolor','r','markeredgecolor','r','marker','s');
        set(E2,'color','b','linewidth',lw,'markersize',markersize,'linestyle','-',...
            'markerfacecolor','b','markeredgecolor','b','marker','^');
        set(E3,'color','k','linewidth',lw,'markersize',markersize,'linestyle','-',...
            'markerfacecolor','k','markeredgecolor','k','marker','o');
        
        % Set labels
        xlabel('Disparity (deg)');
        if dd == 1; % Shared ylabel
            ylabel('Spikes per second');
        end
        set(gca,'xtick',[-0.9:0.45:0.9],'ytick',0:50:200);
        ylim([0,225])
        
        if dd == 2; % Remove yticklabels for right-hand plot
            set(gca,'yticklabel',[]);
            leg = legend('Correlated','Half-matched','Anticorrelated');
        end
        title(sprintf('%i%% dot density',dds(dd)));
    end
    set_plot_params(fig);
    
    
    for dd = 1:2
        currentCell = Base(densityIndex(dd)).Cells(cellIndex(dd));
        c = currentCell.correlatedResponse*2.5;
        hm = currentCell.halfmatchedResponse*2.5;
        ac = currentCell.anticorrelatedResponse*2.5;
        
        dx = currentCell.Dx;
        
        subplot(2,2,dd+2); hold on;
        
        
        hmSlope = currentCell.regHmRegular(2); hmOffset = currentCell.regHmRegular(3);
        acSlope = currentCell.regAcRegular(2); acOffset = currentCell.regAcRegular(3);
        myC = [min(c),max(c)]; % Without *2.5
        
        hmLine = myC*hmSlope + hmOffset*2.5; 
        acLine = myC*acSlope + acOffset*2.5;
        
        
        plotHm=plot(c,hm,'b o','markerfacecolor','b','markersize',8,'linewidth',2);
        plotAc=plot(c,ac,'k o','markerfacecolor','w','markersize',8,'linewidth',2);
        legend(sprintf('Half-matched (slope: %.2f)',hmSlope),sprintf('Anticorrelated (slope: %.2f)',acSlope));
        
        plot(myC,hmLine,'b -','linewidth',2);
        plot(myC,acLine,'k --','linewidth',2);
        xlabel('Correlated response (spikes/s)');
        
        if dd == 1;
            ylabel('Response (spikes/s)');
        end
        
        set(gca,'xtick',0:50:200,'ytick',0:50:200);
        
        if dd == 2;
            set(gca,'yticklabel',[]);
            %legend([plotHm,plotAc],'Half-matched','Anticorrelated');
        end
        
        
        
        
        
        
        xlim([0,225]); ylim([0,225])
        
    end
    
    set_plot_params(fig);
    
    dd=2;
    figure(fig); b = subplot(2,2,dd);

    currentCell = Base(densityIndex(dd)).Cells(cellIndex(dd));
    dx = currentCell.Dx;
    hm = currentCell.halfmatchedResponse*2.5;


    a=figure();
    hold on;
    E2 = errorbar(dx,hm,currentCell.halfmatchedSEM*2.5);
    set(E2,'color','b','linewidth',2,'markersize',5,'linestyle','-',...
    'markerfacecolor','b','markeredgecolor','b','marker','^');


    
        
    
    set(gca,'ytick',[55,70,85],'ylim',[55,85]);



    xlim([-0.9,0.9])
    set(gca,'xtick',[-0.9,0,0.9],'fontsize',16);


    inset_size = 0.125;
    inset_handle = gca;
    new_handle = copyobj(inset_handle,fig);
    close(a);
    ax=get(b,'Position');
    set(new_handle,'Position',[1.1*ax(1)+ax(3)-inset_size 0.95*ax(2)+ax(4)-inset_size inset_size inset_size])
    
    set(gcf,'Position',[200,200,1200,800]);
    savefig('fig2.fig');
    

end