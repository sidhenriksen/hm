function make_auc_plot(myFig,evt);

    load('/sid/Ephys/HalfMatched/hm/CuratedCells.mat');


    myFig=figure(); hold on;
    

    
    dds = [1,4];
    for dd = 1:length(dds);
        figure(myFig);
        subplot(1,2,dd); hold on;
        
        if dd == 2;
            % This has one entry for each observer
            psychData = load('hm_psych_data.mat');
            psychData = psychData.hmPerformance;

            initials = {'AD','DS','MS','SH','SR'};

            nSubs = length(psychData);

            cols = rand([nSubs,3]);

            for sub = 1:nSubs;
                plot([psychData(sub),psychData(sub)],[0,0.3], ...
                    'color',cols(sub,:),'linewidth',3,'linestyle',':');
                text(1,0.125+0.03*sub,initials{sub},'color',cols(sub,:));
            end
        end
        
        
        currentCells = Base(dds(dd)).Cells;
        
        AUCs = cat(1,currentCells.HMauc);
        hmReg = cat(1,currentCells.regHm);
        hmSlope = hmReg(:,2);
        
        
        plot(AUCs,hmSlope,'k o','markerfacecolor',[0.1,0.4,0.1],'markersize',10);
        
        xlabel('ROC AUC','fontsize',18);
        ylabel('Half-matched slope','fontsize',18)
        set(gca,'fontsize',16,'xtick',0:0.25:1,'ytick',-.3:0.15:0.3);
        xlim([0,1]); ylim([-0.3,0.3]);
        
        hold on;
        plot([0.5,0.5],[-1,1],'r --','linewidth',2);
        plot([-1,1],[0,0],'r --','linewidth',2);
        

        
        
    end
    
    
    
    
    set(myFig,'color','white')
    

    
    

end