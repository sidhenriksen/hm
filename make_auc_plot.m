function make_auc_plot(myFig,evt);

    load('/sid/Ephys/HalfMatched/hm/CuratedCells.mat');


    myFig=figure(); hold on;
    

    
    dds = [1,2];
    for dd = 1:length(dds);
        figure(myFig);
        subplot(1,2,dd); hold on;
        
        if dd == 2;
            % This has one entry for each observer
            psychData = load('hm_psych_data.mat');
            psychData = psychData.hmPerformance;

            initials = {'AD (H)','DS (H)','MS (H)','SH (H)','SR (H)','Lem (M)'};
            psychData = [psychData;0.83];

            nSubs = length(psychData);

            cols = gen_diverse_colors(nSubs);

            for sub = 1:nSubs;
                plot([psychData(sub),psychData(sub)],[0,0.3], ...
                    'color',cols(sub,:),'linewidth',2,'linestyle',':');
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

function cols = gen_diverse_colors(N,varargin);
    % Usage: cols = gen_diverse_colors(N,[minDiff]);
    % Function to return a diverse set of colors that differ by
    % at least some absolute distance (can be specified as varargin)
    % N is number of colors
    % Optional argument minDiff gives the minimum
    % squared distance between colors (increasing this will give larger
    % differences between colors).
    if nargin == 1
        minDiff = 0.3;
    else
        minDiff = varargin{1};
    end

    cols = rand(N,3);

    for c = 2:N;
       previousCols = cols(1:c-1,:);
       currentCol = repmat(cols(c,:),[size(previousCols,1),1]);

       allDiffs = sum((previousCols-currentCol).^2,2);
       counter = 0;
       while any(allDiffs < minDiff);
           counter = counter+1;
           cols(c,:) = rand(1,3);
           currentCol = repmat(cols(c,:),[size(previousCols,1),1]);

           allDiffs = sum((previousCols-currentCol).^2,2);
           
           if counter > 1000;
               error('Too many iterations; increment allowable iterations, change minDiff or use smaller N');
           end
               
       end

    fprintf('Counter: %i\n',counter);
    end

end