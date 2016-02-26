function fig7();

    load('../CuratedCells.mat');
    
    mypath = mfilename('fullpath');
    basepath = mypath(1:(end-12));
    addpath(basepath);
    humanPsychData = psych_analysis();
    %monkeyPsychData = monkey_analysis()


    initials = {'AD (H)','DS (H)','MS (H)','SH (H)','SR (H)','Lem_{fov} (M)','Lem_{ecc} (M)'};
    psychData = [humanPsychData;0.63;0.64]; % This is for lem

    %psychData = sqrt(2)*norminv(psychData,0,1); % Convert to d'

    nSubs = length(psychData);

    cols = gen_diverse_colors(nSubs);


    myFig=figure(); hold on;
    

    
    dds = [5,24];
    for dd = 1:length(dds);
        
        figure(myFig);
        subplot(1,2,dd); hold on;
        
        if dd == 2;

            for sub = 1:nSubs;
                plot([psychData(sub),psychData(sub)],[0,0.3], ...
                    'color',cols(sub,:),'linewidth',2,'linestyle',':');
                
                %if sub < 6
                text(0.9,0+0.05*sub,sprintf(initials{sub}),'color',cols(sub,:),'fontsize',15);
                %else
                %    text(psychData(sub)*0.95,0.35,sprintf(initials{sub}),'color',cols(sub,:),'fontsize',12);
                %end
            end
        end
        
        
        currentCells = Base(dd).Cells;
        
        
        
        AUCs = cat(1,currentCells.HMauc);
        dprime = aROC2dprime(AUCs);
        
        hmReg = cat(1,currentCells.regHm);
        hmSlope = hmReg(:,2);
        fnames = {currentCells.filename};
        
        hold on;
        plot([0.5,0.5],[-1,1],'k --','linewidth',2);
        plot([-1,100],[0,0],'k --','linewidth',2);
        
        
        
        for j = 1:length(hmSlope);
            % Map the recording type (foveal vs eccentric) to the
            % appropriate colour
            rec_number = strcmp(fnames{j},Base(dd).exptlist);
            pen_info = Base(dd).penetrationlist(:,rec_number);
            dist_origin = sqrt(sum((pen_info).^2));
            
            col = cols(6+(dist_origin > 8),:);
            
            plot(AUCs(j),hmSlope(j),'k o','markerfacecolor',col,'markersize',10);
        end
        
        xlabel('AUROC','fontsize',18);
        ylabel('Half-matched slope','fontsize',18)
        set(gca,'xtick',0:0.2:1,'ytick',-1:0.1:1);
        xlim([0.4,1]); ylim([-0.1,0.4]);
        

        

        title(sprintf('%i%% dot density',dds(dd)));
        
    end
    
    
    
    
    set_plot_params(myFig);
    

    savefig('fig7.fig');
    

end

function cols = gen_diverse_colors(N,minDiff);
    % Usage: cols = gen_diverse_colors(N,<minDiff>);
    % Function to return a diverse set of colors that differ by
    % at least some absolute distance (can be specified as varargin)
    % N is number of colors
    % Optional argument minDiff gives the minimum
    % squared distance between colors (increasing this will give larger
    % differences between colors).
    if nargin == 1
        minDiff = 0.3;
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

    %fprintf('Counter: %i\n',counter);
    end

end