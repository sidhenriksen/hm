function PlotCurated
    plotDots = 1;
    plotSlope = 1;


    load('CuratedCells.mat');

    densities = cat(1,Base.density);
    
    % This is included just for the sake of colors
    fileNames = {};
    myCols = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;rand(10,3)];
    monkeys = {'jbe','lem'};
    if plotDots
        figure();
        %densities = densities(1);
        count = 0;
        for dd = [1,4];
            count = count+1;
            
            currentData = Base(dd).Cells;
            nCells = length(currentData);
            
            cols = rand(nCells,3);
            subplot(1,round(length(densities)/2),count); hold on;
            hold on;
            for cell = 1:nCells
                fileNameExists = strcmp(currentData(cell).filename, fileNames);
                if any(fileNameExists);
                    currentColor = myCols(fileNameExists,:);
                else
                    fileNames{length(fileNames)+1} = currentData(cell).filename;
                    currentColor = myCols(length(fileNames),:);
                end

                % Skip if we have NaNs...
                if ~isnan(currentData(cell).DDI);
                    currentMonkey = currentData(cell).filename(9:11);

                    % Only plot monkeys in monkeys
                    if strcmp(currentMonkey,monkeys{1});
                        plot(currentData(cell).regHm(1),currentData(cell).regAc(1),'k','marker','s','markersize',round(1+currentData(cell).DDI*20),'markerfacecolor',currentColor);
                    elseif strcmp(currentMonkey,monkeys{2});
                        plot(currentData(cell).regHm(1),currentData(cell).regAc(1),'k','marker','o','markersize',round(1+currentData(cell).DDI*20),'markerfacecolor',currentColor);
                    end


                end
            end
            
            
            hm = cat(1,currentData.regHm);
            ac = cat(1,currentData.regAc);
            XYData.x = hm(:,1);
            XYData.y = ac(:,1);
            setappdata(gca,'XYData',XYData);
            setappdata(gca,'Cells',currentData);
            setappdata(gca,'subplot',count);

            xlabel('Correlated-halfmatched r','fontsize',24)
            ylabel('Correlated-AC r', 'fontsize', 24);
            xlim([-1,1]); ylim([-1,1]);
            set(gca,'XTick',-1:0.5:1,'fontsize',20);
            set(gca,'YTick',-1:0.5:1);

            title(sprintf('Density: %i%', densities(dd)),'fontsize', 25);
            
            
            % Put some text in the window showing how many cells
            jbeCount = 0;
            lemCount = 0;
            for j = 1:nCells
                jbeCount = jbeCount + strcmp(currentData(j).filename(9:11),'jbe');
                lemCount = lemCount + strcmp(currentData(j).filename(9:11),'lem');
            end
            text(-0.75,0.8,['jbe: ',num2str(jbeCount),' cells']);
            text(-0.75,0.7,['lem: ',num2str(lemCount),' cells']);
                

            % Draw some lines to split into quadrants
            line([0,0],[-1,1],'linewidth',1,'color','red','linestyle','--');
            line([-1,1],[0,0],'linewidth',1,'color','red','linestyle','--');
        end
        set(gcf,'Color','White','windowbuttondownfcn',@TC_callback)
    end

    if plotSlope
        % This is included just for the sake of colors
        fileNames = {};
        if ~exist('myCols','var');
            myCols = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;rand(10,3).^2];
        end

        monkeys = {'jbe','lem'};
        a=figure();
        
        %densities = densities(1);
        count = 0;
        for dd = [1,4];
            count = count+1;
            
            currentData = Base(dd).Cells;
            nCells = length(currentData);

            cols = rand(nCells,3);
            subplot(1,round(length(densities)/2),count); hold on;

            for cell = 1:nCells
                fileNameExists = strcmp(currentData(cell).filename, fileNames);
                if any(fileNameExists);
                    currentColor = myCols(fileNameExists,:);
                else
                    fileNames{length(fileNames)+1} = currentData(cell).filename;
                    currentColor = myCols(length(fileNames),:);
                end

                % Skip if we have NaNs...
                if ~isnan(currentData(cell).DDI);
                    currentMonkey = currentData(cell).filename(9:11);
                    if strcmp(currentMonkey,monkeys{1});
                        plot(currentData(cell).regHm(2),currentData(cell).regAc(2),'k','marker','s','markersize',...
                            round(1+currentData(cell).DDI*20),'markerfacecolor',currentColor);
                    elseif strcmp(currentMonkey,monkeys{2});
                        plot(currentData(cell).regHm(2),currentData(cell).regAc(2),'k','marker','o','markersize', ...
                            round(1+currentData(cell).DDI*20),'markerfacecolor',currentColor);
                    end


                end
            end
            
            hm = cat(1,currentData.regHm);
            ac = cat(1,currentData.regAc);
            XYData.x = hm(:,2);
            XYData.y = ac(:,2);
            setappdata(gca,'XYData',XYData);
            setappdata(gca,'Cells',currentData);
            setappdata(gca,'subplot',count);

            xlabel('Correlated-halfmatched slope','fontsize',24)
            ylabel('Correlated-AC slope', 'fontsize', 24);
            ylim([-1,1]);
            %set(gca,'XTick',-1:0.5:1,'fontsize',20);
            set(gca,'YTick',-1:0.5:1,'fontsize',20);

            title(sprintf('Density: %i%', densities(dd)),'fontsize', 30);
            if dd==1;
                currentXlim = [-0.35,0.35];
            end
            
            
                
            
            % Put some text in the window showing how many cells
            jbeCount = 0;
            lemCount = 0;
            for j = 1:nCells
                jbeCount = jbeCount + strcmp(currentData(j).filename(9:11),'jbe');
                lemCount = lemCount + strcmp(currentData(j).filename(9:11),'lem');
            end
            text(-0.25,0.8,['jbe: ',num2str(jbeCount),' cells']);
            text(-0.25,0.7,['lem: ',num2str(lemCount),' cells']);
            
            % Draw some lines to split into quadrants
            plot([0,0],[-1,1],'linewidth',1,'color','red','linestyle','--');
            plot([-1,1],[0,0],'linewidth',1,'color','red','linestyle','--');
            xlim(currentXlim);
            
        end
        set(gcf,'Color','White','windowbuttondownfcn',{@TC_callback});

        
    end
    

end

function TC_callback(myFig,evt,varargin);

    currentPoint = get(gca,'currentpoint');
    
    axisData = getappdata(gca);
    figData = getappdata(gcf);
    
    % Get data
    Cells = axisData.Cells;
    XYData = axisData.XYData;
    subplotNumber = axisData.subplot;
    
    % Find point with closest Euclidean distance
    distance = (XYData.x-currentPoint(1,1)).^2 + (XYData.y - currentPoint(2,2)).^2;
    
    if nargin == 2
        whichCell = find(distance == min(distance)); 
        whichCell = whichCell(1); % If equidistant, just choose first one
    else
        otherID = varargin{1};
        for j = 1:length(Cells);
            
            uniqueID = quickhash([Cells(j).filename, num2str(Cells(j).cellnumber)]);
            if uniqueID == otherID;
                whichCell = j;
                break
            else
                whichCell = 0;
            end
        end
   
    end
    

    
    % If this hasn't been called before, assign a new figure handle
    if ~isfield(figData,'TCHandle');
        newFig = figure();
        set(newFig,'color','white');
        figData.TCHandle = newFig;
        setappdata(myFig,'TCHandle',newFig);
    end
    
    figure(figData.TCHandle);
    
    if subplotNumber > 2;
        maxSubplots = subplotNumber;
    else
        maxSubplots = 2;
    end
    
    
    subplot(1,maxSubplots,subplotNumber);
    if whichCell > 0;
        cellName = [Cells(whichCell).filename(18:24),'-cell',num2str(Cells(whichCell).cellnumber)];
        
        set(gca,'fontsize',18);
        Dx = Cells(whichCell).Dx;
        plot(Dx,Cells(whichCell).halfmatchedResponse*2.5,'b -- o', ...
            Dx,Cells(whichCell).correlatedResponse*2.5,'r -- o', ...
            Dx,Cells(whichCell).anticorrelatedResponse*2.5,'k -- o', ...
            'linewidth',3,'markersize',5,'markerfacecolor','k');
    elseif nargin > 2
        cla
        cellName = '';
    else
        error('Error: This shouldn''t happen...');
    end
    
    xlabel('Disparity (deg)','fontsize',20);
    ylabel('Spikes per second','fontsize',20);
    density = (subplotNumber==1)*5 + (subplotNumber==2)*24;
    title(['Density: ', num2str(density),'; ',cellName],'fontsize',22);

        

    selectionType = get(myFig,'selectionType');
    shiftPressed = strcmp(selectionType,'extend');
    ctrlPressed = strcmp(selectionType,'alt');
    
    % This will only get executed if uniqueID is not passed to TC_callback
    % (which we only do when we call TC_callback recursively)

    
    if shiftPressed || ctrlPressed

        % If this hasn't been called before, assign a new figure handle
        if ~isfield(figData,'scatterHandle');
            scatterFig = figure();
            set(scatterFig,'color','white');
            figData.scatterHandle = scatterFig;
            setappdata(myFig,'scatterHandle',scatterFig);
        end
    
        figure(figData.scatterHandle);
    
        subplot(1,maxSubplots,subplotNumber)
        if whichCell > 0;
            set(gca,'fontsize',18);
            if shiftPressed
                plot(Cells(whichCell).correlatedResponse*2.5,Cells(whichCell).halfmatchedResponse*2.5,'o k',...
                    'markersize',4,'markerfacecolor','b');
                xlabel('Correlated response','fontsize',20);
                ylabel('Halfmatched response','fontsize',20);
                b = Cells(whichCell).regHm(2); m = Cells(whichCell).regHm(3);
                corrMin = min(Cells(whichCell).correlatedResponse); corrMax = max(Cells(whichCell).correlatedResponse);
                delta = corrMax*0.25; corrMin = corrMin-delta; corrMax = corrMax+delta;
                line([corrMin,corrMax]*2.5,[m+corrMin*b,m+corrMax*b]*2.5,'linewidth',2,'linestyle','-','color','r');
                title(['r= ', num2str(Cells(whichCell).regHm(1))])
                
            else
                plot(Cells(whichCell).correlatedResponse*2.5,Cells(whichCell).anticorrelatedResponse*2.5,'o k',...
                    'markersize',4,'markerfacecolor','r');
                xlabel('Correlated response','fontsize',20);
                ylabel('Anticorrelated response','fontsize',20);
                b = Cells(whichCell).regAc(2); m = Cells(whichCell).regAc(3);
                corrMin = min(Cells(whichCell).correlatedResponse); corrMax = max(Cells(whichCell).correlatedResponse);
                delta = corrMax*0.25; corrMin = corrMin-delta; corrMax = corrMax+delta;
                line([corrMin,corrMax]*2.5,[m+corrMin*b,m+corrMax*b]*2.5,'linewidth',2,'linestyle','-','color','k');
                title(['r= ', num2str(Cells(whichCell).regAc(1))])
            end
            %reg = 
        elseif nargin > 2
            cla
        else
            error('Error: This shouldn''t happen...');
        end
    end

    if nargin == 2
        figure(myFig); 
        if subplotNumber == 1;
            subplot(1,2,2);
        elseif subplotNumber == 2;
            subplot(1,2,1);
        end

        uniqueID = quickhash([Cells(whichCell).filename, num2str(Cells(whichCell).cellnumber)]);

        TC_callback(myFig,evt,uniqueID);

    end
    
end

function myHash = quickhash(string)
    string = double(string);
    myHash = mod(sum(string(:).^2.3),2^32-1);
end