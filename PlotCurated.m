function PlotCurated2(varargin)

    % Create figure
    myFig=figure();
    
    mh = uimenu(myFig,'Label','Plot');
    mainx = uimenu(mh,'Label','Main: x axis');
    mainy = uimenu(mh,'Label','Main: y axis');
    
    % Takes optional for which curated cells file
    if nargin > 0
        load(varargin{1});
        ddmax = 2;
    else
        load('CuratedCells.mat');
        ddmax = 4;
    end
    

    densities = cat(1,Base.density);
    
    % This is included just for the sake of colors
    fileNames = {};
    myCols = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;rand(10,3)];
    monkeys = {'jbe','lem'};

    setappdata(myFig,'Base',Base);
    
    %%% This is a bit too much work right now. Might do this later.
    
    % X menu
    menux_CHMr = uimenu(mainx,'Label','Correlated-HM r','Checked', 'on','Callback',{@make_subplots,Base,'CHMr'});
    menux_CACr = uimenu(mainx,'Label','Correlated-AC r','Checked', 'off','Callback',{@make_subplots,Base,'CACr'});
    menux_CMeanr = uimenu(mainx,'Label','Correlated-Mean CAC r','Checked', 'off','Callback',{@make_subplots,Base,'CMeanr'});
    
    menux_CHMslope = uimenu(mainx,'Label','Correlated-HM slope','Checked', 'off','Callback',{@make_subplots,Base,'CHMslope'});
    menux_CACslope = uimenu(mainx,'Label','Correlated-AC slope','Checked', 'off','Callback',{@make_subplots,Base,'CACslope'});
    menux_CMeanslope = uimenu(mainx,'Label','Correlated-Mean CAC slope','Checked', 'off','Callback',{@make_subplots,Base,'CMeanslope'});
    
    
    % Y menu
    menuy_CHMr = uimenu(mainy,'Label','Correlated-HM r','Checked', 'on','Callback',{@make_subplots,Base,'CHMr'});
    menuy_CACr = uimenu(mainy,'Label','Correlated-AC r','Checked', 'off','Callback',{@make_subplots,Base,'CACr'});
    menuy_CMeanr = uimenu(mainy,'Label','Correlated-Mean CAC r','Checked', 'off','Callback',{@make_subplots,Base,'CMeanr'});
    
    menuy_CHMslope = uimenu(mainy,'Label','Correlated-HM slope','Checked', 'off','Callback',{@make_subplots,Base,'CHMslope'});
    menuy_CACslope = uimenu(mainy,'Label','Correlated-AC slope','Checked', 'off','Callback',{@make_subplots,Base,'CACslope'});
    menuy_CMeanslope = uimenu(mainy,'Label','Correlated-Mean CAC slope','Checked', 'off','Callback',{@make_subplots,Base,'CMeanslope'});
    

    menuPenetrations = uimenu(mh,'Label','Penetrations','Callback',{@plot_penetrations,Base});
    menuRank = uimenu(mh,'Label','Rank','Callback',{@PlotRank});
    menuSmile = uimenu(mh,'Label','Smile','Callback',{@smile});
    
    
    count = 0;
    for dd = [1,ddmax];
        count = count+1;
            
        currentData = Base(dd).Cells;
            
            
        subplot(1,2,count); 
        setappdata(gca,'subplot',count);
        setappdata(gca,'density',Base(dd).density);
        
        % Set the x and y axis data
        setappdata(myFig,'xdata','CHMr');
        setappdata(myFig,'ydata','CACr');
        
        make_subplots(myFig,NaN,Base,'Correlation');
    end
    set(gcf,'Color','White','windowbuttondownfcn',{@TC_callback},'windowkeypressfcn',{@keydown_callback});

end

function make_subplots(myFig,evt,Base,type)
    
    count = 0;
    ddmax = length(Base);
    ch = get(gcf,'Children');
    
    % Need to redo this thing...
    plotMenu = get(ch(1),'Children');
    xmainMenu = plotMenu(5);
    ymainMenu = plotMenu(4);
    xmainChildren = get(xmainMenu,'Children');
    ymainChildren = get(ymainMenu,'Children');
    
    axchange = '';
    
    if ~isempty(events(evt));
        rootMenu = evt.Source.Parent.Label;
        if strcmp(rootMenu,'Main: x axis');

            setappdata(gcf,'xdata',type)
            axchange = 'x';

        elseif strcmp(rootMenu,'Main: y axis');
            setappdata(gcf,'ydata',type);
            axchange = 'y';
        else
            error('Error: This isn''t supposed to happen');
        end
    end
    
    
    for dd = [1,ddmax];
        count = count+1;
            
        currentData = Base(dd).Cells;
            
            
        subplot(1,2,count); 
        setappdata(gca,'subplot',count);
        setappdata(gca,'density',Base(dd).density);
        
        plot_data(currentData);
        
        % Adjust the ticks appropriately
        if strcmp(axchange,'x');
            mainChildren = xmainChildren;
        elseif strcmp(axchange,'y');
            mainChildren = ymainChildren;
        end
        
        if ~strcmp(axchange,'');
            for k = 1:length(mainChildren);
                current = mainChildren(k);

                if strcmp(current.Label,evt.Source.Label);
                    current.Checked = 'on';
                else
                    current.Checked = 'off';
                end
            end
        end        
           
    end

end





function plot_data(currentData);
    monkeys = {'jbe','lem'};
    figAppdata = getappdata(gcf);
    axAppdata = getappdata(gca);
    
    % If colours are not defined, define them. Otherwise just get them
    % from appdata
    if ~isfield(figAppdata,'colors');
        myCols = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;rand(10,3).^2];
        setappdata(gcf,'colors',myCols);
    else
        myCols = figAppdata.colors;
    end
    
    nCells = length(currentData);
    
    fileNames = {};
    cla;
    hold on;
    X = zeros(1,nCells);
    Y = zeros(1,nCells);
    
    xtype = getappdata(gcf,'xdata');
    ytype = getappdata(gcf,'ydata');
    
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
            
            % We choose from a range of data types to plot
            switch xtype
                case 'CHMr'
                    x = currentData(cell).regHm(1);
                    
                case 'CACr'
                    x = currentData(cell).regAc(1);
                    
                case 'CMeanr'
                    [x,~,~] = regression(currentData(cell).correlatedResponse, ...
                            (currentData(cell).correlatedResponse+ ...
                            currentData(cell).anticorrelatedResponse)/2);
                    
                case 'CHMslope'
                    x = currentData(cell).regHm(2);
                    
                case 'CACslope'
                    x = currentData(cell).regAc(2);
                    
                case 'CMeanslope'
                    [~,x,~] = regression(currentData(cell).correlatedResponse, ...
                            (currentData(cell).correlatedResponse+ ...
                            currentData(cell).anticorrelatedResponse)/2);
            end
            
            switch ytype
                case 'CHMr'
                    y = currentData(cell).regHm(1);
                    
                case 'CACr'
                    y = currentData(cell).regAc(1);
                    
                case 'CMeanr'
                    [y,~,~] = regression(currentData(cell).correlatedResponse, ...
                            (currentData(cell).correlatedResponse+ ...
                            currentData(cell).anticorrelatedResponse)/2);
                    
                case 'CHMslope'
                    y = currentData(cell).regHm(2);
                    
                case 'CACslope'
                    y = currentData(cell).regAc(2);
                    
                case 'CMeanslope'
                    [~,y,~] = regression(currentData(cell).correlatedResponse, ...
                            (currentData(cell).correlatedResponse+ ...
                            currentData(cell).anticorrelatedResponse)/2);
                    
            end
                
            
            % Different shapes for jbe and lem
            if strcmp(currentMonkey,monkeys{1});
                plot(x,y,'k','marker','s','markersize',...
                    round(1+currentData(cell).DDI*20),'markerfacecolor',currentColor);
            elseif strcmp(currentMonkey,monkeys{2});
                plot(x,y,'k','marker','o','markersize', ...
                    round(1+currentData(cell).DDI*20),'markerfacecolor',currentColor);
            end
            
            X(cell) = x;
            Y(cell) = y;


        end
    end
    
    switch xtype
        case 'CHMr'
            xlab = 'Correlated-HM r';
            xlims = [-1,1];

        case 'CACr'
            xlab = 'Correlated-AC r';
            xlims = [-1,1];            
            
        case 'CMeanr'
            xlab = 'Correlated-Mean CAC r';
            xlims = [-0.2,1];
            
        case 'CHMslope'
            xlab = 'Correlated-HM slope';
            xlims = [-0.35,0.35];
            
        case 'CACslope'
            xlab = 'Correlated-AC slope';
            xlims = [-1,1];
           
        case 'CMeanslope'
            xlab = 'Correlated-Mean CAC slope';
            xlims = [-1,1];
    end
    
    switch ytype
        case 'CHMr'
            ylab = 'Correlated-HM r';
            ylims = [-1,1];

        case 'CACr'
            ylab = 'Correlated-AC r';
            ylims = [-1,1];            
            
        case 'CMeanr'
            ylab = 'Correlated-Mean CAC r';
            ylims = [-0.2,1];
            
        case 'CHMslope'
            ylab = 'Correlated-HM slope';
            ylims = [-0.35,0.35];
            
        case 'CACslope'
            ylab = 'Correlated-AC slope';
            ylims = [-1,1];
            
        case 'CMeanslope'
            ylab = 'Correlated-Mean CAC slope';
            ylims = [-1,1];
    end
    
        
    XYData.x = X;
    XYData.y = Y;
    setappdata(gca,'XYData',XYData);
    setappdata(gca,'Cells',currentData);

    xlabel(xlab,'fontsize',20)
    ylabel(ylab, 'fontsize', 20);
    ylim(ylims);
    xlim(xlims);
    %set(gca,'XTick',-1:0.5:1,'fontsize',18);
    set(gca,'YTick',-1:0.5:1,'fontsize',18);

    if isfield(axAppdata,'density');
        title(sprintf('Density: %i%', axAppdata.density),'fontsize', 24);
    else
        title('Density unknown','fontsize',18);
    end
    

    % Put some text in the window showing how many cells
    jbeCount = 0;
    lemCount = 0;
    for j = 1:nCells
        jbeCount = jbeCount + strcmp(currentData(j).filename(9:11),'jbe');
        lemCount = lemCount + strcmp(currentData(j).filename(9:11),'lem');
    end
    text(0.4,0.8,['jbe: ',num2str(jbeCount),' cells']);
    text(0.4,0.7,['lem: ',num2str(lemCount),' cells']);
    

    % Draw some lines to split into quadrants
    plot([0,0],[-1,1],'linewidth',1,'color','red','linestyle','--');
    plot([-1,1],[0,0],'linewidth',1,'color','red','linestyle','--');
    
    plot([-1,1],[-1,1],'k -');

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
    if isfield(figData,'TCHandle');
        tcFigDeleted = ~ishandle(figData.TCHandle);
    end
    
    if ~isfield(figData,'TCHandle') || tcFigDeleted;
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
    
    xlabel('Disparity (deg)','fontsize',18);
    ylabel('Spikes per second','fontsize',18);
    density = (subplotNumber==1)*5 + (subplotNumber==2)*24;
    title(['Density: ', num2str(density),'; ',cellName],'fontsize',22);

        

    selectionType = get(myFig,'selectionType');
    shiftPressed = strcmp(selectionType,'extend');
    ctrlPressed = strcmp(selectionType,'alt');
    
    % This will only get executed if uniqueID is not passed to TC_callback
    % (which we only do when we call TC_callback recursively)

    
    if shiftPressed || ctrlPressed

        % If this hasn't been called before, assign a new figure handle
        if isfield(figData,'scatterHandle');
            figDeleted = ~ishandle(figData.scatterHandle);
        end
        
        if ~isfield(figData,'scatterHandle') || figDeleted;
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
                xlabel('Correlated response','fontsize',18);
                ylabel('Halfmatched response','fontsize',18);
                b = Cells(whichCell).regHm(2); m = Cells(whichCell).regHm(3);
                corrMin = min(Cells(whichCell).correlatedResponse); corrMax = max(Cells(whichCell).correlatedResponse);
                delta = corrMax*0.25; corrMin = corrMin-delta; corrMax = corrMax+delta;
                line([corrMin,corrMax]*2.5,[m+corrMin*b,m+corrMax*b]*2.5,'linewidth',2,'linestyle','-','color','r');
                title(['r= ', num2str(Cells(whichCell).regHm(1))])
                
            else
                plot(Cells(whichCell).correlatedResponse*2.5,Cells(whichCell).anticorrelatedResponse*2.5,'o k',...
                    'markersize',4,'markerfacecolor','r');
                xlabel('Correlated response','fontsize',18);
                ylabel('Anticorrelated response','fontsize',18);
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
    figure(myFig);
    
end

function keydown_callback(myFig,evt)
    
    keydown = evt.Key;
    Base = getappdata(myFig,'Base');
    currentType = getappdata(myFig,'type');
    
    if strcmp(keydown,'f');
        count = 0;
        for dd = [1,4];
            count = count+1;
            
            currentData = Base(dd).Cells;
            
            subplot(1,2,count); 
            setappdata(gca,'subplot',count);
            setappdata(gca,'density',Base(dd).density);
            
            % Check what's currently plotted and plot whatever is not
            % currently on the screen
            if strcmp(currentType,'slope');
                plot_ACHM(currentData)
            elseif strcmp(currentType,'ACHM');
                plot_slopes(currentData);
            end
        end
        
    end

end

function myHash = quickhash(string)
    string = double(string);
    myHash = mod(sum(string(:).^2.3),2^32-1);
end


function plot_penetrations(~,~,Base);
    figure(); hold on;
    set(gcf,'color','white');
    SubBase = Base(length(Base));
    XYs = SubBase.penetrationlist;
    exptlist = SubBase.exptlist;
    
    
    for j = 1:length(exptlist);
        fname = exptlist{j};
        monkey = fname(9:11);
        session = fname(18:24);
        
        % Find how many cells come from that session:
        count = 0;
        
        for c = 1:length(SubBase.Cells);
            currentFilename = SubBase.Cells(c).filename;
            count = count + strcmp(currentFilename(18:24),session);
        end
        
        if strcmp(monkey,'jbe'); % colour coded by monkey
            col = 'red';
        else
            col = 'blue';
        end
        
        
        XY = XYs(:,j);
        plot(XY(1),XY(2),'k o','markerfacecolor',col);
        
        text(XY(1)-0.1,XY(2)+0.1,[session, ' (',num2str(count),')']);
        
    end
   
end