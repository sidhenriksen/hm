function PlotCurated(varargin)

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
    myCols = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;rand(15,3)];
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
    menux_HMCACslope = uimenu(mainx,'Label','Mean CAC slope/HM slope','Checked', 'off','Callback',{@make_subplots,Base,'HMCACslope'});
    
    
    % Y menu
    menuy_CHMr = uimenu(mainy,'Label','Correlated-HM r','Checked', 'on','Callback',{@make_subplots,Base,'CHMr'});
    menuy_CACr = uimenu(mainy,'Label','Correlated-AC r','Checked', 'off','Callback',{@make_subplots,Base,'CACr'});
    menuy_CMeanr = uimenu(mainy,'Label','Correlated-Mean CAC r','Checked', 'off','Callback',{@make_subplots,Base,'CMeanr'});
    
    menuy_CHMslope = uimenu(mainy,'Label','Correlated-HM slope','Checked', 'off','Callback',{@make_subplots,Base,'CHMslope'});
    menuy_CACslope = uimenu(mainy,'Label','Correlated-AC slope','Checked', 'off','Callback',{@make_subplots,Base,'CACslope'});
    menuy_CMeanslope = uimenu(mainy,'Label','Correlated-Mean CAC slope','Checked', 'off','Callback',{@make_subplots,Base,'CMeanslope'});
    menuy_HMCACslope = uimenu(mainy,'Label','Mean CAC slope/HM slope','Checked', 'off','Callback',{@make_subplots,Base,'HMCACslope'});
    

    menuStats = uimenu(mh,'Label','Stats','Callback',{@show_stats,Base});
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
    xmainMenu = plotMenu(length(plotMenu));
    ymainMenu = plotMenu(length(plotMenu)-1);
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
    allExptNames = figAppdata.Base(1).exptlist;
    
    % If colours are not defined, define them. Otherwise just get them
    % from appdata
    if ~isfield(figAppdata,'colors');
        myCols = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;rand(20,3).^2];
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
    
    lambda = 1/sqrt(10);
    
    for cell = 1:nCells
        
        whichFile = strcmp(currentData(cell).filename,allExptNames);
        currentColor = myCols(whichFile,:);
        
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
                    corrResp = currentData(cell).correlatedResponse';
                    meanCACResp = (currentData(cell).correlatedResponse+ ...
                            currentData(cell).anticorrelatedResponse)'/2;
                    [x,~,~] = regression2(meanCACResp,corrResp);
                    %[~,x,~] = fit_bothsubj2error(corrResp,meanCACResp,lambda);
                    
                case 'CHMslope'
                    x = currentData(cell).regHm(2);
                    
                case 'CACslope'
                    x = currentData(cell).regAc(2);
                    
                case 'CMeanslope'
                    corrResp = currentData(cell).correlatedResponse';
                    meanCACResp = (currentData(cell).correlatedResponse+ ...
                            currentData(cell).anticorrelatedResponse)'/2;
                    [~,x,~] = regression2(meanCACResp,corrResp);
                    
                case 'HMCACslope'
                    corrResp = currentData(cell).correlatedResponse';
                    hmResp = currentData(cell).halfmatchedResponse';
                    meanCACResp = (currentData(cell).correlatedResponse+ ...
                            currentData(cell).anticorrelatedResponse)'/2;
                    [~,hmslope,~] = regression2(hmResp,corrResp); % Change this to fit_bothsubj2error?
                    [~,meanCACslope,~] = regression2(meanCACResp,corrResp);
                    x = hmslope/meanCACslope;
            end
            
            switch ytype
                case 'CHMr'
                    y = currentData(cell).regHm(1);
                    
                case 'CACr'
                    y = currentData(cell).regAc(1);
                    
                case 'CMeanr'
                    corrResp = currentData(cell).correlatedResponse';
                    meanCACResp = (currentData(cell).correlatedResponse+ ...
                            currentData(cell).anticorrelatedResponse)'/2;
                    [y,~,~] = regression2(meanCACResp,corrResp);
                    
                case 'CHMslope'
                    y = currentData(cell).regHm(2);
                    
                case 'CACslope'
                    y = currentData(cell).regAc(2);
                    
                case 'CMeanslope'
                    corrResp = currentData(cell).correlatedResponse';
                    meanCACResp = (currentData(cell).correlatedResponse+ ...
                            currentData(cell).anticorrelatedResponse)'/2;
                    [~,y,~] = regression2(meanCACResp,corrResp);
                    
                case 'HMCACslope'
                    corrResp = currentData(cell).correlatedResponse';
                    hmResp = currentData(cell).halfmatchedResponse';
                    meanCACResp = (currentData(cell).correlatedResponse+ ...
                            currentData(cell).anticorrelatedResponse)'/2;
                    [~,hmslope,~] = regression2(hmResp,corrResp); % Change this to fit_bothsubj2error?
                    [~,meanCACslope,~] = regression2(meanCACResp,corrResp);
                    y = hmslope/meanCACslope; % Change this to fit_bothsubj2error?
                    
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
            
        case  'HMCACslope'
            xlab = 'HM slope/mean CAC slope';
            xlims = [-0.25,1.25];
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
            
        case  'HMCACslope'
            ylab = 'HM slope/mean CAC slope';
            ylims = [-0.25,1.25];
    end
    
        
    XYData.x = X;
    XYData.y = Y;
    setappdata(gca,'XYData',XYData);
    setappdata(gca,'Cells',currentData);

    xlabel(xlab,'fontsize',20)
    ylabel(ylab, 'fontsize', 20);
    ylim(ylims);
    xlim(xlims);
    set(gca,'XTick',-1:0.25:1,'fontsize',18);
    set(gca,'YTick',-1:0.25:1,'fontsize',18);

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
        newFig = figure(); setappdata(newFig,'plotCAC','Off');
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
    
    %cacMenu = uimenu(myFig,'Label','Plot CAC','Checked','Off','Callback',{@plot_CAC_in_TC_plot,);
    
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

function show_stats(myFig,evt,Base);
    regHm5 = cat(1,Base(1).Cells.regHm);
    regHm24 = cat(1,Base(length(Base)).Cells.regHm);
    
    [Hr5,Pr5,CIr5,Statsr5] = ttest(regHm5(:,1),0);
    [Hm5,Pm5,CIm5,Statsm5] = ttest(regHm5(:,2),0);
    
    [Hr24,Pr24,CIr24,Statsr24] = ttest(regHm24(:,1),0);
    [Hm24,Pm24,CIm24,Statsm24] = ttest(regHm24(:,2),0);
    
    figure();
    xlim([0,1]);
    ylim([0,1]);
    set(gca,'visible','off')
    x = -0.1;
    
    % 5% density
    fs = 13;
    text(x,0.925,sprintf(['one-sample t-tests testing whether slopes and \ncorrelation coefficients', ...
        ' are drawn from sample with mean 0']));
    text(x,0.7,'5% density','fontsize',16);
    text(x,0.6,sprintf('r: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsr5.df,Statsr5.tstat,Pr5,CIr5(1),CIr5(2)),'fontsize',fs)
    text(x,0.525,sprintf('slope: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm5.df,Statsm5.tstat,Pm5,CIm5(1),CIm5(2)),'fontsize',fs);
    
    % 5% density
    text(x,0.3,'24% density','fontsize',16);
    text(x,0.2,sprintf('r: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsr24.df,Statsr24.tstat,Pr24,CIr24(1),CIr24(2)),'fontsize',fs)
    text(x,0.125,sprintf('slope: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm24.df,Statsm24.tstat,Pm24,CIm24(1),CIm24(2)),'fontsize',fs);
    set(gcf,'color','white','position',[100,200,550,400])
end