function PlotCurated(varargin)
    %
    % Author: Sid Henriksen (2015). Email: sid.henriksen@gmail.com.
    %
    % Note, the code for this browser is not well-documented, nor is it
    % intended to be. The README file in the repo documents the structure
    % of the CuratedCells.mat file, which this script uses. 
    %
    
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
        ddmax = 2;
    end
    

    densities = cat(1,Base.density);
    
    % This is included just for the sake of colors
    fileNames = {};
    myCols = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;rand(15,3)];
    monkeys = {'jbe','lem'};

    setappdata(myFig,'Base',Base);
    
    % This first segment defines the menus; all call make_subplots to
    % populate the subplots, together with Base (consisting of all the
    % data) and a string which specifies what we're plotting.
    
    % X menu
    menux_CHMr = uimenu(mainx,'Label','Half-matched r','Checked', 'on','Callback',{@make_subplots,Base,'CHMr'});
    menux_CACr = uimenu(mainx,'Label','Anticorrelated r','Checked', 'off','Callback',{@make_subplots,Base,'CACr'});
    menux_CMeanr = uimenu(mainx,'Label','Mean CAC r','Checked', 'off','Callback',{@make_subplots,Base,'CMeanr'});
    
    menux_CHMslope = uimenu(mainx,'Label','Half-matched slope','Checked', 'off','Callback',{@make_subplots,Base,'CHMslope'});
    menux_CACslope = uimenu(mainx,'Label','Anticorrelated slope','Checked', 'off','Callback',{@make_subplots,Base,'CACslope'});
    menux_CMeanslope = uimenu(mainx,'Label','Mean CAC slope','Checked', 'off','Callback',{@make_subplots,Base,'CMeanslope'});
    menux_HMCACslope = uimenu(mainx,'Label','Mean CAC slope/HM slope','Checked', 'off','Callback',{@make_subplots,Base,'HMCACslope'});
    
    %%% ROC menu
    menux_ROC = uimenu(mainx,'Label','ROC/d''');
        menux_HMauc = uimenu(menux_ROC,'Label','Half-matched AUC','Checked', 'off','Callback',{@make_subplots,Base,'HM AUC'});
        menux_HMaucZ = uimenu(menux_ROC,'Label','Half-matched AUC (Z-score)','Checked', 'off','Callback',{@make_subplots,Base,'HM AUC Z'});
        menux_HMdprime = uimenu(menux_ROC,'Label','Half-matched d''','Checked', 'off','Callback',{@make_subplots,Base,'HM d'''});
        menux_Cauc = uimenu(menux_ROC,'Label','Correlated AUC','Checked', 'off','Callback',{@make_subplots,Base,'C AUC'});
        menux_CaucZ = uimenu(menux_ROC,'Label','Correlated AUC (Z-score)','Checked', 'off','Callback',{@make_subplots,Base,'C AUC Z'});
        menux_Cdprime = uimenu(menux_ROC,'Label','Correlated d''','Checked', 'off','Callback',{@make_subplots,Base,'C d'''});
    

    % Y menu
    menuy_CHMr = uimenu(mainy,'Label','Half-matched r','Checked', 'on','Callback',{@make_subplots,Base,'CHMr'});
    menuy_CACr = uimenu(mainy,'Label','Anticorrelated r','Checked', 'off','Callback',{@make_subplots,Base,'CACr'});
    menuy_CMeanr = uimenu(mainy,'Label','Mean CAC r','Checked', 'off','Callback',{@make_subplots,Base,'CMeanr'});

    menuy_CHMslope = uimenu(mainy,'Label','Half-matched slope','Checked', 'off','Callback',{@make_subplots,Base,'CHMslope'});
    menuy_CACslope = uimenu(mainy,'Label','Anticorrelated slope','Checked', 'off','Callback',{@make_subplots,Base,'CACslope'});
    menuy_CMeanslope = uimenu(mainy,'Label','Mean CAC slope','Checked', 'off','Callback',{@make_subplots,Base,'CMeanslope'});
    menuy_HMCACslope = uimenu(mainy,'Label','Mean CAC slope/HM slope','Checked', 'off','Callback',{@make_subplots,Base,'HMCACslope'});

    %%% ROC menu
    menuy_ROC = uimenu(mainy,'Label','ROC/d''');
        menuy_HMauc = uimenu(menuy_ROC,'Label','Half-matched AUC','Checked', 'off','Callback',{@make_subplots,Base,'HM AUC'});
        menuy_HMaucZ = uimenu(menuy_ROC,'Label','Half-matched AUC (Z-score)','Checked', 'off','Callback',{@make_subplots,Base,'HM AUC Z'});
        menuy_HMdprime = uimenu(menuy_ROC,'Label','Half-matched d''','Checked', 'off','Callback',{@make_subplots,Base,'HM d'''});

        menuy_Cauc = uimenu(menuy_ROC,'Label','Correlated AUC','Checked', 'off','Callback',{@make_subplots,Base,'C AUC'});
        menuy_CaucZ = uimenu(menuy_ROC,'Label','Correlated AUC (Z-score)','Checked', 'off','Callback',{@make_subplots,Base,'C AUC Z'});
        menuy_Cdprime = uimenu(menuy_ROC,'Label','Correlated d''','Checked', 'off','Callback',{@make_subplots,Base,'C d'''});
    
    % Misc
    menuPlotDDIROC = uimenu(mh,'Label','DDI-ROC','Callback',{@make_auc_plot});
    menuStats = uimenu(mh,'Label','Stats','Callback',{@show_stats,Base});
    menuPenetrations = uimenu(mh,'Label','Penetrations','Callback',{@plot_penetrations,Base});
    menuRank = uimenu(mh,'Label','Rank','Callback',{@PlotRank});
    menuSmile = uimenu(mh,'Label','Smile','Callback',{@smile});
    
    
    
    % This section initializes the plot and figure data
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
        setappdata(myFig,'xy_errorbars',0);
        
        make_subplots(myFig,NaN,Base,'Correlation');
    end
    set(gcf,'Color','White','windowbuttondownfcn',{@TC_callback},'windowkeypressfcn',{@keydown_callback});

end

function make_subplots(myFig,evt,Base,type)
    % Function to populate the subplots of the main figure
    % This will take Base, which is the structure generated by
    % CurateCells.m, and a string as argument. The string should specify
    % what we're plotting on the relevant axis.
    count = 0;
    ddmax = length(Base);
    ch = get(gcf,'Children');
    
    % This doesn't work for older version of Matlab; I don't have easy
    % access to an older version so debugging becomes a bit tricky.
    plotMenu = get(ch(1),'Children');
    xmainMenu = plotMenu(length(plotMenu));
    ymainMenu = plotMenu(length(plotMenu)-1);
    xmainChildren = get(xmainMenu,'Children');
    ymainChildren = get(ymainMenu,'Children');
    
    axchange = '';
    
    if ~isempty(events(evt));
        rootMenu = evt.Source.Parent;
        
        
        if strcmp(rootMenu.Label,'Main: x axis');

            setappdata(gcf,'xdata',type)
            axchange = 'x';

        elseif strcmp(rootMenu.Label,'Main: y axis');
            setappdata(gcf,'ydata',type);
            axchange = 'y';
        else
            rootMenu = rootMenu.Parent;
                    
            if strcmp(rootMenu.Label,'Main: x axis');

                setappdata(gcf,'xdata',type)
                axchange = 'x';

            elseif strcmp(rootMenu.Label,'Main: y axis');
                setappdata(gcf,'ydata',type);
                axchange = 'y';
            end
            
            %error('Error: This isn''t supposed to happen');
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
    % This is a function that actually plots the data in the highlighted
    % subplot. This function will take the figure data that has been set
    % elsewhere and will plot the correct thing. This is really lengthy.
    
    monkeys = {'jbe','lem'};
    myFig = gcf;
    figAppdata = getappdata(myFig);
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
            
            % We choose from a range of data types to plot;
            % This is really long, but simply chooses the correct values
            % to plot on the x and y axes.
            switch xtype
                case 'CHMr'
                    x = currentData(cell).regHmRegular(1);
                    
                case 'CACr'
                    x = currentData(cell).regAcRegular(1);
                    
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
                    
                case 'HM AUC';
                    x = currentData(cell).HMauc;
                    
                case 'HM AUC Z';
                    x = currentData(cell).HMaucZ;
                    
                case 'HM d''';
                    x = currentData(cell).HMdprime;
                    
                case 'C AUC';
                    x = currentData(cell).Cauc;
                    
                case 'C AUC Z';
                    x = currentData(cell).CaucZ;
                    
                case 'C d''';
                    x = currentData(cell).Cdprime;
            end
            
            switch ytype
                case 'HM AUC';
                    y = currentData(cell).HMauc;
                    
                case 'HM d''';
                    y = currentData(cell).HMdprime;
                case 'CHMr'
                    y = currentData(cell).regHmRegular(1);
                    
                case 'CACr'
                    y = currentData(cell).regAcRegular(1);
                    
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
                    
                case 'HM AUC';
                    y = currentData(cell).HMauc;
                case 'HM AUC Z';
                    y = currentData(cell).HMaucZ;
                    
                case 'HM d''';
                    y = currentData(cell).HMdprime;
                    
                case 'C AUC';
                    y = currentData(cell).Cauc;
                case 'C AUC Z';
                    y = currentData(cell).CaucZ;
                    
                case 'C d''';
                    y = currentData(cell).Cdprime;
                    
            end
                
            
            % Get errorbar data
            xbar = []; ybar = [];
            errorbars_on=getappdata(myFig,'xy_errorbars');
            
            if errorbars_on
                switch xtype
                    case 'CHMr'
                        xbar = currentData(cell).CHm_r_CI;

                    case 'CACr'
                        xbar = currentData(cell).CAc_r_CI;

                    case 'CHMslope'
                        xbar = currentData(cell).CHm_slope_CI;

                    case 'CACslope'
                        xbar = currentData(cell).CAc_slope_CI;
                end


                switch ytype
                    case 'CHMr'
                        ybar = currentData(cell).CHm_r_CI;

                    case 'CACr'
                        ybar = currentData(cell).CAc_r_CI;

                    case 'CHMslope'
                        ybar = currentData(cell).CHm_slope_CI;

                    case 'CACslope'
                        ybar = currentData(cell).CAc_slope_CI;
                end
            end

            %% This is where we actually plot the data
            % Different shapes for jbe and lem
            if strcmp(currentMonkey,monkeys{1});
                plot(x,y,'k','marker','s','markersize',...
                    round(1+currentData(cell).DDI*20),'markerfacecolor',currentColor);
                
                if ~isempty(xbar);
                    xbar = abs(xbar-x);
                    E=herrorbar(x,y,xbar(1),xbar(2));
                    set(E,'linewidth',2,'color',currentColor);
                end
                
                if ~isempty(ybar);
                    ybar = abs(ybar-y);
                   E = errorbar(x,y,ybar(1),ybar(2));
                   set(E,'linewidth',2,'color',currentColor);
                end
                        
                    
            elseif strcmp(currentMonkey,monkeys{2});
                plot(x,y,'k','marker','o','markersize', ...
                    round(1+currentData(cell).DDI*20),'markerfacecolor',currentColor);
                
                
                if ~isempty(xbar);
                    xbar = abs(xbar-x);
                    E=herrorbar(x,y,xbar(1),xbar(2));
                    set(E,'linewidth',2,'color',currentColor);
                end
                
                if ~isempty(ybar);
                    ybar = abs(ybar-y);
                   E = errorbar(x,y,ybar(1),ybar(2));
                   set(E,'linewidth',2,'color',currentColor);
                end
            end
            
            X(cell) = x;
            Y(cell) = y;


        end
    end
    
    % Again, very lengty, but simply chooses appropriate xlims/ylims and 
    % axis labels
    switch xtype
        case 'CHMr'
            xlab = 'HM r';
            xlims = [-1,1];

        case 'CACr'
            xlab = 'AC r';
            xlims = [-1,1];            
            
        case 'CMeanr'
            xlab = 'Mean CAC r';
            xlims = [-0.2,1];
            
        case 'CHMslope'
            xlab = 'HM slope';
            xlims = [-0.35,0.35];
            
        case 'CACslope'
            xlab = 'AC slope';
            xlims = [-1,1];
           
        case 'CMeanslope'
            xlab = 'Mean CAC slope';
            xlims = [-1,1];
            
        case  'HMCACslope'
            xlab = 'HM slope/mean CAC slope';
            xlims = [-0.25,1.25];
            
        case 'HM AUC'
            xlab = 'Half-matched AUC';
            xlims = [0,1];
            
        case 'HM AUC Z'
            xlab = 'Half-matched AUC (Z-score)';
            xlims = [0,1];
            
        case 'HM d'''
            xlab = 'Half-matched d''';
            xlims = [-0.5,5];
            
        case 'C AUC'
            xlab = 'Correlated AUC';
            xlims = [0,1];
            
        case 'C AUC Z'
            xlab = 'Correlated AUC (Z-score)';
            xlims = [0,1];
            
        case 'C d'''
            xlab = 'Correlated d''';
            xlims = [-0.5,10];
            
           
    end
    
    switch ytype
        case 'CHMr'
            ylab = 'HM r';
            ylims = [-1,1];

        case 'CACr'
            ylab = 'Anticorrelated r';
            ylims = [-1,1];            
            
        case 'CMeanr'
            ylab = 'Mean CAC r';
            ylims = [-0.2,1];
            
        case 'CHMslope'
            ylab = 'HM slope';
            ylims = [-0.35,0.35];
            
        case 'CACslope'
            ylab = 'AC slope';
            ylims = [-1,1];
            
        case 'CMeanslope' %for half-matched and anticorrelated
            ylab = 'Mean CAC slope';
            ylims = [-1,1];
            
        case  'HMCACslope'
            ylab = 'HM slope/mean CAC slope';
            ylims = [-0.25,1.25];
            
        case 'HM AUC'
            ylab = 'Half-matched AUC';
            ylims = [0,1];
            
        case 'HM AUC Z'
            ylab = 'Half-matched AUC (Z-score)';
            ylims = [0,1];
            
        case 'HM d'''
            ylab = 'Half-matched d''';
            ylims = [-0.5,5];
                
        case 'C AUC'
            ylab = 'Correlated AUC';
            ylims = [0,1];
            
        case 'C AUC Z'
            ylab = 'Correlated AUC (Z-score)';
            ylims = [0,1];
            
        case 'C d'''
            ylab = 'Correlated d''';
            ylims = [-0.5,10];
    end
    
    % Select appropriate ticks for purteh lookin' plot
    if ~isempty(strfind(xtype,'d'''));
        xticks = -1:2:10;
    else
        xticks = -1:0.5:1;
    end
    
    if ~isempty(strfind(ytype,'d'''));
        yticks = -1:2:10;
    else
        yticks = -1:0.5:1;
    end
    
    % Set XYData!
    XYData.x = X;
    XYData.y = Y;
    setappdata(gca,'XYData',XYData);
    setappdata(gca,'Cells',currentData);

    xlabel(xlab,'fontsize',20)
    ylabel(ylab, 'fontsize', 20);
    ylim(ylims);
    xlim(xlims);
    
    set(gca,'XTick',xticks,'fontsize',18);
    set(gca,'YTick',yticks,'fontsize',18);

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
    

    % Draw some lines to split into quadrants; where to draw this depends
    % on what we're plotting
    if ~isempty(strfind(ytype,'AUC'));
        plot([-20,20],[0.5,0.5],'linewidth',1,'color','red','linestyle','--');
       
    else
        plot([-20,20],[0,0],'linewidth',1,'color','red','linestyle','--'); 
    end
    
    if ~isempty(strfind(xtype,'AUC'));
        plot([0.5,0.5],[-20,20],'linewidth',1,'color','red','linestyle','--');
    else
        plot([0,0],[-20,20],'linewidth',1,'color','red','linestyle','--');
    end
    
    

end


function TC_callback(myFig,evt,varargin);
    % Callback function to plot tuning curves
    
    figData = getappdata(myFig);
    % If this hasn't been called before, assign a new figure handle
    if isfield(figData,'TCHandle');
        tcFigDeleted = ~ishandle(figData.TCHandle);
    end
    
    if ~isfield(figData,'TCHandle') || tcFigDeleted;
        newFig = figure(); setappdata(newFig,'plotCAC','Off');
        set(newFig,'color','white');
        figData.TCHandle = newFig;
        setappdata(myFig,'TCHandle',newFig);
        setappdata(myFig,'Errorbars','Off');
        setappdata(figData.TCHandle,'PseudoParent',myFig);
        toggleMenu = uimenu(figData.TCHandle,'Label','Toggle');
        toggleErrorbarsSEM = uimenu(toggleMenu,'Label','SEM error bars (95% CIs)','Checked','off','Callback',{@set_errorbars,myFig,evt,'SEM'});
        toggleErrorbarsSD = uimenu(toggleMenu,'Label','Error bars (+/- 1 SD)','Checked','off','Callback',{@set_errorbars,myFig,evt,'SD'});
        toggleErrorbarsBootstrap = uimenu(toggleMenu,'Label','Bootstrap error bars (95% CIs)','Checked','off','Callback',{@set_errorbars,myFig,evt,'bootstrap'});
        setappdata(figData.TCHandle,'ErrorbarsToggled','None');
    end
    
    errorbarsToggled = getappdata(figData.TCHandle,'ErrorbarsToggled');
    currentAxis = myFig.CurrentAxes;
    axisData = getappdata(currentAxis);
    

    currentPoint = get(currentAxis,'currentpoint');    

    
            
        
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
    

    

    
    
    figure(figData.TCHandle);
    
    if subplotNumber > 2;
        maxSubplots = subplotNumber;
    else
        maxSubplots = 2;
    end
    
    
    subplot(1,maxSubplots,subplotNumber); hold off;
    if whichCell > 0;
        cellName = [Cells(whichCell).filename(18:24),'-cell',num2str(Cells(whichCell).cellnumber)];
        
        set(gca,'fontsize',18);
        Dx = Cells(whichCell).Dx;
        menus = get(figData.TCHandle,'Children');
        
        toggleMenu = menus(1); errorbarMenu = toggleMenu.Children;
        
        if ~any(strcmpi({errorbarMenu.Checked},'On'));
            plot(Dx,Cells(whichCell).halfmatchedResponse*2.5,'b -- o', ...
            Dx,Cells(whichCell).correlatedResponse*2.5,'r -- o', ...
            Dx,Cells(whichCell).anticorrelatedResponse*2.5,'k -- o', ...
            'linewidth',3,'markersize',5,'markerfacecolor','k');
        else
            
            barType = getappdata(figData.TCHandle,'ErrorbarsToggled');
            
            switch barType
                case 'SEM'
                    errorbarHM1 = Cells(whichCell).halfmatchedSEM*1.96;
                    errorbarHM2 = Cells(whichCell).halfmatchedSEM*1.96;
                    
                    errorbarC1 = Cells(whichCell).correlatedSEM*1.96;
                    errorbarC2 = Cells(whichCell).correlatedSEM*1.96;
                    
                    errorbarAC1 = Cells(whichCell).anticorrelatedSEM*1.96;
                    errorbarAC2 = Cells(whichCell).anticorrelatedSEM*1.96;
                    
                case 'SD'
                    errorbarHM1 = Cells(whichCell).halfmatchedSD;
                    errorbarHM2 = Cells(whichCell).halfmatchedSD;
                    
                    errorbarC1 = Cells(whichCell).correlatedSD;
                    errorbarC2 = Cells(whichCell).correlatedSD;
                    
                    errorbarAC1 = Cells(whichCell).anticorrelatedSD;
                    errorbarAC2 = Cells(whichCell).anticorrelatedSD;
                    
                case 'bootstrap'
                    %errorbarC1 = Cells(whichCell).ciLowC;
                    %errorbarC2 = Cells(whichCell).ciHighC;
                    
                    %errorbarHM1 = Cells(whichCell).ciLowHm;
                    %errorbarHM2 = Cells(whichCell).ciHighHm;
                    
                    %errorbarAC1 = Cells(whichCell).ciLowAc;
                    %errorbarAC2 = Cells(whichCell).ciHighAc;
                    
                    errorbarC1 = Cells(whichCell).correlatedResponse-Cells(whichCell).ciLowC;
                    errorbarC2 = Cells(whichCell).ciHighC - Cells(whichCell).correlatedResponse;
                    errorbarHM1 = Cells(whichCell).halfmatchedResponse-Cells(whichCell).ciLowHm;
                    errorbarHM2 = Cells(whichCell).ciHighHm - Cells(whichCell).halfmatchedResponse;
                    errorbarAC1 = Cells(whichCell).anticorrelatedResponse-Cells(whichCell).ciLowAc;
                    errorbarAC2 = Cells(whichCell).ciHighAc - Cells(whichCell).anticorrelatedResponse;
            end
            
            E1=errorbar(Dx,Cells(whichCell).halfmatchedResponse*2.5,...
                errorbarHM1*2.5,errorbarHM2*2.5);
            set(E1,'linewidth',3,'color','b','marker','o','linestyle','--',...
            'markersize',5,'markerfacecolor','b');
            hold on;
            E2=errorbar(Dx,Cells(whichCell).correlatedResponse*2.5,...
                errorbarC1*2.5,errorbarC2*2.5);
            set(E2,'linewidth',3,'color','r','marker','o','linestyle','--',...
            'markersize',5,'markerfacecolor','r')
            
            E3=errorbar(Dx,Cells(whichCell).anticorrelatedResponse*2.5,...
                errorbarAC1*2.5,errorbarAC2*2.5);
            set(E3,'linewidth',3,'color','k','marker','o','linestyle','--',...
            'markersize',5,'markerfacecolor','k')
            hold off;
        end
        
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
                %b = Cells(whichCell).regHm(2); m = Cells(whichCell).regHm(3);
                [r,b,m] = regression2(Cells(whichCell).halfmatchedResponse*2.5,Cells(whichCell).correlatedResponse*2.5);
                corrMin = min(Cells(whichCell).correlatedResponse)*2.5; corrMax = max(Cells(whichCell).correlatedResponse)*2.5;
                delta = corrMax*0.25; corrMin = corrMin-delta; corrMax = corrMax+delta;
                line([corrMin,corrMax],[m+corrMin*b,m+corrMax*b],'linewidth',2,'linestyle','-','color','r');
                title(['r= ', num2str(r)])
                
            else
                plot(Cells(whichCell).correlatedResponse*2.5,Cells(whichCell).anticorrelatedResponse*2.5,'o k',...
                    'markersize',4,'markerfacecolor','r');
                xlabel('Correlated response','fontsize',18);
                ylabel('Anticorrelated response','fontsize',18);
                
                [r,b,m] = regression2(Cells(whichCell).anticorrelatedResponse*2.5,Cells(whichCell).correlatedResponse*2.5);
                corrMin = min(Cells(whichCell).correlatedResponse)*2.5; corrMax = max(Cells(whichCell).correlatedResponse)*2.5;
                delta = corrMax*0.25; corrMin = corrMin-delta; corrMax = corrMax+delta;
                line([corrMin,corrMax],[m+corrMin*b,m+corrMax*b],'linewidth',2,'linestyle','-','color','k');
                title(['r= ', num2str(r)])
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
       
    elseif strcmp(keydown,'e');
        errorbars_on = getappdata(myFig,'xy_errorbars');
        setappdata(myFig,'xy_errorbars',~errorbars_on);
        make_subplots(myFig,NaN,Base,'Correlation');
    end

end

function myHash = quickhash(string);
    string = double(string);
    myHash = mod(sum(string(:).^2.3),2^32-1);
end


function plot_penetrations(~,~,Base);
    figure(); hold on;
    set(gcf,'color','white');
    SubBase = Base(1);
    XYs = SubBase.penetrationlist;
    exptlist = SubBase.exptlist;
    
    
    for j = 1:length(exptlist);
        fname = exptlist{j};
        monkey = fname(9:11);
        session = fname(18:24);
        
        % Find how many cells come from that session:
        count = 0;
        dws = [];
        for c = 1:length(SubBase.Cells);
            currentFilename = SubBase.Cells(c).filename;
            
            
            if strcmp(currentFilename(18:24),session);
                count = count + 1;
                dws = [dws,SubBase.Cells(c).dw];
            end
        end
        dws(dws==1) = 0.15;
        if strcmp(monkey,'jbe'); % colour coded by monkey
            col = 'red';
        else
            col = 'blue';
        end
        
        
        XY = XYs(:,j);
        plot(XY(1),XY(2),'k o','markerfacecolor',col,'markersize',5*exp(mean(dws)*4));
        
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
    
    % Do this for each monkey separately
    lem_idx = zeros(1,length(Base(1).Cells));
    for j = 1:length(Base(1).Cells);
        currentCell = Base(1).Cells(j);
        lem_idx(j) = any(strfind(currentCell.filename,'lem'));
    end
    lem_idx = logical(lem_idx);
    jbe_idx = ~lem_idx;    
    
    [~,Pm5_lem,CIm5_lem,Statsm5_lem] = ttest(regHm5(lem_idx,2),0);    
    [~,Pm24_lem,CIm24_lem,Statsm24_lem] = ttest(regHm24(lem_idx,2),0);
    
    [~,Pm5_jbe,CIm5_jbe,Statsm5_jbe] = ttest(regHm5(jbe_idx,2),0);
    [~,Pm24_jbe,CIm24_jbe,Statsm24_jbe] = ttest(regHm24(jbe_idx,2),0);
    
    [~,Pm5v24,CIm5v24,Statsm5v24] = ttest(regHm5(:,2),regHm24(:,2));
    [~,Pm5v24_jbe,CIm5v24_jbe,Statsm5v24_jbe] = ttest(regHm5(jbe_idx,2),regHm24(jbe_idx,2));
    [~,Pm5v24_lem,CIm5v24_lem,Statsm5v24_lem] = ttest(regHm5(lem_idx,2),regHm24(lem_idx,2));
    
    
    
    figure();
    xlim([0,1]);
    ylim([0,1]);
    set(gca,'visible','off')
    x = -0.1;
    
    % 5% density
    fs = 13;
    text(x,1.0,sprintf(['one-sample t-tests testing whether slopes and \ncorrelation coefficients', ...
        ' are drawn from sample with mean 0']));
    
    y1 = 0.9;
    dy=0.04;
    text(x,y1,'5% density','fontsize',16);
    text(x,y1-dy,sprintf('r: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsr5.df,Statsr5.tstat,Pr5,CIr5(1),CIr5(2)),'fontsize',fs)
    text(x,y1-dy*2,sprintf('slope: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm5.df,Statsm5.tstat,Pm5,CIm5(1),CIm5(2)),'fontsize',fs);
    
    y2 = y1-dy*3.5;
    text(x,y2,'By monkey:','fontsize',14);
    text(x,y2-dy,sprintf('lem slope: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm5_lem.df,Statsm5_lem.tstat,Pm5_lem,CIm5_lem(1),CIm5_lem(2)),'fontsize',fs);
    text(x,y2-dy*2,sprintf('jbe slope: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm5_jbe.df,Statsm5_jbe.tstat,Pm5_jbe,CIm5_jbe(1),CIm5_jbe(2)),'fontsize',fs);
    
    
    y3 = y2-dy*5;
    % 5% density
    text(x,y3,'24% density','fontsize',16);
    text(x,y3-dy,sprintf('r: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsr24.df,Statsr24.tstat,Pr24,CIr24(1),CIr24(2)),'fontsize',fs)
    text(x,y3-dy*2,sprintf('slope: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm24.df,Statsm24.tstat,Pm24,CIm24(1),CIm24(2)),'fontsize',fs);
    
    y4 = y3-dy*3.5;
    text(x,y4,'By monkey:','fontsize',14);
    text(x,y4-dy,sprintf('lem slope: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm24_lem.df,Statsm24_lem.tstat,Pm24_lem,CIm24_lem(1),CIm24_lem(2)),'fontsize',fs);
    text(x,y4-dy*2,sprintf('jbe slope: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm24_jbe.df,Statsm24_jbe.tstat,Pm24_jbe,CIm24_jbe(1),CIm24_jbe(2)),'fontsize',fs);
    
    y5 = y4-dy*5;
    text(x,y5,'5%-24% slope comparisons','fontsize',16);
    text(x,y5-dy,sprintf('combined: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm5v24.df,Statsm5v24.tstat,Pm5v24,CIm5v24(1),CIm5v24(2)),'fontsize',fs);
    text(x,y5-dy*2,sprintf('lem: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm5v24_lem.df,Statsm5v24_lem.tstat,Pm5v24_lem,CIm5v24_lem(1),CIm5v24_lem(2)),'fontsize',fs);
    text(x,y5-dy*3,sprintf('jbe: t(%i)=%.2f, P=%.2d [95%% CIs: (%.2f, %.2f)]',Statsm5v24_jbe.df,Statsm5v24_jbe.tstat,Pm5v24_jbe,CIm5v24_jbe(1),CIm5v24_jbe(2)),'fontsize',fs);
    
    set(gcf,'color','white','position',[100,200,800,700])
end

function set_errorbars(menu,menuEvt,myFig,evt,barType)
    

    % We only want to be able to tick one at a time..
    parentMenu = menu.Parent;
    siblingMenus = parentMenu.Children;
    for m = 1:length(siblingMenus);
        currentMenu = siblingMenus(m);
        
        if menu ~= currentMenu;
            set(currentMenu,'Checked','Off');
        end
        
    end

    if strcmpi(menu.Checked,'On');
        set(menu,'Checked','Off');
    else
        set(menu,'Checked','On');
    end
    
    
    
    
    
    figData = getappdata(myFig);
    
    setappdata(figData.TCHandle,'ErrorbarsToggled',barType);
    
    TC_callback(myFig,evt);
    
end