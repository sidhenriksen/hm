function PlotRank(varargin)

    load('CuratedCells.mat');
    
    myFig = figure();
    set(gcf,'color','white')
    plotNormal(Base)
    
    mh = uimenu(myFig,'Label','Plot');
    
    % For the main plot
    main = uimenu(mh,'Label','Normal rank','checked','off','callback',{@keydown_callback,Base,'normal'});
    
    % For the alternative plot
    alt = uimenu(mh,'Label','By AC slope');
    corr = uimenu(alt,'Label','Correlated','checked','off','callback',{@keydown_callback,Base,'alt','corr'});
    hm = uimenu(alt,'Label','Half-matched','checked','off','callback',{@keydown_callback,Base,'alt','hm'});
    ac = uimenu(alt,'Label','Anticorrelated','checked','off','callback',{@keydown_callback,Base,'alt','ac'});
    
end

function plotNormal(Base)
    for d = [1,4];
        
        nCells = length(Base(d).Cells);
        dxLengths = zeros(1,nCells);
        for j = 1:nCells;
            dxLengths(j) = length(Base(d).Cells(j).Dx);
        end

        dxL = min(dxLengths);
        
        if d == 1;
            rankData.lowDensity = zeros(nCells,dxL);
        elseif d == 4;
            rankData.highDensity = zeros(nCells,dxL);
        end
        
        correlatedRank = zeros(nCells,dxL);
        halfmatchedRank = zeros(nCells,dxL);
        anticorrelatedRank = zeros(nCells,dxL);
        
        for cell = 1:nCells;
            %% First we compute the rank order for each cell
            
            % Extract and normalize data
            
            correlatedResponse = Base(d).Cells(cell).correlatedResponse;
            halfmatchedResponse = Base(d).Cells(cell).halfmatchedResponse;
            anticorrelatedResponse = Base(d).Cells(cell).anticorrelatedResponse;
            
            [cRank,cIndex] = sort(correlatedResponse,'descend');
            hmRank = halfmatchedResponse(cIndex);
            acRank = anticorrelatedResponse(cIndex);
            
            Dx = Base(d).Cells(cell).Dx;
            
            
            % Resample disparity values
            newDx = linspace(min(Dx),max(Dx),dxL);
            
            cRank2 = interp1(Dx,cRank,newDx);
            hmRank2 = interp1(Dx,hmRank,newDx);
            acRank2 = interp1(Dx,acRank,newDx);
            
            maxC = max([cRank2,hmRank2,acRank2]); minC = min([cRank2,hmRank2,acRank2]);
            correlatedRank(cell,:) = (cRank2-minC)/(maxC-minC);
            halfmatchedRank(cell,:) = (hmRank2-minC)/(maxC-minC);
            anticorrelatedRank(cell,:) = (acRank2-minC)/(maxC-minC);
            

        end
        
        rank = 1:9;
        
        if d == 1;
            subCount = 1;
            maxSub = 2;
        elseif d == 4; 
            subCount = 2;
            maxSub = 2;
        else
            subCount = 3;
            maxSub = 3;
        end
        subplot(1,maxSub,subCount);

        corrRank = mean(correlatedRank); hmRank = mean(halfmatchedRank); acRank = mean(anticorrelatedRank);
        
        N = size(correlatedRank,1);
        corrSEM = std(correlatedRank)/sqrt(N); hmSEM = std(halfmatchedRank)/sqrt(N); acSEM = std(anticorrelatedRank)/sqrt(N);
        
    
        hold off;
        E1=errorbar(rank,corrRank,corrSEM*1.96);
        set(E1,'color',[0.8,0.1,0.1],'marker','o','markerfacecolor',[0.8,0.1,0.1], ...
            'linestyle',':','markersize',10,'linewidth',3);
        
        hold on;
        E2=errorbar(rank,hmRank,hmSEM*1.96);
        set(E2,'color',[0.1,0.1,0.8],'marker','o','markerfacecolor',[0.1,0.1,0.8], ...
            'linestyle',':','markersize',10,'linewidth',3);
        E3=errorbar(rank,acRank,acSEM*1.96);
        set(E3,'color',[0,0,0],'marker','o','markerfacecolor',[0,0,0], ...
            'linestyle',':','markersize',10,'linewidth',3);
        
        xlabel('Rank','fontsize',20);
        ylabel('Mean normalized response','fontsize',20);
        leg = legend('Correlated','Half-matched','Anti-correlated');
        set(leg,'fontsize',14);
        set(gca,'fontsize',18,'ytick',0:0.25:1);
        title(['Density: ',num2str(Base(d).density),' %'],'fontsize',22);
        ylim([0,1]);
        
    end
end


function plotAlt(Base,type)

    for d = [1,4];
        
        nCells = length(Base(d).Cells);
        dxLengths = zeros(1,nCells);
        for j = 1:nCells;
            dxLengths(j) = length(Base(d).Cells(j).Dx);
        end

        dxL = min(dxLengths);
        
        if d == 1;
            rankData.lowDensity = zeros(nCells,dxL);
        elseif d == 4;
            rankData.highDensity = zeros(nCells,dxL);
        end
        
        
        acSlopes = cat(1,Base(d).Cells.regAc);
        acSlopes = acSlopes(:,2);
        q50 = quantile(acSlopes,0.5);
        
        nLower = sum(acSlopes < q50);
        nUpper = sum(acSlopes > q50);

        
        lowerPlotRank = zeros(nLower,dxL);
        upperPlotRank = zeros(nUpper,dxL);
        lowerCorrelatedRank = zeros(nLower,dxL);
        lowerHalfmatchedRank = zeros(nLower,dxL);
        lowerAnticorrelatedRank = zeros(nLower,dxL);
        
        upperCorrelatedRank = zeros(nUpper,dxL);
        upperHalfmatchedRank = zeros(nUpper,dxL);
        upperAnticorrelatedRank = zeros(nUpper,dxL);
        
        lowerCell = 0; upperCell = 0;
        for cell = 1:nCells;
            %% First we compute the rank order for each cell
            
            % Extract and normalize data
            
            correlatedResponse = Base(d).Cells(cell).correlatedResponse;
            halfmatchedResponse = Base(d).Cells(cell).halfmatchedResponse;
            anticorrelatedResponse = Base(d).Cells(cell).anticorrelatedResponse;
            
            [cRank,cIndex] = sort(correlatedResponse,'descend');
            hmRank = halfmatchedResponse(cIndex);
            acRank = anticorrelatedResponse(cIndex);
            
            Dx = Base(d).Cells(cell).Dx;
            
            
            % Resample disparity values
            newDx = linspace(min(Dx),max(Dx),dxL);
            
            cRank2 = interp1(Dx,cRank,newDx);
            hmRank2 = interp1(Dx,hmRank,newDx);
            acRank2 = interp1(Dx,acRank,newDx);
            
            maxC = max([cRank2,hmRank2,acRank2]); minC = min([cRank2,hmRank2,acRank2]);
            
            acSlope = Base(d).Cells(cell).regAc(2);
            
            if acSlope > q50;
                upperCell = upperCell+1;
                upperCorrelatedRank(upperCell,:) = (cRank2-minC)/(maxC-minC);
                upperHalfmatchedRank(upperCell,:) = (hmRank2-minC)/(maxC-minC);
                upperAnticorrelatedRank(upperCell,:) = (acRank2-minC)/(maxC-minC);
                
                switch type
                    case 'corr'
                        upperPlotRank(upperCell,:) = (cRank2-minC)/(maxC-minC);
                    case 'hm'
                        upperPlotRank(upperCell,:) = (hmRank2-minC)/(maxC-minC);
                    case 'ac'
                        upperPlotRank(upperCell,:) = (acRank2-minC)/(maxC-minC);
                end
                
                % Remove this after:
                upperHalfmatchedRank(upperCell,:) = (acRank2-minC)/(maxC-minC);
            else
                lowerCell = lowerCell+1;
                switch type
                    case 'corr'
                        lowerPlotRank(lowerCell,:) = (cRank2-minC)/(maxC-minC);
                    case 'hm'
                        lowerPlotRank(lowerCell,:) = (hmRank2-minC)/(maxC-minC);
                    case 'ac'
                        lowerPlotRank(lowerCell,:) = (acRank2-minC)/(maxC-minC);
                end
                
            end
            

        end
        
        rank = 1:9;
        
        if d == 1;
            subCount = 1;
            maxSub = 2;
        elseif d == 4; 
            subCount = 2;
            maxSub = 2;
        else
            subCount = 3;
            maxSub = 3;
        end
        subplot(1,maxSub,subCount); hold on;
        
        upperRank = mean(upperPlotRank);
        lowerRank = mean(lowerPlotRank);
        
        N_lower = size(lowerPlotRank,1);
        N_upper = size(upperPlotRank,1);
        
        lowerSEM = std(lowerPlotRank)/sqrt(N_lower);
        upperSEM = std(upperPlotRank)/sqrt(N_upper);
        hold off;
        E1=errorbar(rank,lowerRank,lowerSEM*1.96);
        set(E1,'color',[0.1,0.1,0.8],'marker','o','markerfacecolor',[0.1,0.1,0.8], ...
            'linestyle',':','markersize',10,'linewidth',3);
        hold on;
        E2=errorbar(rank,upperRank,upperSEM*1.96);
        set(E2,'color',[0.8,0.1,0.8],'marker','o','markerfacecolor',[0.8,0.1,0.8], ...
            'linestyle',':','markersize',10,'linewidth',3);
        
        leg = legend([E1,E2],'Lower','Upper');
                
        xlabel('Rank','fontsize',20);
        ylabel('Mean normalized response','fontsize',20);
        
        set(leg,'fontsize',14);
        set(gca,'fontsize',18,'ytick',0:0.25:1);
        title(['Density: ',num2str(Base(d).density),' %'],'fontsize',22);
        ylim([0,1]);
        
    end
end

function [P,Q] = resampleCell(Cell,K,varargin)

% So what we want to do here:
% Figure out what integers make it so that the length of the 
% disparity tuning curve for Cell ends up being K following resampling.
% We're just doing that with a nested loop as the space is so small
    N = length(Cell.Dx);
    if nargin == 2
        multiplier = 1;
    else
        multiplier = varargin{1};
    end
    
    Ps = 1:(25*multiplier);
    Qs = 1:(25*multiplier);
        
    
    for Q = Qs;
        for P = Ps;
            cost = intCost(N,P,Q,K);
            if cost < 0.001;
                break
            end
        end
        if cost < 0.001;
            break
        end
    end
    
    if cost >= 0.001
        warning('Failed to find a solution; extending the range');
        [P,Q] = resampleCell(Cell,K,multiplier+1);
    end
end

function cost = intCost(N,P,Q,K)
    % Usage: cost = intCost(N,P,Q,K);
    % N: Length of array
    % P: Numerator to be optimized
    % Q: Denominator to be optimized
    % K: Desired length
    % This function is a stepping stone for figuring out
    % how much you have to resample your array in order to
    % get the desired size.
        
    cost = sum((K-N*P/Q).^2);
    %sprintf('Current cost: %2.f',cost)
end

function keydown_callback(myFig,evt,Base,type,varargin)
    if nargin > 4
        plotType = varargin{1};
    else
        plotType = 'none';
    end
    
    if strcmp(type,'normal');
        plotNormal(Base);
        
    elseif strcmp(type,'alt');
        plotAlt(Base,plotType);
    end

end