function PlotRank(varargin)

    load('CuratedCells.mat');
    
    figure();
    set(gcf,'color','white')
    for d = [1,4];
        
        nCells = length(Base(d).Cells);
        dxLengths = zeros(1,nCells);
        for j = 1:nCells;
            dxLengths(j) = length(Base(d).Cells(j).Dx);
        end

        dxL = min(dxLengths);
        costs = zeros(1,100);
        if d == 1;
            rankData.lowDensity = zeros(nCells,dxL);
        elseif d == 4;
            rankData.highDensity = zeros(nCells,dxL);
        end
        
        for cell = 1:nCells;
            %% First we need to resample the tuning curves to make them align
            currentDx = Base(d).Cells(cell).Dx;
            
            % Get resampling coefficients
            [P,Q] = resampleCell(Base(d).Cells(cell),dxL);
            
            % Resample disparity values
            Base(d).Cells(cell).Dx = linspace(currentDx(1),currentDx(end),dxL);

            % Resample responses
            Base(d).Cells(cell).correlatedResponse = resample(Base(d).Cells(cell).correlatedResponse,P,Q);
            Base(d).Cells(cell).anticorrelatedResponse = resample(Base(d).Cells(cell).anticorrelatedResponse,P,Q);
            Base(d).Cells(cell).halfmatchedResponse = resample(Base(d).Cells(cell).halfmatchedResponse,P,Q);
            
            % And SEMs
            Base(d).Cells(cell).correlatedSEM= resample(Base(d).Cells(cell).correlatedSEM,P,Q);
            Base(d).Cells(cell).anticorrelatedSEM = resample(Base(d).Cells(cell).anticorrelatedSEM,P,Q);
            Base(d).Cells(cell).halfmatchedSEM = resample(Base(d).Cells(cell).halfmatchedSEM,P,Q);
            

        end
        
        %% Now we can start doing the ranks
        
        % Extract the tuning curve data
        correlatedData = cat(1,Base(d).Cells.correlatedResponse);
        anticorrelatedData = cat(1,Base(d).Cells.anticorrelatedResponse);
        halfmatchedData = cat(1,Base(d).Cells.halfmatchedResponse);
        
        % Normalize it
        normMatrix = repmat(max(correlatedData,[],2),[1,dxL]);
        
        correlatedData = correlatedData./normMatrix;
        anticorrelatedData = anticorrelatedData./normMatrix;
        halfmatchedData = halfmatchedData./normMatrix;
        
        % Create the ranks
        [correlatedRank,columnIndices] = sort(correlatedData,2,'descend');
        
        % Need to initialize these into their appropriate shapes
        anticorrelatedRank = zeros(size(correlatedRank));
        halfmatchedRank = zeros(size(correlatedRank));
        
        %columnIndices = columnIndices;
        
        rowIndices = repmat(1:nCells,[1,dxL]);
        
        
        indices = sub2ind([nCells,dxL],rowIndices(:),columnIndices(:));
        
        %anticorrelatedRank = anticorrelatedData(rowIndices',columnIndices(:));
        anticorrelatedRank(:) = anticorrelatedData(indices);
        halfmatchedRank(:) = halfmatchedData(indices);
        
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
        
        corrRank = mean(correlatedRank); hmRank = mean(halfmatchedRank); acRank = mean(anticorrelatedRank);
        plot(rank,corrRank,'color',[0.8,0.1,0.1],'marker','o','markerfacecolor',[0.8,0.1,0.1], ...
            'linestyle',':','markersize',10,'linewidth',3);
        plot(rank,hmRank,'color',[0.1,0.1,0.8],'marker','o','markerfacecolor',[0.1,0.1,0.8], ...
            'linestyle',':','markersize',10,'linewidth',3);
        plot(rank,acRank,'color',[0,0,0],'marker','o','markerfacecolor',[0,0,0], ...
            'linestyle',':','markersize',10,'linewidth',3);
        
        xlabel('Rank','fontsize',20);
        ylabel('Mean normalized response','fontsize',20);
        leg = legend('Correlated','Half-matched','Anti-correlated');
        set(leg,'fontsize',14);
        set(gca,'fontsize',18,'ytick',0:0.25:1);
        title(['Density: ',num2str(Base(d).density),' %'],'fontsize',22);
        ylim([0.25,0.65]);
        
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

function cost = intCost(N,P,Q,K);
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

