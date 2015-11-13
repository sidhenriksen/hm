function out = plot_TCs_base(cellStructure,varargin)
    % To use this, load CuratedCells.mat and feed it a cell structure.
    % This will plot tuning curves for the specified cells
    % Example:
    % X = load('CuratedCells.mat');
    % plot_TCs_base(X.Base(1).Cells(15:20))
    
    plotErrorbars = 0; slopeText = 0;
    if nargin > 1; % Still not implemented...
        for j = 1:length(varargin);
            switch varargin{j};
                
                case 'errorbar';
                    plotErrorbars = varargin{j+1};
                    
                case 'slope'
                    if strcmpi(varargin{j+1},'on');
                        slopeText = 1;
                    else
                        slopeText = 0;
                    end
            
            end
                        
        end
            
    end
        
        
    nCells = size(cellStructure,2);
    nRows = floor(sqrt(nCells));
    nCols = ceil(sqrt(nCells));
    while nRows * nCols < nCells;
        nCols = nCols+1;
    end
    
    out = figure(); 
    for iCell = 1:nCells;
        subplot(nRows,nCols,iCell); hold on;
        Dx = cellStructure(iCell).Dx;
        corr = cellStructure(iCell).correlatedResponse;
        corrSEM = cellStructure(iCell).correlatedSEM;
        
        anticorr = cellStructure(iCell).anticorrelatedResponse;
        acSEM = cellStructure(iCell).anticorrelatedSEM;
        
        hm = cellStructure(iCell).halfmatchedResponse;
        hmSEM = cellStructure(iCell).halfmatchedSEM;

        if plotErrorbars
            E=errorbar(Dx,corr,corrSEM*1.96);
            set(E,'linewidth',4,'markerfacecolor','r','linestyle', '--','color','r', ...
                'marker','o','markersize',12,'markeredgecolor','k');

            E = errorbar(Dx,anticorr,acSEM*1.96);
            set(E,'linewidth',4,'markerfacecolor','k','linestyle','--','color','k',...
                'marker','o','markersize',12,'markeredgecolor','k');

            E = errorbar(Dx,hm,hmSEM*1.96);
            set(E,'linewidth',4,'markerfacecolor','b','linestyle','--','color','b',...
                'marker','o','markersize',12,'markeredgecolor','k');
            
        else
            plot(Dx,corr,'linewidth',4,'markerfacecolor','r','linestyle', '--','color','r', ...
                'marker','o','markersize',12,'markeredgecolor','k');
            plot(Dx,anticorr,'linewidth',4,'markerfacecolor','k','linestyle','--','color','k',...
                'marker','o','markersize',12,'markeredgecolor','k');
            plot(Dx,hm,'linewidth',4,'markerfacecolor','b','linestyle','--','color','b',...
                'marker','o','markersize',12,'markeredgecolor','k');
        end

        xlabel('Disparity (deg)','fontsize',25);
        ylabel('Mean spikes per trial','fontsize',25);
        set(gca,'fontsize',20);
        
        % Plot uncorrelated if we have that data...
        if ~isempty(cellStructure(iCell).uncorrelatedResponse)
            plot(Dx,Dx*0 + cellStructure(iCell).uncorrelatedResponse, 'color',[0.15,0.8,0.1], ...
                'linewidth',2,'linestyle','--')
        end
        
        cellname = [cellStructure(iCell).filename(18:24),'-',num2str(cellStructure(iCell).cellnumber)];
        title(cellname,'fontsize',20)
        
        if slopeText
            xLim = get(gca,'xlim');
            yLim = get(gca,'yLim');
            
            hmSlope = cellStructure(iCell).regHm(2);
            text(xLim(2)-0.5,yLim(2)-yLim(2)*0.05,sprintf('C-HM slope: %.2f',hmSlope),'fontsize',16)
           
        end
            
    end
    set(out,'color','white');
end