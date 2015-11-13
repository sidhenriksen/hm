function run_slope_CIs();
    % Function to get the confidence intervals for the slopes

    densities = [1,2];
    N = 100000;
    load('CuratedCells.mat');

    
    for dd = 1:length(densities);
        
        currentBase = Base(densities(dd));
        
        % List of all files
        fileList = unique({currentBase.Cells.filename});
        
        % List to track which files have been processed
        % These are ordered in Base so we only have to load each once
        processedFiles = {};
        
        for cell = 1:length(currentBase.Cells);
            currentCell = currentBase.Cells(cell);
            currentFile = currentCell.filename;
            
            matchCurrentFile = strcmp(processedFiles,currentFile);
            
            % If match is empty or currentFile is not in processedFiles,
            % then we run load the data...
            if isempty(matchCurrentFile);
                matchCurrentFile = 0;
            end
            
            if ~any(matchCurrentFile);
                loadedData = load(currentFile);
                if isfield(loadedData,'cExpt');
                    AllExpt = loadedData.cExpt;
                else
                    AllExpt = loadedData.AllExpt;
                end
            end
            
            correlatedTCs = get_bootstrap_TCs(currentCell,0,N,AllExpt);
            halfmatchedTCs = get_bootstrap_TCs(currentCell,0.5,N,AllExpt);
            anticorrelatedTCs = get_bootstrap_TCs(currentCell,1,N,AllExpt);
            
            CHm_r = zeros(1,N);
            CAc_r = zeros(1,N);
            CHm_slope = zeros(1,N);
            CAc_slope = zeros(1,N);
            
            for j = 1:N;
                currentC = correlatedTCs(:,j)';
                currentHm = halfmatchedTCs(:,j)';
                currentAc = anticorrelatedTCs(:,j)';
                
                lambda = sqrt(10); % Ratio of variances
                [b_CHm,m_CHm] = fit_bothsubj2error(currentHm,currentC,lambda);
                m_CHm = 1./m_CHm; b_CHm = -b_CHm/m_CHm;
                rHm = corr(currentHm',b_CHm + m_CHm*currentC');
            
                [b_Ac,m_Ac] = fit_bothsubj2error(currentAc,currentC,1);
                m_Ac = 1./m_Ac; b_Ac = -b_Ac/m_Ac;
                rAc = corr(currentAc',b_Ac + m_Ac*currentC');
                
                CHm_r(j) = rHm;
                CHm_slope(j) = m_CHm;
                
                CAc_r(j) = rAc;
                CAc_slope(j) = m_Ac;
                
                
            end
            
            Base(densities(dd)).Cells(cell).CHm_r_CI = [quantile(CHm_r,0.025),quantile(CHm_r,0.975)];
            Base(densities(dd)).Cells(cell).CHm_slope_CI = [quantile(CHm_slope,0.025),quantile(CHm_slope,0.975)];
            
            Base(densities(dd)).Cells(cell).CAc_r_CI = [quantile(CAc_r,0.025),quantile(CAc_r,0.95)];
            Base(densities(dd)).Cells(cell).CAc_slope_CI = [quantile(CAc_slope,0.025),quantile(CAc_slope,0.975)];
            
            nTC = 25; indices = randperm(N,nTC);
            if ~exist('a','var');
                a = figure();
            end
            hold off;
            for j = 1:nTC;
                plot(correlatedTCs(:,j),'color','red','linewidth',2);
                plot(halfmatchedTCs(:,j),'color','blue','linewidth',2);
                plot(anticorrelatedTCs(:,j),'color','black','linewidth',2);
                title(sprintf('Slope: M(m_CHm) = %.2f, cHm slope: %.2f, CI: [%.2f, %.2f]', mean(CHm_slope), currentCell.regHm(2),quantile(CHm_slope,0.025),quantile(CHm_slope,0.975)))
                if j == 1;
                    hold on;
                end
            end
            
            
        end
        
        
    end

save('CuratedCells.mat','Base');
end



function bootstrapMeans = get_bootstrap_TCs(currentCell,mixac,N,varargin);

    if nargin == 2;
        counts = GetCounts(currentCell,mixac);
    else
        AllExpt = varargin{1};
        counts = GetCounts(currentCell,mixac,AllExpt);
    end
    
    
    allCounts = [counts{:}];
    nCounts = length(allCounts);
    
    allDxIndices = zeros(1,nCounts); % These are not actual disparities, but disparity indices
    
    assert(length(counts) == length(currentCell.Dx),'Error: Number of spike count categories must match number of disparities');
    
    start = 1;
    means = zeros(1,length(counts));
    allResiduals = zeros(1,length(allCounts));
    
    Ns = zeros(1,length(counts));
    
    for dx = 1:length(counts);
        stop = start+length(counts{dx})-1;
        allDxIndices(start:stop) = dx;
        
        currentCounts = sqrt(counts{dx});
        currentMean = mean(currentCounts);
        
        means(dx) = mean(counts{dx});
        
        allResiduals(start:stop) = currentCounts-currentMean;
        
        
        start = stop+1;
    
    end
    
    
    %N = 5000; % Number of repeats    
    
    bootstrapMeans = zeros(length(counts),N);

    for dx = 1:length(counts);
        
        
        K = length(counts{dx});
        
        indices = randi(length(allResiduals),[1,N*K]);
        
        currentResiduals = reshape(allResiduals(indices),[K,N]);
        
        bootstrapMeans(dx,:) = (mean(currentResiduals)+sqrt(means(dx))).^2;
        
    end
       
    
    
end