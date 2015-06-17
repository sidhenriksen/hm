function [ciLow,ciHigh] = bootstrap_hm(currentCell,mixac,varargin);
    
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
    
    
    N = 5000; % Number of repeats    
    
    bootstrapMeans = zeros(length(counts),N);

    for dx = 1:length(counts);
        
        
        K = length(counts{dx});
        
        indices = randi(length(allResiduals),[1,N*K]);
        
        currentResiduals = reshape(allResiduals(indices),[K,N]);
        
        bootstrapMeans(dx,:) = (mean(currentResiduals)+sqrt(means(dx))).^2;
        
    end
       
    
    ciLow = quantile(bootstrapMeans',0.05);
    ciHigh = quantile(bootstrapMeans',0.95);
    
end