function [AUCs,data,dprimes] = ROCanalysis(AllExpt,mixac,dd,varargin);

    zNorm = 'none';
    for j = 1:length(varargin);
        switch varargin{j};
            case 'znorm'
                zNorm = 'znorm';
        end
    end
    
    % This is to deal with lemM116
    if isfield(AllExpt,'Spikes');
        P = run_ANOVA(AllExpt);
        [Dx,MeanResponse] = AnalyzeNewLaminar(AllExpt);
        
    
    else
        P = 0.001;
        [Dx,MeanResponse] = AnalyzeLaminar(AllExpt);
    
    end
    
    whichCells = find(P < 0.01);
    
    dd_values = get_stimvalues(AllExpt,'dd');
    mixac_values = get_stimvalues(AllExpt,'mixac');
    
    dd_idx = (dd_values == dd); 
    
    
    assert(sum(dd_idx)>0,'Error: Invalid density value');
    assert(sum(mixac_values == mixac)>0,'Error: Invalid mixac value');
    
    AUCs = zeros(1,length(whichCells));
    data = struct([]);
    
    dprimes = zeros(1,length(whichCells));
    %fprintf(whichCells)
    
    for k = 1:length(whichCells);
        
        whichCell = whichCells(k);

        
        
        
        if isfield(AllExpt,'Spikes');
            correlatedTC = MeanResponse(:,dd_idx,1,1,whichCell);
            [~,maxDx] = max(correlatedTC);
            [~,minDx] = min(correlatedTC);
            prefSpikeCounts = get_spikecounts(AllExpt,whichCell,mixac,Dx(maxDx),dd,zNorm);
            nullSpikeCounts = get_spikecounts(AllExpt,whichCell,mixac,Dx(minDx),dd,zNorm);
        else
            correlatedTC = MeanResponse(:,1,1);
            [~,maxDx] = max(correlatedTC);
            [~,minDx] = min(correlatedTC);
            prefSpikeCounts = get_spikecounts2(AllExpt,mixac,Dx(maxDx),dd,zNorm);
            nullSpikeCounts = get_spikecounts2(AllExpt,mixac,Dx(minDx),dd,zNorm);
        end
        
        n = length(prefSpikeCounts);
        m = length(nullSpikeCounts);
        
        binaryVector = [ones(1,n),zeros(1,m)];
        
        [sortedSpikeCounts,sortedIndices] = sort([prefSpikeCounts,nullSpikeCounts]);
        
        nonRepeats = find(diff(sortedSpikeCounts) > 0);
        
        sortedBinaryVector = binaryVector(sortedIndices);
        cumsumTrue = cumsum(sortedBinaryVector)/n;
        cumsumFalse = cumsum(~sortedBinaryVector)/m;
        
        TP = 1-cumsumTrue(nonRepeats); TP = [1,TP,0];
        FP = 1-cumsumFalse(nonRepeats); FP = [1,FP,0];
        TP = TP(length(TP):-1:1); FP = FP(length(FP):-1:1);
        
        pooledVariance = (var(prefSpikeCounts)*(n-1) + var(nullSpikeCounts)*(m-1))/(n+m-2);
        
        dprimes(k) = (mean(prefSpikeCounts)-mean(nullSpikeCounts))/sqrt(pooledVariance);
        
        data(k).TP = TP;
        data(k).FP = FP;
        
                
        AUCs(k) = trapz(FP,TP);
                
    end
        
end

function spikeCounts = get_spikecounts(AllExpt,neuron,mixac,dx,dd,varargin);
    % Usage: spikeCounts = get_spikecounts(AllExpt,neuron,dx,dd);
    % AllExpt: AllExpt structure
    % neuron: cell index
    % dx: disparity value
    % dd: dot density value
    
    znorm = 0;
    for j = 1:length(varargin);
        switch varargin{j};
            case 'znorm'
                znorm = 1;
        end
        
    end

    if ~isfield(AllExpt,'Expt');
        AllExpt.Expt.Header = AllExpt.Header;
        AllExpt.Expt.Stimvals = AllExpt.Stimvals;
        AllExpt.Expt.Trials = AllExpt.Trials;
        AllExpt.Expt.Comments = AllExpt.Comments;
        
    end
    
    trialid = AllExpt.Spikes{neuron}.trialid;
    allTrials = cat(1,AllExpt.Expt.Trials.id);
    
    % Find indices of trials where we have both spike data
    % and trial information
    [~,validTrialList] = intersect(allTrials,trialid);
   
    blockStartids = AllExpt.Expt.Header.BlockStartid;
    blockStartids = [blockStartids,10^6];
    nBlocks = length(blockStartids)-1;
    
    blockMeans = zeros(1,nBlocks);
    blockSDs = zeros(1,nBlocks);
    trialBlockid = zeros(1,length(trialid));
    
    for block = 1:nBlocks
        
        currentTrials = find(...
            (trialid >= blockStartids(block)) .* ...
            (trialid < blockStartids(block+1)));
        
        currentSpikes = {AllExpt.Spikes{neuron}.Spikes{currentTrials}};
        currentCounts = zeros(size(currentSpikes));
        
        for trial = 1:length(currentSpikes);
            currentCounts(trial) = length(currentSpikes{trial} > 0);
        end
        
        blockMeans(block) = mean(currentCounts);
        blockSDs(block) = std(currentCounts);
        
        trialBlockid(currentTrials) = block;
        
    end
    
    %%% All the stimulus dimensions you want to look at
    disparities = cat(1,AllExpt.Expt.Trials.dx);
    disparities = disparities(validTrialList);
    
    disparities = round(disparities,3); % Deal with the 0.001 issue
    
    mixacs = cat(1,AllExpt.Expt.Trials.mixac);
    mixacs = mixacs(validTrialList);
    
    densities = cat(1,AllExpt.Expt.Trials.dd);
    densities = densities(validTrialList);
    
    try
        dotSizes = cat(1,AllExpt.Expt.Trials.dw);
        dotSizes = dotSizes(validTrialList);
    catch
        dotSizes = ones(length(validTrialList),1);
    end
       
            
        
    % Take care of the uniqueness issue:
    dx_unique = unique(disparities);
    [dx_values,idx] = quasi_unique(dx_unique,0.001);
    diff = setdiff(dx_unique,dx_values);
    idxDiff = setdiff(1:length(dx_unique),idx);
    for k = 1:length(diff);
        disparities(disparities == diff(k)) = dx_unique(idxDiff(k)+1);
    end
    
    dw_values = unique(dotSizes);
    
   
    % Select current trials based on input (picks largest available dot
    % size)
    currentTrials = find(...
        (disparities == dx) .* ...
        (mixacs == mixac) .* ...
        (dotSizes == dw_values(length(dw_values))) .* ...
        (densities == dd));
        %%%% 13927


    currentSpikes = cat(1,AllExpt.Spikes{neuron}.Spikes(currentTrials));
    
    currentTrialBlocks = trialBlockid(currentTrials);
    trialBlockMeans = blockMeans(currentTrialBlocks);
    trialBlockSDs = blockSDs(currentTrialBlocks);

    spikeCounts = zeros(1,length(currentSpikes));

    for t = 1:length(currentSpikes);

        spikes = currentSpikes{t}>0;

        spikeCounts(t) = sum(spikes);
    end    

    
    if znorm
        spikeCounts = (spikeCounts-trialBlockMeans)./(trialBlockSDs);
    end
end

function spikeCounts = get_spikecounts2(cExpt,mixac,dx,dd,varargin);
    % Usage: spikeCounts = get_spikecounts2(cExpt,mixac,dx,dd)
    % Special function for getting spikeCounts for lemM116
    
    znorm = 0;
    for j = 1:length(varargin);
        switch varargin{j};
            case 'znorm'
                znorm = 1;
        end
                
    end
    
    dd_vals = cat(1,cExpt.Trials.dd);
    dx_vals = cat(1,cExpt.Trials.dx);
    mixac_vals = cat(1,cExpt.Trials.mixac);
    
    trialid = cat(1,cExpt.Trials.id);
    
    % This is for the z-score normalization
    blockStartids = cExpt.Header.BlockStart;
    blockStartids = [blockStartids,10^6];
    nBlocks = length(blockStartids)-1;
    
    blockMeans = zeros(1,nBlocks);
    blockSDs = zeros(1,nBlocks);
    trialBlockid = zeros(1,length(trialid));
    
    for block = 1:nBlocks
        
        currentTrials = find(...
            (trialid >= blockStartids(block)) .* ...
            (trialid < blockStartids(block+1)));
        
        
        currentSpikes = {cExpt.Trials(currentTrials).Spikes};
        currentCounts = zeros(size(currentSpikes));
        
        for trial = 1:length(currentSpikes);
            currentCounts(trial) = length(currentSpikes{trial} > 0);
        end
        
        blockMeans(block) = mean(currentCounts);
        blockSDs(block) = std(currentCounts);
        
        trialBlockid(currentTrials) = block;
        
    end
    
    
    trials = find((dd_vals == dd) .* (dx_vals == dx) ...
        .* (mixac_vals == mixac));
    
    allSpikes = {cExpt.Trials.Spikes};
    
    spikeCounts = zeros(1,length(trials));
    
    for j = 1:length(trials);
        
        spikeCounts(j) = sum(allSpikes{trials(j)} > 0);
        
    end
    
    if znorm
        trialByBlock = trialBlockid(trials);
        spikeCounts = (spikeCounts-blockMeans(trialByBlock))./blockSDs(trialByBlock);
    end
end