function spikeCounts = GetCounts(Cell,mixac,varargin)
    % We probably want to handle this differently for lemM116
    
    
    if nargin == 2;
        load(Cell.filename)
    else
        AllExpt = varargin{1};
    end
    
    isM116 = any(strfind(Cell.filename,'M116'));
    
    if isM116
        if nargin == 2;
            AllExpt = cExpt;
        end
        allDxs = [AllExpt.Trials.dx];
        allDds = [AllExpt.Trials.dd];
        allMixac = [AllExpt.Trials.mixac];
        allCes = [AllExpt.Trials.ce];
        allMes = [AllExpt.Trials.me];
        
        dx_values = unique(allDxs);
        
    else    
        whichCell = find(cat(1,AllExpt.Header.cellnumber) == Cell.cellnumber);

        nTrials = length(AllExpt.Spikes{whichCell}.Spikes);

        allDxs = round(cat(1,AllExpt.Expt.Trials.dx),3);


        allDds = cat(1,AllExpt.Expt.Trials.dd);
        allMixac = cat(1,AllExpt.Expt.Trials.mixac);
        
        if isfield(AllExpt.Expt.Trials,'ce');
            allCes = cat(1,AllExpt.Expt.Trials.ce);
        else
            allCes = ones(length(allMixac),1);
        end
        
        if isfield(AllExpt.Expt.Trials,'me');
            allMes = cat(1,AllExpt.Expt.Trials.me);
        else
            allMes = zeros(length(allMixac),1);
        end
        
        % So we need to account for the fact that not all neurons are going to
        % be present for the entire experiment...
        
        % Find all trial ids
        allTrialIds = [AllExpt.Expt.Trials.id];
        
        % And the trial ids that the neuron was available for
        thisCellTrialIds = [AllExpt.Spikes{whichCell}.trialid];
        
        % Find the indices 
        [~,ordinalIndex] =  intersect(allTrialIds,thisCellTrialIds);
        
        % Only consider the parameter values for which this neuron was
        % present
        allDxs = allDxs(ordinalIndex); allDds = allDds(ordinalIndex);
        allMixac = allMixac(ordinalIndex); allCes = allCes(ordinalIndex);
        allMes = allMes(ordinalIndex);
    


        dx_values = unique(allDxs);

        spikeCounts = cell(1,length(dx_values));
    end
    
    

    
    for dx = 1:length(dx_values);        
    
        currentTrials = find((allDxs == dx_values(dx)) .* (allDds == Cell.density) ...
                        .* (allMixac == mixac) .* (allCes == 1) .* (allMes == 0));


        currentSpikeCounts = zeros(1,length(currentTrials));

        for trial = 1:length(currentTrials);

            if isM116;
                currentSpikes = AllExpt.Trials(currentTrials(trial)).Spikes;
            else
                currentSpikes = AllExpt.Spikes{whichCell}.Spikes{currentTrials(trial)};
            end
            currentSpikeCounts(trial) = sum(currentSpikes > 0);

        end
        spikeCounts{dx} = currentSpikeCounts;
    end
    
    
    
end