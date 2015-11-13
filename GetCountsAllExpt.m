function spikeCounts = GetCountsAllExpt(AllExpt,Cell,mixac)
    % We probably want to handle this differently for lemM116

    whichCell = cat(1,AllExpt.Header.cellnumber) == Cell.cellnumber;
    
    nTrials = length(AllExpt.Spikes{whichCell}.Spikes);
    
    allDxs = round(cat(1,AllExpt.Expt.Trials.dx),3);
    
    
    allDds = cat(1,AllExpt.Expt.Trials.dd);
    allMixac = cat(1,AllExpt.Expt.Trials.mixac);
    
    allDxs = allDxs(1:nTrials); allDds = allDds(1:nTrials); allMixac = allMixac(1:nTrials);
    
    dx_values = unique(allDxs);
    
    spikeCounts = cell(1,length(dx_values));
    for dx = 1:length(dx_values);
        
        currentTrials = find((allDxs == dx_values(dx)) .* (allDds == Cell.density) ...
                        .* (allMixac == mixac));
                    
        currentSpikeCounts = zeros(1,length(currentTrials));
        
        for trial = 1:length(currentTrials);
            
            currentSpikes = AllExpt.Spikes{whichCell}.Spikes{currentTrials(trial)};
            
            currentSpikeCounts(trial) = length(currentSpikes > 0);
        end
        spikeCounts{dx} = currentSpikeCounts;
    end
    
end