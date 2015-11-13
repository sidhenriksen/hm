function RunStats()
    spikeCounts = 0;
    ttests = 1;
    load('CuratedCells.mat');
    
    % We prob
    for dd = [1,4];
        currentData = Base(dd).Cells;
        nCells = length(currentData);
        
        lastFileName = '';
        if spikeCounts
            for cell = 2%:nCells;

                % Only load if we're encountering a new file
                % The cells are ordered so this works fine
                if ~strcmp(lastFileName,currentData(cell).filename);
                    load(currentData(cell).filename);
                end
                currentData(cell).density = Base(dd).density;

                hmSpikeCounts = GetCounts(AllExpt,currentData(cell),0.5);
                corrSpikeCounts = GetCounts(AllExpt,currentData(cell),0);

                % Okay, so we now have extracted the half-matched and
                % correlated spike counts... Next we want to randomly sample
                % from the bins to construct 95% confidence intervals.


            end
        end
        
        if ttests
            rs = currentData.
            
            
        end
    end

end

function spikeCounts = GetCounts(AllExpt,Cell,mixac)
    % We probably want to handle this differently for lemM116

    whichCell = cat(1,AllExpt.Header.cellnumber) == Cell.cellnumber;
    
    nTrials = length(AllExpt.Spikes{whichCell}.Spikes);
    
    allDxs = cat(1,AllExpt.Expt.Trials.dx);
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