function [Dx,MeanResponse,SDResponse,SEM,RMS,rootSpikeCount] = AnalyzeLaminar(cExpt)
    % Extract data
    NumNeurons = 1;
    
    allDisparities = [cExpt.Trials.dx];
    allMixacs = [cExpt.Trials.mixac];
    allDensities = [cExpt.Trials.dd];
    allSpikeCounts = [cExpt.Trials.count];
    allMes = [cExpt.Trials.me];
    allCes = [cExpt.Trials.ce];

    dx_values = unique(allDisparities);
    corr_values = unique(allMixacs);
    dd_values = unique(allDensities);

    % Ignore corr_values that are negative. Negative values of % dots
    % anti-correlated is not valid.
    corr_values = corr_values(corr_values >= 0);

    spikeCountM = zeros(length(dx_values),length(dd_values),length(corr_values),1,NumNeurons);
    rootSpikeCount = zeros(length(dx_values),length(dd_values),length(corr_values),1,NumNeurons);
    spikeCountSD = zeros(length(dx_values),length(dd_values),length(corr_values),1,NumNeurons);
    Ns = zeros(length(dx_values),length(dd_values),length(corr_values),1,NumNeurons);

    for dx = 1:length(dx_values);
        for corr = 1:length(corr_values);
            for dd = 1:length(dd_values);
                currentTrials = ...
                    find((allDisparities == dx_values(dx)) .* ...
                    (allMixacs == corr_values(corr)) .* ...
                    (allDensities == dd_values(dd)) .* ...
                    (allMes==0) .* (allCes==1)); % This last one is to exclude monocular and uncorrelated stimuli.

                
                
                currentSpikesCell= {cExpt.Trials(currentTrials).Spikes};
                currentSpikeCounts = cellfun(@sum_more_than_zero,currentSpikesCell);
                
                currentSpikeCountMean = mean(currentSpikeCounts);
                currentRootSpikeCount = mean(sqrt(currentSpikeCounts));
                CurrentSpikeCountSD = std(currentSpikeCounts);
                
                

                spikeCountM(dx,dd,corr,1,:) = currentSpikeCountMean;
                rootSpikeCount(dx,dd,corr,1,:) = currentRootSpikeCount;
                spikeCountSD(dx,dd,corr,1,:) = CurrentSpikeCountSD;
                
                Ns(dx,dd,corr,1,:) = length(currentTrials);
            end
        end
    end
    
    RMS = sqrt(mean(spikeCountSD.^2));
    Dx = dx_values;
    MeanResponse = spikeCountM;
    SDResponse = spikeCountSD;

    SEM = sqrt(SDResponse./Ns)*1.96;
    
end

function mt0 = sum_more_than_zero(a);
    
    mt0 = sum(a > 0);

end