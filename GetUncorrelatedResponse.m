function [MeanResponse,SDResponse,SEM] = GetUncorrelatedResponse(AllExpt);
    % Ignore data from probes that aren't cells:



    validTrialList = find(cat(1,AllExpt.Expt.Trials.excluded)==0);
    trialId = cat(1,AllExpt.Expt.Trials.Trial);
    trialId = trialId(validTrialList);
    
    %%% All the stimulus dimensions you want to look at
    disparities = cat(1,AllExpt.Expt.Trials.dx);
    disparities = disparities(validTrialList);
    disparities = round(disparities,2);
    
    
    mixacs = cat(1,AllExpt.Expt.Trials.mixac);
    mixacs = mixacs(validTrialList);
    
    try
        densities = cat(1,AllExpt.Expt.Trials.dd);
        densities = densities(validTrialList);
    catch
        densities = ones(length(validTrialList),1);
    end
    
    try
        dotSizes = cat(1,AllExpt.Expt.Trials.dw);
        dotSizes = dotSizes(validTrialList);
    catch
        dotSizes = ones(length(validTrialList),1);
    end
       
    if isfield(AllExpt.Expt.Trials,'ce');
        ces = cat(1,AllExpt.Expt.Trials.ce);
    else
        ces = ones(length(validTrialList),1);
    end        
        
    dd_values = unique(densities);
    dw_values = unique(dotSizes);
    %ce_values = unique(ces);
    
    nNeurons = sum(cat(1,AllExpt.Header.cellnumber) > 0);

    % If we want to exclude 0 disparity
    %dx_values = dx_values(dx_values ~= 0);

    spikeCountM = zeros(length(dd_values),length(dw_values),nNeurons);
    spikeCountSD = zeros(length(dd_values),length(dw_values),nNeurons);
    Ns = zeros(length(dd_values),length(dw_values),nNeurons);
    
    for neuron = 1:nNeurons;
        allCurrentTrials = AllExpt.Spikes{neuron}.Trial;
        [~,currentIndex,~] = intersect(trialId,allCurrentTrials);
        currentTrialsMask = zeros(size(trialId)); currentTrialsMask(currentIndex)=1;

        for dd = 1:length(dd_values);
            for dw = 1:length(dw_values);
                currentTrials = ...
                    (densities == dd_values(dd)) .* ...
                    (ces == 1);


                [~,~,ordinalIndex] = intersect(find(currentTrials),currentIndex);

                currentSpikes = cat(1,AllExpt.Spikes{neuron}.Spikes(ordinalIndex));

                miniSpikeCount = zeros(1,length(currentSpikes));

                for t = 1:length(currentSpikes);

                    spikes = currentSpikes{t}>0;

                    miniSpikeCount(t) = sum(spikes);
                end



                currentMean = mean(miniSpikeCount);
                currentSD = std(miniSpikeCount);

                spikeCountM(dd,dw,neuron) = currentMean;
                spikeCountSD(dd,dw,neuron) = currentSD;
                Ns(dd,dw,neuron) = length(miniSpikeCount);
            end
        end

    end
        
    MeanResponse = spikeCountM;
    SDResponse = spikeCountSD;
    SEM = sqrt(SDResponse./Ns);
    
end