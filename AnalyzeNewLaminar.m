function [Dx,MeanResponse,SDResponse,SEM,RMS,rootSpikeCount] = AnalyzeNewLaminar(AllExpt);
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
    
    if isfield(AllExpt.Expt.Trials,'me');
        mes = cat(1,AllExpt.Expt.Trials.me);
    else
        mes = zeros(length(validTrialList),1);
    end
    
    dx_values = unique(disparities);
    corr_values = unique(mixacs);
    dd_values = unique(densities);
    dw_values = unique(dotSizes);
    %ce_values = unique(ces);
    
    nNeurons = sum(cat(1,AllExpt.Header.cellnumber) > 0);

    % If we want to exclude 0 disparity
    %dx_values = dx_values(dx_values ~= 0);

    spikeCountM = zeros(length(dx_values),length(dd_values),length(corr_values),length(dw_values),nNeurons);
    rootSpikeCount = zeros(length(dx_values),length(dd_values),length(corr_values),length(dw_values),nNeurons);
    spikeCountSD = zeros(length(dx_values),length(dd_values),length(corr_values),length(dw_values),nNeurons);
    Ns = zeros(length(dx_values),length(dd_values),length(corr_values),length(dw_values));
    
    for neuron = 1:nNeurons;
        allCurrentTrials = AllExpt.Spikes{neuron}.Trial;
        [~,currentIndex,~] = intersect(trialId,allCurrentTrials);
        currentTrialsMask = zeros(size(trialId)); currentTrialsMask(currentIndex)=1;
        
    
        
        
        for dx = 1:length(dx_values);
            for dd = 1:length(dd_values);
                for corr = 1:length(corr_values);
                    for dw = 1:length(dw_values);
                    currentTrials = ...
                        (disparities == dx_values(dx)) .* ...
                        (mixacs == corr_values(corr)) .* ...
                        (dotSizes == dw_values(dw)) .* ...
                        (densities == dd_values(dd)) .* ...
                        currentTrialsMask .* ...
                        (ces == 1) .* ...
                        (mes == 0);
                        %%%% 13927
                        
                        [~,~,ordinalIndex] = intersect(find(currentTrials),currentIndex);
                        
                        currentSpikesCell = AllExpt.Spikes{neuron}.Spikes(ordinalIndex);
                        
                        miniSpikeCount = cellfun(@sum_more_than_zero,currentSpikesCell);                        
                        
                        currentMean = mean(miniSpikeCount);
                        currentRootMean = mean(sqrt(miniSpikeCount));
                        currentSD = std(miniSpikeCount);
                        
                        %if isnan(currentMean);
                        %    currentMean = 0;
                        %    currentSD = 0;
                        %end
                        
                        spikeCountM(dx,dd,corr,dw,neuron) = currentMean;
                        rootSpikeCount(dx,dd,corr,dw,neuron) = currentRootMean;
                        spikeCountSD(dx,dd,corr,dw,neuron) = currentSD;
                        Ns(dx,dd,corr,dw,neuron) = length(miniSpikeCount);
                        

                        
                    end
                end
            end
        end
        
    end
    % I think this comes out to be the right answer...
    
    RMS = sqrt(mean(spikeCountSD.^2));
    Dx = dx_values;
    MeanResponse = spikeCountM;
    SDResponse = spikeCountSD;

    SEM = SDResponse./sqrt(Ns);
    
end

function mt0 = sum_more_than_zero(a);

    mt0 = sum(a > 0);
end