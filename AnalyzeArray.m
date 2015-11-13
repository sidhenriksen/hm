function [Dx,MeanResponse,SDResponse,SEM] = AnalyzeArray(AllExpt);
    NumTotalTrials = length(AllExpt.Expt.Trials);
    count = 0;
    for j = 1:NumTotalTrials;
        if AllExpt.Expt.Trials(j).excluded == 0
            count = count+1;
            ValidTrialList(count) = j;
        end
    end

    NumTrials = length(ValidTrialList);
    NumNeurons = length(AllExpt.Spikes);

    disparity = zeros(1,NumTrials);
    mixac = zeros(1,NumTrials);
    density = zeros(1,NumTrials);

    SpikeCount = zeros(NumNeurons,NumTrials);

    for neuron = 1:NumNeurons;  
        for j = 1:NumTrials;
            trial = ValidTrialList(j);
            if AllExpt.Expt.Trials(trial).excluded == 0;

                disparity(j) = AllExpt.Expt.Trials(trial).dx;
                mixac(j) = AllExpt.Expt.Trials(trial).mixac;
                density(j) = AllExpt.Expt.Trials(trial).dd;

            
                NumSpikes = length(AllExpt.Spikes{neuron}.Spikes{j});
                % Alternatively:
                % NumSpikes = ...
                %          sum(AllExpt.Spikes{neuron}.Spikes{trial} > 0);
                SpikeCount(neuron,j) = NumSpikes;
            end
        end
    end
    dx_values = unique(disparity);
    corr_values = unique(mixac);
    dd_values = unique(density);

    % If we want to exclude 0 disparity
    dx_values = dx_values(dx_values ~= 0);


    TuningCurves_AllNeurons = zeros(length(dx_values),length(corr_values),length(dd_values),NumNeurons);
    TuningCurves_AllNeurons_SD = zeros(length(dx_values),length(corr_values),length(dd_values),NumNeurons);

    for dx = 1:length(dx_values);
        for corr = 1:length(corr_values);
            for dd = 1:length(dd_values);
                CurrentTrials = ...
                    (disparity == dx_values(dx)) .* ...
                    (mixac == corr_values(corr)) .* ...
                    (density == dd_values(dd));

                CurrentSpikeCountMean = mean(SpikeCount(:,logical(CurrentTrials)),2);
                CurrentSpikeCountSD = std(SpikeCount(:,logical(CurrentTrials)),[],2);

                TuningCurves_AllNeurons_SD(dx,corr,dd,:) = CurrentSpikeCountSD;
                TuningCurves_AllNeurons(dx,corr,dd,:) = CurrentSpikeCountMean;
            end
        end
    end
    Dx = dx_values;
    MeanResponse = TuningCurves_AllNeurons;
    SDResponse = TuningCurves_AllNeurons_SD;

    N = NumTrials/(length(dx_values)*length(dd_values)*length(corr_values));

    SEM = sqrt(SDResponse/N)*1.96;
    
end