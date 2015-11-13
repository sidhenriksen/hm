function stimvalue = get_stimvalues(Expt,name);
    % Function to return the stimulus values used in an experiment
    if isfield(Expt,'Stimvals');
        whereStimVal = 0;
    else 
        whereStimVal = 1;
    end
    
    
    
    if whereStimVal
        
        AllExpt = Expt;
        
        
        
        count = 0;
        
        NumTotalTrials = length(AllExpt.Expt.Trials);
        
        
        for j = 1:NumTotalTrials;
            if AllExpt.Expt.Trials(j).excluded == 0
                count = count+1;
                ValidTrialList(count) = j;
            end
        end

        NumTrials = length(ValidTrialList);
        
        all_stimvalues = zeros(1,NumTrials);
        
        if isfield(AllExpt.Expt.Trials,name);
            for j = 1:NumTrials;
                trial = ValidTrialList(j);
                if AllExpt.Expt.Trials(trial).excluded == 0;
                    all_stimvalues(j) = eval(['AllExpt.Expt.Trials(trial).',name]);
                end
            end
        else
            if isfield(AllExpt.Expt.Stimvals,'dw');
                all_stimvalues = AllExpt.Expt.Stimvals.dw;
            else
                all_stimvalues = [];
            end
        end
    else
        NumTrials = length(Expt.Trials);
        
        all_stimvalues = zeros(1,NumTrials);
        
        % Might have to add a thing that excludes bad trials here?
        for trial = 1:NumTrials;
            all_stimvalues(trial) = eval(['Expt.Trials(trial).',name]);
        end
        
    end
    
    stimvalue = unique(all_stimvalues);
end