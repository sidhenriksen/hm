function CurateCells(fileName,varargin);
    % Usage: CurateCells(input)
    % Curates input file and puts it into CuratedCells.mat
    % Usage: CurateCells(input,'filter',0)
    % Curated without filter
    % Usage: CurateCells(input,'output','otherfile.mat')
    % Stores data in otherfile.mat instead of default CuratedCells.mat
    
    inDir = '/sid/Ephys/HalfMatched/hm/';
    
    if nargin > 1
        for j = 1:length(varargin)
            switch varargin{j}
                case 'output'
                    outFile = varargin{j+1};
                case 'filter'
                    filter = varargin{j+1};
            end
        end
    end
    
    if ~exist('filter','var')
        filter = 1;
    end
    
    if ~exist('outFile','var');
        outFile = 'CuratedCells.mat';
    end
    
    allFiles = what(inDir);
    allMatFiles = allFiles.mat;
    
    if any(strcmp(outFile,allMatFiles));
        load([inDir,outFile]);
    else
        %Base = struct();
    end
    
    load(fileName);
    
    if exist('cExpt','var');
        newLaminar = 0;
    else
        newLaminar = 1;
    end

    % We need different procedures for laminar recordings and multi-electrode
    % arrays
    if newLaminar
        [Dx,MeanResponse,SDResponse,SEM,RMS,rootSpikeCount] = AnalyzeNewLaminar(AllExpt);
    
        Expt = AllExpt;
            
        if isfield(AllExpt.Expt.Trials,'ce');
            [ucResp,ucSD,ucSEM] = GetUncorrelatedResponse(AllExpt);
        end
    else
        [Dx,MeanResponse,SDResponse,SEM,RMS,rootSpikeCount] = AnalyzeLaminar(cExpt);
        Expt = cExpt;
    end


    dd_values = get_stimvalues(Expt,'dd');
    
    dd_indices = 1:length(dd_values);
        
    corr_values = get_stimvalues(Expt,'mixac');

    try
        dw_values = get_stimvalues(Expt,'dw');
    catch
        dw_values = 1;
    end

    corr_values = corr_values(corr_values >= 0);

    if newLaminar;
        Ps = {};
        
        for j = 1:length(dd_values);
            
            P = run_ANOVA(Expt,dd_values(j));
            % Only include this if it's one of our target densities
            if (dd_values(j) == 24)||(dd_values(j)==5);
                Ps{j} = P;
            else
                Ps{j} = zeros(size(P));
            end
        end
        P_thresh = 0.01;
        Ps = cat(1,Ps{:}) < P_thresh;
        Ps = prod(Ps,1);
        
        dims = size(MeanResponse);
        dw = dims(4); % Always select the largest dot size

    else
        Ps = 0.001;
        dw = 1;
    end

    neuronsAfterANOVA = find(Ps);
    
    

    
    MaxResp = max(squeeze(MeanResponse(:,1,1,dw,neuronsAfterANOVA)));
    %if any(isnan(MaxResp));
    %    dwqwe
    %end

    neuronsAfterMax = neuronsAfterANOVA(MaxResp >= 2);

    if filter
        PlotNeurons = neuronsAfterMax;
    else
        PlotNeurons = 1:size(MeanResponse,5);
    end
    
    
    meanResp = squeeze(MeanResponse(:,:,:,dw,:));
    meanRootResp = squeeze(rootSpikeCount(:,:,:,dw,:));
    sdResp = squeeze(SDResponse(:,:,:,dw,:));
    
    % We're going to split the recordings into densities
    % And then split the densities into cells
    for dd = 1:length(dd_values);
        
        % In this segment we sort out the densities
        
        % If Base exists (i.e. we're not making a new file)
        % then append to existing densities if density matches
        % Otherwise, start from scratch.
        if exist('Base','var')
            dd_count = length(Base); % Number of densities
            existingDensities = cat(1,Base.density);

        else
            dd_count = 0; % Number of densities
            existingDensities = [];
        end
        
        
        existingDensity = find(dd_values(dd) == existingDensities);
        
        % If the density isn't already added, we start a new density
        if isempty(existingDensity);
            existingDensity = dd_count+1;  
            Base(existingDensity).density = dd_values(dd);
            count = 0;
        else
            count = length(Base(existingDensity).Cells);
        end
        
        
        [HMaucs,~,HMdprimes] = ROCanalysis(Expt,0.5,dd_values(dd));
        [Caucs,~,Cdprimes] = ROCanalysis(Expt,0,dd_values(dd));
        HMaucsZ = ROCanalysis(Expt,0.5,dd_values(dd),'znorm');
        CaucsZ =  ROCanalysis(Expt,0,dd_values(dd),'znorm');
        
        
        for cell = 1:length(PlotNeurons);
            
            currentCell = PlotNeurons(cell);
            if newLaminar
                cellnumber = AllExpt.Header(currentCell).cellnumber;
            else
                cellnumber = 1;
            end
            
            count = count+1;
            
            
            % Extract current responses for the cell/density combo
            corrResp = meanResp(:,dd_indices(dd),1,currentCell)';
            halfmatchedResp = meanResp(:,dd_indices(dd),2,currentCell)';
            anticorrResp = meanResp(:,dd_indices(dd),3,currentCell)';
            
            corrRootResp = meanRootResp(:,dd_indices(dd),1,currentCell)';
            halfmatchedRootResp = meanRootResp(:,dd_indices(dd),2,currentCell)';
            anticorrRootResp = meanRootResp(:,dd_indices(dd),3,currentCell)';
            
            corrSD = sdResp(:,dd_indices(dd),1,currentCell)';
            hmSD = sdResp(:,dd_indices(dd),2,currentCell)';
            acSD = sdResp(:,dd_indices(dd),3,currentCell)';
            
            SEM2 = squeeze(SEM);
            corrSEM = SEM2(:,dd_indices(dd),1,currentCell);
            hmSEM = SEM2(:,dd_indices(dd),2,currentCell);
            acSEM = SEM2(:,dd_indices(dd),3,currentCell);
            
            
            
            % Get the RMS error
            currentRMS = squeeze(RMS(1,dd_indices(dd),:,dw,currentCell));
            
            % Now compute DDI
            sigma = currentRMS(1); % For correlated
            rMax = max(corrResp); 
            rMin = min(corrResp);
            DDI = (rMax-rMin)/(rMax-rMin+2*sigma);

            
            lambda = sqrt(10); % Ratio of the variances
            
            % Carry out regression between the two tuning curves
            [rHm,m1,b1] = regression2(halfmatchedResp,corrResp); % Correlated - half-matched
            
            [rAc,m2,b2] = regression2(anticorrResp,corrResp); % Correlated - anti-correlated
            
            [b3,m3] = fit_bothsubj2error(halfmatchedResp,corrResp,lambda);
            m3 = 1./m3; b3 = -b3/m3;
            rHm2 = corr(halfmatchedResp',b3 + m3*corrResp');
            
            [b4,m4] = fit_bothsubj2error(anticorrResp,corrResp,1);
            m4 = 1./m4; b4 = -b4/m4;
            rAc2 = corr(anticorrResp',b4 + m4*corrResp');
            
            % Now we store everything in our big structure
            Base(existingDensity).Cells(count).cellnumber = cellnumber;
            Base(existingDensity).Cells(count).filename = fileName;
            
            Base(existingDensity).Cells(count).regHm = [rHm2,m3,b3];
            Base(existingDensity).Cells(count).regAc = [rAc2,m4,b4];
            
            Base(existingDensity).Cells(count).regHmRegular = [rHm,m1,b1];
            Base(existingDensity).Cells(count).regAcRegular = [rAc,m2,b2];
            Base(existingDensity).Cells(count).Dx = Dx;
            
            Base(existingDensity).Cells(count).correlatedResponse = corrResp;
            Base(existingDensity).Cells(count).halfmatchedResponse = halfmatchedResp;
            Base(existingDensity).Cells(count).anticorrelatedResponse = anticorrResp;
            
            Base(existingDensity).Cells(count).correlatedSEM = corrSEM;
            Base(existingDensity).Cells(count).halfmatchedSEM = hmSEM;
            Base(existingDensity).Cells(count).anticorrelatedSEM = acSEM;
            
            Base(existingDensity).Cells(count).correlatedSD = corrSD;
            Base(existingDensity).Cells(count).halfmatchedSD = hmSD;
            Base(existingDensity).Cells(count).anticorrelatedSD = acSD;
            
            Base(existingDensity).Cells(count).RMS = currentRMS; % This gives for corr/hm/ac
            Base(existingDensity).Cells(count).DDI = DDI; % For corr only
            
            Base(existingDensity).Cells(count).HMauc = HMaucs(cell); % HM ROC AUC values
            Base(existingDensity).Cells(count).HMdprime = HMdprimes(cell); 
            
            Base(existingDensity).Cells(count).Cauc = Caucs(cell); 
            Base(existingDensity).Cells(count).Cdprime = Cdprimes(cell); 
            
            Base(existingDensity).Cells(count).HMaucZ = HMaucsZ(cell); % HM ROC AUC values Z score
            Base(existingDensity).Cells(count).CaucZ = CaucsZ(cell); 
            
            Base(existingDensity).Cells(count).density = dd_values(dd);
            
            Base(existingDensity).Cells(count).dw = dw_values(dw);
            
            
            
            
            
            if newLaminar
                if isfield(AllExpt.Expt.Trials,'ce');
                    
                    Base(existingDensity).Cells(count).uncorrelatedResponse = ucResp(dd_indices(dd),dw,currentCell);
                    Base(existingDensity).Cells(count).uncorrelatedSEM = ucSEM(dd_indices(dd),dw,currentCell);
                end
                
            end
            
            % Get bootstrap CIs for correlated, half-matched and anticorrelated
            [ciLowC,ciHighC] = bootstrap_hm(Base(existingDensity).Cells(count),0,Expt);
            [ciLowHm,ciHighHm] = bootstrap_hm(Base(existingDensity).Cells(count),0.5,Expt);
            [ciLowAc,ciHighAc] = bootstrap_hm(Base(existingDensity).Cells(count),1,Expt);
                
            Base(existingDensity).Cells(count).ciLowHm = ciLowHm;
            Base(existingDensity).Cells(count).ciHighHm = ciHighHm;

            Base(existingDensity).Cells(count).ciLowC = ciLowC;
            Base(existingDensity).Cells(count).ciHighC = ciHighC;

            Base(existingDensity).Cells(count).ciLowAc = ciLowAc;
            Base(existingDensity).Cells(count).ciHighAc = ciHighAc;
            
        end
        
        %%% Add penetration information 
        %newBase(existingDensity) = AddPenetrationInformation(Base(existingDensity));
 
    end
    save([inDir,outFile],'Base');
end

    
