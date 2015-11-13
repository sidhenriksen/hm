function compute_supplementary()

    load('CuratedCells.mat');

    for mixac = [0,0.5];
        for dd = [1,2]
            Ps = zeros(1,length(Base(dd).Cells));
            fnames = {};
            for j = 1:length(Base(dd).Cells);

                currentCell = Base(dd).Cells(j);

                if ~any(strcmp(fnames,currentCell.filename));
                    data=load(currentCell.filename);
                    if isfield(data,'cExpt');
                        Expt = data.cExpt;
                    else
                        Expt = data.AllExpt;
                    end
                end
                counts = GetCounts(currentCell,mixac,Expt);
                allCounts = [counts{:}];
                for k = 1:length(counts);
                    counts{k} = sqrt(counts{k});
                end

                spikeNum = cellfun(@length,counts);
                nDisp = length(spikeNum);

                

                conditions = zeros(1,sum(spikeNum));
                allResiduals = zeros(1,sum(spikeNum));
                start = 1;

                tc = cellfun(@mean,counts);

                for k = 1:nDisp;
                    % This is to compute the stuff for the ANOVA
                    stop = start+spikeNum(k)-1;
                    conditions(start:stop) = k;
                    

                    % This is for the DDI
                    residuals = (mean(counts{k})-counts{k});
                    allResiduals(start:stop) = residuals;
                    
                    start = stop+1;
                end

                RMSerror = sqrt(mean(allResiduals.^2));
                Rmax = max(tc); Rmin = min(tc);

                DDI = (Rmax-Rmin)/((Rmax-Rmin) + 2*RMSerror);

                P = anova1(allCounts,conditions,'off');
                Ps(j) = P;

                if mixac == 0;
                    Base(dd).Cells(j).Pc = P;
                    Base(dd).Cells(j).DDI = DDI;

                elseif mixac == 0.5;
                    Base(dd).Cells(j).Phm = P;
                    Base(dd).Cells(j).DDIhm = DDI;

                end

            end
        
        end
        
        
        
    end

    save('CuratedCells.mat','Base');
end
