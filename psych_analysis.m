function out = psych_analysis();
    
    mypath = mfilename('fullpath');
    basepath = mypath(1:end-14);
    inDir = [basepath,'halfmatched_psych_data/'];

    MyNames = dir([inDir,'*.csv']);
    for fname = 1:length(MyNames);
        MyFiles{fname} = MyNames(fname).name;
    end
    NumFiles = length(MyNames);% NumFiles = NumFiles(1);
    CustomDisp = [-0.48, -0.03, 0.03, 0.48];
    MergeFiles = {
        };

    Exclude = [];
    Merged = [];
    Count = NumFiles - length(Exclude) - length(MergeFiles);
    IncludeCount = 0;

    AnswerMatrix = zeros(4,4,3);

    Colours = {'Red';'Green';'Blue';'Black';[1,0,1]; 'Cyan'; [1,1,0]};

    for file = 1:NumFiles;
        if ~any(Merged == file) && ~any(Exclude == file);
            merge = 0;
            for merger = 1:length(MergeFiles);
                if any(file == MergeFiles{merger});
                    merge = 1;
                    mergewith = setdiff(MergeFiles{merger},file);
                end
            end

            if merge
                CurrentData = [];
                FilesToMerge = [file,mergewith];
                for mergefile = FilesToMerge;
                    CurrentFile = fopen([inDir,MyFiles{mergefile}],'r');
                    ReadData = cell2mat(textscan(CurrentFile, '%f%f%f%f%f', 'Delimiter',',', 'HeaderLines',2));
                    CurrentData = [CurrentData; ReadData];
                    Merged = [Merged,mergefile];
                end

            else

                CurrentFile = fopen([inDir,MyFiles{file}],'r');
                CurrentData = cell2mat(textscan(CurrentFile, '%f%f%f%f%f', 'Delimiter',',', 'HeaderLines',2));
            end


            Disparities = unique(CurrentData(:,1));
            RefreshRates = unique(CurrentData(:,4));
            DotMatchLevels = unique(CurrentData(:,5));


            for disp = 1:length(Disparities);
                for RF = 1:length(RefreshRates);
                    for DM = 1:length(DotMatchLevels);
                        CurrentIndices = ...
                            (CurrentData(:,1) == Disparities(disp)) .* ...
                            (CurrentData(:,4) == RefreshRates(RF)) .* ...
                            (CurrentData(:,5) == DotMatchLevels(DM));

                        CurrentIndices = logical(CurrentIndices);

                        CorrectAns = 1 + (Disparities(disp) > 0)*2;

                        Percent = sum(CurrentData(CurrentIndices,2) == CorrectAns)/sum(CurrentIndices);

                        AnswerMatrix(disp,RF,DM,IncludeCount+1) = Percent;
                    end
                end
            end
            IncludeCount = IncludeCount +1;

            
            initials = MyFiles{file}; initials = initials(1:2);
            fprintf('Initials: %s\n',initials)
            

        end
    end

  
    out = squeeze(mean(mean(AnswerMatrix(:,:,2,:))));
    
    
    
    
    
end