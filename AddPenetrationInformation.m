function newBase = AddPenetrationInformation(SubBase);
    count = 0;
    lemExpts = {};
    for j = 1:length(SubBase.Cells);
        fname = SubBase.Cells(j).filename;
        %if strcmp(fname(9:11),'lem');
        count = count+1;
        lemExpts{count} = fname;
        %end
    end
    
    uniqueExpts = unique(lemExpts);

    % Skips lemM116
    XYs = zeros(2,length(uniqueExpts)-1);
    
    for expt = 1:length(uniqueExpts);
        fname = uniqueExpts{expt};
        prefix = fname(1:25);
        uflName = [prefix,'10.ufl'];

        
        file= fopen(uflName,'r');
        
        % This is so that it will handle M116 which doesn't have separate
        % channels
        if file == -1;
            uflName = [prefix,'ufl'];
            file = fopen(uflName,'r');
        end
        
        contents = fscanf(file,'%s');
        
        x = str2double(sandwich_string(contents,'rf',','));
        y = str2double(sandwich_string(contents,',',':'));
        fx = str2double(sandwich_string(contents,'fx=',','));
        fy = str2double(sandwich_string(contents,'fy=',','));
        
        fx(isnan(fx)) = 0;
        fy(isnan(fy)) = 0;
        
        fclose(file);
        
        XYs(1,expt) = x-fx;
        XYs(2,expt) = y-fy;
        
    end
    
    
    SubBase.exptlist = uniqueExpts;
    SubBase.penetrationlist = XYs;
    
    
    newBase = SubBase;
end