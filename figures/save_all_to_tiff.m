

c = dir;

for j = 1:length(c);
    
    fname = c(j).name;    
    if any(strfind(fname,'.fig'));
        
        fprintf('Found figure: %s\n',fname);
        
        uiopen(fname,1);
        save_to_tiff(gcf);
    end
end