function string = sandwich_string(input,start,stop);
    % Usage: myString = sandwich_string(input,start,stop)
    % This will return the first case of a string sandwiched
    % between start and stop, where start and stop are both strings
    % E.g., 
    % S = 'the quick brown fox jumped over the lazy dog';
    % myString = sandwich_string(S,'fox ', ' dog');
    % Will return 'jumped over the lazy'
    
    nStart = length(start);
    nStop = length(stop);
    
    for j = 1:length(input)-nStop;
        current = input(j:j+nStart-1);
        
        if strcmp(current,start);
            
            for k = j:length(input)-nStop;
                current = input(k:k+nStop-1);
                
                if strcmp(current,stop);
                    break
                end
            end
            
            break
        end
        
    end
    
    if j ~= length(input);
        string = input(j+nStart:k-1);
        
    else
        string = '';
    end

end