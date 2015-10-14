function [ MASCc, sepvalc ] = cutoffcorrelation( MASC, sepval )
% This function returns an abridged correlation function and associated
% distnace vector that are truncated after the correlation (MASC) crosses
% zero the nth time

count = 0; 
n = 2;

% initialize cutoff just incase MASC never crosses zero
cutoff = length(MASC);

for i=1:length(MASC)
    if(MASC(i) < 0 ) 
        count = count + 1;
    end;
    if(count >= n)
        cutoff = i; 
        break;
    end;

end;

MASCc = MASC(1: cutoff); % the cut off structure function
sepvalc = sepval(1:cutoff); %the cut off sepval

end

