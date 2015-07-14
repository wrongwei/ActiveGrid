% make time sequence with logarithmically increasing spacing
%makelogtime.m
% tlog = makelogtime(begin, end, Npoints);
%Greg Bewley 01/08/04
 
function tlog = makelogtime(abegin, aend, Npoints)
 
% all time points, from 0 to 1:
scl =  ( (0:(Npoints-1)) / (Npoints-1) );
 
% scale the time points to cover the desired range on the log scale:
scl = (log10(aend)-log10(abegin)) * scl + log10(abegin);
 
% recover linear coordinates:
tlog = 10.^scl;