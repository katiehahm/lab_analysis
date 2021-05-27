function [index] = findTindex(time,timearr)
%5/26/21
% returns the index number of the first value in time array larger than given time value

index = (find(timearr > time,1));

end

