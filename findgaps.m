function [broken_ind,start_ind,end_ind]=findgaps(input_ind)
% finds gaps (greater than 0) in between indeces in a vector. If there are
% no gaps then findgaps returns results as if input was the only gap found
% i.e. broken_ind = input_ind & start_ind/end_ind are the first and last
% indeces of the input_ind. 
% rechecked for bugs by SDK on 1/5/2017
% 
% Inputs:
%   1) input_ind: indexes (not logicals) of desired condition
%   
% Outputs:
%   1) broken_ind: a matrix of all indexes for each break by row. Since
%   rectangular index, non-values indexes are 0 (sorry not NaN for
%   backwards comptability)
%   2) start_ind: start indeces of each break
%   3) end_ind: end indeces of each break


if isempty(input_ind)
    broken_ind = [];
    start_ind = [];
    end_ind = [];
    return
end

gaps =find(abs(diff(input_ind)) > 1);
broken_ind = zeros(length(gaps),50);
start_ind = NaN(1,length(gaps));
end_ind = NaN(1,length(gaps));
if ~isempty(gaps)
    for gapind = 1:length(gaps)+1
        if gapind == 1
            temp = input_ind(1:gaps(gapind));
        elseif gapind == length(gaps)+1
            temp = input_ind(gaps(gapind-1)+1:end);
        else
            temp = input_ind(gaps(gapind-1)+1:gaps(gapind));
        end
        broken_ind(gapind,1:length(temp)) = temp;
        start_ind(gapind) = temp(1);
        end_ind(gapind) = temp(end);
    end
else
    start_ind = input_ind(1);
    end_ind = input_ind(end);
    if size(input_ind,1) > 1%row vector
        broken_ind = input_ind';
    else
        broken_ind = input_ind;
    end
end

end