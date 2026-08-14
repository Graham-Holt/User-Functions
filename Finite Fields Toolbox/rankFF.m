function [r,isFullRank] = rankFF(A,modulo)

r = size(A,1) - sum(ismember(rrefFF(A,modulo),zeros(1,size(A,2)),'rows'));

isFullRank = r == min(size(A));