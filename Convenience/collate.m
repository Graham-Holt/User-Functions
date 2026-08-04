function V = collate(varargin)
% collate(v1,v2,...) collates same-size vectors.
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% collate(v1,v2,...)
% collate(A1,A2,...)
% V = collate(___)
% 
%% Description
% collate(v1,v2,...) returns a row vector of the collated elements of the
% input vectors
%
% collate(A1,A2,...) takes in arrays as input and returns a new array with
% the columns collated

n = length(varargin); m = size(varargin{1},2);
for k = 1:m
    for j = 1:n
        V(:,n*(k-1) + j) = varargin{j}(:,k);
    end
end