function V = vanderd(v,nCols,d)
% vanderd(v) computes Vandermonde matrices.
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% vanderd(v)
% vanderd(v,nCols)
% V = vanderd(___)
% 
%% Description
% vanderd(v) returns a square Vandermonde matrix with the values in "v"
%
% vanderd(v,nCols) returns a rectangular matrix with "nCol" columnds
% 
% vanderd(___,d) applies the given order of derivative to the Vandermonde
% matrix.

if nargin < 2
    nCols = length(v);
end
if isempty(nCols)
    nCols = length(v);
end
if nargin < 3
    d = zeros(1,length(v));
end

del = nCols - (1:nCols);
for k = 1:length(v)
    V(k,:) = [v(k).^del((1+d(k)):nCols) zeros(1,d(k))];

    for j = 1:d(k)
        V(k,:) = V(k,:).*circshift(del,1-j);
    end
end