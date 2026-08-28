function V = cvander(v,nCols,d)
% cvander(v) computes confluent Vandermonde matrices.
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% cvander(v)
% cvander(v,nCols)
% V = cvander(___)
% 
%% Description
% cvander(v) returns a square Vandermonde matrix with the values in "v"
%
% cvander(v,nCols) returns a rectangular matrix with "nCol" columns
% 
% cvander(___,d) applies the given order of derivative to the corresponding
% row of the rectangular Vandermonde matrix.

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

end