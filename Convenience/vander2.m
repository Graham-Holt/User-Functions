function V = vander2(v,nCols,d)
% vander2(v) computes Vandermonde matrices.
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% vander2(v)
% vander2(v,nCols)
% V = vander2(___)
% 
%% Description
% vander2(v) returns a square Vandermonde matrix with the values in "v"
%
% vander2(v,nCols) returns a rectangular matrix with "nCol" columnds
% 
% vander2(___,d) applies the given order of derivative to the Vandermonde
% matrix.

if nargin<2
    nCols = length(v);
end
if isempty(nCols)
    nCols = length(v);
end
if nargin<3
    d = zeros(1,nCols);
end

del = nCols - (1:nCols);
for k = 1:nCols
    % Create shifted Vandemonde matrix
    V(k,:) = [v(k).^del((1+d(k)):nCols) zeros(1,d(k))];

    % Multiply coefficients using power rule
    for j = 1:d(k)
        V(k,:) = V(k,:).*circshift(del,1-j);
    end
end