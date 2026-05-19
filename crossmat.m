function A = crossmat(a)
% crossmat(a) converts between vectors and skew-symmetric matrices
% 
% Graham Holt, Mar. 2026. Updated May 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% crossmat(a)
% A = crossmat(a)
% 
%% Description
% crossmat(a) takes in a 3D array of column vectors or skew-symmetric 
% matrices and returns either a 3D array of skew-symmetric cross-product 
% matrices or a 2D array of column vectors


if size(a,1)~=3
    error('Input must be an 3x1xd or 3x3xd array');
end

if size(a,2) == 1
    for k = 1:size(a,3)
        A(:,:,k) = [0 -a(3,k) a(2,k); a(3,k) 0 -a(1,k); -a(2,k) a(1,k) 0];
    end
elseif size(a,2) == 3
    for k = 1:size(a,3)
        A(:,k) = [a(3,2,k); a(1,3,k); a(2,1,k)];
    end
else
    error('Input must be an 3x1xd or 3x3xd array');
end