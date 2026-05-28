function D = pagediag(V)
% pagediag(V) performs the diag() function along array pages
% 
% Graham Holt, May 2026. Updated May 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% pagediag(V)
% D = crossmat(V)
% 
%% Description
% pagediag(V) takes in a 2D array of column vectors or 3D array of matrices 
% and returns either a 3D array of diagonal matrices according to the
% column vectors or a 2D matrices of column vectors with the diagonal 
% elements of the input matrices

switch ndims(V)
    case 2
        for k = 1:size(V,2)
            D(:,:,k) = diag(V(:,k));
        end
    case 3
        for k = 1:size(V,3)
            D(:,k) = diag(V(:,:,k));
        end
    otherwise
    error('Improper size for input, must be 2D array of column vectors or 3D array of matrices.');
end