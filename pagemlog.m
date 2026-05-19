function L = pagemlog(M)
% pagemlog(M) performs matrix logarithm across array pages
% 
% Graham Holt, March 2026. Updated May 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% pagemlog(M)
% L = pagemlog(M)
% 
%% Description
% pagemlog(M) takes in a 3D array of matrices and performs the matrix
% logarithm across each page of the input

L = zeros(size(M));
for k = 1:size(M,3)
    L(:,:,k) = logm(M(:,:,k));
end