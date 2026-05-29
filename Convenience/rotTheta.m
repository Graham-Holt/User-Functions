function R = rotTheta(theta,u)
% rotTheta(theta,u) generates a rotation matrix according to the Principle 
% Axis Theorem
% 
% Graham Holt, March 2026. Updated May 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% rotTheta(theta,u)
% R = rotTheta(theta,u)
% 
%% Description
% rotTheta(theta,u) takes in angles "theta" (in radians) and 3D vectors "u" 
% and returns the matrix which rotates a 3D vector about "u" by "theta" in
% sequence.

if size(u,1)~=3
    error('Input must be columnm vector, or array of column vectors, with dimension 3');
end

R = eye(3); u = u./vecnorm(u,2,1);
for k = 1:length(theta)
    R = expm(theta(k)*crossmat(u(:,k)))*R;
end