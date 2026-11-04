function P = numIntPartition(N)
% numIntPartition(n) computes the number of integer partitions of "n"
% 
% Graham Holt, October 2026. Updated October 2026.
% Embry-Riddle Aeronautical University
% 
%% Syntax
% numIntPartition(n)
% P = numIntPartition(n)
% 
%% Description
% numIntPartition(n) returns the number of partitions of the integer (N)
% without determining each partition individually

k = 1; p = 1;
while length(p) < N
    n = ((-1)^k)*ceil(0.5*(k+1));
    m = ((-1)^(k-1))*ceil(0.5*k);
    
    n = 0.5*n*(3*n - 1);
    m = 0.5*m*(3*m - 1);

    p = [p; zeros(n-m-1,1); (-1)^(0.5*k*(k+3))];

    k = k + 1;
end

P = 1;
for k = 1:N
    P = [P*p(1:k), P];
end
P = P(1);