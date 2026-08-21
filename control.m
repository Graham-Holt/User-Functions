function varargout = control(varargin)

if nargin < 6
    error('Too few input arguments.');
elseif nargin > 6
    error('Too many input arguments.');
end

A = varargin{1};
B = varargin{2};
C = varargin{3};
D = varargin{4};

[n,m] = size(B);
p = size(C,1);

eig = varargin{5};
method = lower(char(varargin{6}));

if isscalar(method)
    method = [method, 'x'];
else
    method = method(1:2);
end

switch method(1)
case 'x'
    K = placeX(A,B,eig(1:n));
    L = placeX(A',C',eig((n+1):end)).';
case 'y'
    K = placeY(A,B,C,D,eig);
    L = [];

    Kinf = K(1:(n-p+1),1:(n-m+1)) - K(1:(n-p+1),(n-m+2):end)/K((n-p+2):end,(n-m+2):end)*K((n-p+2):end,1:(n-m+1));
otherwise
    error('Improper method chosen.');
end

if nargout > 2
    
switch method
case 'xx'
    Kr = -(A - B*K)\B;
case 'xy'
    Kr = D - (C + D*K)/(A - B*K)*B;
case 'yx'
    Kr = -(A - B*Kinf/(eye(p) + D*Kinf)*C)\B*(eye(m) + Kinf/(eye(p) + D*Kinf)*D);
case 'yy'
    Kr = -(A - B*Kinf/(eye(p) + D*Kinf)*C)\B*(eye(m) + Kinf/(eye(p) + D*Kinf)*D);
    Kr = (eye(p) + D*Kinf)\(D + C*Kr);
otherwise
    error('Improper method chosen.');
end
Kr = pinv(Kr);

end

switch nargout
case {0,1}
    varargout = {K};
case 2
    varargout = {K,L};
case 3
    varargout = {K,L,Kr};
otherwise
    error('Invalid number of outputs.');
end