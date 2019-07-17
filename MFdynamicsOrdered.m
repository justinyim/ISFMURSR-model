function xdot = MFdynamicsOrdered(t,X,M,F,v,varargin)
% Justin Yim - 2018

if nargin >= 6
    f_concatenated = varargin{1};
    concatenated = f_concatenated(t,X);
else
    concatenated = [];
end

Mm = M(t,X);
Fm = F(t,X);

acc = Mm\Fm;
vel = v(t,X);

xdot = [vel;acc(1:(length(X)-length(vel)-length(concatenated)));concatenated];

end