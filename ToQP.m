function [H,G,A_cons,b] = ToQP(varargin)
%% 转化为带约束二次型并求解
% [H,G,A_cons,b] = ToQP(F,Phi,Phid,Q,R,nu,Nc,Delta_dis,Xf,Yref,umin,umax,u0,deltaUmin,deltaUmax,Ymin,Ymax)
% min J=1/2*DeltaU'*H*DeltaU+DeltaU'*G; A_cons*DeltaU≤b;
% narginchk(12,16);
F = varargin{1};
Phi = varargin{2};
Phid = varargin{3};
Q = varargin{4};
R = varargin{5};
nu = varargin{6};
Nc = varargin{7};
Delta_dis = varargin{8};
Xf = varargin{9};
Yref = varargin{10};
umin = varargin{11};
umax = varargin{12};
u0 = varargin{13};
% H,G
H=Phi'*Q*Phi + R;
M = F*Xf+Phid*Delta_dis-Yref;
G = Phi'*Q*M;
%%约束矩阵
C1 = zeros(nu*Nc,nu);
C2 = zeros(nu*Nc,nu*Nc);
for i = 1:Nc
    C1(((i-1)*nu+1):(i*nu),:) = eye(nu,nu);
end
for i = 1:Nc
    C2(((i-1)*nu+1):(nu*Nc),((i-1)*nu+1):(i*nu)) = C1(((i-1)*nu+1):(nu*Nc),:);
end
A_cons=[-C2;C2];
Umax=zeros(nu*Nc,1);
Umin=zeros(nu*Nc,1);
for i = 1:Nc
    Umin(((i-1)*nu+1):(i*nu),:) = umin;
    Umax(((i-1)*nu+1):(i*nu),:) = umax;
end
b = [-Umin+C1*u0;Umax-C1*u0];
end