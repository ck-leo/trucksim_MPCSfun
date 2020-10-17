%% 程序说明
%%问题描述：车辆从坐标原点出发，以期望速度v=1m/s跟踪一条直线轨迹y=2。采样时间50ms。
%%①基于二自由度运动学模型，下述跟踪目标设计及车速约束为低车速状况
%%②高车速时应改变跟踪目标设计及车速约束，但跟踪效果降低，应采用动力学模型
%%③车辆系统为非线性系统，采用线性化求解
% 状态空间：[dX/dt,dY/dt,dphi/dt]'=[v*cos(phi),v*sin(phi),v/L*tan(kesif)]';
% discretization:[dX/dt  ] [0  0 -vr*sin(phir)] [X  ] 
%   (Akesi*kesi) [dY/dt  ]=[0  0  vr*。(phir)]*[Y  ]+         
%                [dphi/dt] [0  0      0       ] [phi] 
%                [cos(phir)              0            ]         
%    (Bkesi*u)   [sin(phir)              0            ]*[v    ]+
%                [tan(kesifr)/L vr/(L*(cos(kesifr))^2)] [kesif] 
%                [vr*phir*sin(phir)           ]
%    (dis)       [-vr*phir*cos(phir)          ]                 (d)
%                [-vr*phir/(L*(cos(kesifr))^2)]
%% 初始值设定(轨迹点，采样周期，Np，Nc)
clc;clear;close;
tic;
N_sim=100;                     %% 参考/仿真轨迹点数量
Delta_t=0.05;                  %% 采样周期
Np=5;                          %% 预测时域
Nc=4;                          %% 控制时域
%%设定跟踪目标
Y_goal=zeros(3,N_sim+Np);
% for j=1:N_sim+Np
%     Y_goal(1,j)=j*Delta_t*1;   %% X方向目标值
%     Y_goal(2,j)=4;             %% Y方向目标值
%     Y_goal(3,j)=0;             %% 航向角/横摆角目标值
% end 
for j=1:N_sim+Np
    Y_goal(1,j)=2*cos(j*Delta_t-pi/2);   %% X方向目标值
    Y_goal(2,j)=2*sin(j*Delta_t-pi/2)+2;             %% Y方向目标值
    Y_goal(3,j)=-cot(j*Delta_t-pi/2);             %% 航向角/横摆角目标值
end
% 线性化状态空间矩阵的初始值
l = 2.5;                         %% 轴距
kesi0 = [0;0;pi/3];                 %% 状态变量初始值
u0 = [1;0];                       %% 控制变量初始值:车速和转向角
% 干扰因素
dis0 = kesi0;                        %% 干扰变量初始值
dis0(1,1) = kesi0(3,1)*sin(kesi0(3,1))*u0(1,1);
dis0(2,1) = -kesi0(3,1)*cos(kesi0(3,1))*u0(1,1);
dis0(3,1) = -u0(1,1)*u0(2,1)/l/(cos(u0(2,1)/l))^2;
Delta_dis = zeros(3,1);          %% 干扰量的增量(k时刻与k-1时刻的差)向量
y0 = kesi0;                          %% 输出变量初始值
Xf = zeros(6,1);                 %% 增量型矩阵的状态变量的初始值(delta_x与y)
phi0 = kesi0(3,1);                  %% 航向角初始值
v_r0 = u0(1,1);                   %% 车速的初始值
delta_f0 = u0(2,1);               %% 前轮偏角(转向角)的初始值
uu = zeros(2,N_sim);             %% 全时范围内的控制量增量
u1 = zeros(2,N_sim);             %% 全时范围内的控制量
y1 = zeros(3,N_sim);             %% 全时范围内的输出量

%% 开始仿真
for k=1:N_sim
%% 状态空间矩阵
% 线性化_连续时间系统-泰勒展开(每一时刻更新计算)
% dkesi/dt = Akesi*kesi+Bkesi*u+d
% y = Ckesi*kesi
Akesi=[0 0 -v_r0*sin(phi0);0 0 v_r0*cos(phi0);0 0 0];                              %% 系统矩阵
Bkesi=[cos(phi0) 0;sin(phi0) 0;tan(delta_f0)/l v_r0/l/((cos(delta_f0))^2)];      %% 控制矩阵
Bd=eye(3,3);                                                                     %% 干扰矩阵
Ckesi=[1 0 0;0 1 0;0 0 1];                                                         %% 输出矩阵
%  近似离散化
Ackesi = eye(3,3)+Akesi*Delta_t;
Bckesi=Bkesi*Delta_t;
Bdkesi=Bd*Delta_t;
Cckesi=Ckesi;
% 采用matlab的离散化函数
% [Ackesi,Bckesi,Cckesi,Dckesi] = c2dm(Akesi,Bkesi,Ckesi,Dkesi,Delta_t);
% 增量矩阵形式
[ny,nkesi]=size(Cckesi);               %% ny为输出变量个数，nkesi为状态变量个数
[nkesi,nu]=size(Bckesi);            %% nkesi为状态变量个数，nu为控制变量个数
[nkesi,nd]=size(Bdkesi);           %% nkesi为状态变量个数，nd为干扰变量个数
%% 预测模型
%%输出预测模型：Y=F*X+Phi*DeltaU+Phid*DeltaDis；
%%代价函数：    J=(Y_goal-Y)'*Q*(Y_goal-Y)+DeltaU'*R*DeltaU
%%求解F、Phi、Phid；
[F,Phi,Phid,A_e, B_e,C_e] = mpcgain(Ackesi,Bckesi,Bdkesi,Cckesi,Nc,Np);
%计算目标跟踪权重Q
q=[1 0 0;0 5 0;0 0 0.5];
Q=zeros(Np*ny,Np*ny);
for j=1:Np
    Q((j-1)*ny+1:j*ny,(j-1)*ny+1:j*ny)=q;
end
 %计算控制量权重R
r=[1 0;0 1]*0.1;
R=eye(Nc*nu,Nc*nu);   
for j=1:Nc
    R((j-1)*nu+1:j*nu,(j-1)*nu+1:j*nu)=r;
end
%计算期望输出目标：转化为当前时刻预测时域内的期望目标
Yref=zeros(Np*ny,1);
for j=1:Np       
    Yref((j-1)*ny+1:j*ny,1)=Y_goal(:,k+j-1);
end
%% 转化为带约束二次型并求解
% J=1/2*DeltaU'*H*DeltaU+DeltaU'*G; A_cons*DeltaU≤b;
% 输入的约束
umin=[0.5;-0.64];     %% 最小车速为0.5m/s，最小转向角为-0.64
umax=[1.5;0.64];      %% 最大车速为1.5m/s，最大转向角为0.64
[H,G,A_cons,b] = ToQP(F,Phi,Phid,Q,R,nu,Nc,Delta_dis,Xf,Yref,umin,umax,u0);
% 求解
DeltaU = QPhild(H,G,A_cons,b);
%% 滚动优化
deltau=DeltaU(1:nu);    %% 滚动优化：取第一个控制量增量为k时刻的控制增量
u0=u0+deltau;       %% 滚动优化：计算k时刻的控制量
uu(:,k)=deltau;   %% 保存k时刻的控制量增量
u1(:,k)=u0;        %% 保存k时刻的控制量
y1(:,k)=y0;        %% 保存k时刻的输出量
kesi0_old=kesi0;        %% 更新k时刻状态变量
dis_old=dis0;      %% 更新上一时刻干扰量
kesi0=Ackesi*kesi0+Bckesi*u0+Bdkesi*dis0;   %% 计算k+1时刻状态量
y0=Cckesi*kesi0;                  %% 计算k+1时刻输出量
Xf=[kesi0-kesi0_old;y0];         %% 计算k+1时刻状态量(增量矩阵形式)
% 状态空间矩阵迭代更新
% 计算k+1时刻的状态量
phi0=kesi0(3,1);
v_r0=u0(1,1);
delta_f0=u0(2,1);
dis0(1,1)=kesi0(3,1)*sin(kesi0(3,1))*u0(1,1);
dis0(2,1)=-kesi0(3,1)*cos(kesi0(3,1))*u0(1,1);
dis0(3,1)=-u0(1,1)*u0(2,1)/l/(cos(u0(2,1)/l))^2;
Delta_dis=dis0-dis_old;             % 计算k+1时刻的Delata_dis
end
toc;
%% 作图
k=0:(N_sim-1);
figure
subplot(411)
plot(y1(1,:),y1(2,:),'k-','LineWidth',1)
hold on
plot(Y_goal(1,:),Y_goal(2,:),'r*','LineWidth',0.5)
xlabel('横向位置X(m)')
ylabel('纵向位置Y(m)')
legend('汽车行驶轨迹','汽车参考轨迹')
subplot(412)
plot(k*Delta_t,y1(3,:),'k-','LineWidth',1)
ylabel('横摆角')
xlabel('时间(s)')
subplot(413)
plot(k*Delta_t,u1(1,:),'k-','LineWidth',1)
xlabel('时间(s)')
ylabel('车速(m/s)')
subplot(414)
plot(k*Delta_t,u1(2,:),'k-','LineWidth',1)
xlabel('时间(s)')
ylabel('前轮偏角')

