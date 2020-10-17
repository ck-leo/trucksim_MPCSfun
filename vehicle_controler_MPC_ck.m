%% ����˵��
%%��������������������ԭ��������������ٶ�v=1m/s����һ��ֱ�߹켣y=2������ʱ��50ms��
%%�ٻ��ڶ����ɶ��˶�ѧģ�ͣ���������Ŀ����Ƽ�����Լ��Ϊ�ͳ���״��
%%�ڸ߳���ʱӦ�ı����Ŀ����Ƽ�����Լ����������Ч�����ͣ�Ӧ���ö���ѧģ��
%%�۳���ϵͳΪ������ϵͳ���������Ի����
% ״̬�ռ䣺[dX/dt,dY/dt,dphi/dt]'=[v*cos(phi),v*sin(phi),v/L*tan(kesif)]';
% discretization:[dX/dt  ] [0  0 -vr*sin(phir)] [X  ] 
%   (Akesi*kesi) [dY/dt  ]=[0  0  vr*��(phir)]*[Y  ]+         
%                [dphi/dt] [0  0      0       ] [phi] 
%                [cos(phir)              0            ]         
%    (Bkesi*u)   [sin(phir)              0            ]*[v    ]+
%                [tan(kesifr)/L vr/(L*(cos(kesifr))^2)] [kesif] 
%                [vr*phir*sin(phir)           ]
%    (dis)       [-vr*phir*cos(phir)          ]                 (d)
%                [-vr*phir/(L*(cos(kesifr))^2)]
%% ��ʼֵ�趨(�켣�㣬�������ڣ�Np��Nc)
clc;clear;close;
tic;
N_sim=100;                     %% �ο�/����켣������
Delta_t=0.05;                  %% ��������
Np=5;                          %% Ԥ��ʱ��
Nc=4;                          %% ����ʱ��
%%�趨����Ŀ��
Y_goal=zeros(3,N_sim+Np);
% for j=1:N_sim+Np
%     Y_goal(1,j)=j*Delta_t*1;   %% X����Ŀ��ֵ
%     Y_goal(2,j)=4;             %% Y����Ŀ��ֵ
%     Y_goal(3,j)=0;             %% �����/��ڽ�Ŀ��ֵ
% end 
for j=1:N_sim+Np
    Y_goal(1,j)=2*cos(j*Delta_t-pi/2);   %% X����Ŀ��ֵ
    Y_goal(2,j)=2*sin(j*Delta_t-pi/2)+2;             %% Y����Ŀ��ֵ
    Y_goal(3,j)=-cot(j*Delta_t-pi/2);             %% �����/��ڽ�Ŀ��ֵ
end
% ���Ի�״̬�ռ����ĳ�ʼֵ
l = 2.5;                         %% ���
kesi0 = [0;0;pi/3];                 %% ״̬������ʼֵ
u0 = [1;0];                       %% ���Ʊ�����ʼֵ:���ٺ�ת���
% ��������
dis0 = kesi0;                        %% ���ű�����ʼֵ
dis0(1,1) = kesi0(3,1)*sin(kesi0(3,1))*u0(1,1);
dis0(2,1) = -kesi0(3,1)*cos(kesi0(3,1))*u0(1,1);
dis0(3,1) = -u0(1,1)*u0(2,1)/l/(cos(u0(2,1)/l))^2;
Delta_dis = zeros(3,1);          %% ������������(kʱ����k-1ʱ�̵Ĳ�)����
y0 = kesi0;                          %% ���������ʼֵ
Xf = zeros(6,1);                 %% �����;����״̬�����ĳ�ʼֵ(delta_x��y)
phi0 = kesi0(3,1);                  %% ����ǳ�ʼֵ
v_r0 = u0(1,1);                   %% ���ٵĳ�ʼֵ
delta_f0 = u0(2,1);               %% ǰ��ƫ��(ת���)�ĳ�ʼֵ
uu = zeros(2,N_sim);             %% ȫʱ��Χ�ڵĿ���������
u1 = zeros(2,N_sim);             %% ȫʱ��Χ�ڵĿ�����
y1 = zeros(3,N_sim);             %% ȫʱ��Χ�ڵ������

%% ��ʼ����
for k=1:N_sim
%% ״̬�ռ����
% ���Ի�_����ʱ��ϵͳ-̩��չ��(ÿһʱ�̸��¼���)
% dkesi/dt = Akesi*kesi+Bkesi*u+d
% y = Ckesi*kesi
Akesi=[0 0 -v_r0*sin(phi0);0 0 v_r0*cos(phi0);0 0 0];                              %% ϵͳ����
Bkesi=[cos(phi0) 0;sin(phi0) 0;tan(delta_f0)/l v_r0/l/((cos(delta_f0))^2)];      %% ���ƾ���
Bd=eye(3,3);                                                                     %% ���ž���
Ckesi=[1 0 0;0 1 0;0 0 1];                                                         %% �������
%  ������ɢ��
Ackesi = eye(3,3)+Akesi*Delta_t;
Bckesi=Bkesi*Delta_t;
Bdkesi=Bd*Delta_t;
Cckesi=Ckesi;
% ����matlab����ɢ������
% [Ackesi,Bckesi,Cckesi,Dckesi] = c2dm(Akesi,Bkesi,Ckesi,Dkesi,Delta_t);
% ����������ʽ
[ny,nkesi]=size(Cckesi);               %% nyΪ�������������nkesiΪ״̬��������
[nkesi,nu]=size(Bckesi);            %% nkesiΪ״̬����������nuΪ���Ʊ�������
[nkesi,nd]=size(Bdkesi);           %% nkesiΪ״̬����������ndΪ���ű�������
%% Ԥ��ģ��
%%���Ԥ��ģ�ͣ�Y=F*X+Phi*DeltaU+Phid*DeltaDis��
%%���ۺ�����    J=(Y_goal-Y)'*Q*(Y_goal-Y)+DeltaU'*R*DeltaU
%%���F��Phi��Phid��
[F,Phi,Phid,A_e, B_e,C_e] = mpcgain(Ackesi,Bckesi,Bdkesi,Cckesi,Nc,Np);
%����Ŀ�����Ȩ��Q
q=[1 0 0;0 5 0;0 0 0.5];
Q=zeros(Np*ny,Np*ny);
for j=1:Np
    Q((j-1)*ny+1:j*ny,(j-1)*ny+1:j*ny)=q;
end
 %���������Ȩ��R
r=[1 0;0 1]*0.1;
R=eye(Nc*nu,Nc*nu);   
for j=1:Nc
    R((j-1)*nu+1:j*nu,(j-1)*nu+1:j*nu)=r;
end
%�����������Ŀ�꣺ת��Ϊ��ǰʱ��Ԥ��ʱ���ڵ�����Ŀ��
Yref=zeros(Np*ny,1);
for j=1:Np       
    Yref((j-1)*ny+1:j*ny,1)=Y_goal(:,k+j-1);
end
%% ת��Ϊ��Լ�������Ͳ����
% J=1/2*DeltaU'*H*DeltaU+DeltaU'*G; A_cons*DeltaU��b;
% �����Լ��
umin=[0.5;-0.64];     %% ��С����Ϊ0.5m/s����Сת���Ϊ-0.64
umax=[1.5;0.64];      %% �����Ϊ1.5m/s�����ת���Ϊ0.64
[H,G,A_cons,b] = ToQP(F,Phi,Phid,Q,R,nu,Nc,Delta_dis,Xf,Yref,umin,umax,u0);
% ���
DeltaU = QPhild(H,G,A_cons,b);
%% �����Ż�
deltau=DeltaU(1:nu);    %% �����Ż���ȡ��һ������������Ϊkʱ�̵Ŀ�������
u0=u0+deltau;       %% �����Ż�������kʱ�̵Ŀ�����
uu(:,k)=deltau;   %% ����kʱ�̵Ŀ���������
u1(:,k)=u0;        %% ����kʱ�̵Ŀ�����
y1(:,k)=y0;        %% ����kʱ�̵������
kesi0_old=kesi0;        %% ����kʱ��״̬����
dis_old=dis0;      %% ������һʱ�̸�����
kesi0=Ackesi*kesi0+Bckesi*u0+Bdkesi*dis0;   %% ����k+1ʱ��״̬��
y0=Cckesi*kesi0;                  %% ����k+1ʱ�������
Xf=[kesi0-kesi0_old;y0];         %% ����k+1ʱ��״̬��(����������ʽ)
% ״̬�ռ�����������
% ����k+1ʱ�̵�״̬��
phi0=kesi0(3,1);
v_r0=u0(1,1);
delta_f0=u0(2,1);
dis0(1,1)=kesi0(3,1)*sin(kesi0(3,1))*u0(1,1);
dis0(2,1)=-kesi0(3,1)*cos(kesi0(3,1))*u0(1,1);
dis0(3,1)=-u0(1,1)*u0(2,1)/l/(cos(u0(2,1)/l))^2;
Delta_dis=dis0-dis_old;             % ����k+1ʱ�̵�Delata_dis
end
toc;
%% ��ͼ
k=0:(N_sim-1);
figure
subplot(411)
plot(y1(1,:),y1(2,:),'k-','LineWidth',1)
hold on
plot(Y_goal(1,:),Y_goal(2,:),'r*','LineWidth',0.5)
xlabel('����λ��X(m)')
ylabel('����λ��Y(m)')
legend('������ʻ�켣','�����ο��켣')
subplot(412)
plot(k*Delta_t,y1(3,:),'k-','LineWidth',1)
ylabel('��ڽ�')
xlabel('ʱ��(s)')
subplot(413)
plot(k*Delta_t,u1(1,:),'k-','LineWidth',1)
xlabel('ʱ��(s)')
ylabel('����(m/s)')
subplot(414)
plot(k*Delta_t,u1(2,:),'k-','LineWidth',1)
xlabel('ʱ��(s)')
ylabel('ǰ��ƫ��')

