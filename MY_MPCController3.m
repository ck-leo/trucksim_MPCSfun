function [sys,x0,str,ts] = MY_MPCController3(t,x,u,flag)
% �ú�����д�ĵ�3��S������������MATLAB�汾��R2011a��
% �޶��ڳ����˶�ѧģ�ͣ�������Ϊ�ٶȺ�ǰ��ƫ�ǣ�ʹ�õ�QPΪ�°汾��QP�ⷨ
% [sys,x0,str,ts] = MY_MPCController3(t,x,u,flag)
% is an S-function implementing the MPC controller intended for use 
% with Simulink. The argument md, which is the only user supplied
% argument, contains the data structures needed by the controller. The
% input to the S-function block is a vector signal consisting of the 
% measured outputs and the reference values for the controlled 
% outputs. The output of the S-function block is a vector signal
% consisting of the control variables and the estimated state vector,
% potentially including estimated disturbance states.

% SFUNTMPL General MATLAB S-Function Template
%   With MATLAB S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an MATLAB S-function syntax is:
%       [SYS,X0,STR,TS,SIMSTATECOMPLIANCE] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG   RESULT             DESCRIPTION
%   -----  ------             --------------------------------------------
%   0      [SIZES,X0,STR,TS]  Initialization, return system sizes in SYS,
%                             initial state in X0, state ordering strings
%                             in STR, and sample times in TS.
%   1      DX                 Return continuous state derivatives in SYS.
%   2      DS                 Update discrete states SYS = X(n+1)
%   3      Y                  Return outputs in SYS.
%   4      TNEXT              Return next time hit for variable step sample
%                             time in SYS.
%   5                         Reserved for future (root finding).
%   9      []                 Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%      SYS(1) = Number of continuous states.
%      SYS(2) = Number of discrete states.
%      SYS(3) = Number of outputs.
%      SYS(4) = Number of inputs.
%               Any of the first four elements in SYS can be specified
%               as -1 indicating that they are dynamically sized. The
%               actual length for all other flags will be equal to the
%               length of the input, U.
%      SYS(5) = Reserved for root finding. Must be zero.
%      SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%               has direct feedthrough if U is used during the FLAG=3
%               call. Setting this to 0 is akin to making a promise that
%               U will not be used during FLAG=3. If you break the promise
%               then unpredictable results will occur.
%      SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%      X0     = Initial state conditions or [] if no states.
%
%      STR    = State ordering strings which is generally specified as [].
%
%      TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0      0,      : Continuous sample time.
%                     0      1,      : Continuous, but fixed in minor step
%                                      sample time.
%                     PERIOD OFFSET, : Discrete sample time where
%                                      PERIOD > 0 & OFFSET < PERIOD.
%                     -2     0];     : Variable step discrete sample time
%                                      where FLAG=4 is used to get time of
%                                      next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying more than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].
%
%      SIMSTATECOMPLIANCE = Specifices how to handle this block when saving and
%                           restoring the complete simulation state of the
%                           model. The allowed values are: 'DefaultSimState',
%                           'HasNoSimState' or 'DisallowSimState'. If this value
%                           is not speficified, then the block's compliance with
%                           simState feature is set to 'UknownSimState'.


%   Copyright 1990-2010 The MathWorks, Inc.

%
% The following outlines the general structure of an S-function.
%
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
%   case 1,
%     sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
%   case 4,
%     sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
%   case 9,
%     sys=mdlTerminate(t,x,u);
  %%%%%%%%%%%%%%%%
  % Unused flags %
  %%%%%%%%%%%%%%%%
  case {1,4,9},
    sys=[];
    
  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 3;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [0;0;0];

%
% str is always an empty matrix
%
str = [];
global U
U = [0;0];
%
% initialize the array of sample times
%
ts  = [0.1 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
% function sys=mdlDerivatives(t,x,u)
% 
% sys = [];

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u) 

sys = x;

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

global a b u_piao;
global X_ref Y_ref;
clc;
% 5km��ֱ��·
if isempty(X_ref)
    N_S = 1500;     % ��·����
    ds = 0.4;       % ·�߼��
    X_ref = (0:0.4:N_S)';
    N = length(X_ref);
    Y_ref = X_ref;
end

global U;
global kesi;
tic
Nx = 3;% ״̬������
Nu = 2;% �������ĸ���
Np = 60;% Ԥ�ⲽ��
Nc = 30;% ���Ʋ���
Row = 10;% �ɳ�����
fprintf('Update start, t=%6.3f\n',t)
% Trucksim���Ϊǰ�����ĵ�λ�ˣ���Ҫת��Ϊ��������

t_d = u(3)*pi/180; % Trucksim�����Ϊ�Ƕȣ��Ƕ�ת��Ϊ����

% %ֱ��·��
r = findproj(u,X_ref,Y_ref);
% r(1) = 5*t/sqrt(2);
% r(2) = 5*t/sqrt(2);
% r(3) = 0;
vd1 = 5;
vd2 = 0;
% �뾶Ϊ25m��Բ�ι켣���ٶ�Ϊ5m/s
% r(1) = 25*sin(0.2*t);
% r(2) = 25+10-25*cos(0.2*t);
% r(3) = 0.2*t;
% vd1 = 5;
% vd2 = 0.104;                % �ο�deltaf
% % �뾶Ϊ25m��Բ�ι켣���ٶ�Ϊ3m/s
% r(1) = 25*sin(0.12*t);
% r(2) = 25+10-25*cos(0.12*t);
% r(3) = 0.2*t;
% vd1 = 3;
% vd2 = 0.104;
% % �뾶Ϊ25m��Բ�ι켣���ٶ�Ϊ10m/s
% r(1) = 25*sin(0.4*t);
% r(2) = 25+10-25*cos(0.4*t);
% r(3) = 0.4*t;
% vd1 = 10;
% vd2 = 0.104;
% % �뾶Ϊ25m��Բ�ι켣���ٶ�Ϊ4m/s
% r(1) = 25*sin(0.16*t);
% r(2) = 25+10-25*cos(0.16*t);
% r(3) = 0.16*t;
% vd1 = 4;
% vd2 = 0.104;
kesi = zeros(Nx+Nu,1);
kesi(1) = u(1)-r(1); % u(1) == X(1)
kesi(2) = u(2)-r(2); % u(2) == X(2)
kesi(3) = t_d-r(3);  % u(3) == X(3)
kesi(4) = U(1);
kesi(5) = U(2);
fprintf('Update start,u(1)=%4.2f\n',U(1))
fprintf('Update start,u(2)=%4.2f\n',U(2))

T = 0.1;
T_all = 40; % ��ʱ�趨���ܵķ���ʱ�䣬��Ҫ�����Ƿ�ֹ���������켣Խ��
% Mobile Robot Parameters
L = 3.5;            % ���
% Mobile Robot variables

% �����ʼ��
u_piao = zeros(Nx,Nu);
Q = 100*eye(Nx*Np,Nx*Np);
R = 5*eye(Nu*Nc);
a  = [1 , 0 , -vd1*sin(t_d)*T;
      0 , 1 ,  vd1*cos(t_d)*T;
      0 , 0 ,         1;];
b = [cos(t_d)*T ,            0;
     sin(t_d)*T ,            0;
     tan(vd2)*T/L , vd1*T/((cos(vd2))^2);];
A_cell = cell(2,2);
B_cell = cell(2,1);
A_cell{1,1} = a;
A_cell{1,2} = b;
A_cell{2,1} = zeros(Nu,Nx);
A_cell{2,2} = eye(Nu);
B_cell{1,1} = b;
B_cell{2,1} = eye(Nu);
A = cell2mat(A_cell);
B = cell2mat(B_cell);
C = [1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;];
PHI_cell = cell(Np,1);
THETA_cell = cell(Np,Nc);
for j = 1:1:Np
    PHI_cell{j,1} = C*A^j;
    for k = 1:1:Nc
        if k<=j
            THETA_cell{j,k} = C*A^(j-k)*B;
        else
            THETA_cell{j,k} = zeros(Nx,Nu);
        end
    end
end
PHI = cell2mat(PHI_cell); % size(PHI) = [Nx*Np Nx+Np]
THETA = cell2mat(THETA_cell); % size(THETA) = [Nx*Np Nu*(Nc+1)]
H_cell = cell(2,2);
H_cell{1,1} = THETA'*Q*THETA + R;
H_cell{1,2} = zeros(Nu*Nc,1);
H_cell{2,1} = zeros(1,Nu*Nc);
H_cell{2,2} = Row;
H = cell2mat(H_cell);
H = (H+H')/2;               % H��������ǶԳƣ��Գƻ�

error = PHI*kesi;
f_cell = cell(1,2);
f_cell{1,1} = 2*error'*Q*THETA;
f_cell{1,2} = 0;
f = cell2mat(f_cell);

%% ����ΪԼ����������
%����ʽԼ��
A_t = zeros(Nc,Nc);
for p = 1:1:Nc
    for q = 1:1:Nc
        if q<=p
            A_t(p,q) = 1;
        else
            A_t(p,q) = 0;
        end
    end
end
A_I = kron(A_t,eye(Nu)); % ������ڿ˻�
Ut = kron(ones(Nc,1),U);        % UΪ��һʱ�̵�Uֵ
umin = [-0.2;-0.02;];    % ά������Ʊ���������ͬ,��ʾu-ur����Сֵ�����ֵ
umax = [0.2;0.02;];    
delta_umin = [-0.05;-0.0082;];
delta_umax = [0.05;0.0082;];
Umin = kron(ones(Nc,1),umin);
Umax = kron(ones(Nc,1),umax);
A_cons_cell = {A_I , zeros(Nu*Nc,1);-A_I , zeros(Nu*Nc,1)};
b_cons_cell = {Umax-Ut;-Umin+Ut};
A_cons = cell2mat(A_cons_cell);
b_cons = cell2mat(b_cons_cell);
% ״̬��Լ��
M = 10;
delta_Umin = kron(ones(Nc,1),delta_umin);
delta_Umax = kron(ones(Nc,1),delta_umax);
lb = [delta_Umin;0]; % ״̬���½磬��������ʱ���ڿ����������ɳ�����
ub = [delta_Umax;M]; % ״̬���Ͻ磬��������ʱ���ڿ����������ɳ�����

%% ��ʼ������
options = optimset('Algorithm','interior-point-convex');
[X,fval,exitflag] = quadprog(H,f,A_cons,b_cons,[],[],lb,ub,[],options);
%% �������
u_piao(1) = X(1);               % deltaU
u_piao(2) = X(2);
U(1) = kesi(4)+u_piao(1); % ���ڴ洢��һ��ʱ�̵Ŀ���������ΪdeltaU+��u(k-1)-ur(k-1)��=u(k)-ur(k)
U(2) = kesi(5)+u_piao(2);
u_real(1) = U(1)+vd1;       % ��Ϊu(k)
u_real(2) = U(2)+vd2;
sys = u_real;
toc
% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
% function sys=mdlGetTimeOfNextVarHit(t,x,u)
% 
% sampleTime = 1;    %  Example, set the next hit to be one second later.
% sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
% function sys=mdlTerminate(t,x,u)
% 
% sys = [];

% end mdlTerminate
