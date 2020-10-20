function [sys,x0,str,ts] = MPC_controller(t,x,u,flag)
% (MATLAB version:R2019a)
% vehicle dynamics and control
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
%% INITIALIZATION    
    [sys,x0,str,ts]=mdlInitializeSizes;
	basic_state_size_ = 4;% number of state:lateral error,lateral error rate,heading error, heading error rate
    controls_ = 1;% number of controls:delta_f
    vertical_ = 5;% Nc
    M_SU = 4455;M_US1 = 570;M_US2 = 735;
    mass_ = M_SU+M_US1+M_US2;
    lf_ = 1110/1000;
    lr_ = 2790/1000;
    cf_ = 2 * (20164.4-15677.2)/(2*pi/180);
    cr_ = cf_;
    iz_ = 34802.6;
    global matrix_a_ matrix_a_coeff_ matrix_b_ matrix_d_ pre_matrix_d_ matrix_state_ pre_matrix_state_ matrix_q_ matrix_r_ pre_control_;
    global trajectory_analyzer_;
    matrix_a_ = zeros(basic_state_size_,basic_state_size_);
    matrix_a_coeff_ = zeros(basic_state_size_,basic_state_size_);
    matrix_b_ = zeros(basic_state_size_,controls_);
    matrix_d_ = zeros(basic_state_size_,1);
    pre_matrix_d_ = zeros(basic_state_size_,1);
    matrix_state_ = zeros(basic_state_size_,1);
    pre_matrix_state_ = zeros(basic_state_size_,1);
    matrix_a_(1, 2) = 1.0;
    matrix_a_(2, 3) = (cf_ + cr_) / mass_;
    matrix_a_(3, 4) = 1.0;
    matrix_a_(4, 3) = (lf_ * cf_ - lr_ * cr_) / iz_;
    matrix_a_coeff_(2, 2) = -(cf_ + cr_) / mass_;
    matrix_a_coeff_(2, 4) = (lr_ * cr_ - lf_ * cf_) / mass_;
    matrix_a_coeff_(3, 4) = 1.0;
    matrix_a_coeff_(4, 2) = (lr_ * cr_ - lf_ * cf_) / iz_;
    matrix_a_coeff_(4, 4) = -1.0 * (lf_ * lf_ * cf_ + lr_ * lr_ * cr_) / iz_;
    matrix_b_(2, 1) = cf_ / mass_;
    matrix_b_(4, 1) = lf_ * cf_ / iz_;
    matrix_q_ = eye(basic_state_size_);
    matrix_r_ = eye(controls_);
    pre_control_ = zeros(vertical_ * controls_ ,1);
%% generate trajectory
% 5km straight road
    N_S = 1500;     % 道路长度
    ds = 0.01;       % 路线间距
    X_ref = (0:ds:N_S)';
    N = length(X_ref);
    Y_ref = X_ref;
    yaw_ref = pi/4*ones(N,1);
    kappa_ref = zeros(N,1);
    trajectory_analyzer_ = [X_ref,Y_ref,yaw_ref,kappa_ref];
% 1km sine road
% load('D:\TRUNK\E_matlab-ros-simulation\map\test_S\trajectory.mat');
% trajectory_analyzer_ = trajectory;
% clear trajectory;
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
sizes.NumDiscStates  = 5;    % lateral error, lateral error rate , heading error, heading error rate , steer angle
sizes.NumOutputs     = 5;    % lateral error, lateral error rate , heading error, heading error rate , steer angle
sizes.NumInputs      = 10;   % x,y,theta,vx,w = dtheta / dt
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [0;0;0;0;0];

%
% str is always an empty matrix
%
str = [];
%
% initialize the array of sample times
%
ts  = [0.01 0];

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
%% 设置参数
global matrix_a_ matrix_a_coeff_ matrix_b_ matrix_d_ pre_matrix_d_ matrix_state_ pre_matrix_state_ matrix_q_ matrix_r_ pre_control_;
global trajectory_analyzer_;
basic_state_size_ = size(matrix_a_,1);% number of state:lateral error,lateral error rate,heading error, heading error rate
controls_ = size(matrix_b_,2);% number of controls:delta_f
horizon_ = 20;% Np
vertical_ = size(pre_control_,1)/controls_;% Nc
% Row = 10;% 松弛因子
% weight param
weight_lateral_error = 1;
weight_lateral_error_rate = 0;
weight_heading_error = 1;
weight_heading_error_rate = 0;
weight_steer = 2;
ts_ = 0.01;
wheel_single_direction_max_degree_ = 30;            % max wheel degree(deg)
curr_u = u(1:5,1);      % current input: x,y,phi,linear_v,anglar_v
pre_u = u(6:10,1);      % pre input: x,y,phi,linear_v,anglar_v;
curr_u(3) = curr_u(3)*pi/180; % deg -> rad
curr_u(5) = curr_u(5)*pi/180; % deg -> rad
curr_u(4) = curr_u(4)/3.6;    % km/h -> m/s
pre_u(3) = pre_u(3)*pi/180; % deg -> rad
pre_u(5) = pre_u(5)*pi/180; % deg -> rad
pre_u(4) = pre_u(4)/3.6;    % km/h -> m/s

global matrix_a_ matrix_a_coeff_ matrix_b_ matrix_d_ pre_matrix_d_ matrix_state_ pre_matrix_state_ matrix_q_ matrix_r_ pre_control_;
global trajectory_analyzer_;

% The path points are more dense than the vehicle location points
%% start MPC calculation

fprintf('Update start, t=%6.3f\n',t)
% Trucksim输出为前轴中心的位姿，需要转化为后轴中心
% debug: lateral_error;lateral_error_rate;heading_error;heading_error_rate
matrix_state_ = ComputeLateralErrors(curr_u,trajectory_analyzer_);
% matrix_state_(2) = 0;matrix_state_(4) = 0;
pre_matrix_state_ = x(1:4);
% pre_matrix_state_(2) = 0;pre_matrix_state_(4) = 0;
% pre_matrix_state_ = ComputeLateralErrors(pre_u,trajectory_analyzer_);
%% update matrix
% update matrix_ad_,matrix_bd_,matrix_dd_,matrix_q_,matrix_r_
linear_v = curr_u(4);
pre_linear_v = pre_u(4);
matrix_a_(2, 2) = matrix_a_coeff_(2, 2) / linear_v;
matrix_a_(2, 4) = matrix_a_coeff_(2, 4) / linear_v;
matrix_a_(4, 2) = matrix_a_coeff_(4, 2) / linear_v;
matrix_a_(4, 4) = matrix_a_coeff_(4, 4) / linear_v;
matrix_ad_ = inv(eye(basic_state_size_) - ts_ * 0.5 * matrix_a_) * (eye(basic_state_size_) + ts_ * 0.5 * matrix_a_);
matrix_bd_ = matrix_b_ * ts_;
matrix_d_(2,1) = matrix_a_(2,4) - linear_v;
matrix_d_(4,1) = matrix_a_(4,4);
pre_matrix_d_(2,1) = matrix_a_coeff_(2,4) / pre_linear_v - pre_linear_v;
pre_matrix_d_(4,1) = matrix_a_coeff_(4,4) / pre_linear_v;
matrix_dd_ = matrix_d_ * matrix_state_(4) * ts_;
pre_matrix_dd_ = pre_matrix_d_ * pre_matrix_state_(4) * ts_;
matrix_q_ = diag([weight_lateral_error,weight_lateral_error_rate,weight_heading_error,weight_heading_error_rate]);
matrix_r_ = weight_steer;

% update  matrix_lower_,matrix_upper_,matrix_state_,pre_matrix_state_
matrix_lower_ = -wheel_single_direction_max_degree_ * pi/180;
matrix_upper_ = wheel_single_direction_max_degree_ * pi/180;

% update reference,pre_reference
reference = zeros(basic_state_size_*horizon_,1);
pre_reference = reference;

%% calculate control cmd
[control] = SolveLinearMPC(matrix_ad_,matrix_bd_,matrix_dd_,pre_matrix_dd_,matrix_q_,matrix_r_,...
    matrix_lower_,matrix_upper_,matrix_state_,pre_matrix_state_,reference,pre_reference,pre_control_);
pre_control_ = control;
sys = [matrix_state_;pre_control_(1)];

% end mdlUpdate

%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
sys = x;    % matrix_state;control

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
