function [sys,x0,str,ts] = T_302(t,x,u,flag)
  
switch flag,

  %%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes();

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case { 2, 4, 9 },
    sys = [];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end
% end csfunc

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes()
sizes = simsizes;
sizes.NumContStates  = 2;   % Two continuous states: V and C_X_tank
sizes.NumDiscStates  = 0;   % No discrete states
sizes.NumOutputs     = 2;   % Three outputs: V, C_X_tank, and Ea
sizes.NumInputs      = 3;   % Three inputs: Q_T, C_X_T, Q_OUT
sizes.DirFeedthrough = 1;   % Direct feedthrough is enabled
sizes.NumSampleTimes = 1;   % Single sample time

sys = simsizes(sizes);
x0  = [2.253 0.01832];  % Initial conditions for V and C_X_tank
str = [];
ts  = [0 0];    % Continuous sample time

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)

 % States
 V = x(1);             % Volume
 C_linpDNA_tank = x(2);      % Concentration in the tank

 % Inputs
 Q_T = u(1);                 % Input flow rate
 C_linpDNA_T = u(2);         % Input concentration
 Q_OUT = u(3);               % Output flow rate
 
%%


 % Parameter values
  dVdt=Q_T-Q_OUT;
  dC_linpDNA_dt=(C_linpDNA_T*Q_T-C_linpDNA_tank*Q_OUT-C_linpDNA_tank*(Q_T-Q_OUT))/V;

sys = [dVdt; dC_linpDNA_dt];

% end mdlDerivatives

%
%=============================================================================
% mdlOutputs
% Return the outputs of the system.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

 % States
 V = x(1);                   % Volume
 C_linpDNA_tank = x(2);      % Concentration in the tank



 % Outputs: V, C_X_tank, and Ea
 sys = [V; C_linpDNA_tank];

% end mdlOutputs
