function [sys,x0,str,ts] = M_501(t,x,u,flag)
  
switch flag,

  %%%%%%%%%%%%%%%%%%
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
sizes.NumContStates  = 6;   % There are 6 continuous states: Volume and 5 concentrations
sizes.NumDiscStates  = 0;   % No discrete states
sizes.NumOutputs     = 6;   % 6 outputs: V, C_LNP_tank, C_ethanol_tank, C_sucrose_tank, C_imp1_tank, C_imp2_tank
sizes.NumInputs      = 7;   % 7 inputs: Q_T, C_LNP_T, C_ethanol_T, C_sucrose_T, C_imp1_T, C_imp2_T, Q_OUT
sizes.DirFeedthrough = 1;   % Direct feedthrough is enabled
sizes.NumSampleTimes = 1;   % Single sample time

sys = simsizes(sizes);

% Initial conditions (x0) - these values define the initial state of the system:
x0  = [4.827 0.1650 0.001595 0.03772 0.0001 0.0001]; % Initial conditions for V, C_LNP_tank, C_ethanol_tank, C_sucrose_tank, C_imp1_tank, C_imp2_tank
str = [];
ts  = [0 0];                % Continuous sample time

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)

 % States - these are the current values of the states at time t
 V = x(1);                  % Volume in the tank (current state)
 C_LNP_tank = x(2);         % Concentration of LNP in the tank (current state)
 C_ethanol_tank = x(3);     % Concentration of Ethanol in the tank (current state)
 C_sucrose_tank = x(4);     % Concentration of Sucrose in the tank (current state)
 C_imp1_tank = x(5);        % Concentration of Impurity 1 in the tank (current state)
 C_imp2_tank = x(6);        % Concentration of Impurity 2 in the tank (current state)
      
 % Inputs - these are the values of the inputs at time t
 Q_T = u(1);                % Input flow rate
 C_LNP_T = u(2);            % Input concentration of LNP
 C_ethanol_T = u(3);        % Input concentration of Ethanol
 C_sucrose_T = u(4);        % Input concentration of Sucrose
 C_imp1_T = u(5);           % Input concentration of Impurity 1
 C_imp2_T = u(6);           % Input concentration of Impurity 2
 Q_OUT = u(7);              % Output flow rate
 
 % Parameter values - these represent the change rates of the system states
 dVdt = Q_T - Q_OUT;   % Change in volume
 dC_LNP_dt = (C_LNP_T * Q_T - C_LNP_tank * Q_OUT - C_LNP_tank * (Q_T - Q_OUT)) / V; % Change in concentration of LNP
 dC_ethanol_dt = (C_ethanol_T * Q_T - C_ethanol_tank * Q_OUT - C_ethanol_tank * (Q_T - Q_OUT)) / V; % Change in concentration of Ethanol
 dC_sucrose_dt = (C_sucrose_T * Q_T - C_sucrose_tank * Q_OUT - C_sucrose_tank * (Q_T - Q_OUT)) / V; % Change in concentration of Sucrose
 dC_imp1_dt = (C_imp1_T * Q_T - C_imp1_tank * Q_OUT - C_imp1_tank * (Q_T - Q_OUT)) / V; % Change in concentration of Impurity 1
 dC_imp2_dt = (C_imp2_T * Q_T - C_imp2_tank * Q_OUT - C_imp2_tank * (Q_T - Q_OUT)) / V; % Change in concentration of Impurity 2

sys = [dVdt; dC_LNP_dt; dC_ethanol_dt; dC_sucrose_dt; dC_imp1_dt; dC_imp2_dt];

% end mdlDerivatives

%
%=============================================================================
% mdlOutputs
% Return the outputs of the system.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

% Outputs - these are the values of the system states that will be returned as outputs
sys = [x(1);  % Volume
       x(2);  % Concentration of LNP in the tank
       x(3);  % Concentration of Ethanol in the tank
       x(4);  % Concentration of Sucrose in the tank
       x(5);  % Concentration of Impurity 1 in the tank
       x(6)]; % Concentration of Impurity 2 in the tank

% end mdlOutputs
