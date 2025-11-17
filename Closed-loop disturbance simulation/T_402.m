function [sys,x0,str,ts] = T_402(t,x,u,flag)
  
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization 
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
sizes.NumContStates  = 3;   % Two continuous states: V and C_X_tank
sizes.NumDiscStates  = 0;   % No discrete states
sizes.NumOutputs     = 3;   % Three outputs: V, C_X_tank, and Ea
sizes.NumInputs      = 4;   % Three inputs: Q_T, C_X_T, Q_OUT
sizes.DirFeedthrough = 1;   % Direct feedthrough is enabled
sizes.NumSampleTimes = 1;   % Single sample time

sys = simsizes(sizes);
x0  = [0.4507 0.44 0.0091];   % Initial conditions for V and C_X_tank
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

 % States
 V = x(1);                  % Volume
 C_mRNA_tank = x(2);        % Concentration in the tank
 C_impurities_tank = x(3);  % Concentration in the tank

 % Inputs
 Q_T = u(1);                % Input flow rate
 C_mRNA_T = u(2);           % Input concentration
 Q_OUT = u(3);              % Output flow rate
 C_impurities_T = u(4);     % Input concentration
 

 % Parameter values
 dVdt = Q_T - Q_OUT;   % Change in volume
 dC_mRNA_dt = (C_mRNA_T * Q_T - C_mRNA_tank * Q_OUT - C_mRNA_tank * (Q_T - Q_OUT)) / V; % Change in concentration
 dC_impurities_dt = (C_impurities_T * Q_T - C_impurities_tank * Q_OUT - C_impurities_tank * (Q_T - Q_OUT)) / V; % Change in concentration

sys = [dVdt; dC_mRNA_dt; dC_impurities_dt];

% end mdlDerivatives

%
%=============================================================================
% mdlOutputs
% Return the outputs of the system.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

 % States
 V = x(1);                      % Volume
 C_mRNA_tank= x(2);             % Concentration in the tank
 C_impurities_tank= x(3);       % Concentration in the tank


 % Outputs: V, C_X_tank, and Ea
 sys = [V; C_mRNA_tank;  C_impurities_tank];

% end mdlOutputs
