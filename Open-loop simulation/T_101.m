function [sys,x0,str,ts] = T_101(t,x,u,flag)
% T_101 - S-function for a dynamic system in Simulink
% Inputs:
%   t    - Current simulation time
%   x    - State vector (contains system states)
%   u    - Input vector (external inputs to the system)
%   flag - Action flag indicating the operation to perform
% Outputs:
%   sys  - Depending on flag, it can represent derivatives, outputs, etc.
%   x0   - Initial state values
%   str  - Reserved for future use (always empty here)
%   ts   - Sample times

switch flag
    %%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%
    case 0
        % Initialize sizes, states, and sample times
        [sys,x0,str,ts] = mdlInitializeSizes();

    %%%%%%%%%%%%%%%
    % Derivatives %
    %%%%%%%%%%%%%%%
    case 1
        % Compute the derivatives of the continuous states
        sys = mdlDerivatives(t,x,u);

    %%%%%%%%%%%
    % Outputs %
    %%%%%%%%%%%
    case 3
        % Compute the system's outputs
        sys = mdlOutputs(t,x,u);

    %%%%%%%%%%%%%%%%%%%
    % Unhandled flags %
    %%%%%%%%%%%%%%%%%%%
    case {2, 4, 9}
        % Return empty for unhandled flags
        sys = [];

    %%%%%%%%%%%%%%%%%%%%
    % Unexpected flags %
    %%%%%%%%%%%%%%%%%%%%
    otherwise
        % Error handling for unexpected flags
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
% end of function T_101

%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
function [sys,x0,str,ts] = mdlInitializeSizes()
% Define the sizes of states, inputs, outputs, and sample times

sizes = simsizes;                   % Initialize the sizes structure
sizes.NumContStates  = 2;           % Number of continuous states
sizes.NumDiscStates  = 0;           % Number of discrete states
sizes.NumOutputs     = 2;           % Number of outputs
sizes.NumInputs      = 3;           % Number of inputs
sizes.DirFeedthrough = 1;           % Direct feedthrough flag
sizes.NumSampleTimes = 1;           % Number of sample times

sys = simsizes(sizes);              % Apply sizes to the system
x0  = [2.567 11.02];                % Initial state values
% State 1: Initial Volume (2.567)
% State 2: Initial Concentration in the tank (11.02)

str = [];                           % Reserved for future use
ts  = [0 0];                        % Sample time definition

% end of function mdlInitializeSizes

%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
function sys = mdlDerivatives(t,x,u)
% Compute the time derivatives of the states

% States
V = x(1);                           % State 1: Current Volume of the tank
C_X_tank = x(2);                    % State 2: Current Concentration in the tank

% Inputs
Q_T = u(1);                         % Input 1: Flow rate into the tank
C_X_T = u(2);                       % Input 2: Concentration of inflow
Q_OUT = u(3);                       % Input 3: Flow rate out of the tank

% Rate of change calculations
dVdt = Q_T - Q_OUT;                 % Rate of change of volume
dC_X_dt = (C_X_T*Q_T - C_X_tank*Q_OUT - ...
           C_X_tank*(Q_T - Q_OUT)) / V; % Rate of change of concentration

sys = [dVdt; dC_X_dt];              % Return the derivatives as a vector
% State 1 Derivative: Change in Volume (dVdt)
% State 2 Derivative: Change in Concentration (dC_X_dt)

% end of function mdlDerivatives

%=============================================================================
% mdlOutputs
% Return the outputs for the system.
%=============================================================================
function sys = mdlOutputs(t,x,u)
% Define the system outputs based on the current states

sys = [x(1); x(2)];                 % Output the current states
% Output 1: Current Volume (x(1))
% Output 2: Current Concentration in the tank (x(2))

% end of function mdlOutputs
