function [sys,x0,str,ts] = T_501(t,x,u,flag)
% T_501 - S-function for a tank model with concentration of a substance.
% Inputs:
% t - Time
% x - State variables (Volume and Concentration)
% u - Input variables (Flow rates and Concentration)
% flag - Operation flag determining function behavior

switch flag
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0
        % Initialize system sizes, initial states, and sample times
        [sys,x0,str,ts] = mdlInitializeSizes();

    %%%%%%%%%%%%%%%
    % Derivatives %
    %%%%%%%%%%%%%%%
    case 1
        % Compute state derivatives
        sys = mdlDerivatives(t,x,u);

    %%%%%%%%%%%
    % Outputs %
    %%%%%%%%%%%
    case 3
        % Compute system outputs
        sys = mdlOutputs(t,x,u);

    %%%%%%%%%%%%%%%%%%%
    % Unhandled flags %
    %%%%%%%%%%%%%%%%%%%
    case {2, 4, 9}
        % No action for these flags
        sys = [];

    %%%%%%%%%%%%%%%%%%%%
    % Unexpected flags %
    %%%%%%%%%%%%%%%%%%%%
    otherwise
        % Handle unexpected flags with an error
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
% End T_501 function

%=============================================================================
% mdlInitializeSizes
% Define system sizes, initial conditions, and sample times
%=============================================================================
function [sys,x0,str,ts] = mdlInitializeSizes()
    sizes = simsizes;
    sizes.NumContStates = 2;    % Two continuous states: V and C_X_tank
    sizes.NumDiscStates = 0;    % No discrete states
    sizes.NumOutputs = 2;       % Two outputs: V and C_X_tank
    sizes.NumInputs = 3;        % Three inputs: Q_T, C_X_T, Q_OUT
    sizes.DirFeedthrough = 1;   % Enable direct feedthrough
    sizes.NumSampleTimes = 1;   % Single sample time

    sys = simsizes(sizes);
    x0 = [0.4507, 0.4223];  % Initial conditions for V and C_X_tank
    str = [];                % Reserved for future use
    ts = [0, 0];             % Continuous sample time

% End mdlInitializeSizes

%=============================================================================
% mdlDerivatives
% Compute the derivatives for the continuous states
%=============================================================================
function sys = mdlDerivatives(t,x,u)
    % States
    V = x(1);                      % State 1: Volume of the tank (L)
    C_X_tank = x(2);                % State 2: Concentration of substance X in the tank (g/L)

    % Inputs
    Q_T = u(1);                     % Input 1: Input flow rate (L/h)
    C_X_T = u(2);                   % Input 2: Input concentration of substance X (g/L)
    Q_OUT = u(3);                   % Input 3: Output flow rate (L/h)

    % Derivatives
    dVdt = Q_T - Q_OUT;  % Rate of change of volume (L/h)
    dC_X_dt = (C_X_T * Q_T - C_X_tank * Q_OUT - C_X_tank * (Q_T - Q_OUT)) / V;  % Rate of change of concentration of X (g/L/h)

    % Return the state derivatives
    sys = [dVdt; dC_X_dt];

% End mdlDerivatives

%=============================================================================
% mdlOutputs
% Compute the outputs of the system
%=============================================================================
function sys = mdlOutputs(t,x,u)
    % Outputs:
    % x(1) - Volume of the tank
    % x(2) - Concentration of substance X in the tank

    sys = [x(1); x(2)];  % Return outputs: Volume (V) and concentration (C_X_tank)

% End mdlOutputs
