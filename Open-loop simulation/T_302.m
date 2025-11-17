function [sys,x0,str,ts] = T_302(t,x,u,flag)
% T_302 - S-function for a tank model simulating volume and concentration dynamics.
% Inputs:
% t - Time
% x - State variables
% u - Input variables
% flag - Operation flag determining function behavior

switch flag
    %%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%
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
% End T_302 function

%=============================================================================
% mdlInitializeSizes
% Define system sizes, initial conditions, and sample times
%=============================================================================
function [sys,x0,str,ts] = mdlInitializeSizes()
    sizes = simsizes;
    sizes.NumContStates = 2;   % Two continuous states: V and C_linpDNA_tank
    sizes.NumDiscStates = 0;   % No discrete states
    sizes.NumOutputs = 2;      % Two outputs: V and C_linpDNA_tank
    sizes.NumInputs = 3;       % Three inputs: Q_T, C_linpDNA_T, Q_OUT
    sizes.DirFeedthrough = 1;  % Enable direct feedthrough
    sizes.NumSampleTimes = 1;  % Single sample time

    sys = simsizes(sizes);
    x0 = [2.253 0.01832];  % Initial conditions for V and C_linpDNA_tank
    str = [];              % Reserved for future use
    ts = [0 0];            % Continuous sample time

% End mdlInitializeSizes

%=============================================================================
% mdlDerivatives
% Compute the derivatives for the continuous states
%=============================================================================
function sys = mdlDerivatives(t,x,u)
    % States
    V = x(1);                   % State 1: Volume of the tank (L)
    C_linpDNA_tank = x(2);      % State 2: Concentration in the tank (g/L)

    % Inputs
    Q_T = u(1);                 % Input 1: Input flow rate (L/h)
    C_linpDNA_T = u(2);         % Input 2: Input concentration (g/L)
    Q_OUT = u(3);               % Input 3: Output flow rate (L/h)

    % Derivatives
    dVdt = Q_T - Q_OUT;  % Rate of change of tank volume (L/h)
    dC_linpDNA_dt = (C_linpDNA_T * Q_T - C_linpDNA_tank * Q_OUT - ...
                     C_linpDNA_tank * (Q_T - Q_OUT)) / V; % Rate of change of concentration (g/L/h)

    % Return the state derivatives
    sys = [dVdt; dC_linpDNA_dt];

% End mdlDerivatives

%=============================================================================
% mdlOutputs
% Compute the outputs of the system
%=============================================================================
function sys = mdlOutputs(t,x,u)
    % States
    V = x(1);                   % State 1: Volume of the tank (L)
    C_linpDNA_tank = x(2);      % State 2: Concentration in the tank (g/L)

    % Outputs:
    % V: Volume of the tank
    % C_linpDNA_tank: Concentration of linear plasmid DNA in the tank
    sys = [V; C_linpDNA_tank];  % Return outputs: Volume and concentration

% End mdlOutputs
