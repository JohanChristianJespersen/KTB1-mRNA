function [sys,x0,str,ts] = T_402(t,x,u,flag)
% T_402 - S-function for a tank model simulating volume and concentrations of mRNA and impurities.
% Inputs:
% t - Time
% x - State variables
% u - Input variables
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
% End T_402 function

%=============================================================================
% mdlInitializeSizes
% Define system sizes, initial conditions, and sample times
%=============================================================================
function [sys,x0,str,ts] = mdlInitializeSizes()
    sizes = simsizes;
    sizes.NumContStates = 3;    % Three continuous states: V, C_mRNA_tank, C_impurities_tank
    sizes.NumDiscStates = 0;    % No discrete states
    sizes.NumOutputs = 3;       % Three outputs: V, C_mRNA_tank, and C_impurities_tank
    sizes.NumInputs = 4;        % Four inputs: Q_T, C_mRNA_T, Q_OUT, C_impurities_T
    sizes.DirFeedthrough = 1;   % Enable direct feedthrough
    sizes.NumSampleTimes = 1;   % Single sample time

    sys = simsizes(sizes);
    x0 = [0.4507, 0.44, 0.0091];  % Initial conditions for V, C_mRNA_tank, C_impurities_tank
    str = [];                     % Reserved for future use
    ts = [0, 0];                  % Continuous sample time

% End mdlInitializeSizes

%=============================================================================
% mdlDerivatives
% Compute the derivatives for the continuous states
%=============================================================================
function sys = mdlDerivatives(t,x,u)
    % States
    V = x(1);                      % State 1: Volume of the tank (L)
    C_mRNA_tank = x(2);            % State 2: Concentration of mRNA in the tank (g/L)
    C_impurities_tank = x(3);      % State 3: Concentration of impurities in the tank (g/L)

    % Inputs
    Q_T = u(1);                    % Input 1: Input flow rate (L/h)
    C_mRNA_T = u(2);               % Input 2: Input concentration of mRNA (g/L)
    Q_OUT = u(3);                  % Input 3: Output flow rate (L/h)
    C_impurities_T = u(4);         % Input 4: Input concentration of impurities (g/L)

    % Derivatives
    dVdt = Q_T - Q_OUT;  % Rate of change of tank volume (L/h)
    dC_mRNA_dt = (C_mRNA_T * Q_T - C_mRNA_tank * Q_OUT - ...
                  C_mRNA_tank * (Q_T - Q_OUT)) / V; % Rate of change of mRNA concentration (g/L/h)
    dC_impurities_dt = (C_impurities_T * Q_T - C_impurities_tank * Q_OUT - ...
                        C_impurities_tank * (Q_T - Q_OUT)) / V; % Rate of change of impurities concentration (g/L/h)

    % Return the state derivatives
    sys = [dVdt; dC_mRNA_dt; dC_impurities_dt];

% End mdlDerivatives

%=============================================================================
% mdlOutputs
% Compute the outputs of the system
%=============================================================================
function sys = mdlOutputs(t,x,u)
    % States
    V = x(1);                      % State 1: Volume of the tank (L)
    C_mRNA_tank = x(2);            % State 2: Concentration of mRNA in the tank (g/L)
    C_impurities_tank = x(3);      % State 3: Concentration of impurities in the tank (g/L)

    % Outputs:
    % V: Volume of the tank
    % C_mRNA_tank: Concentration of mRNA in the tank
    % C_impurities_tank: Concentration of impurities in the tank
    sys = [V; C_mRNA_tank; C_impurities_tank];  % Return outputs: Volume, mRNA concentration, impurities concentration

% End mdlOutputs
