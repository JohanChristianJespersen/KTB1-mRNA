function [sys,X_0,str,ts] = CSTR(t,x,u,flag,X_0)
% CSTR Function simulating the behavior of a Continuous Stirred Tank Reactor (CSTR)
% Inputs:
% t - time
% x - state variables
% u - input variables
% flag - operation flag determining what the function should do
% X_0 - initial conditions of the system

switch flag
    case 0 % initialization;
        % Define system sizes (states, inputs, outputs, etc.)
        sizes = simsizes;
        
        % Set the number of continuous states (5 state variables)
        sizes.NumContStates = 4;
        
        % No discrete states in this system
        sizes.NumDiscStates = 0;
        
        % Define the number of outputs (6 outputs)
        sizes.NumOutputs = 5;
        
        % Number of inputs (6 input variables)
        sizes.NumInputs = 5;
        
        % No direct feedthrough from input to output
        sizes.DirFeedthrough = 0;
        
        % One sample time for the system
        sizes.NumSampleTimes = 1;
        
        % Initialize system sizes
        sys = simsizes(sizes);
        
        % Initial conditions of the states (volume, concentrations of various compounds)
        % x(1): Volume of the reactor
        % x(2): Biomass concentration
        % x(3): Glucose concentration
        % x(4): Nitrogen concentration
        X_0 = [100 8.51 1.118 3.94];
     
     
        % Not used here
        str = [];
        
        % Set sample time as continuous
        ts = [0 0];
    
    case 1 % derivatives
        % This block computes the state derivatives (rate of change of state variables)

        % Extract current states
        V = x(1);          % Volume of the reactor
        C_X_S = x(2);      % Concentration of biomass
        C_GLU_S = x(3);    % Concentration of glucose
        C_N_S = x(4);      % Concentration of nitrogen
        
        % Extract input variables
        C_X_F = u(1);      % Concentration of biomass in feed
        C_GLU_F = u(2);    % Concentration of glucose in feed
        C_N_F = u(3);      % Concentration of nitrogen in feed
        Q_IN = u(4);       % Feed flow rate
        Q_OUT = u(5);      % Out flow rate

        % Define system parameters
        mu_max = 0.57;        % Maximum growth rate of biomass (1/h)
        Y_X_GLU = 0.6;        % Yield of biomass from glucose (g_X/g_GLU)
        Y_X_N = 7.85;         % Yield of biomass from nitrogen (g_X/g_N)
        K_GLU = 1.8;          % Constant for glucose (g_GLY/g_X)
        K_N = 0.0884;         % Constant for nitrogen (g_N/g_X)
        
        % Reactor geometry (changeable if necessary)
        d = 0.458;             % Diameter of the reactor (m)
        A = pi/4 * d^2;        % Cross-sectional area of the reactor (m^2)      
    
        % Growth and production rates
        mu_1 = (mu_max * C_GLU_S / (C_GLU_S + K_GLU * C_X_S) * C_N_S / (C_N_S + K_N * C_X_S));               % Biomass growth
        mu_2 = (-1 / Y_X_GLU*mu_max * C_GLU_S / (C_GLU_S + K_GLU * C_X_S) * C_N_S / (C_N_S + K_N * C_X_S));  % Glucose consumption
        mu_3 = (-1 / Y_X_N*mu_max *C_GLU_S / (C_GLU_S + K_GLU * C_X_S) * C_N_S / (C_N_S + K_N * C_X_S));     % Nitrogen consumption

        % Calculate the height of the fluid in the reactor
        h = (V / 1000) / A;     % Height (m)

        % State derivatives (rate of change of states)
        dVdt = (Q_IN - Q_OUT); % Change in volume
        dC_X_Sdt = (mu_1 * C_X_S + (Q_IN * C_X_F - Q_OUT * C_X_S-C_X_S*(Q_IN-Q_OUT))/ V);                    % Change in biomass concentration
        dC_GLU_Sdt = (mu_2 * C_X_S + (Q_IN * C_GLU_F - Q_OUT * C_GLU_S-C_GLU_S*(Q_IN-Q_OUT)) / V);           % Change in glucose concentration
        dC_N_Sdt = ( mu_3 * C_X_S + (Q_IN * C_N_F - Q_OUT * C_N_S-C_N_S*(Q_IN-Q_OUT)) / V);                  % Change in nitrogen concentration
       

        % Return the state derivatives as the output
        sys = [dVdt; dC_X_Sdt; dC_GLU_Sdt; dC_N_Sdt];
    
    case 3 % outputs
        % This block computes and returns the current output variables (e.g., reactor states)

        % Reactor geometry and calculation of height (as in case 1)
        d = 0.458;             % Reactor diameter (m)
        A = pi / 4 * d^2;      % Reactor cross-sectional area (m^2)
        h = (x(1) / 1000) / A; % Height of the fluid in the reactor (m)

        % Return the current state values and fluid height as the output
         % sys = [Volume; Biomass Concentration; Glucose Concentration; Nitrogen Concentration; Height]
        sys = [x(1); x(2); x(3); x(4); h];

    case {2, 4, 9}
        % For flags 2, 4, and 9, no action is taken (not used here)
        sys = [];

    otherwise
        % If an unhandled flag is passed, throw an error
        error(['unhandled flag = ', num2str(flag)]);
end
