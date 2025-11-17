function [sys, X_0, str, ts] = Linearization(t, x, u, flag, X_0)
% CSTR Function simulating the behavior of a Continuous Stirred Tank
% Reactor (CSTR) for linearization

% Persistent variable to store C_Plasmid_F
persistent C_Plasmid_F_stored;

switch flag
    case 0 % Initialization
        % Define system sizes
        sizes = simsizes;
        
        % Number of continuous states - Here, we have 3 continuous states
        sizes.NumContStates = 3;  % Volume, Plasmid, Linearized Plasmid
        
        % No discrete states
        sizes.NumDiscStates = 0;
        
        % Number of outputs, including Ea (Enzyme activity), impurities, and height
        sizes.NumOutputs = 6; % Volume, Plasmid, Linearized Plasmid, Impurities, Ea, Height
        
        % Number of inputs
        sizes.NumInputs = 3;   % C_Plasmid_F (Plasmid concentration), Q_IN (flow rate), Q_OUT (outflow rate)
        
        % No direct feedthrough from input to output
        sizes.DirFeedthrough = 0;
        
        % One sample time (continuous time)
        sizes.NumSampleTimes = 1;
        
        % Initialize sizes
        sys = simsizes(sizes);
        
        % Initial conditions (X_0) - these are the starting values for the system states
        X_0 = [13.2; 0.002027; 0.00857]; % Initial conditions for Volume, Plasmid, and Linearized Plasmid
        
        % Not used here, set as empty
        str = [];
        
        % Continuous-time system (no discrete-time sample periods)
        ts = [0 0];
        
        % Initialize the persistent variable to store C_Plasmid_F
        C_Plasmid_F_stored = 0;
    
    case 1 % Derivatives
        % Extract states - These represent the current values of the system states at time t
        V = x(1);              % Volume (L) in the reactor (continuous state)
        C_Plasmid_S = x(2);    % Concentration of plasmid (g/L) in the reactor (continuous state)
        C_Lin_S = x(3);        % Concentration of Linearized plasmid (g/L) in the reactor (continuous state)
       
        % Extract inputs - These are the current values of the inputs at time t
        C_Plasmid_F = u(1);    % Feed concentration of plasmid (g/L)
        Q_IN = u(2);           % Feed flow rate (L/h)
        Q_OUT = u(3);          % Outflow rate (L/h)
        
        % Store the current value of C_Plasmid_F in the persistent variable for later use
        C_Plasmid_F_stored = C_Plasmid_F;
        
        % System parameters
        k_cat = 792;           % enzyme turnover number (h^-1)
        k_m = (2.2 * 10^-8) * 660;   % michelias menten constant mol/L * g/mol = g/L 
     
        % Linearization using the Michealis-Menten equation for enzyme kinetics
        Ea = (C_Plasmid_S * 0.8)*0.9;                             % Enzyme cocentration (g/L)
        v_max = k_cat * Ea;                                 % Maxium reaction rate (g/L*h)
        v = (v_max * C_Plasmid_S) / (k_m + C_Plasmid_S);    % Linearization rate (g/L*h)
        
        % Reactor geometry
        d = 0.212;                                            % Diameter of the reactor (m)
        A = pi / 4 * d^2;                                   % Cross-sectional area of the reactor (m^2)
        
        % Fluid height - calculates the height of fluid in the reactor based on volume
        h = (V / 1000) / A;                                 % Height (m)
        
        % State derivatives - these represent the rates of change for each state
        dVdt = Q_IN - Q_OUT;                                % Volume change rate (L/h)
        dC_Plasmid_Sdt = ((Q_IN * C_Plasmid_F - Q_OUT * C_Plasmid_S) / V) - v * C_Plasmid_S;  % Plasmid concentration change rate
        dC_Lin_Sdt = v * C_Plasmid_S - (Q_OUT * C_Lin_S / V); % Linearized plasmid concentration change rate

        % Return the state derivatives (rate of change for each state)
        sys = [dVdt; dC_Plasmid_Sdt; dC_Lin_Sdt];
    
    case 3 % Outputs
        % Extract current states - these are the current values of the states at time t
        V = x(1); % Volume (L)
        
        % Reactor geometry
        d = 0.212;               % Diameter of the reactor (m)
        A = pi / 4 * d^2;      % Cross-sectional area of the reactor (m^2)
        
        % Calculate the fluid height (h) based on the current volume (V)
        h = (V / 1000) / A;    % Height (m)

        % Calculate enzyme activity (Ea) using the stored persistent value of C_Plasmid_F
        C_Plasmid_S = x(2);
        Ea = (C_Plasmid_S * 0.8)*0.9; % Enzyme concentration (g/L)
        
        % Calculate impurities based on the concentration of Linearized plasmid
        C_Lin_S = x(3);
        impurities = C_Lin_S * 0.05; % Impurities in the reactor
        
        % Return the outputs: Volume, Plasmid concentration, Linearized Plasmid concentration,
        % Impurities, Enzyme concentration (Ea), Fluid Height (h)
        sys = [x(1); x(2); x(3); impurities; Ea; h];
    
    case {2, 4, 9}
        % No action for unused flags
        sys = [];
    
    otherwise
        % Handle unhandled flags
        error(['Unhandled flag = ', num2str(flag)]);
end
