function [sys, X_0, str, ts] = IVT(t, x, u, flag, X_0)
% CSTR Function simulating the behavior of a Continuous Stirred Tank
% Reactor (CSTR) for in vitro transcription

% Persistent variables to store input feed concentrations
persistent C_Plasmid_F_stored;
persistent C_Cleancap_F_stored;
persistent C_MgAcetate_F_stored;
persistent C_Nulceotides_F_stored;

switch flag
    case 0 % Initialization
        % Define system sizes
        sizes = simsizes;
        
        % Number of continuous states - Here, we have 2 continuous states
        sizes.NumContStates = 2;  % Volume, mRNA concentration
        
        % No discrete states
        sizes.NumDiscStates = 0;
        
        % Number of outputs, including impurities, enzyme activity, and height
        sizes.NumOutputs = 9; % Volume, mRNA, Impurities, Ea, Plasmid, Cleancap, MgAcetate, Nucleotides, Height
        
        % Number of inputs
        sizes.NumInputs = 3;   % Feed flow rate, Outflow rate, Feed concentration of plasmid
        
        % Direct feedthrough from input to output
        sizes.DirFeedthrough = 1;
        
        % One sample time (continuous-time system)
        sizes.NumSampleTimes = 1;
        
        % Initialize sizes
        sys = simsizes(sizes);
        
        % Initial conditions (X_0) - These are the initial states of the system
        X_0 = [20; 0.09135]; % Initial Volume, mRNA concentration
        
        % Not used, set as empty
        str = [];
        
        % Continuous-time system
        ts = [0 0];
        
        % Initialize the persistent variables with initial values
        C_Plasmid_F_stored = 0.0202; % Feed concentration of plasmid (g/L)
        C_Cleancap_F_stored = 0.0248; % Concentration of Cleancap (M)
        C_MgAcetate_F_stored = 0.050; % Concentration of MgAcetate (M)
        C_Nulceotides_F_stored = 0.031; % Concentration of Nucleotides (M)
    
    case 1 % Derivatives
        % Extract states - These are the current values of the system's continuous states
        V = x(1);               % Volume (L)
        C_mRNA_S = x(2);        % mRNA concentration (g/L)
        
        % Extract inputs - These are the current values of the inputs at time t
        Q_IN = u(1);            % Feed flow rate (L/h)
        Q_OUT = u(2);           % Outflow rate (L/h)
        C_Plasmid_F = u(3);     % Feed concentration of plasmid (g/L)
        
        % Store feed concentration in the persistent variable for later use
        C_Plasmid_F_stored = C_Plasmid_F;

        % IVT reaction using Michealis-Menten kinetics        
        v_max = 5.63;            % Maximum reaction rate (g/L*h)
        k_m = (3.16 * 10^-6);    % Michaelis-Menten constant (g/L)
        v = ((v_max * (C_Plasmid_F * 0.1)) / (k_m + (C_Plasmid_F * 0.1))); % Reaction rate

        % Reactor geometry
        d = 0.266;               % Diameter of the reactor (m)
        A = pi / 4 * d^2;        % Cross-sectional area of the reactor (m^2)
        
        % Fluid height
        h = (V / 1000) / A;      % Height (m)
        
        % State derivatives - these represent the rates of change for each state
        dVdt = Q_IN - Q_OUT;     % Volume change (L/h)
        dC_mRNA_Sdt = (v * C_Plasmid_F * 0.1) - (Q_OUT * C_mRNA_S / V);  % mRNA concentration change rate

        % Return the state derivatives (rate of change for each state)
        sys = [dVdt; dC_mRNA_Sdt];
    
    case 3 % Outputs
        % Extract current states
        V = x(1);               % Volume (L)
        C_mRNA_S = x(2);        % mRNA concentration (g/L)

        % Reactor geometry
        d = 0.266;               % Diameter of the reactor (m)
        A = pi / 4 * d^2;        % Cross-sectional area of the reactor (m^2)
        
        % Calculate the fluid height (h) based on the current volume (V)
        h = (V / 1000) / A;      % Height (m)
        
        % Calculate impurities based on mRNA concentration
        Y_impurities = 0.05;      % Impurity yield factor
        C_impurities = Y_impurities * C_mRNA_S;  % Impurity concentration
        
        % Calculate enzyme concentration (Ea) using the stored persistent value of plasmid feed concentration
        Q_IN = u(1);              % Feed flow rate (L/h)
        Ea_T7 = Q_IN * 0.063;     % Enzyme concentration (g/L)
        
        % Access stored feed concentrations
        plasmid = C_Plasmid_F_stored;
        Cleancap = C_Cleancap_F_stored;
        MgAcetate = C_MgAcetate_F_stored;
        Nulceotides = C_Nulceotides_F_stored;
        
        % Return the outputs: Volume, mRNA concentration, impurities, enzyme concentration (Ea), and feed concentrations
        sys = [x(1); x(2); C_impurities; Ea_T7; plasmid; Cleancap; MgAcetate; Nulceotides; h];
    
    case {2, 4, 9}
        % No action for unused flags
        sys = [];
    
    otherwise
        % Handle unhandled flags
        error(['Unhandled flag = ', num2str(flag)]);
end
