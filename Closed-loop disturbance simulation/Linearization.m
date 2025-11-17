function [sys, X_0, str, ts] = Linearization(t, x, u, flag, X_0)
% CSTR Function simulating the behavior of a Continuous Stirred Tank Reactor (CSTR)

% Persistent variable to store C_Plasmid_F
persistent C_Plasmid_F_stored;

switch flag
    case 0 % Initialization
        % Define system sizes
        sizes = simsizes;
        
        % Number of continuous states
        sizes.NumContStates = 3; 
        
        % No discrete states
        sizes.NumDiscStates = 0;
        
        % Number of outputs, including Ea
        sizes.NumOutputs = 6; % Volume, Plasmid, Linearized Plasmid, Impurities, Ea, Height
        
        % Number of inputs
        sizes.NumInputs = 3;
        
        % No direct feedthrough from input to output
        sizes.DirFeedthrough = 0;
        
        % One sample time
        sizes.NumSampleTimes = 1;
        
        % Initialize sizes
        sys = simsizes(sizes);
        
        % Initial conditions
        X_0 = [13.2; 0.002027; 0.00857]; % Volume, Plasmid, Linearized Plasmid, Impurities
        
        % Not used
        str = [];
        
        % Continuous-time system
        ts = [0 0];
        
        % Initialize the persistent variable
        C_Plasmid_F_stored = 0;
    
    case 1 % Derivatives
        % Extract states
        V = x(1);              % Volume (L)
        C_Plasmid_S = x(2);    % Concentration of plasmid (g/L)
        C_Lin_S = x(3);        % Concentration of Linearized plasmid (g/L)
       
        % Extract inputs
        C_Plasmid_F = u(1);    % Feed concentration of plasmid (g/L)
        Q_IN = u(2);           % Feed flow rate (L/h)
        Q_OUT = u(3);          % Outflow rate (L/h)
        
        % Store C_Plasmid_F in the persistent variable
        C_Plasmid_F_stored = C_Plasmid_F;
        
        % System parameters
        k_cat = 792;           % Linearization rate constant (h^-1)
        k_m=(2.2*10^-8)*660;   % mol/L * g/mol(molar weight DNA) = g/L
     
        
        % Linearization by Michealis-Menten
        Ea = C_Plasmid_S * 0.8;                             % g/L
        v_max = k_cat * Ea;                                 % g/L*h
        v = (v_max * C_Plasmid_S) / (k_m + C_Plasmid_S);    % Linearization rate g/L*h
        
        % Reactor geometry
        d = 0.212;                                            % Diameter of the reactor (m)
        A = pi / 4 * d^2;                                   % Cross-sectional area of the reactor (m^2)
        
        % Fluid height
        h = (V / 1000) / A;                                 % Height (m)
        
        % State derivatives
        dVdt = Q_IN - Q_OUT;                                % Volume change (L/h)
        dC_Plasmid_Sdt = ((Q_IN * C_Plasmid_F - Q_OUT * C_Plasmid_S) / V) - v *C_Plasmid_S ;
        dC_Lin_Sdt = v *C_Plasmid_S  - (Q_OUT * C_Lin_S / V);   % Lin-pDNA change rate

        % Return the state derivatives
        sys = [dVdt; dC_Plasmid_Sdt; dC_Lin_Sdt];
    
    case 3 % Outputs
        % Extract current states
        V = x(1); % Volume (L)
        
        % Reactor geometry
        d = 0.212;               % Diameter of the reactor (m)
        A = pi / 4 * d^2;      % Cross-sectional area of the reactor (m^2)
        
        % Calculate the fluid height
        h = (V / 1000) / A;    % Height (m)

        % Calculate Ea using the stored persistent value of C_Plasmid_F
        C_Plasmid_S = x(2);
        Ea = C_Plasmid_S * 0.8; % Enzyme activity

        C_Lin_S = x(3);
        impurities = C_Lin_S * 0.05; % Enzyme activity
        
        % Return the outputs: Volume, Plasmid, Linearized Plasmid, Impurities, Ea, Height
        sys = [x(1); x(2); x(3); impurities;  Ea; h];
    
    case {2, 4, 9}
        % No action for unused flags
        sys = [];
    
    otherwise
        % Handle unhandled flags
        error(['Unhandled flag = ', num2str(flag)]);
end
