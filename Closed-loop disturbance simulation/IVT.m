function [sys, X_0, str, ts] = IVT(t, x, u, flag, X_0)
% CSTR Function simulating the behavior of a Continuous Stirred Tank Reactor (CSTR)

% Persistent variable to store C_Plasmid_F
persistent C_Plasmid_F_stored;
% Persistent variable to store Cleancap
persistent C_Cleancap_F_stored;
% Persistent variable to store MgAcetate
persistent C_MgAcetate_F_stored;
% Persistent variable to store Nucleotides
persistent C_Nulceotides_F_stored;


switch flag
    case 0 % Initialization
        % Define system sizes
        sizes = simsizes;
        
        % Number of continuous states
        sizes.NumContStates = 2; 
        
        % No discrete states
        sizes.NumDiscStates = 0;
        
        % Number of outputs, including Ea
        sizes.NumOutputs = 9; % Volume, Plasmid, Linearized Plasmid, Impurities, Ea, Height
        
        % Number of inputs
        sizes.NumInputs = 3;
        
        % No direct feedthrough from input to output
        sizes.DirFeedthrough = 1;
        
        % One sample time
        sizes.NumSampleTimes = 1;
        
        % Initialize sizes
        sys = simsizes(sizes);
     

        % Initial conditions
        X_0 = [20; 0.09135]; % Volume, Plasmid, Linearized Plasmid, Impurities
        
        % Not used
        str = [];
        
        % Continuous-time system
        ts = [0 0];
        
        % Initialize the persistent variable
        C_Plasmid_F_stored = 0.0202; %g/L
          % Initialize the persistent variable
        C_Cleancap_F_stored = 0.0248; %M
          % Initialize the persistent variable
        C_MgAcetate_F_stored = 0.050; %M
          % Initialize the persistent variable
        C_Nulceotides_F_stored = 0.031; %M
        

    case 1 % Derivatives
        % Extract states
        V = x(1);               % Volume (L)
        C_mRNA_S=x(2);          % mRNA in the system
     
     
       
        % Extract inputs
        Q_IN = u(1);           % Feed flow rate (L/h)
        Q_OUT = u(2);          % Outflow rate (L/h)
        C_Plasmid_F = u(3);    % Feed concentration of plasmid (g/L)
      
        % Store constants in the persistent variable
        C_Plasmid_F_stored = C_Plasmid_F;

        
        % IVT by Michealis-Menten        
        v_max=5.63;

      
        k_m=(3.16*10^-6); %g/L
        v=((v_max * (C_Plasmid_F*0.1)) / (k_m + (C_Plasmid_F*0.1)));

        % Reactor geometry
        d = 0.266;               % Diameter of the reactor (m)
        A = pi / 4 * d^2;      % Cross-sectional area of the reactor (m^2)
        
        % Fluid height
        h = (V / 1000) / A;    % Height (m)
        
        % State derivatives
        dVdt = Q_IN - Q_OUT;   % Volume change (L/h)
        dC_mRNA_Sdt = (v*C_Plasmid_F*0.1)  - (Q_OUT * C_mRNA_S / V);

        % Return the state derivatives
       
        sys = [dVdt; dC_mRNA_Sdt];
    


    case 3 % Outputs
        % Extract current states
        V = x(1);              % Volume (L)
        C_mRNA_S = x(2);        % mRNA concentration (g/L)

        % Reactor geometry
        d = 0.266;               % Diameter of the reactor (m)
        A = pi / 4 * d^2;      % Cross-sectional area of the reactor (m^2)
        
        % Calculate the fluid height
        h = (V / 1000) / A;    % Height (m)
        
        % yield impurities
         Y_impurities = 0.05;
         C_impurities=Y_impurities*C_mRNA_S;

        % Calculate Ea using the stored persistent value of C_Plasmid_F
        Q_IN=u(1);
        Ea_T7 = Q_IN * 0.063;  % Enzyme conc g/L
        plasmid=C_Plasmid_F_stored;
        Cleancap=C_Cleancap_F_stored;
        MgAcetate=C_MgAcetate_F_stored;
        Nulceotides=C_Nulceotides_F_stored;
 
        % Return the outputs: Volume, Plasmid, Linearized Plasmid, Impurities, Ea, Height
        sys = [x(1); x(2); C_impurities; Ea_T7; plasmid; Cleancap; MgAcetate; Nulceotides; h];
    
    case {2, 4, 9}
        % No action for unused flags
        sys = [];
    
    otherwise
        % Handle unhandled flags
        error(['Unhandled flag = ', num2str(flag)]);
end