function [sys, X_0, str, ts] = Degradation_reactor(t, x, u, flag, X_0)
% Degradation Reactor Function simulating the behavior of a Continuous Stirred Tank Reactor (CSTR)
% In this reactor, degradation of plasmidDNA is modeled.

switch flag
    case 0 % Initialization
        % Define system sizes
        sizes = simsizes;
        
        % Number of continuous states: Added one extra state for degradation product
        sizes.NumContStates = 9; % Volume, mRNA, Impurities, T7 Ea, Linearized DNA, Cleancap, Magnesiumacetate, Nucleotides, Degradation
        
        % No discrete states
        sizes.NumDiscStates = 0;
        
        % Number of outputs, including degradation product
        sizes.NumOutputs = 11; % Volume, mRNA, Impurities, T7 Ea, Linearized DNA, Cleancap, Magnesiumacetate, Nucleotides, Degradation, Enzyme concentration
        
        % Number of inputs: Includes feed concentrations for all substances
        sizes.NumInputs = 9;
        
        % No direct feedthrough from input to output
        sizes.DirFeedthrough = 1;
        
        % One sample time (continuous-time system)
        sizes.NumSampleTimes = 1;
        
        % Initialize sizes
        sys = simsizes(sizes);
        
        % Initial conditions (X_0) for each of the 9 states
        X_0 = [1.13; 0.0069; 0.02184; 0.007060; 0.1; 0.0248; 0.031; 0.050; 0.00001]; 
        % [Volume; mRNA concentration; Impurities; T7 enzyme concentration; Linearized plasmid DNA concentration; 
        %  Cleancap concentration; MgAcetate concentration; Nucleotide concentration; Degradation product concentration]
        
        % Not used, set as empty
        str = [];
        
        % Continuous-time system (no sample time)
        ts = [0 0];
    
    case 1 % Derivatives
        % Extract states: These represent the current values of the system's continuous states
        V = x(1);                % Volume (L)
        C_mRNA_Sy = x(2);        % mRNA concentration (g/L)
        C_impurities_Sy = x(3);  % Impurities concentration (g/L)
        C_T7_Sy = x(4);          % T7 Enzyme concentration (g/L)
        C_linpDNA_Sy = x(5);     % Linearized plasmid DNA concentration (g/L)
        C_cleancap_Sy = x(6);    % Cleancap concentration (g/L)
        C_magnesium_Sy = x(7);   % Magnesium acetate concentration (M)
        C_nucleo_Sy = x(8);      % Nucleotide concentration (M)
        C_degrade_Sy = x(9);     % Degradation product concentration (abortive transcripts) (g/L)
     
        % Extract inputs: These are the current input concentrations and flow rates
        Q_T = u(1);              % Input flow rate (L/h)
        Q_OUT = u(2);            % Output flow rate (L/h)
        C_mRNA_T = u(3);         % mRNA input concentration (g/L)
        C_impurities_T = u(4);   % Impurities input concentration (g/L)
        C_T7_T = u(5);           % T7 Ea input concentration (g/L)
        C_linpDNA_T = u(6);      % Linearized DNA input concentration (g/L)
        C_cleancap_T = u(7);     % Cleancap input concentration (g/L)
        C_magnesium_T = u(8);    % Magnesium acetate input concentration (M)
        C_nucleo_T = u(9);       % Nucleotide input concentration (M)
        
        % System parameters for degradation
        k_m = 0.27;  % Michaelis-Menten constant for degradation (g/L)
        
        % Degradation reaction using Michaelis-Menten kinetics
        Ea_DNase = C_linpDNA_Sy * 1;  % enzyme concentration (based on linearized DNA)
        v_max = 1188000 * Ea_DNase;  % Maximum degradation rate (g/L*h)
        v = (v_max * C_linpDNA_Sy) / (k_m + C_linpDNA_Sy); % Degradation rate
        
        % Reactor geometry
        d = 0.1084;               % Diameter of the reactor (m)
        A = pi / 4 * d^2;         % Cross-sectional area of the reactor (m^2)
        
        % Fluid height based on volume and reactor cross-sectional area
        h = (V / 1000) / A;       % Fluid height (m)
        
        % State derivatives: These represent the rates of change of each state
        dVdt = Q_T - Q_OUT;  % Rate of change of volume (L/h)
        dC_mRNA_dt = (C_mRNA_T * Q_T - C_mRNA_Sy * Q_OUT - C_mRNA_Sy * (Q_T - Q_OUT)) / V; % mRNA concentration rate change
        dC_impurities_dt = (C_impurities_T * Q_T - C_impurities_Sy * Q_OUT - C_impurities_Sy * (Q_T - Q_OUT)) / V; % Impurities concentration rate change
        dC_T7_dt = (C_T7_T * Q_T - C_T7_Sy * Q_OUT - C_T7_Sy * (Q_T - Q_OUT)) / V; % T7 enzyme concentration rate change
        dC_linpDNA_dt = (Q_T * C_linpDNA_T - Q_OUT * C_linpDNA_Sy) / V - v * C_linpDNA_Sy; % Linearized DNA concentration rate change
        dC_cleancap_dt = (C_cleancap_T * Q_T - C_cleancap_Sy * Q_OUT - C_cleancap_Sy * (Q_T - Q_OUT)) / V; % Cleancap concentration rate change
        dC_magnesium_dt = (C_magnesium_T * Q_T - C_magnesium_Sy * Q_OUT - C_magnesium_Sy * (Q_T - Q_OUT)) / V; % Magnesium acetate concentration rate change
        dC_nucleo_dt = (C_nucleo_T * Q_T - C_nucleo_Sy * Q_OUT - C_nucleo_Sy * (Q_T - Q_OUT)) / V; % Nucleotide concentration rate change
        dC_degrade_dt = (v * C_linpDNA_Sy) - (Q_OUT * C_degrade_Sy / V); % Degradation product concentration rate change
        
        % Return the state derivatives
        sys = [dVdt; dC_mRNA_dt; dC_impurities_dt; dC_T7_dt; dC_linpDNA_dt; dC_cleancap_dt; dC_magnesium_dt; dC_nucleo_dt; dC_degrade_dt];
    
    case 3 % Outputs
        % Extract the reactor geometry and calculate fluid height
        d = 0.1084;              % Reactor diameter (m)
        A = pi / 4 * d^2;       % Reactor cross-sectional area (m^2)
        h = (x(1) / 1000) / A;  % Fluid height (m) based on the current volume

        % Extract the enzyme activity from the input flow
        C_linpDNA_T = u(6);     % Linearized plasmid DNA input concentration (g/L)
        Ea_DNase = C_linpDNA_T * 1;  % Enzyme concentration
        
        % Return the outputs
        sys = [x(1); x(2); x(3); x(4); x(5); x(6); x(7); x(8); x(9); Ea_DNase; h];
    
    case {2, 4, 9}
        % No action for unused flags
        sys = [];
    
    otherwise
        % Handle unhandled flags
        error(['Unhandled flag = ', num2str(flag)]);
end
