function [sys, X_0, str, ts] = Degradation_reactor(t, x, u, flag, X_0)
% CSTR Function simulating the behavior of a Continuous Stirred Tank Reactor (CSTR)

switch flag
    case 0 % Initialization
        % Define system sizes
        sizes = simsizes;
        
        % Number of continuous states
        sizes.NumContStates = 9; % Added one more state for "degrade"
        
        % No discrete states
        sizes.NumDiscStates = 0;
        
        % Number of outputs, including Ea
        sizes.NumOutputs = 11; % Outputs now include "degrade"
        
        % Number of inputs
        sizes.NumInputs = 9;
        
        % No direct feedthrough from input to output
        sizes.DirFeedthrough = 1;
        
        % One sample time
        sizes.NumSampleTimes = 1;
        
        % Initialize sizes
        sys = simsizes(sizes);
        
        % Initial conditions
        X_0 = [1.13; 0.0069; 0.02184; 0.007060; 0.1; 0.0248; 0.031; 0.050; 0.00001]; % Added initial condition for "degrade"
        
        % Not used
        str = [];
        
        % Continuous-time system
        ts = [0 0];
        
    case 1 % Derivatives
        % Extract states
        V = x(1);                % Volume
        C_mRNA_Sy = x(2);        % mRNA concentration
        C_impurities_Sy = x(3);  % Impurities concentration
        C_T7_Sy = x(4);          % Ea T7 concentration
        C_linpDNA_Sy = x(5);     % LinpDNA concentration
        C_cleancap_Sy = x(6);    % cleancap concentration
        C_magnesium_Sy = x(7);   % Magnesium ion concentration
        C_nucleo_Sy = x(8);      % Nucleotide concentration
        C_degrade_Sy = x(9);     % Degradation product concentration (abortive transcripts)
     
        % Extract inputs
        Q_T = u(1);              % Input flow rate
        Q_OUT = u(2);            % Output flow rate
        C_mRNA_T = u(3);         % mRNA input concentration
        C_impurities_T = u(4);   % Impuriteis input concentration
        C_T7_T = u(5);           % T7 Ea input concentration
        C_linpDNA_T = u(6);      % LinpDNA input concentration
        C_cleancap_T = u(7);     % cleancap input concentration
        C_magnesium_T = u(8);    % Magnesium ion input concentration
        C_nucleo_T = u(9);       % Nucleotide input concentration
        

        % System parameters
        %k_m = 0.27; 
        k_m = 0.27;  

        % Degradation by Michaelis-Menten (Modified calculation)
        %Ea_DNase = C_linpDNA_Sy * 1;  % Simplified enzyme activity calculation
        %v_max = k_p * Ea_DNase;  % Maximum rate (this is still required for degradation rate)
        Ea_DNase = C_linpDNA_Sy * 1;
        v_max=1188000*Ea_DNase; %g/L*h
        v = (v_max * C_linpDNA_Sy) / (k_m + C_linpDNA_Sy); % Degradation rate

       
        % Reactor geometry
        d = 0.1084;               % Diameter of the reactor (m)
        A = pi / 4 * d^2;      % Cross-sectional area of the reactor (m^2)
        
        % Fluid height
        h = (V / 1000) / A;    % Height (m)
        
        % State derivatives
        dVdt = Q_T - Q_OUT; % Change in volume
        dC_mRNA_dt = (C_mRNA_T * Q_T - C_mRNA_Sy * Q_OUT - C_mRNA_Sy * (Q_T - Q_OUT)) / V; % Nucleotide change
        dC_impurities_dt = (C_impurities_T * Q_T - C_impurities_Sy * Q_OUT - C_impurities_Sy * (Q_T - Q_OUT)) / V; % Magnesium change
        dC_T7_dt = (C_T7_T * Q_T - C_T7_Sy * Q_OUT - C_T7_Sy * (Q_T - Q_OUT)) / V; % mRNA change
        dC_linpDNA_dt = (Q_T * C_linpDNA_T - Q_OUT * C_linpDNA_Sy)  / V - v*C_linpDNA_Sy ; % Linearized DNA change
        dC_cleancap_dt = (C_cleancap_T * Q_T - C_cleancap_Sy * Q_OUT - C_cleancap_Sy * (Q_T - Q_OUT)) / V; % cleancap change
        dC_magnesium_dt = (C_magnesium_T * Q_T - C_magnesium_Sy * Q_OUT - C_magnesium_Sy * (Q_T - Q_OUT)) / V; % Magnesium ion change
        dC_nucleo_dt = (C_nucleo_T * Q_T - C_nucleo_Sy * Q_OUT - C_nucleo_Sy * (Q_T - Q_OUT)) / V; % Nucleotide change
        dC_degrade_dt =  (v*C_linpDNA_Sy) - (Q_OUT * C_degrade_Sy / V); % Degradation product change
        


        % Return derivatives
        sys = [dVdt; dC_mRNA_dt; dC_impurities_dt; dC_T7_dt; dC_linpDNA_dt; dC_cleancap_dt; dC_magnesium_dt; dC_nucleo_dt; dC_degrade_dt];

    
    case 3 % Outputs
        % Return the outputs, including the new "degrade" state
     
        d = 0.1084;             % Reactor diameter (m)
        A = pi / 4 * d^2;      % Reactor cross-sectional area (m^2)
        h = (x(1) / 1000) / A; % Height of the fluid in the reactor (m)

        % Extract state
        %C_linpDNA_Sy = x(5); % Linearized plasmid DNA concentration
      C_linpDNA_T = u(6);
        % Outputs: Enzyme activity is proportional to plasmid concentration
        Ea_DNase = C_linpDNA_T * 1; % Enzyme activity (same as in tankEE_test)
        
        sys = [x(1); x(2); x(3); x(4); x(5); x(6); x(7); x(8); x(9); Ea_DNase; h];
    
    case {2, 4, 9}
        % No action for unused flags
        sys = [];
    
    otherwise
        % Handle unhandled flags
        error(['Unhandled flag = ', num2str(flag)]);
end
