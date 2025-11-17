function [sys, X_0, str, ts] = Alkalinelysis(t, x, u, flag, X_0)
% PFR Function simulating the behavior of a Plug Flow Reactor (PFR)
% Inputs:
% t - time
% x - state variables
% u - input variables
% flag - operation flag determining what the function should do
% X_0 - initial conditions of the system

persistent C_NaOH_F_constant; % Persistent variable to store C_NaOH_F

switch flag
    case 0 % initialization;
        % Initialize the system sizes and state variables
        
        sizes = simsizes; 
        
        sizes.NumContStates = 6;        % Number of continuous states: Volume, C_X_S, C_DIS_S, C_plasmid_S, C_de_S, C_endo_S
        sizes.NumDiscStates = 0;        % No discrete states
        sizes.NumOutputs = 8;           % Number of outputs: 6 states + Alkaline_buffer and reactor length
        sizes.NumInputs = 4;            % Number of inputs: C_X_F, Q_IN, Q_OUT, C_NaOH_F
        sizes.DirFeedthrough = 0;       % No direct feedthrough from input to output
        sizes.NumSampleTimes = 1;       % Only one sample time
        
        % Set system sizes
        sys = simsizes(sizes);
        
        % Set initial conditions for all state variables
        X_0 = [2.70; 0.52; 0.099; 0.00001; 0.00001; 0.00001]; % Initial values for [V, C_X_S, C_DIS_S, C_plasmid_S, C_de_S, C_endo_S]
        
        % Not used (string for state descriptions)
        str = [];
        
        % Sample time [0, 0] (continuous system)
        ts = [0 0];
        
        % Initialize the persistent variable with a default value if not already set
        if isempty(C_NaOH_F_constant)
            C_NaOH_F_constant = 0.2; % Default initial value for NaOH concentration
        end
    
    case 1 % derivatives (calculate rate of change of states)
        % Extract current state variables
        V = x(1);           % Volume of the reactor (L)
        C_X_S = x(2);       % Concentration of biomass (g/L)
        C_DIS_S = x(3);     % Concentration of disruoted biomass (g/L)
        
        % Extract input variables
        C_X_F = u(1);       % Concentration of biomass in the feed (g/L)
        Q_IN = u(2);        % Feed flow rate (L/h)
        Q_OUT = u(3);       % Outflow rate (L/h)

        % Set C_NaOH_F from input only during the first call or at t = 0
        if isempty(C_NaOH_F_constant) || t == 0
            C_NaOH_F_constant = u(4); % Update NaOH concentration constant from the input
        end
        
        % Constants
        k_lysis = 20000 * C_NaOH_F_constant;    % Lysis rate constant (depends on NaOH concentration)
        Y_plasmid = 0.0011;                     % Plasmid yield (g plasmid / g biomass)
        Y_de = 0.0003;                          % Denatured plasmid yield (g plasmid / g biomass)
        Y_endo = 0.0275;                        % Endotoxin yield (g endotoxin / g plasmid)

        % State derivatives (rate of change of states)
        dVdt = (Q_IN - Q_OUT);                  % Change in volume (based on flow rates)
        
        % Cell lysis concentration rate (including input and output flow)
        dC_cell_Sdt = -k_lysis * C_X_S + (Q_IN * C_X_F - Q_OUT * C_X_S - C_X_S * (Q_IN - Q_OUT) / V);
        
        % Cell disrupted concentration rate (depends on biomass lysis)
        dC_dis_Sdt = k_lysis * C_X_S - (Q_OUT * C_DIS_S) / V;
        
        % Plasmid production rate (based on plasmid concentration and yield)
        dC_plasmid_Sdt = dC_dis_Sdt * Y_plasmid;
        
        % Denatured plasmid production rate (based on plasmid concentration and yield)
        dC_de_Sdt = dC_dis_Sdt * Y_de;
        
        % Endotoxin production rate (based on plasmid concentration and yield)
        dC_endo_Sdt = dC_plasmid_Sdt * Y_endo;

        % Return the derivatives of the states
        sys = [dVdt; dC_cell_Sdt; dC_dis_Sdt; dC_plasmid_Sdt; dC_de_Sdt; dC_endo_Sdt];
    
    case 3 % outputs (return the system outputs)
        % Reactor geometry and calculation of reactor length
        d = 0.099;          % Reactor diameter (m)
        A = pi * (d/2)^2;   % Cross-sectional area of the reactor (m^2)
        L = x(1)/1000 / A;  % Reactor length (m), calculated based on volume and cross-sectional area

        % Return outputs: Volume, concentrations, NaOH buffer concentration, and reactor length
        sys = [x(1); x(2); x(3); x(4); x(5); x(6); C_NaOH_F_constant; L];

    case {2, 4, 9}
        % No action for unused flags
        sys = [];

    otherwise
        % Handle any unhandled flags by returning an error message
        error(['unhandled flag = ', num2str(flag)]);
end
