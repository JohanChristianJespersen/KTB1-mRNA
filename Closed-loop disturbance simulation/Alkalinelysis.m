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
        sizes = simsizes;
        
        sizes.NumContStates = 6;        % Volume, C_X_S, C_DIS_S, C_plasmid_S, C_de_S, C_endo_S
        sizes.NumDiscStates = 0;
        sizes.NumOutputs = 8;           % 6 states + Alkaline_buffer + reactor length
        sizes.NumInputs = 4;            % C_X_F, Q_IN, Q_OUT, C_NaOH_F
        sizes.DirFeedthrough = 0;
        sizes.NumSampleTimes = 1;
        
        sys = simsizes(sizes);
        X_0 = [2.70; 0.52; 0.099; 0.00001; 0.00001; 0.00001]; % Initial conditions
        
        str = [];
        ts = [0 0];
        
        % Initialize persistent constant
        if isempty(C_NaOH_F_constant)
            C_NaOH_F_constant = 0.2;                % Default initial value
        end
    
    case 1 % derivatives
        % Extract current states
        V = x(1);           % Volume of the reactor
        C_X_S = x(2);       % Concentration of cell
        C_DIS_S = x(3);     % Concentration of plasmid  
        
        % Extract input variables
        C_X_F = u(1);       % Concentration of biomass in feed
        Q_IN = u(2);        % Feed flow rate
        Q_OUT = u(3);       % Outflow rate

        % Set C_NaOH_F from input only during the first call
        if isempty(C_NaOH_F_constant) || t == 0
            C_NaOH_F_constant = u(4); % Assign C_NaOH_F as a constant
        end
        
        % Constants
        k_lysis = 20000*C_NaOH_F_constant;    % Lysis rate constant
        Y_plasmid = 0.0011;                 % g pDNA/g biomass (wcw) wet cell weight
        Y_de = 0.0003;                      % Denatured plasmid release rate constant g/g
        Y_endo = 0.0275;                    % Endotoxin release rate constant g/g

        % State derivatives (rate of change of states)
        dVdt = (Q_IN - Q_OUT);              % Change in volume
        
        % Reaction rates incorporating residence time (tau)
        dC_cell_Sdt = -k_lysis * C_X_S + (Q_IN * C_X_F - Q_OUT * C_X_S-C_X_S*(Q_IN-Q_OUT)/ V); % Change in biomass concentration

        dC_dis_Sdt = k_lysis * C_X_S - ((Q_OUT * C_DIS_S)/ V);
        
        dC_plasmid_Sdt = dC_dis_Sdt * Y_plasmid; % Change in plasmid concentration

        dC_de_Sdt = dC_dis_Sdt * Y_de;           % Change in denatured plasmid concentration

        dC_endo_Sdt = dC_plasmid_Sdt * Y_endo;   % Change in endotoxin concentration

        sys = [dVdt; dC_cell_Sdt; dC_dis_Sdt; dC_plasmid_Sdt; dC_de_Sdt; dC_endo_Sdt];
    
    case 3 % outputs
        % Reactor geometry and calculation of reactor length
        d = 0.099; % Reactor diameter (m)
        A = pi * (d/2)^2 ; % Reactor cross-sectional area (m^2)
        L = x(1)/1000 / A; % Length of the reactor (m)

        sys = [x(1); x(2); x(3); x(4); x(5); x(6); C_NaOH_F_constant; L]; % Volume, concentrations, Alkaline_buffer, and length

    case {2, 4, 9}
        sys = [];

    otherwise
        error(['unhandled flag = ', num2str(flag)]);
end