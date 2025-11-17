function [sys,x0,str,ts] = M_501(t,x,u,flag)
  
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes();

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case { 2, 4, 9 },
    sys = [];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end
% end csfunc

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes()
sizes = simsizes;
sizes.NumContStates  = 6;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 7;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [4.827 0.1650 0.001595 0.03772 0.0001 0.0001 ];
str = [];
ts  = [0 0];

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)

 %states
 V = x(1);
 C_LNP_tank = x(2);
 C_ethanol_tank=x(3);
 C_sucrose_tank=x(4);
 C_imp1_tank=x(5);
 C_imp2_tank=x(6);
      
 %inputs
 Q_T = u(1);
 C_LNP_T = u(2);
 C_ethanol_T=u(3);
 C_sucrose_T=u(4);
 C_imp1_T=u(5);
 C_imp2_T=u(6);
 Q_OUT = u(7);
 

 % parameter values
 dVdt=Q_T-Q_OUT;
 dC_LNP_dt=(C_LNP_T*Q_T-C_LNP_tank*Q_OUT-C_LNP_tank*(Q_T-Q_OUT))/V;
 dC_ethanol_dt=(C_ethanol_T*Q_T-C_ethanol_tank*Q_OUT-C_ethanol_tank*(Q_T-Q_OUT))/V;
 dC_sucrose_dt=(C_sucrose_T*Q_T-C_sucrose_tank*Q_OUT-C_sucrose_tank*(Q_T-Q_OUT))/V;
 dC_imp1_dt=(C_imp1_T*Q_T-C_imp1_tank*Q_OUT-C_imp1_tank*(Q_T-Q_OUT))/V;
 dC_imp2_dt=(C_imp2_T*Q_T-C_imp2_tank*Q_OUT-C_imp2_tank*(Q_T-Q_OUT))/V;
 

sys = [dVdt;dC_LNP_dt;dC_ethanol_dt; dC_sucrose_dt; dC_imp1_dt;dC_imp2_dt];

% end mdlDerivatives

function sys=mdlOutputs(t,x,u)

sys = [x(1);x(2);x(3);x(4);x(5); x(6)];

% end mdlOutputs