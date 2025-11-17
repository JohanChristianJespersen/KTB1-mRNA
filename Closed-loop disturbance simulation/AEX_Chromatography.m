function [sys,x0,str,ts,simStateCompliance] = AEX_Chromatography(t,x,u,flag,gridsize)

switch flag,
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0,
        [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes(gridsize);
    %%%%%%%%%%%%%%%
    % Derivatives %
    %%%%%%%%%%%%%%%
    case 1,
        sys = mdlDerivatives(t,x,u,gridsize);
    %%%%%%%%%%
    % Update %
    %%%%%%%%%%
    case 2,
        sys = mdlUpdate(t,x,u);
    %%%%%%%%%%%%
    % Outputs % 
    %%%%%%%%%%%%
    case 3,
        sys = mdlOutputs(t,x,u,gridsize);
    %%%%%%%%%%%%%%%%%%%%%%%
    % GetTimeOfNextVarHit %
    %%%%%%%%%%%%%%%%%%%%%%%
    case 4,
        sys = mdlGetTimeOfNextVarHit(t,x,u);
    %%%%%%%%%%%%%
    % Terminate %
    %%%%%%%%%%%%%
    case 9,
        sys = mdlTerminate(t,x,u);
    %%%%%%%%%%%%%%%%%%%%
    % Unexpected flags %
    %%%%%%%%%%%%%%%%%%%%
    otherwise
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes(gridsize)

par = struct(...
    'components',3,...          % Number of components
    'nconc',3,...               % Number of concentrations
    'flow',1);                  % Flow rate as input

sizes = simsizes;
sizes.NumContStates  = par.nconc * par.components * gridsize;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = par.nconc * par.components * gridsize;
sizes.NumInputs      = par.components + par.flow; % Flow + components
sizes.DirFeedthrough = 1;                         % Inputs affect outputs
sizes.NumSampleTimes = 1;                         % At least one sample time is needed

sys = simsizes(sizes);

% Initialize the initial conditions
x0  = zeros(par.nconc * par.components * gridsize ,1); 
str = [];                       % str is always an empty matrix
ts  = [0 0];                    % Initialize sample timee

% Specify the block simStateCompliance
simStateCompliance = 'UnknownSimState';

%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
function sys = mdlDerivatives(t,x,u,gridsize)

% Parameters
par = struct(...
    'components',3,...          
    'length', 199E-2,...                            % [dm]           
    'epsbed',  0.35,...                             % no unit
    'epspar',  0.8,...                              % no unit
    'rad', 450E-6,...                               % Converteret til [dm]
    'keff',  [0.00409, 0.03223, 0.00409],...        % from cm/s to [dm/h] 
   ...
     'ads', [3.6E-10, 3.6E-7, 4.0],...             % Langmuir adsorption coefficient [dm^3/g/h]
     'des', [1.44E-13, 3.6E-10, 0.5],...           % Langmuir desorption coefficient [1/h]
     'qmax', [1.44E-6, 1.02E-6, 1.6]);              % Langmuir maximum concentration [g/dm^3]

nconc = 3;                                          % Number of concentrations
h = par.length / (gridsize - 1);                    % Length of each cell

% Extract flow and components
flow = u(1);                                        % Flow is the first input
conc = u(2:end);                                    % Remaining inputs are concentrations

% Velocity based on flow
par.velocity = flow / (50.26 * 0.35) * 100;         % Ensure correct units
par.dax=(par.rad*par.velocity*0.35)/0.09;           % Unit [dm^2/h]

% Prepare sparse matrices
M = h/6 * diag(ones(gridsize-1,1),1) + ...
    2*h/3 * diag(ones(gridsize,1)) + ...
    h/6 * diag(ones(gridsize-1,1),-1);    
M(gridsize,gridsize) = h/3; 
M(1,1) = 1; 
M(1,2) = 0;
M = sparse(M);

C = (1/2) * diag(ones(gridsize-1,1), 1) + ...
    - (1/2) * diag(ones(gridsize-1,1),-1); 
C(1,1) = -1/2; 
C(gridsize,gridsize) =  1/2; 
C = sparse(C);

A = -1/h * diag(ones(gridsize-1,1), 1) + ...
    2/h * diag(ones(gridsize,1)) + ...
    -1/h * diag(ones(gridsize-1,1),-1);     
A(1,1) = 1/h; 
A(gridsize,gridsize) = 1/h; 
A = sparse(A);

% Reshape states and inputs for vectorized calculations
c = reshape(x, gridsize, par.components, nconc);

% Calculate derivatives
dcdt = zeros(size(c)); 
for j = 1:par.components
    % Mass balance for the first concentration
    dcdt(:,j,1) = -(par.velocity/par.epsbed * C + par.dax * A) * c(:,j,1) - ...
                  (1-par.epsbed)/par.epsbed * par.keff(j) * 3.0 / par.rad * M * (c(:,j,1) - c(:,j,2));

    % Boundary condition
    dcdt(1,j,1) = dcdt(1,j,1) - par.velocity/par.epsbed * (c(1,j,1) - conc(j));
    
    % Langmuir model for the third concentration
    sum = max(0,1 - (c(:,3,1) / par.qmax(j))); 
    dcdt(:,j,3) = par.ads(j) * par.qmax(j) * sum .* c(:,j,2) - ...
                  par.des(j) * c(:,j,3);
    
    % Pore concentration
    dcdt(:,j,2) = -(1-par.epspar)/par.epspar * dcdt(:,j,3) + ...
                   par.keff(j) * 3 / par.rad / par.epspar * (c(:,j,1) - c(:,j,2));     
end

sys = reshape(dcdt, [], 1); % Output as a column vector

%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
function sys = mdlUpdate(t,x,u)

sys = []; % No updates needed

%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
function sys = mdlOutputs(t,x,u,gridsize)

ncomp = 3; % Number of components
nconc = 3; % Number of concentrations

% Reshape the state vector to ensure it matches expected output size
sys = reshape(x, gridsize, ncomp, nconc); % Shape it based on gridsize

sys = sys(:); % Ensure output is a column vector

%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.
%=============================================================================
function sys = mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1; % Set the next hit to be one second later.
sys = t + sampleTime;

%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
function sys = mdlTerminate(t,x,u)

sys = []; % No termination tasks needed