%DATAFILE for FS1 - Solid

% Source term
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(0.0 + 0.*x.*y);

% Dirichlet
data.bcDir{1} = @(x, y, t, param)(0.*x.*y); 
data.bcDir{2} = @(x, y, t, param)(0.*x.*y); 

% Neumann
data.bcNeu{1} = @(x, y, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, t, param)(0.*x.*y);

% Normal Pressure
data.bcPrex   = @(x, y, t, param)(0.*x.*y);

% BC flag
data.flag_dirichlet{1}    =  [6];
data.flag_neumann{1}      =  [];
data.flag_FSinterface{1}  =  [5];
data.flag_pressure{1}     =  [];

data.flag_dirichlet{2}    =  [6];
data.flag_neumann{2}      =  [];
data.flag_FSinterface{2}  =  [5];
data.flag_pressure{2}     =  [];

data.u0{1} = @(x, y, t, param)(0.*x.*y);
data.u0{2} = @(x, y, t, param)(0.*x.*y);
data.du0{1} = @(x, y, t, param)(0.*x.*y);
data.du0{2} = @(x, y, t, param)(0.*x.*y);

% material parameters 
data.Young   = 2.0*10^6; %2.5*10^6
data.Poisson = 0.35;
data.Density = 2.0; %0.1
data.Material_Model   = 'StVenantKirchhoff';%'StVenantKirchhoff', Linear, % Attention: Neohookean not available in 2D
data.model   = 'CSM';

% NonLinear Solver
data.NonLinearSolver.tol               = 1e-6;
data.NonLinearSolver.maxit             = 35;

% Time options
data.time.t0         = 0;
data.time.dt         = 0.001; 
data.time.tf         = 15;
data.time.gamma      = 1/2;
data.time.beta       = 1/4;