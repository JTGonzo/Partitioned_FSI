%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

%DATAFILE Fluid
U_bar = 51.3;

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y);
data.force{2} = @(x, y, z, t, param)(0.*x.*y);
data.force{3} = @(x, y, z, t, param)(0.*x.*y);

% Dirichlet
data.bcDir_t  = @(t)( 0.5*(1-cos(pi/0.2*t)).*(t<0.0) + 1 * (t>=0.0));

data.bcDir{1} = @(x, y, z, t, param)(data.bcDir_t(t) * U_bar.*(x==0) + 0.*x.*y); 
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y); 
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y); 

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{3} = @(x, y, z, t, param)(0.*x.*y);

% initial condition
data.u0{1}    = @(x, y, z, t, param)(0.*x.*y);
data.u0{2}    = @(x, y, z, t, param)(0.*x.*y);
data.u0{3}    = @(x, y, z, t, param)(0.*x.*y);

% BC flag
data.flag_dirichlet{1}    =  [1 5];
data.flag_neumann{1}      =  [2 3 4];
data.flag_FSinterface{1}  =  [6];
data.flag_ring{1}         =  []; 
data.flag_ALE_fixed{1}    =  [1 5 2];


data.flag_dirichlet{2}    =  [1 5 4];
data.flag_neumann{2}      =  [2 3];
data.flag_FSinterface{2}  =  [6];
data.flag_ring{2}         =  []; 
data.flag_ALE_fixed{2}    =  [1 5 4];

data.flag_dirichlet{3}    =  [1 5 3];
data.flag_neumann{3}      =  [2 4];
data.flag_FSinterface{3}  =  [6];
data.flag_ring{3}         =  []; 
data.flag_ALE_fixed{3}    =  [1 5 3];

% Model parameters
data.dynamic_viscosity  =   1.82*1e-4;
data.density              =   1.18*1e-3;
%data.Stabilization        =   'SUPG';

% Nonlinear solver
data.NonLinearSolver.tol         = 1e-6; 
data.NonLinearSolver.maxit       = 30;
 
% Linear solver
%   If parallel pool available, use gmres with AdditiveSchwarz preconditioner
%   Otherwise use direct solver
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end

% if poolsize > 0
%     
%     % Linear Solver
%     data.LinearSolver.type              = 'gmres'; % MUMPS, backslash, gmres
%     data.LinearSolver.tol               = 1e-8;
%     data.LinearSolver.maxit             = 500;
%     data.LinearSolver.gmres_verbosity   = 5;
%     data.LinearSolver.mumps_reordering  = 4;
%     
%     % Preconditioner
%     data.Preconditioner.type              = 'AdditiveSchwarz'; % AdditiveSchwarz, None, ILU
%     data.Preconditioner.local_solver      = 'MUMPS'; % matlab_lu, MUMPS
%     data.Preconditioner.overlap_level     = 2;
%     data.Preconditioner.mumps_reordering  = 7;
%     data.Preconditioner.num_subdomains    = poolsize; %poolsize, number of subdomains
%     data.Preconditioner.coarse_level             = 'None'; % None, Aggregation, SmoothedAggregation
%     data.Preconditioner.coarse_num_aggregates    = 200;
%     data.Preconditioner.coarse_smoother_iter     = 1;
%     data.Preconditioner.coarse_smoother_dumping  = 1;
% 
% 
% else
    
    % Linear Solver
    data.LinearSolver.type              = 'backslash'; % MUMPS, backslash, gmres
    data.LinearSolver.mumps_reordering  = 7;

    % Preconditioner
    data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU
% end


%% Time Setting
data.time.BDF_order  = 2;
data.time.t0         = 0;
data.time.dt         = 0.001; 
data.time.tf         = 5;
data.time.nonlinearity  = 'semi-implicit';

%% Output options
data.Output.DragLift.computeDragLift = 1;
data.Output.DragLift.factor          = 2/(1.18*1e-3*51.3^2*1);
data.Output.DragLift.flag            = [5 6];
data.Output.DragLift.filename        = 'Results/AeroForces.txt';