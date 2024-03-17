clear all; close all; clc;
[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');
addpath('C_Files/')
dim      =  2;

%% Load P1 mesh
load('beam.mat')
fem        = {'P2', 'P1'};
t      = [];
param = [];
vtk_filename = 'Figures/CFD_'; %[];
data_file = 'CFD_data';
use_SUPG = false;
quad_order   = 4;

%%
eval(data_file);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)
    
    eval(['DATA.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
    
end

DATA.param = param;

[ MESH, FE_SPACE_v, FE_SPACE_p] = build( dim, elements, vertices, boundaries, fem{1}, quad_order, DATA, 'CFD' );

MESH.internal_dof_c{MESH.dim+1} = 1:FE_SPACE_p.numDof;

totSize = FE_SPACE_v.numDof + FE_SPACE_p.numDof;

%% Gather Time Setting
BDF_order = DATA.time.BDF_order;
t0        = DATA.time.t0;
dt        = DATA.time.dt;
tf        = DATA.time.tf;
t         = DATA.time.t0;
k_t       = 0;

BDFhandler = BDF_TimeAdvance( BDF_order);

v0  = [];
for k = 1 : FE_SPACE_v.numComponents
  v0  = [v0; DATA.u0{k}(  MESH.nodes(1,:), MESH.nodes(2,:), t0, param )'];
end

u = [v0; zeros(FE_SPACE_p.numDof,1)];

if ~isempty(vtk_filename)
    CFD_export_solution(MESH.dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, 0);
end

BDFhandler.Initialize( v0 );
for bd = 2 : BDF_order
    BDFhandler.Append( v0 );
end

%% Create Fluid Assembler Object
FluidModel = CFD_Assembler( MESH, DATA, FE_SPACE_v, FE_SPACE_p );

%% Assemble Constant Terms
[A_Stokes] = FluidModel.compute_Stokes_matrix();

Mv = FluidModel.compute_mass_velocity();
Mp = FluidModel.compute_mass_pressure();
M  = blkdiag(DATA.density * Mv, 0*Mp);

%% PreProcessing for Drag and Lift Computation
compute_AerodynamicForces = 0 ;
if isfield(DATA, 'Output') && isfield(DATA.Output, 'DragLift')
    if DATA.Output.DragLift.computeDragLift == 1
        compute_AerodynamicForces = true;
    end
end

if compute_AerodynamicForces
    AeroF_x(k_t+1)  = 0;
    AeroF_y(k_t+1)  = 0;
    AeroF_z(k_t+1)  = 0;
    dofs_drag    = [];
    
    for j = 1 : length(DATA.Output.DragLift.flag)
        Dirichlet_side         = find(MESH.boundaries(MESH.bc_flag_row,:) == DATA.Output.DragLift.flag(j));
        Dirichlet_side         = unique(Dirichlet_side);
        Dirichlet_dof          = MESH.boundaries(1:MESH.numBoundaryDof,Dirichlet_side);
        dofs_drag              = [dofs_drag; Dirichlet_dof(:)];
    end
    dofs_drag = unique(dofs_drag);
    
    fileDragLift = fopen(DATA.Output.DragLift.filename, 'w+');
    fprintf(fileDragLift, 'Time          F_x          F_y          F_z');
    fprintf(fileDragLift, '\n%1.4e  %1.4e  %1.4e  %1.4e', t, AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
end

%% Time Loop
while ( t < tf )
    
    iter_time = tic;
    
    t       = t   + dt;
    k_t     = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n',t0,t,tf);
    
    v_BDF = BDFhandler.RhsContribute( );
    u_BDF = [v_BDF; zeros(FE_SPACE_p.numDof,1)];
    alpha = BDFhandler.GetCoefficientDerivative();
    
    v_extrapolated = BDFhandler.Extrapolate();
    U_k            = zeros(totSize,1);
            
    % Assemble matrix and right-hand side
    [C1] = FluidModel.compute_convective_Oseen_matrix( v_extrapolated );
            
    F_NS = 1/dt * M * u_BDF;
    C_NS = alpha/dt * M + A_Stokes + C1;
                       
    % Apply boundary conditions
    [A, b, u_D]   =  CFD_ApplyBC(C_NS, F_NS, FE_SPACE_v, FE_SPACE_p, MESH, DATA, t, 0, u);
            
    % Solve
%    fprintf('\n -- Solve A x = b ... ');
    U_k(MESH.internal_dof) = A \ b;
    U_k(MESH.Dirichlet_dof) = u_D;
            
    fprintf('\n -- Norm(U_np1 - U_n) / Norm( U_n ) = %1.2e \n', norm(U_k - u) / norm(u));
                
    u = U_k;
    
    %% Update BDF
    BDFhandler.Append( u(1:FE_SPACE_v.numDof) );
   
    %% Export to VTK
    if ~isempty(vtk_filename)
        CFD_export_solution(MESH.dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, k_t);
    end
    
        %% Compute_DragLift
    if compute_AerodynamicForces
        
        Z              = zeros(FE_SPACE_v.numDofScalar,1);
        Z(dofs_drag)   = 1;
        
        W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
        W(1:FE_SPACE_v.numDofScalar)        = Z;
        AeroF_x(k_t+1) = DATA.Output.DragLift.factor*(W'*(-C_NS*u + F_NS));
        
        W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
        W(FE_SPACE_v.numDofScalar+[1:FE_SPACE_v.numDofScalar])  = Z;
        AeroF_y(k_t+1)  = DATA.Output.DragLift.factor*(W'*(-C_NS*u  + F_NS));
        
        AeroF_z(k_t+1) = 0.0;
        
        fprintf('\n *** F_x = %e, F_y = %e, F_z = %e *** \n',  AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
        fprintf(fileDragLift, '\n%1.4e  %1.4e  %1.4e  %1.4e', t, AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
    end
end

if compute_AerodynamicForces
    fclose(fileDragLift);
end
%     
% % Postprocessing
% [t, Drag, Lift] = importAerodynamicForces( DATA.Output.DragLift.filename );
%              
% handle = figure;
% subplot(2,1,1)
% plot(t,Drag,'-b','LineWidth',2)
% xlabel('Time [s]')
% grid on
%         
% subplot(2,1,2)
% plot(t,Lift,'-r','LineWidth',2)
% xlabel('Time [s]')
% grid on
