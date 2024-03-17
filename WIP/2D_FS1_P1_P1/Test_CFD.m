% This solver was largely built starting from the redbkit monolithic FSI code which was
% deconstructed to its baseline elements, reassembled and added to (as needed), and used 
% to solve specific FSI benchmarks in a Partitioned fashion for subsequent coupling analysis

clear all; close all; clc;

%% Initialize problem and solver data 
init_data;

use_CML = false; %true - to confirm CML Navier-Stokes integration

%% Build finite element meshes and initialize boundary conditions
[MESH.Fluid, FE_SPACE_v, FE_SPACE_p] = build_F( dim, meshFluid.elements, meshFluid.vertices, meshFluid.boundaries, fem_F{1}, quad_order, DATA.Fluid, 'CFD', meshFluid.rings );
[MESH.Solid, FE_SPACE_s ] = build_S( dim, meshSolid.elements, meshSolid.vertices, meshSolid.boundaries, fem_S, quad_order, DATA.Solid, 'CSM', meshSolid.rings );
[FE_SPACE_g] = build_G( MESH.Fluid, fem_F{1}, dim, quad_order );

MESH.Fluid.internal_dof_c{MESH.dim+1} = 1:FE_SPACE_p.numDof;

%% Definte interface mapping and segregate domains (internal/interface)
[MESH] = interface(DATA, MESH, FE_SPACE_v, FE_SPACE_p, FE_SPACE_s, dim);

% undisturbed/orginal node coordinates for mesh updating
Fluid_ReferenceNodes = MESH.Fluid.nodes(1:dim,:);

%% Initialize fluid and solid domain variables
init_fluid;
init_solid;

%% Initialize coupling and acceleration scheme variables
init_couple;

%% Initialize variables if CML solver used
if use_CML
    CML_init;
end
%% Initalize Geometry Time Advance
[MESH, DATA, Solid_Extension] = meshjac(DATA, MESH, FE_SPACE_v, FE_SPACE_g, dim);

%% Initialize Linear/Non-Linear Solver Settings
tol        = DATA.Solid.NonLinearSolver.tol;
maxIter    = DATA.Solid.NonLinearSolver.maxit;

uS_n = u0; 

%% Time Loop
while ( t < tf )
    
    iter_time = tic;
    cellfun(@increase_time,objects); % update current step info for all sub-functions 
        
    t   = t   + dt;
    k_t = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );
    
    if use_CML
        % Predict the solution
        Sol.u = Sol.uPrev ;
        Sol.uDot = (pmc.gamma - 1)/pmc.gamma * Sol.uDotPrev ;
        Sol.p = Sol.pPrev ;
        Sol.aleDisp = Sol.aleDispPrev ;
        Sol.aleVel = Sol.aleVelPrev ;
    end

    % Enter the non-linear FP iteration loop    
    for nLiter = 1:Couple.imax
               
        %clear MESH.Fluid.vertices MESH.Fluid.nodes MESH.Fluid.jac MESH.Fluid.invjac MESH.Fluid.h
        % move the ALE mesh per the solid displacement and compute the ALE velocity
        if nLiter > 1     
            [MESH, d_Fn_t, ALE_velocity] = movemesh(MESH, Fluid_ReferenceNodes, Solid_Extension, TimeAdvanceS, uS, d_Fn, dim, dt, nLiter, ALE_velocity);
        else
            [~, ~, ALE_velocity] = movemesh(MESH, Fluid_ReferenceNodes, Solid_Extension, TimeAdvanceS, uS, d_Fn, dim, dt, nLiter, ALE_velocity);
        end
        
        % solve the Navier-Stokes equations (whatever formulation)
        if use_CML
                if nLiter > 1
                    Sol.aleDisp(:,1) = d_Fn_t(1:MESH.Fluid.numNodes);
                    Sol.aleDisp(:,2) = d_Fn_t(MESH.Fluid.numNodes+1:end);
                    
                    crdNew = MESH.Fluid.nodes';
                    crdNew = [crdNew zeros(MESH.Fluid.numNodes,1)];
                end
                
                Sol.aleVel(:,1) = ALE_velocity(1:MESH.Fluid.numNodes);
                Sol.aleVel(:,2) = ALE_velocity(MESH.Fluid.numNodes+1:end);
                
                [Sol, NSnormIndicator] = navierStokes(solver, fluid, pmc, Sol, cnn, crdNew, ndof, nen, nElem, MESH);
         
                [~, trac] = IntegratedOutput(Sol, crdNew, fluid, cnn, nen, ndof, MESH);
                                
        else
            [uF, C_NS, F_NS, trac] =  solve_fluid(MESH, DATA, FE_SPACE_v, FE_SPACE_p, TimeAdvanceF, TimeAdvanceS, t, dt, ALE_velocity, uS);
        end
        
        if nLiter > 4
            break
        end
    end
    
    if use_CML
        % Copy current variables to previous variables
        Sol.uPrev = Sol.u ;
        Sol.uDotPrev = Sol.uDot ;
        Sol.pPrev = Sol.p ;
        Sol.aleDispPrev = Sol.aleDisp ;
        Sol.aleVelPrev = Sol.aleVel ;
        
        uF = [Sol.u(:,1);Sol.u(:,2);Sol.p];
    end
    
%% Store current time-step data and update time-integration variables
    d_Fn = d_Fn_t;
    uFprev = uF;
    
    % update time advance variables for integration 
    TimeAdvanceS.Update( uS );
    TimeAdvanceF.Append( uF(1:FE_SPACE_v.numDof) );
     
%% Export Domain and Iteration Data 
   
    if (mod(k_t,Couple.probeOutFreq)==0)       
        % Export watchpoint physical data
        fprintf(fileId3,'%d %d %d\n',t,uF(WidF,1),uF(MESH.Fluid.numNodes+WidF,1));
        
        % Export lift and drag coefficients
        [~] = AeroForces(MESH, DATA, k_t, t, uF, C_NS, F_NS, FE_SPACE_v, FE_SPACE_p, dofs_drag, fileDragLift);
        
        % Fluid velocity and pressure on deformed mesh    
        if ~isempty(vtk_filename)
            CFD_export_solution(dim, uF(1:FE_SPACE_v.numDof), uF(1+FE_SPACE_v.numDof:end), ...
                MESH.Fluid.vertices, MESH.Fluid.elements, MESH.Fluid.numNodes, [vtk_filename,'Fluid'], k_t);
        end

    end
    
    % view time step duration for monitoring 
    iter_time = toc(iter_time);
    fprintf('\n-------------- iter time: %3.2f s -----------------',iter_time);
end
fclose(fileDragLift);
fclose(fileId1);
fclose(fileId2);
fclose(fileId3);