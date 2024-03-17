function [MESH, d_Fn_t, ALE_velocity] = movemesh(MESH, Fluid_ReferenceNodes, Solid_Extension, TimeAdvanceS, uS, d_Fn, dim, dt, nLiter, ALE_velocity)
    % first iteration of time step the mesh has not additionally deformed
    if nLiter == 1     
        u_gamma= [];    
        duS = TimeAdvanceS.velocity(uS); % apply interface velocity  
        for k = 1 : MESH.Solid.dim
            tmp  = duS(MESH.Solid.numNodes*(k-1)+MESH.Solid.dof_interface{k});
            u_gamma = [u_gamma; tmp(MESH.Interface_SFmap{k})];
        end
        ALE_velocity(MESH.Fluid.Gamma_global) = u_gamma;
        d_Fn_t = d_Fn;
    else
        % deform Fluid mesh by Solid-Extension Mesh Motion technique
        d_F = FSI_SolidExtension(MESH, uS, Solid_Extension);

        % update the fluid node/vertice coordinates
        Fluid_def_nodes    = Fluid_ReferenceNodes + d_F;
        Fluid_def_vertices = Fluid_def_nodes(1:dim, 1:MESH.Fluid.numVertices);

        d_F = reshape(d_F',dim*MESH.Fluid.numNodes,1);

        % Compute Fluid mesh velocity: w = 1/dt * ( d_f^(n+1) - d_f^n )
        % This part should be checked! inconsistent with BDF integrator
        ALE_velocity  =  1/dt * ( d_F - d_Fn );
        
        % Set the ALE mesh velocity at the interface (solid boundary vel)
        duS = TimeAdvanceS.velocity(uS);      
        for k = 1 : MESH.Solid.dim
            tmp  = duS(MESH.Solid.numNodes*(k-1)+MESH.Solid.dof_interface{k});
            ALE_velocity(MESH.Fluid.numNodes*(k-1)+MESH.Fluid.dof_interface{k}) = tmp(MESH.Interface_SFmap{k});
        end
        
        % output mesh displacement array 
        d_Fn_t          =  d_F;

        % update Fluid MESH (for NS solver)
        MESH.Fluid.vertices = Fluid_def_vertices;
        MESH.Fluid.nodes    = Fluid_def_nodes;

        % update mesh jacobian, determinant and inverse
        [MESH.Fluid.jac, MESH.Fluid.invjac, MESH.Fluid.h] = geotrasf(dim, MESH.Fluid.vertices, MESH.Fluid.elements);
    end
end