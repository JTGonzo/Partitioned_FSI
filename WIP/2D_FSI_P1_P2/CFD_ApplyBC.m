function [A_in, F_in, u_D, u_gamma] =  CFD_ApplyBC(A, F, FE_SPACE, FE_SPACE_p, MESH, DATA, t, TimeAdvanceS, uS, zero_Dirichlet)
%% Initialize some function parameters
if nargin < 7
    t = [];
end

if isempty(A)
    A = sparse(FE_SPACE.numDof + FE_SPACE_p.numDof, FE_SPACE.numDof + FE_SPACE_p.numDof);
end

if isempty(F)
    F = sparse(FE_SPACE.numDof + FE_SPACE_p.numDof, 1);
end

if nargin < 10
    zero_Dirichlet = 0;
end

param = DATA.param;
u_D = [];
u_gamma = [];

%% Pressure condition
for k = 1 : MESH.Fluid.dim
        if ~isempty(MESH.Fluid.Pressure_side{k}) 
            % evaluate the surface quadrature nodes and basis function
            [csi,wi,phi]  =  bound_quad(FE_SPACE.quad_order, 0, 1);
            eta            =  1 - csi;
            nqn            =  length(csi);
                
            nof         = length(MESH.Fluid.Pressure_side{k});
            nbn         = MESH.Fluid.numBoundaryDof;
                
            Rrows       = zeros(nbn*nof,1);
            Rcoef       = Rrows;
            
            % compute global coordinates of quadrature nodes 
            xlt = zeros(nof,nqn); ylt = xlt;
            coord_ref = [eta; csi];
            for j = 1 : 2
                dof = MESH.Fluid.boundaries(j,MESH.Fluid.Pressure_side{k});
                vtemp = MESH.Fluid.vertices(1,dof);
                xlt = xlt + vtemp'*coord_ref(j,:);
                vtemp = MESH.Fluid.vertices(2,dof);
                ylt = ylt + vtemp'*coord_ref(j,:);
            end
            
            % Evaluate/import applied pressures at the quadrature nodes
            pressure = DATA.bcPrex(xlt,ylt,t,param);
            one       = ones(nof,nqn);
            pressure = pressure.*one;
            
            % compute the boundary elemet length    
            x    =  MESH.Fluid.vertices(1,MESH.Fluid.boundaries(1:2, MESH.Fluid.Pressure_side{k}));
            y    =  MESH.Fluid.vertices(2,MESH.Fluid.boundaries(1:2, MESH.Fluid.Pressure_side{k}));    
            side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
            
            % compute the boundary element pressure force acting on
            % the face normal
            for l = 1 : nof
                face = MESH.Fluid.Pressure_side{k}(l);
                    
                 pressure_loc  = pressure(l,:).*wi;
                 pressure_loc  = pressure_loc(1,:)';
                    
                 Rrows(1+(l-1)*nbn:l*nbn) = MESH.Fluid.boundaries(1:nbn,face);
                 Rcoef(1+(l-1)*nbn:l*nbn) = MESH.Fluid.Normal_Faces(k,face)*side_length(l)*phi*pressure_loc;
            end
             %  sparse vector of boundary pressure forces
                F = F + sparse(Rrows+(k-1)*MESH.Fluid.numNodes,1,Rcoef,FE_SPACE.numDof + FE_SPACE_p.numDof,1);                
        end
end

%% Dirichlet condition
% extract coordinates of dirichlet nodes and evaluate the BC
for k = 1 : MESH.Fluid.dim
      x  = MESH.Fluid.nodes(1,MESH.Fluid.Dirichlet_dof_c{k});
      y  = MESH.Fluid.nodes(2,MESH.Fluid.Dirichlet_dof_c{k});
      u_Dirichlet{k} = DATA.bcDir{k}(x,y,t,param);
      u_D = [u_D; u_Dirichlet{k}'];
end

% either a zero or non-zero enforced condition
u_D  = u_D * (1 - zero_Dirichlet);

% extract the solid interface velocity  
duS = TimeAdvanceS.velocity(uS);
for k = 1 : MESH.Fluid.dim
    tmp  = duS(MESH.Solid.numNodes*(k-1)+MESH.Solid.dof_interface{k});
    u_gamma = [u_gamma; tmp(MESH.Interface_SFmap{k})];
end

% boundary and interface Dirichlet conditions
u_BC = [u_D; u_gamma];

%% Apply the BCs to the LHS/RHS matrices
% output the reduced LHS and RHS NS-matrices for the free internal nodes 
F_in = F(MESH.Fluid.II_global) - A(MESH.Fluid.II_global,[MESH.Fluid.Dirichlet_dof; MESH.Fluid.Gamma_global])*u_BC;
A_in = A(MESH.Fluid.II_global,MESH.Fluid.II_global);

end
