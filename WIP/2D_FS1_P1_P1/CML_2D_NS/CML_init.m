%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Integrating CML Navier Stokes                  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nonlinear iteration data:
solver.nLIterMin = Couple.imin ;
solver.nLIterMax = Couple.imax ;
solver.dt = DATA.Fluid.time.dt ;
solver.maxSteps = Couple.maxSteps ;
solver.rhoinfty = 0.0 ;
solver.nLTol = Couple.rtol;

% Fluid properties:
fluid.dens = DATA.Fluid.density;
fluid.visc = DATA.Fluid.dynamic_viscosity;
fluid.gravFrc = [0, 0];

% Initial boundary conditions:
fluid.vel0 = [-U_bar, 0];  
fluid.pres0 = 0.0 ;
            
% Gen-alpha time integration parameters
pmc.alphaM = 0.5*(3-solver.rhoinfty)/(1+solver.rhoinfty) ;
pmc.alpha = 1/(1+solver.rhoinfty) ;
pmc.gamma = 0.5 + pmc.alphaM - pmc.alpha ;

% Initialize other variables
crd = MESH.Fluid.nodes';
crd = [crd zeros(MESH.Fluid.numNodes,1)];
cnn = MESH.Fluid.elements(1:3,:)';
nen = MESH.Fluid.numElemDof;
nElem = MESH.Fluid.numElem;
ndof = MESH.Fluid.numNodes;

%% Initialize Fluid Boundary conditions
fluid.DirichletU = MESH.Fluid.Dirichlet_dof_c{1};
fluid.DirichletUval = v0(MESH.Fluid.Dirichlet_dof_c{1});

fluid.DirichletV = MESH.Fluid.Dirichlet_dof_c{2};
fluid.DirichletVval = v0(MESH.Fluid.numNodes+MESH.Fluid.Dirichlet_dof_c{2});

% Fluid variables
Sol.u = zeros(ndof,2,1);
Sol.uDot = zeros(ndof,2,1);
Sol.u(:,1,:) = fluid.vel0(1) ;
Sol.u(:,2,:) = fluid.vel0(2) ;
Sol.uAlpha = zeros(ndof,2,1) ;
Sol.uDotAlpha = zeros(ndof,2,1) ;
Sol.uPrev = Sol.u ;
Sol.uDotPrev = Sol.uDot ;
Sol.p = fluid.pres0*ones(ndof,1);
Sol.pPrev = Sol.p ;

% ALE mesh variables
Sol.aleDisp = zeros(ndof,2);
Sol.aleDispPrev = zeros(ndof,2);
Sol.aleVel = zeros(ndof,2);
Sol.aleVelPrev = zeros(ndof,2);

crdNew = crd;