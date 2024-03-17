clc
clear all

wrkDir = './' ;
problemString = 'rigidV3' ;
elemType = '4Quad' ;
problemType = '2D' ;

% Coordinate, connectivity and boundary data
dataStr = strcat(wrkDir,'Data_rigid_cylV3.mat') ;
load(dataStr);
ndof = size(crd,1) ;

% Nonlinear iteration data:
solver.nLIterMin = 2 ;
solver.nLIterMax = 10 ;
solver.dt = 0.1 ;
solver.maxSteps = 2000 ;
solver.rhoinfty = 0.0 ;
solver.nLTol = 5e-4 ;
solver.outFreq = 100 ;
solver.intgOutFreq = 1 ;

% Fluid properties:
fluid.dens = 1000.0 ;
fluid.visc = 5.0 ;
fluid.gravFrc = [0, 0];

% Structure properties:
solid.mass = 7853.9816 ;
% damping as [Cxx Cyy Cxy];
solid.damping = [0.0 0.0 0.0];
% stiffness as [Kxx Kyy Kxy];
solid.stiffness = [12402.51062 12402.51062 0.0];
% direction of movement
solid.dom = [1; 1];
solid.gravity = [0; 0];

% Initial boundary conditions:
fluid.vel0 = [1, 0];
fluid.pres0 = 0.0 ;

% Boundary conditions:
fluid.Bc_TopN = size(unique(BCTop),1) ;
fluid.Bc_TopNs = unique(BCTop);
fluid.Bc_TopV = 0;

fluid.Bc_BotN = size(unique(BCBottom),1) ;
fluid.Bc_BotNs = unique(BCBottom);
fluid.Bc_BotV = 0;

fluid.Bc_LeftN = size(unique(BCLeft),1) ;
fluid.Bc_LeftNs = unique(BCLeft);
fluid.Bc_LeftV = [1 0];

fluid.Bc_CylN = size(unique(BCCyl),1) ;
fluid.Bc_CylNs = unique(BCCyl);
fluid.Bc_CylV = [0 0];

fluid.DirichletU = [fluid.Bc_LeftNs; fluid.Bc_CylNs];
fluid.DirichletUval = [fluid.Bc_LeftV(1).*ones(fluid.Bc_LeftN,1); fluid.Bc_CylV(1).*ones(fluid.Bc_CylN,1)];
fluid.DirichletV = [fluid.Bc_LeftNs; fluid.Bc_TopNs; fluid.Bc_BotNs; fluid.Bc_CylNs];
fluid.DirichletVval = [fluid.Bc_LeftV(2).*ones(fluid.Bc_LeftN,1); fluid.Bc_TopV.*ones(fluid.Bc_TopN,1); ...
                 fluid.Bc_BotV.*ones(fluid.Bc_BotN,1); fluid.Bc_CylV(2).*ones(fluid.Bc_CylN,1)];

             
% Gen-alpha time integration parameters
pmc.alphaM = 0.5*(3-solver.rhoinfty)/(1+solver.rhoinfty) ;
pmc.alpha = 1/(1+solver.rhoinfty) ;
pmc.gamma = 0.5 + pmc.alphaM - pmc.alpha ;

% Initialize other variables
nodeId = crd(:,1) ;
crd = crd(:,2:4) ;
cnn = conn(:,1:end) ;
nen = size(cnn,2) ;
nElem = size(cnn,1) ;
ndof = size(crd,1);

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

% Structure variables
ndofS = 1 ;
Sol.uSPrev = zeros(2,ndofS);
Sol.uS = zeros(2,ndofS);
Sol.dispSPrev = zeros(2,ndofS);
Sol.dispS = zeros(2,ndofS);

% ALE mesh variables
Sol.aleDisp = zeros(ndof,2);
Sol.aleDispPrev = zeros(ndof,2);
Sol.aleVel = zeros(ndof,2);
Sol.aleVelPrev = zeros(ndof,2);

crdNew = crd;

%% Time loop starts here
filename1 = sprintf('%s/%s.oisd',wrkDir,problemString);
fileId1 = fopen(filename1,'w');
for timeStep = 1:solver.maxSteps
    fprintf('Time step:%d\n',timeStep);
    
    % Predict the solution
    Sol.u = Sol.uPrev ;
    Sol.uDot = (pmc.gamma - 1)/pmc.gamma * Sol.uDotPrev ;
    Sol.p = Sol.pPrev ;
    Sol.uS = Sol.uSPrev ;
    Sol.dispS = Sol.dispSPrev ;
    Sol.aleDisp = Sol.aleDispPrev ;
    Sol.aleVel = Sol.aleVelPrev ;
    
    % Nonlinear iterations start here
    for nLiter = 1:solver.nLIterMax
        
        % Evaluate integrated values at the surface
        [Length, Force] = IntegratedOutput(Sol, crdNew, BCCyl, fluid, cnn);
        
        % Solve rigid body structural equation
        [Sol] = rigidBody(Sol, solver, solid, Force);
        
        % Solve ALE mesh equation to displace the nodes
        [Sol] = aleMesh(Sol, solver, solid, BCCyl, BCTop, BCBottom,...
                          BCLeft, BCRight, pmc, cnn, crd, elemType, nen,...
                          ndof, nElem);
                      
        clear crdNew              
        crdNew(:,1) = crd(:,1) + Sol.aleDisp(:,1);
        crdNew(:,2) = crd(:,2) + Sol.aleDisp(:,2);
        crdNew(:,3) = crd(:,3);
        
        % Solve Navier-Stokes
        [Sol, NSnormIndicator] = navierStokes(solver, fluid, pmc, Sol, cnn,... 
                                              crdNew, elemType, ndof, nen, nElem,...
                                              BCCyl);
        
        % Check convergence criteria
        if (NSnormIndicator < solver.nLTol)
            break;
        end
    end
 
    fprintf('\n');
    
    % Copy current variables to previous variables
    Sol.uPrev = Sol.u ;
    Sol.uDotPrev = Sol.uDot ;
    Sol.pPrev = Sol.p ;
    Sol.uSPrev = Sol.uS ;
    Sol.dispSPrev = Sol.dispS ;
    Sol.aleDispPrev = Sol.aleDisp ;
    Sol.aleVelPrev = Sol.aleVel ;
    
    % Evaluate integrated values at the surface
    crdNew(:,1) = crd(:,1) + Sol.aleDisp(:,1);
    crdNew(:,2) = crd(:,2) + Sol.aleDisp(:,2);
    crdNew(:,3) = crd(:,3);
    [Length, Force] = IntegratedOutput(Sol, crdNew, BCCyl, fluid, cnn);
    
    
    % Post-process the results
    if (mod(timeStep,solver.intgOutFreq)==0)
       fprintf(fileId1,'====timeStep: %d====\n',timeStep);
       fprintf(fileId1,'Metric:\n');
       fprintf(fileId1,'%e\n',Length);
       fprintf(fileId1,'Traction:\n');
       fprintf(fileId1,'%e %e\n',Force(1),Force(2));
       fprintf(fileId1,'Rigid body Displacement:\n');
       fprintf(fileId1,'%e %e\n',Sol.dispS(1),Sol.dispS(2));
    end
    
    
    if (mod(timeStep,solver.outFreq)==0)
        data = [crd(:,1)+Sol.aleDisp(:,1) crd(:,2)+Sol.aleDisp(:,2) crd(:,3) Sol.u Sol.p Sol.aleDisp];
        filename = sprintf('%s/%s.%d.plt',wrkDir,problemString,timeStep);
        fileId = fopen(filename,'w');

        FileTitle = 'Tecplot data';

        fprintf(fileId,' TITLE = \"%s\"\n',FileTitle);
        fprintf(fileId,' VARIABLES = \"X\", \"Y\", \"Z\" \"U\", \"V\", \"p\", \"dispX\", \"dispY\"\n');
    
        if (strcmp(elemType,'3Tri'))
            fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',(timeStep-1)*solver.dt,size(crd,1),size(conn,1));
        elseif (strcmp(elemType,'4Quad'))
            fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n',(timeStep-1)*solver.dt,size(crd,1),size(conn,1));
        end
    
        fprintf(fileId,'%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n',data');

        if (strcmp(elemType,'3Tri'))
            fprintf(fileId,'%d %d %d\n',conn');
        elseif (strcmp(elemType,'4Quad'))
            fprintf(fileId,'%d %d %d %d\n',conn');
        end
        fclose(fileId);
    end
    
end
fclose(fileId1);