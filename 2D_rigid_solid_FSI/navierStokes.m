%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%        2D Incompressible Navier-Stokes Finite-element Solver            %
%                                                                         % 
%   \rho u_t + \rho u.div u - div(mu*grad u) + grad p = f in \Omega,      %
%                         div u = 0  in \Omega,                           %
%                                                                         %
%       Dirichlet boundary condition        u = g_D  on \Gamma_D,         %
%       Neumann boundary condition du/dn - np = g_N  on \Gamma_N.         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sol, NSnormIndicator] = navierStokes(solver, fluid, pmc, Sol, cnn, crd, ...
                                               elemType, ndof, nen, nElem, BCCyl)

% Quadrature rules for elements
if strcmp(elemType,'3Tri')
    gP = ...
   [1/3,  1/3
    4/3,  1/3
    1/3,  4/3] ;
    gW = ...
   [2/3,  2/3,  2/3] ;
 
    N(:,1) = 0.5.*(2.-gP(:,1)-gP(:,2)) ;
    N(:,2) = 0.5.*(gP(:,1)) ;
    N(:,3) = 0.5.*(gP(:,2)) ;
    
    Nx(:,1) = -0.5.*ones(3,1) ;
    Nx(:,2) =  0.5.*ones(3,1) ;
    Nx(:,3) =  zeros(3,1) ; 
    Ny(:,1) = -0.5.*ones(3,1) ;
    Ny(:,2) =  zeros(3,1) ;
    Ny(:,3) =  0.5.*ones(3,1) ;    
elseif strcmp(elemType,'4Quad')
    gP = ...
   [-5.7735026918962584E-01, -5.7735026918962584E-01
     5.7735026918962584E-01, -5.7735026918962584E-01
    -5.7735026918962584E-01,  5.7735026918962584E-01
     5.7735026918962584E-01,  5.7735026918962584E-01] ;
    gW = [1, 1, 1, 1 ] ;
    
    N(:,1) = 0.25.*(1-gP(:,1)).*(1-gP(:,2)) ;
    N(:,2) = 0.25.*(1+gP(:,1)).*(1-gP(:,2)) ;
    N(:,3) = 0.25.*(1+gP(:,1)).*(1+gP(:,2)) ;
    N(:,4) = 0.25.*(1-gP(:,1)).*(1+gP(:,2)) ;
    
    Nx(:,1) = -0.25.*(1-gP(:,2)) ;
    Nx(:,2) =  0.25.*(1-gP(:,2)) ;
    Nx(:,3) =  0.25.*(1+gP(:,2)) ;
    Nx(:,4) = -0.25.*(1+gP(:,2)) ;
    Ny(:,1) = -0.25.*(1-gP(:,1)) ;
    Ny(:,2) = -0.25.*(1+gP(:,1)) ;
    Ny(:,3) =  0.25.*(1+gP(:,1)) ;
    Ny(:,4) =  0.25.*(1-gP(:,1)) ;
end

Nx = Nx' ;
Ny = Ny' ;
nQuad = length(gW) ;
 
iif = zeros(nen^2*nElem,1); 
jjf = zeros(nen^2*nElem,1);
index = 0;
for i = 1:nen
   for j = 1:nen
      iif(index+1:index+nElem) = double(cnn(:,i)); 
      jjf(index+1:index+nElem) = double(cnn(:,j));  
      index = index + nElem;
   end
end
        
% Satisfy boundary conditions
Sol.u(fluid.DirichletU,1,1) = fluid.DirichletUval ;
Sol.u(fluid.DirichletV,2,1) = fluid.DirichletVval ;
Sol.u(unique(BCCyl(:)),1) = (Sol.aleVel(unique(BCCyl(:)),1)) ;
Sol.u(unique(BCCyl(:)),2) = (Sol.aleVel(unique(BCCyl(:)),2)) ;
 
% Interpolate for alpha values for Gen-alpha
Sol.uAlpha = Sol.uPrev + pmc.alpha.*(Sol.u - Sol.uPrev) ;
Sol.uDotAlpha = Sol.uDotPrev + pmc.alphaM.*(Sol.uDot - Sol.uDotPrev) ;
Sol.aleVelAlpha = Sol.aleVelPrev + pmc.alpha.*(Sol.aleVel - Sol.aleVelPrev) ;
        
% Navier-Stokes equations
xxf = zeros(size(cnn));
yyf = zeros(size(cnn));
ux = zeros(size(cnn));
uy = zeros(size(cnn));

for i=1:nen
   xxf(:,i) = crd(cnn(:,i),1);
   yyf(:,i) = crd(cnn(:,i),2);
   ux(:,i) =  Sol.uAlpha(cnn(:,i),1,1) ;
   uxDot(:,i) = Sol.uDotAlpha(cnn(:,i),1,1) ;
   uy(:,i) =  Sol.uAlpha(cnn(:,i),2,1) ;
   uyDot(:,i) = Sol.uDotAlpha(cnn(:,i),2,1) ;
   pres(:,i) = Sol.p(cnn(:,i),1) ;
   alevelx(:,i) = Sol.aleVelAlpha(cnn(:,i),1);
   alevely(:,i) = Sol.aleVelAlpha(cnn(:,i),2);
end
        
% Form element matrix and assemble Galerkin terms
[LHS, RHS] = GalerkinTerms(Sol, fluid, pmc, solver, iif, jjf, xxf, yyf, ux, uy, ...
                           gW, N, Nx, Ny, nElem, nQuad, nen, ndof, alevelx, alevely);
                                
% Form element matrix and assemble Petrov-Galerkin terms
[LHS, RHS] = PetrovGalerkinTerms1(Sol, fluid, pmc, solver, iif, jjf, ...
                                  xxf, yyf, ux, uy, gW, N, Nx, Ny, ...
                                  nElem, nQuad, nen, ndof, elemType, ...
                                  LHS, RHS, alevelx, alevely);
[LHS, RHS] = PetrovGalerkinTerms2(Sol, fluid, pmc, solver, iif, jjf, ...
                                  xxf, yyf, ux, uy, gW, N, Nx, Ny, ...
                                  nElem, nQuad, nen, ndof, elemType, ...
                                  LHS, RHS, alevelx, alevely);
[LHS, RHS] = PetrovGalerkinTerms3(Sol, fluid, pmc, solver, iif, jjf, ...
                                  xxf, yyf, ux, uy, gW, N, Nx, Ny, ...
                                  nElem, nQuad, nen, ndof, elemType, ...
                                  LHS, RHS, alevelx, alevely);
        
                                      
% Solve the linear system
        
% Select the unknown nodal values
freeNodesU = unique([fluid.DirichletU; unique(BCCyl(:))]);
freeNodesU = setdiff(1:size(crd,1),[freeNodesU]);
freeNodesV = unique([fluid.DirichletV; unique(BCCyl(:))]);
freeNodesV = setdiff(1:size(crd,1),[freeNodesV]);
freeNodesP = setdiff(1:size(crd,1),[]) ;

freeNodes = [freeNodesU';freeNodesV' + size(crd,1); freeNodesP' + 2*size(crd,1)];
        
result = Sol.uAlpha(:,:,1);
result = result(:);
result = [result;Sol.p];
resultDot = Sol.uDotAlpha(:,:,1);
resultDot = resultDot(:);
resultDot = [resultDot; Sol.p];

Increment = LHS(freeNodes,freeNodes)\RHS(freeNodes);
        
% Update the increments
result(freeNodes) = result(freeNodes) + Increment;
resultDot(freeNodes) = resultDot(freeNodes) + (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Increment ;

Sol.uAlpha(:,:,1) = reshape(result(1:2*ndof),[],2);
Sol.uDotAlpha(:,:,1) = reshape(resultDot(1:2*ndof),[],2);
Sol.p = result((2*ndof+1):(3*ndof));

% Update the solution
Sol.u = Sol.uPrev + (1/pmc.alpha)*( Sol.uAlpha - Sol.uPrev ) ;
Sol.uDot = Sol.uDotPrev + (1/pmc.alphaM)*( Sol.uDotAlpha - Sol.uDotPrev ) ;
        
NSnormIndicator =  norm(Increment)/norm(result(freeNodes)) ;
fprintf('NS: %e, ', NSnormIndicator);
clear freeNodes1 freeNodes2 freeNodes3
clear result resultDot
        
end

