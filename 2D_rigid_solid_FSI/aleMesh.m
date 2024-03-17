%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%                      2D ALE Mesh Equation Solver                        %
%                                                                         % 
%              \nabla . \sigma^m = 0,    in \Omega^s,                     %                      
%              \sigma = grad \eta + grad \eta^T + (div \eta)I,            %
%                                                                         %
%    where \sigma is the stress experienced by the ALE mesh due to        %
%    strain induced by structural movement, \eta is the mesh disp-        %
%    lacement for the fluid nodes and I is the identity tensor.           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sol] = aleMesh(Sol, solver, solid, BCCyl, BCTop, BCBottom,...
                          BCLeft, BCRight, pmc, cnn, crd, elemType, nen, ...
                          ndof, nElem)

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
 
% Initialize the displacement
disp = zeros(ndof,2);

% Form the local to global map
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
bcLeft = unique(BCLeft(:));
bcRight = unique(BCRight(:));
bcTop = unique(BCTop(:));
bcBottom = unique(BCBottom(:));
BCBound = [bcLeft; bcRight; bcTop; bcBottom];
BCBound = unique(BCBound(:));
disp(unique(BCCyl(:)),1) = Sol.dispS(1) ;
disp(unique(BCCyl(:)),2) = Sol.dispS(2) ;
disp(BCBound,1) = zeros(size(BCBound,1),1);
disp(BCBound,2) = zeros(size(BCBound,1),1);
        
% ALE mesh equation
xxf = zeros(size(cnn));
yyf = zeros(size(cnn));
aledispx = zeros(size(cnn));
aledispy = zeros(size(cnn));

% Localize the data to each element
for i=1:nen
   xxf(:,i) = crd(cnn(:,i),1);
   yyf(:,i) = crd(cnn(:,i),2);
   aledispx(:,i) =  disp(cnn(:,i),1) ;
   aledispy(:,i) =  disp(cnn(:,i),2) ;
end

% Form element matrix and assemble Galerkin terms
sA1 = zeros(nen^2*nElem,nQuad); 
sA2 = zeros(nen^2*nElem,nQuad);
sA3 = zeros(nen^2*nElem,nQuad);

for p = 1:nQuad  
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
         yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
    if size(J,2)==1
        J = J';
    end
    volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
           
    volume = abs(volume);

    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
        
    index = 0;
    for i = 1:nen
        for j = 1:nen      
            % Galerkin terms for ALE equation
            Aij_1 = gW(p)*(DNDx(:,i).*DNDx(:,j));
            Aij_2 = gW(p)*(DNDy(:,i).*DNDy(:,j));
            Aij_3 = gW(p)*(DNDx(:,i).*DNDy(:,j));
            Aij_1 = Aij_1.*volume;
            Aij_2 = Aij_2.*volume;
            Aij_3 = Aij_3.*volume;
            sA1(index+1:index+nElem,p) = Aij_1;
            sA2(index+1:index+nElem,p) = Aij_2;
            sA3(index+1:index+nElem,p) = Aij_3;

            index = index + nElem;
        end
    end
end
% Summation of all quadrature data
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
 
% Assemble the matrix      
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A3 = sparse(iif,jjf,sA3,ndof,ndof);

% Left-hand side matrix
LHS = [2*A1+A2 A3'; A3 A1+2*A2] + [A1 A3; A3' A2];

% Right-hand side vector
RHS = - (LHS)* disp(:) ;

% Select the unknown nodal values
freeNodes = unique([BCBound; unique(BCCyl(:))]');
freeNodes = setdiff(1:size(crd,1),[freeNodes]);

freeNodes = [freeNodes';freeNodes' + size(crd,1)];
result = disp(:);

% Solve the linear system
result(freeNodes) = LHS(freeNodes,freeNodes)\RHS(freeNodes);
disp = reshape(result(1:2*ndof),[],2);

% Update the ALE displacement and velocity
Sol.aleDisp = disp ;
Sol.aleVel = (Sol.aleDisp - Sol.aleDispPrev)./solver.dt ;

aa = 1;
end