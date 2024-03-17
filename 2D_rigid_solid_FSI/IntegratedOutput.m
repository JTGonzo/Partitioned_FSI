%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Calculate the integrated output values of the forces on the structure  %
%                                                                         %
%  This function evaluated the integral of the stress term on the surfa-  %
%  ce of the structure, i.e.,                                             %
%                                                                         %
%                   \int_\Gamma (\sigma . n) d\Gamma,                     %
%                                                                         %
%  where \sigma is the fluid stress, \sigma = -pI + \mu(grad U + grad U^T)%
%  p being the pressure, I is the Identity matrix, \mu is the fluid visc- %
%  osity and grad U is the gradient of velocity. n is the normal to the   %
%  structural surface.                                                    %
%                                                                         %
%  The integral is ealuated by reduced quadrature integration on the two- %
%  dimensional element just beside the surface.                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Length, Force] = IntegratedOutput(Sol, crd, BCCyl, fluid, cnn)

% Number of element nodes in 2D (4 for fourNodeQuad)
nen = 4 ;

% Swap BCCyl
temp = BCCyl(:,1);
BCCyl(:,1) = BCCyl(:,2);
BCCyl(:,2) = temp;

nElem = size(BCCyl,1);
ndof = size(unique(BCCyl),1) ;

% Get elements corresponding to first layer of boundary layer
[tf,idx2] = ismember(BCCyl(:,1),cnn(:,1));
[tf,idx2(:,2)] = ismember(BCCyl(:,1),cnn(:,2));
[tf,idx2(:,3)] = ismember(BCCyl(:,1),cnn(:,3));
[tf,idx2(:,4)] = ismember(BCCyl(:,1),cnn(:,4));
elemCyl = unique(idx2(:));
elemCyl(elemCyl==0)=[];

cnnCyl = cnn(elemCyl,:);

% Reorder the element points for reduced integration 
for i=1:size(elemCyl,1)
    [tf, idx1(i)] = ismember(BCCyl(i,1),cnnCyl(i,:));
    cnnCylNew(i,1) = cnnCyl(i,idx1(i));
    if (idx1(i) == 3)
        cnnCylNew(i,2) = cnnCyl(i,idx1(i)+1);
        cnnCylNew(i,3) = cnnCyl(i,1);
        cnnCylNew(i,4) = cnnCyl(i,2);
    elseif (idx1(i) == 4)
        cnnCylNew(i,2) = cnnCyl(i,1);
        cnnCylNew(i,3) = cnnCyl(i,2);
        cnnCylNew(i,4) = cnnCyl(i,3);
    end
end

% Shape functions, gauss points and weights for Numerical integration
% Gauss points
gP = [-1/sqrt(3), -1.0
       1/sqrt(3), -1.0] ;
% Gauss weights
gW = [1.0, 1.0];
% Shape functions
N(:,1) = 0.25.*(1-gP(:,1)).*(1-gP(:,2)) ;
N(:,2) = 0.25.*(1+gP(:,1)).*(1-gP(:,2)) ;
N(:,3) = 0.25.*(1+gP(:,1)).*(1+gP(:,2)) ;
N(:,4) = 0.25.*(1-gP(:,1)).*(1+gP(:,2)) ;
     
% Derivative of shape functions
Nx(:,1) = -0.25.*(1-gP(:,2)) ;
Nx(:,2) =  0.25.*(1-gP(:,2)) ;
Nx(:,3) =  0.25.*(1+gP(:,2)) ;
Nx(:,4) = -0.25.*(1+gP(:,2)) ;
Ny(:,1) = -0.25.*(1-gP(:,1)) ;
Ny(:,2) = -0.25.*(1+gP(:,1)) ;
Ny(:,3) =  0.25.*(1+gP(:,1)) ;
Ny(:,4) =  0.25.*(1-gP(:,1)) ;
Nx = Nx' ;
Ny = Ny' ;
% Number of quadrature points
nQuad = size(gW,2);

xxf = zeros(size(cnnCylNew));
yyf = zeros(size(cnnCylNew));
ux = zeros(size(cnnCylNew));
uy = zeros(size(cnnCylNew));

% Localize the data to each element
for i=1:nen
   xxf(:,i) = crd(cnnCylNew(:,i),1);
   yyf(:,i) = crd(cnnCylNew(:,i),2);
   ux(:,i) =  Sol.u(cnnCylNew(:,i),1,1) ;
   uy(:,i) =  Sol.u(cnnCylNew(:,i),2,1) ;
   pres(:,i) = Sol.p(cnnCylNew(:,i),1) ;
end

for p = 1:nQuad  
    % Jacobian for each element
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
         yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
     
    if size(J,2)==1
        J = J';
    end
    
    volume = sqrt(J(:,1).^2 + J(:,3).^2);
    vol = ( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
    vol = abs(vol);
    
    % Normal evaluation
    normal(:,1) = -J(:,3)./volume ;
    normal(:,2) = J(:,1)./volume ;
    
    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(vol,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(vol,1,nen);
    
    % Pressure and gradients evaluation for velocity
    locgradUx = sum(ux.*DNDx,2) ;
    locgradUy = sum(ux.*DNDy,2) ; 
    locgradVx = sum(uy.*DNDx,2) ;
    locgradVy = sum(uy.*DNDy,2) ;
    locP  = sum(repmat(N(p,:),nElem,1).*pres,2);
            
    % Length/ Area of the line/ surface integral in 1D/2D
    A0(:,p) = gW(p).*volume ;
    
    % X- traction term (pressure)
    A1(:,p) = gW(p).*(-locP).*(normal(:,1));
    A1(:,p) = A1(:,p).*volume ;
    
    % Y- traction term (pressure)
    A2(:,p) = gW(p).*(-locP).*(normal(:,2));
    A2(:,p) = A2(:,p).*volume ;
    
    % X- traction term (viscous)
    A3(:,p) = gW(p).*fluid.visc.*(2.*normal(:,1).*locgradUx ...
            + normal(:,2).*(locgradVx + locgradUy)) ;
    A3(:,p) = A3(:,p).*volume ;
    
    % Y- traction term (viscous)
    A4(:,p) = gW(p).*fluid.visc.*(2.*normal(:,2).*locgradVy ...
            + normal(:,1).*(locgradVx + locgradUy)) ;
    A4(:,p) = A4(:,p).*volume ;
    
end
% Summation of all quadrature data
A0 = sum(A0,2);
A1 = sum(A1,2);
A2 = sum(A2,2);
A3 = sum(A3,2);
A4 = sum(A4,2);

Length = [sum(A0,1)];
Force = [sum(A1,1)+sum(A3,1); sum(A2,1)+sum(A4,1)];

end