function [Length, Force, Trac] = IntegratedOutput_old(Sol, crd, fluid, cnn, nen, ndof, MESH)

temp =  MESH.Fluid.boundaries(1:2,MESH.Fluid.surf_elem);
Belem = temp';
nElem = size(Belem,1);
ndof = size(unique(Belem),1) ;
BelemU = unique(Belem);

% % Get elements corresponding to first layer of boundary layer
[tf,idx1] = ismember(BelemU(:,1),cnn(:,1));
[tf,idx1(:,2)] = ismember(BelemU(:,1),cnn(:,2));
[tf,idx1(:,3)] = ismember(BelemU(:,1),cnn(:,3));
Srfelem = unique(idx1(:));
Srfelem(Srfelem==0)=[];

cnnSrf = cnn(Srfelem,:);
idx2 = zeros(size(cnnSrf,1),2);
idx3 = zeros(size(Belem,1),1);

for i = 1:size(Belem,1) 
    for j = 1:size(cnnSrf,1)
        [tf, id_temp] = ismember(Belem(i,:),cnnSrf(j,:));        
        if tf(1,:) == 1
            idx2(j,:) = id_temp;
            idx3(i) = j;
        end
    end
end
cnnFSI = cnnSrf(idx3,:);
idx4 = idx2(idx3,:);

% Reorder the element points for reduced integration
for i = 1:size(cnnFSI,1)
    if idx4(i,1) ~= 1
         cnnFSInew(i,1) = cnnFSI(i,idx4(i,1));
         cnnFSInew(i,2) = cnnFSI(i,idx4(i,2));
         cnnFSInew(i,3) = cnnFSI(i,(6-(idx4(i,1)+idx4(i,2))));       
    else
        cnnFSInew(i,:) = cnnFSI(i,:);
    end
end
 
% Shape functions, gauss points and weights for Numerical integration
% Gauss points
gP = [-1/sqrt(3), -1.0
       1/sqrt(3), -1.0] ;
   
% Gauss weights
gW = [1.0, 1.0];

% Shape functions
N(:,1) = 0.5.*(2.-gP(:,1)-gP(:,2)) ;
N(:,2) = 0.5.*(gP(:,1)) ;
N(:,3) = 0.5.*(gP(:,2)) ;
     
% Derivative of shape functions
Nx(:,1) = -0.5.*ones(3,1) ;
Nx(:,2) =  0.5.*ones(3,1) ;
Nx(:,3) =  zeros(3,1) ; 

Ny(:,1) = -0.5.*ones(3,1) ;
Ny(:,2) =  zeros(3,1) ;
Ny(:,3) =  0.5.*ones(3,1) ;   

Nx = Nx' ;
Ny = Ny' ;

% Number of quadrature points
nQuad = size(gW,2);

xxf = zeros(size(cnnFSInew));
yyf = zeros(size(cnnFSInew));
ux = zeros(size(cnnFSInew));
uy = zeros(size(cnnFSInew));
% pres = zeros(size(cnnFSInew));

% Localize the data to each element
for i=1:nen
   xxf(:,i) = crd(cnnFSInew(:,i),1);
   yyf(:,i) = crd(cnnFSInew(:,i),2);
   ux(:,i) =  Sol.u(cnnFSInew(:,i),1,1) ;
   uy(:,i) =  Sol.u(cnnFSInew(:,i),2,1) ;
   pres(:,i) = Sol.p(cnnFSInew(:,i),1) ;
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

Trac(:,1) = A1+A3; 
Trac(:,2) = A2+A4;

Force = [sum(A1,1)+sum(A3,1); sum(A2,1)+sum(A4,1)];

end