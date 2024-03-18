function [Length, tracSrc] = IntegratedOutput2D_try(Sol, fluid, BCTop, BCBottom, BCLeft, BCRight, crd, cnn, nElem, ndof, nen)

% Get elements corresponding to boundaries
cnnBCtop = getElemBC(BCTop, cnn);
cnnBCbottom = getElemBC(BCBottom, cnn);
cnnBCleft = getElemBC(BCLeft, cnn);
cnnBCright = getElemBC(BCRight, cnn);

% Reorder the element points for reduced integration
[cnnBCrightNew] = reorderElemBC(BCRight, cnnBCright);
[cnnBCtopNew] = reorderElemBC(BCTop, cnnBCtop);
[cnnBCbottomNew] = reorderElemBC(BCBottom, cnnBCbottom);
[cnnBCleftNew] = reorderElemBC(BCLeft, cnnBCleft);

cnnBCnew = [cnnBCtopNew; cnnBCbottomNew; cnnBCleftNew; cnnBCrightNew];
%cnnBCnew = [BCTop; BCBottom; BCLeft; BCRight];

cnnBCnew = [cnnBCnew, ones(size(cnnBCnew, 1), 2)];

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

nElemBC = size(cnnBCnew, 1);
% Form the local to global map
iif = zeros(nen^2*nElemBC,1); 
jjf = zeros(nen^2*nElemBC,1);
index = 0;
for i = 1:nen
   for j = 1:nen
      iif(index+1:index+nElemBC) = double(cnnBCnew(:,i)); 
      jjf(index+1:index+nElemBC) = double(cnnBCnew(:,j));  
      index = index + nElemBC;
   end
end

xxf = zeros(size(cnnBCnew));
yyf = zeros(size(cnnBCnew));
ux = zeros(size(cnnBCnew));
uy = zeros(size(cnnBCnew));

% Localize the data to each element
for i=1:nen
	xxf(:,i) = crd(cnnBCnew(:,i),1);
	yyf(:,i) = crd(cnnBCnew(:,i),2);
	ux(:,i) =  Sol.u(cnnBCnew(:,i),1,1) ;
	uy(:,i) =  Sol.u(cnnBCnew(:,i),2,1) ;
	% pres(:,i) = Sol.p(cnnBCnew(:,i),1) ;
	pres(:,i) = fluid.pInf ;
% 	dens(:,i) = fluid.dens(cnnBCnew(:,i),1) ;
% 	visc(:,i) = fluid.visc(cnnBCnew(:,i),1) ;
end

sA1 = zeros(nen^2*nElemBC,nQuad); 
sA2 = zeros(nen^2*nElemBC,nQuad); 

for p = 1:nQuad  
    % Jacobian for each element
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
		yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
	if size(J,2)==1
		J = J';
    end
    	
	volume = sqrt(J(:,1).^2 + J(:,3).^2); % normalizes the normal vector
	vol =( J(:,1).*J(:,4) -J(:,2).*J(:,3) ); % determinant
    vol = abs(vol);
    
	% Normal evaluation
	normal(:,1) = -J(:,3)./volume ;
	normal(:,2) = J(:,1)./volume ;		
	
    negJacobian = find(vol<0);
    if ~isempty(negJacobian)
       disp('Mesh deformed, Negative Jacobian');
       exit
    end	

	DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
	DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen); 
    
    % Pressure and gradients evaluation for velocity
    locP  = sum(repmat(N(p,:),size(cnnBCnew, 1),1).*pres,2);
% 	localVisc = sum(repmat(N(p,:),nElemBC,1).*visc,2);
    
	locgradUx = sum(ux.*DNDx,2) ;
    locgradUy = sum(ux.*DNDy,2) ; 
    locgradVx = sum(uy.*DNDx,2) ;
    locgradVy = sum(uy.*DNDy,2) ;
	
	sigma11 = -locP + 2*localVisc.*locgradUx;
    sigma12 = localVisc.*(locgradUy + locgradVx);
    sigma21 = localVisc.*(locgradVx + locgradUy);
    sigma22 = -locP + 2*localVisc.*locgradVy;
	            
    % Length/ Area of the line/ surface integral in 1D/2D
     A0(:,p) = gW(p).*volume ;
    
    index = 0;
    for i = 1:nen
        for j = 1:nen 			
			% Galerkin source term (Stress term)			
				Aij_1 = gW(p)*N(p,i)*normal(:, 1).*(sigma11 + sigma21)*N(p,j);
				Aij_2 = gW(p)*N(p,i)*normal(:, 2).*(sigma12 + sigma22)*N(p,j);
				Aij_1 = Aij_1.*volume ;
				Aij_2 = Aij_2.*volume ;
				sA1(index+1:index+nElemBC,p) = Aij_1;
				sA2(index+1:index+nElemBC,p) = Aij_2;
            
            index = index + nElemBC;			
        end   
    end
    
end
% Summation of all quadrature data
A0 = sum(A0,2);

sA1 = sum(sA1,2);
sA2 = sum(sA2,2);

A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);

tracSrc1 = [A1]*[ones(1*ndof,1)];
tracSrc2 = [A2]*[ones(1*ndof,1)];
tracSrc = [tracSrc1; tracSrc2];

tracSrc = [tracSrc; zeros(1*ndof,1)];

Length = [sum(A0,1)];

end
