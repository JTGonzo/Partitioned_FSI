%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%                RBF Mapping Function for interpolation                   %
%                                                                         % 
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% testing
clear all 
Xt = [1,2,0];
Xc = [2,3,0; 4,5,0;6,5,0];
r = 18; 
phi1 = rbf_Function(Xc(1,:),Xc,r)
phi2 = rbf_Function(Xc(2,:),Xc,r)
phi3 = rbf_Function(Xc(3,:),Xc,r)
phi4 = rbf_Function(Xt,Xc,r)
[A_test,C_test] = rbf_Mapping(Xt,Xc,r);
%% function writing

% Inputs:
%   crd_interface_s: the vector of solid coordinates for all the solid nodes
%               (Ns,3) [x,y,z]
%   crd_interface_f: the vector of fluid coordinates for all the fluid interface
%               nodes (Nf,3)
%   r : radius for normalization of the functions
% Outputs:
%   A: sparse mapping function of solid 
function [A,C] = rbf_Mapping(crd_interface_s,crd_interface_f,r)
    % compute size of C and A matrix. 
    % assuming using linear polynomial add on: l + l1 x + l2 y + l3 z
    [n_f,d_f] = size(crd_interface_f);
    [n_s,d_s] = size(crd_interface_s);
    n_C = n_f +4; 
    n_A = n_s;
    d_A = n_f + 4; 
    
    % initialize the empty matrices. 
    A_nonSparse = zeros(n_A,d_A);
    C_nonSparse = zeros(n_C,n_C);
    
    % lets construct C first
    C_nonSparse(5:end,1) = 1; 
    C_nonSparse(5:end, 2) = crd_interface_f(:,1);
    C_nonSparse(5:end, 3) = crd_interface_f(:,2);
    C_nonSparse(5:end, 4) = crd_interface_f(:,3);
    
    C_nonSparse(1,5:end) = 1; 
    C_nonSparse(2,5:end) = transpose(crd_interface_f(:,1));
    C_nonSparse(3,5:end) = transpose(crd_interface_f(:,2));
    C_nonSparse(4,5:end) = transpose(crd_interface_f(:,3));
    
    for i = 5: n_C
        X_target = crd_interface_f(i-4,:);
        X_reference = crd_interface_f; 
        phi = rbf_Function(X_target,X_reference,r);
        C_nonSparse(i,5:end) = transpose(phi);
    end 
    
    % now putting together the A matrix
    A_nonSparse(:,1) = 1; 
    A_nonSparse(:,2) = crd_interface_s(:,1);
    A_nonSparse(:,3) = crd_interface_s(:,2);
    A_nonSparse(:,4) = crd_interface_s(:,3);
    
    for i = 1:n_s
        X_target = crd_interface_s(i,:);
        X_reference = crd_interface_f; 
        phi = rbf_Function(X_target,X_reference,r);
        A_nonSparse(i,5:end) = transpose(phi);
    end 
    
    C = C_nonSparse; 
    A = A_nonSparse;



end 

% Inputs: 
%   Xt: coordinates of the the target point (1,3)
%   Xc : coordinates of the reference point (n,3)
%   r:  radius 
% Outputs:
%   phi: the value of the RBF function (n,1)
function [phi] = rbf_Function(Xt,Xc,r)
    [n,d] = size(Xc);
    
    %compute component wise distance
    distance_comp = repmat(Xt,n,1)- Xc;
    distance_norm_scaled = sum(distance_comp .^2,2) ./r;
    heavside_vector = heaviside(1-distance_norm_scaled);
    heavside_vector(abs(heavside_vector - 0.5) < 1E-10) = 1; % correct equals to  1
    
    % function computation.
    phi_1 = ((1-distance_norm_scaled).* heavside_vector).^4;
    phi_2 = (1+ 4.*distance_norm_scaled) .* heavside_vector;
    phi = phi_1 .* phi_2;
    

end 

