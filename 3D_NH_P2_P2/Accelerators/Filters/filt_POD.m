function [wt, Q, R, count] = filt_POD(V, W, Vprev, Wprev, small, count)
    
%% filter interface data being collected   
    vt = [V Vprev];      %| current and previous residual difference data
    wt = [W Wprev];      %| current and previous solution difference data
    
    if (~isempty(vt))           
         eta = length(vt(1,:));
         sig = (1/eta)*(vt'*vt);     %| form auto-correlation matrix
         [X,lam] = eig(sig);         %| evaluates its eigenvalues & vectors   

         lam_n = diag(sort(diag(lam),'descend'));  %| reorder eigenvalue matrix
         [p, ind] = sort(diag(lam),'descend');     %| store the indices of eigenvalue reordering
         X_n = X(:,ind);                           %| arrange the eigenvectors
         
         % find index where eigenvlaues become tiny
         for i = 1:eta
               if lam_n(i,i)< small
                    c = i;
                    count = count+(eta-c);
                    break
               end

               if i == eta
                    c = eta;
               end
         end
        
         % project the difference matricies onto the eigenvector space
         V_pod = vt*X_n;  %re-verify this multiplication
         W_pod = wt*X_n;
        
         % only keep the most dominant modes
         vt = V_pod(:,1:c);
         wt = W_pod(:,1:c);
    end

%% Output filtered QR-decomposition 
[Q,R] = qr(vt,0); 
end