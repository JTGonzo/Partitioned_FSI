% Anderson Acceleration
classdef AA < handle
%% Set in-function fixed (internal) and dynamic (interface) variable values 
    properties (SetAccess=private)
        r = []; 
        xt = [];
        Vprev = [];
        Wprev = []; 
        Tprev = []; 
        V = [];
        W = [];        
        T  = [];
        n  = 0;
        count = 0;
        finished_primal = false;
        file;
    end
    
    properties (SetAccess = immutable)
        small;
        reuse;
        filter;
    end
    
%% Function initialization and interface-data retention
    methods
        function M = AA(small, reuse, filter, count)
            M.small = small;              %| smallest floating point number
            M.reuse = reuse;              %| number of time steps to reuse
            M.filter = filter;            %| what type of filter to use
            
            M.file = fopen('Results/Filtering.txt','w');
            M.count = count;
        end
        
        function add(M,r,xt)
            if (~isempty(M.r))
                M.V = [r - M.r M.V];     %| Update residual difference matrix
                M.W = [xt - M.xt M.W];   %| Update solution difference matrix
                M.T = [M.n M.T];         %| Update info time identifier    
            end            
            M.r = r;                     %| past residual update
            M.xt = xt;                   %| past intermediary  solution               
        end

%% Secant Jacobian computation
        function [J] = predict(M, r)
            %% Filter the retained information
            if M.filter == 1             %| QR1 -Filter
                [M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, Q, R, M.count] = filt_QR1(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small,M.count);
            elseif M.filter == 2         %| QR2 -Filter
                vt = [M.V M.Vprev];
                if length(vt(1,:))>= 2
                    [M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, Q, R, M.count] = filt_NM(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small,M.count);                   
                else                     
                    vt = [M.V M.Vprev]; 
                    [Q,R] = qr(vt,0);
                end
            else                         %| POD -Filter
               [W, Q, R, M.count] = filt_POD(M.V, M.W, M.Vprev, M.Wprev, M.small, M.count);
            end
            
            %% Solve the least-squares porblem and compute approx. Jacobian
            if ((~isempty(M.V))||(~isempty(M.Vprev)))          
                c = R\(Q'*r);           % least-squares weights
                if M.filter == 3
                    J = W*c;
                else
                    wt = [M.W M.Wprev];
                    J = wt*c;
                end
            else
                if (isempty(M.xt))  %| Min. data from 2 preceeding iterations needed                                          
                    error('Model cannot predict without any data');                  
                else
                    J = zeros(size(M.xt));          
                end
            end
        end
        
%% Function maintenance and monitoring
         function clear(M)
            if M.reuse == 0                 %| no past information is kept
                 M.Vprev = [ ];
                 M.Wprev = [ ];
                 M.Tprev = [ ];
            else 
                 M.Vprev = [M.V M.Vprev];   %| append coupling iteration matrices
                 M.Wprev = [M.W M.Wprev];
                 M.Tprev = [M.T M.Tprev];
                 
                 old = find(M.Tprev < M.n - M.reuse);   %| find outdated diff. information indicies
                 
                 if (~isempty(old))                     %| remove outdated difference information   
                    M.Vprev(:,old) = [];      
                    M.Wprev(:,old) = [];
                    M.Tprev(:,old) = [];
                 end
            end
            
            format = '%4i %4i\n';           %| data format         
            data = [M.n M.count];           %| time-step, coupling iteration, residual 
            fprintf(M.file,format,data);    %| write results to file
                    
            % clear coupling iteration matrices            
            M.r  = [];       
            M.xt = [];
            M.V  = [];
            M.W  = [];
            M.T  = [];
            M.count = 0;
        end
        
        function increase_time(M)             
            M.n = M.n + 1;                  %| increase time step
            M.clear();                      %| update retained info
        end
        
        function closefile(M)
            fclose(M.file);
        end
        
        function LS = ready(M)             
            LS =~ isempty([M.V M.Vprev]);   %| data must be available for Jacobian approx. 
        end
    end
end