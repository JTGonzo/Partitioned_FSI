 % Rank-k Generalized Broyden Update
classdef GB < handle
    properties (SetAccess=private)
        r = []; 
        xt = [];
        Jhat = [];
        Jhat_prev = []; 
        Vprev = [];
        Wprev = []; 
        Tprev = []; 
        V = [];
        W = [];        
        T  = [];
        n  = 0;
        file;
        count = 0;
    end
    
    properties (SetAccess=immutable)
        small;
        reuse;
    end
    
%% Function initialization and interface-data retention
    methods
        function M = GB(small,reuse,count)
            M.small = small;      % Try range of regularization schemes
            M.reuse = reuse;      % 0 for gen. broyden
            
            M.file = fopen('Results/Filtering.txt','w');
            M.count = count;
        end
        
        function add(M,r,xt)
            if (~isempty(M.r))
                M.V = [r - M.r M.V];     %| Update residual difference matrix
                M.W = [xt - M.xt M.W];   %| Update solution difference matrix
                M.T = [M.n M.T];         %| Update info time identifier 
                
                % basic info filtering
                [M.V, M.W, M.T, ~, ~, ~, ~, ~,M.count] = filt_QR1(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small,M.count);
                                 
                if M.n >= 1
                    den = M.V'*M.V;
                    Z = LUinverse(den)*M.V';  
                    Jup = (M.W - M.Jhat_prev*M.V)*Z;
                    M.Jhat = M.Jhat_prev + Jup;     %| Update approx. interface Jacobian
                end              
            end            
            M.r = r;                     %| past residual update
            M.xt = xt;                   %| past intermediary  solutio             
        end
            
%% Secant Jacobian computation        
        function dx = predict(M,r) 
           
            if ((~isempty(M.V))|| (~isempty(M.Jhat)))
                if M.n == 1
                    [Q,R] = qr(M.V,0);
                    c = R\(Q'*r);
                    dx = M.W*c;
                else    
                    dx = M.Jhat*r;        %| Quasi-newton increment
                end
            else
                if (isempty(M.xt))
                    error('Model cannot predict without any data');
                else
                    dx = zeros(size(M.xt));
                end
            end
        end

%% Function maintenance and monitoring
         function clear(M)
            if M.n == 2 %4
                temp = M.V'*M.V; 
                Z = LUinverse(temp)*M.V';            
                M.Jhat = M.W*Z;           
            end
            
            M.Jhat_prev = M.Jhat;        %| save past time-step approx. Jacobian
            
            if M.reuse == 0              %| no past information is kept
                 M.Vprev = [ ];
                 M.Wprev = [ ];
                 M.Tprev = [ ];
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
            M.n = M.n + 1;      %| increase time step
            M.clear();          %| update retained info
        end
        
        function closefile(M)
            fclose(M.file);
        end
        
        function LS = ready(M)
            LS =~isempty(M.Jhat_prev);   %| data must be available for Jacobian approx.
        end
    end
end