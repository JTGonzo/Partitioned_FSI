% Generalized Broyden with Least-Squares Update
classdef MVLS < handle
    properties (SetAccess=private)
        n = 0;
        r = []; 
        xt = [];
        V = [];
        W = [];
        T = [];
        Wprev = [];
        Qprev = [];
        Rprev = [];
        Dummy = [];
        file;
        count = 0;
    end
    
    properties (SetAccess=immutable)
        small;
        limit;
    end

%% Function initialization and interface-data retention    
    methods
        function M = MVLS(small, limit, count)
            M.small = small;         %| smallest floating point number
            M.limit = limit;         %| number of time steps to reuse
            
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
        function dx = predict(M,r)                       
            
            if ((isempty(M.V))&&(isempty(M.Wprev)))                                        
                    error('Model cannot predict without any data');
            end
            
            %% Current time-steps secant contributions        
            if ((~isempty(M.V))||(M.n<=4)) 
                % simple filter to avoid almost linear-dependence 
                vt = [M.V M.Dummy];
                if length(vt(1,:))>= 2
                    smaller = 1e-4;
                    [M.V, M.W, M.T, ~, ~, ~, Q, R, M.count] = filt_NM(M.V, M.W, M.T, M.Dummy, M.Dummy, M.Dummy, smaller, M.count);                   
                else                     
                    vt = [M.V M.Dummy]; 
                    [Q,R] = qr(vt,0);
                end
                b = Q'*r;
                c = R\b;               %| least-squares weights 
                dx = M.W*c;            %| approx. inverse Jacobian
                r = r - Q*b;           %| 
            else
                dx = zeros(length(r),1);
            end
            
            %% Incorporate recursive Jacobian contributions
            if M.n > 4
                if M.n-1 <= M.limit
                    kn = M.n-1;
                else
                    kn = M.limit;
                end
                
                if M.n < 20                                       
                      qq = M.Qprev{kn};
                      b = qq'*r;
                      c = M.Rprev{kn}\b;        %| past step leas-squares weights 
                      dx = dx + M.Wprev{kn}*c;  %| past time step contributions 
                else
                    i = 0;
                    while ((norm(r)>M.small) && (i<min(M.limit,length(M.Wprev))))
                        qq = M.Qprev{kn};
                        b = qq'*r;
                        c = M.Rprev{kn}\b;        %| past step leas-squares weights 
                        dx = dx + M.Wprev{kn}*c;  %| past time step contributions                
                        r = r - qq*b;             %| remove accounted residual info
                        kn = kn - 1;
                        i = i + 1;
                    end
                end
            end
        end

%% Function maintenance and monitoring        
        function clear(M)        
            format = '%4i %4i\n';           %| data format         
            data = [M.n M.count];           %| time-step, coupling iteration, residual 
            fprintf(M.file,format,data);    %| write results to file
            
            M.r  = [];      %} clear current coupling iteration matrices
            M.xt = [];
            M.V  = [];
            M.W  = [];
            M.T  = [];
            M.count = 0;
        end

        function increase_time(M)
            M.n = M.n+1;

            if M.n > 1
                [Q,R] = qr(M.V,0);                     
                %| append coupling iteration matrices 
                if M.n-1 <= M.limit
                    M.Qprev{M.n-1} = Q;
                    M.Rprev{M.n-1} = R;
                    M.Wprev{M.n-1} = M.W;               
                else
                %| remove outdated difference information
                    M.Qprev{1} = [];   
                    M.Rprev{1} = [];
                    M.Wprev{1} = [];
                    
                    for i=2:M.limit
                        M.Qprev{i-1} = M.Qprev{i};
                        M.Rprev{i-1} = M.Rprev{i};
                        M.Wprev{i-1} = M.Wprev{i};
                    end
                    
                    M.Qprev{M.limit} = Q;
                    M.Rprev{M.limit} = R;
                    M.Wprev{M.limit} = M.W;
                end
            end       
            M.clear();
        end
        
        function closefile(M)
            fclose(M.file);
        end
       
        function LS = ready(M)
            if M.n > 4
                LS =~isempty([M.V M.Wprev]); %| data must be available for Jacobian approx. 
            else
                LS =~isempty(M.V);
            end
        end
    end
end