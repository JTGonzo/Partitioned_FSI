classdef Extrapolator < handle
%% Set in-function fixed (internal) and dynamic (interface) variable values 
    properties (SetAccess=private)
        xp;
        x;
        xn;
        x2n;
        
        n  = 0;
        n0 = 0;
        n1 = 1;
        n2 = 2;
    end

%% New time step variable prediction
    methods
        function E = Extrapolator(x)       

            E.initialize(x);                        %| initialize the interface matrices  
            
        end 
        
        function [xe] = predict(E)
            
            switch E.n                
                case E.n0                                 %| No extrapolation for 0 time step 
                    error('Extrapolator cannot predict initial value');                    
                case E.n1
                    xe = E.x;                             %| first time step prediction  (linear)                
                otherwise    
%                case E.n2
                    xe = 2.0*E.x - E.xn;                  %| second time step prediction (linear)                    
%                 otherwise
%                     %xe = 3.0*E.x - 3.0*E.xn + 1.0*E.x2n;  %| all subsequent time step predictions (quadratic)
%                     xe = 2.5*E.x - 2.0*E.xn + 0.5*E.x2n;
            end
            
        end
        
        function add(E,x)
            
            E.xp = x;      % adding new time step information to the predictor array
        
        end
 
%% Dynamically maintaining the data used for extrapolation
        function initialize(E,x)    % initializing the (dispacement) vector size
            E.xp  = x;
            E.x   = x;
            E.xn  = x;
            E.x2n = x;
        end
        
        function shift(E)    % shifting what we define as the previous 3 time-step's data values            
            E.x2n = E.xn;         
            E.xn  = E.x;
            E.x   = E.xp;
        end
        
%% Updating time increment of handle         
        function increase_time(E)            
            E.n = E.n + 1;
            E.shift();            
        end
    end
end