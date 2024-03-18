%GENERALIZEDALPHA_TIMEADVANCE class to Handle the time advancing scheme based 
%on Generalized Alpha scheme

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

classdef GeneralizedAlpha_TimeAdvance < handle
    %% Set global variables
    properties (GetAccess = public, SetAccess = protected)
        M_U;
        M_dU;
        M_d2U;
        M_rhs;
        M_Csi;
        
        M_beta;
        M_gamma;
        M_alpha_f;
        M_alpha_m;
        M_timestep;
        M_stateSize;
    end
    
    methods (Access = public)
        %% Initial solution variables
        function obj = GeneralizedAlpha_TimeAdvance( beta, gamma, alpha_m, alpha_f, timestep )
            obj.M_beta      = beta;
            obj.M_gamma     = gamma;
            obj.M_alpha_m   = alpha_m;
            obj.M_alpha_f   = alpha_f;
            obj.M_timestep  = timestep;
        end
        
        %% Initialize solution space
        function obj = Initialize( obj, U, dU, d2U )
            if size(U,2) > 1
                error('Newmark_TimeAdvance.Initialize: incorrectSize of U')
            end
            obj.M_stateSize = size(U, 1);
            obj.M_U   = U;
            
            if size(dU,1) ~= obj.M_stateSize || size(dU,2) > 1
                error('Newmark_TimeAdvance.Initialize: incorrectSize of dU')
            end
            obj.M_dU  = dU;
            
            if size(d2U,1) ~= obj.M_stateSize || size(d2U,2) > 1
                error('Newmark_TimeAdvance.Initialize: incorrectSize of d2U')
            end
            
            obj.M_d2U = d2U;
            
            obj.M_Csi = 1 / (obj.M_beta * obj.M_timestep^2) * (obj.M_U + obj.M_timestep * obj.M_dU) ...
                + (1 - 2*obj.M_beta)/(2*obj.M_beta) * obj.M_d2U;
            
            obj.M_rhs = (1 - obj.M_alpha_m) / (obj.M_beta * obj.M_timestep^2) * (obj.M_U + obj.M_timestep * obj.M_dU) ...
                + (1 - obj.M_alpha_m - 2*obj.M_beta)/(2*obj.M_beta) * obj.M_d2U;
        end
        
        %% Update solid variables based on current step's displacement 
        function obj = Update( obj, U_np1 )
            
            d2U_np1        = 1 / ( obj.M_beta * obj.M_timestep^2) * U_np1 - obj.M_Csi;
            obj.M_dU       = obj.M_dU + obj.M_timestep * (obj.M_gamma * d2U_np1 + (1 - obj.M_gamma) * obj.M_d2U) ;
            obj.M_d2U      = d2U_np1;
            obj.M_U        = U_np1;
            
            obj.M_Csi = 1 / (obj.M_beta * obj.M_timestep^2) * (obj.M_U + obj.M_timestep * obj.M_dU) ...
                + (1 - 2*obj.M_beta)/(2*obj.M_beta) * obj.M_d2U;
            
            obj.M_rhs = (1 - obj.M_alpha_m) / (obj.M_beta * obj.M_timestep^2) * (obj.M_U + obj.M_timestep * obj.M_dU) ...
                + (1 - obj.M_alpha_m - 2*obj.M_beta)/(2*obj.M_beta) * obj.M_d2U;        
        end
        
        %% Compute current solid velocity/acceleration per domain displacement
        function dU = velocity( obj, uS_t )            
            d2U = 1 / ( obj.M_beta * obj.M_timestep^2) * uS_t - obj.M_Csi;
            dU  = obj.M_dU + obj.M_timestep * (obj.M_gamma * d2U + (1 - obj.M_gamma) * obj.M_d2U) ;
        end
        
        function dU = velocityNL1( obj ) % no current mesh motion             
            dU  = obj.M_dU;
        end
        
        function [ Us, dUs, d2Us, Csi, rhs ] = retrive( obj ) % no current mesh motion             
            Us = obj.M_U;
            dUs = obj.M_dU;
            d2Us = obj.M_d2U;
            Csi = obj.M_Csi;
            rhs = obj.M_rhs;            
        end
        
        function obj = restart( obj, Us, dUs, d2Us, Csi, rhs) % no current mesh motion             
            obj.M_U = Us;
            obj.M_dU = dUs;
            obj.M_d2U = d2Us;
            obj.M_Csi = Csi;
            obj.M_rhs = rhs;            
        end
        
        %% RhsContribute
        function Rhs = RhsContribute( obj )            
            Rhs = obj.M_rhs;            
        end
        
        %% Mass coefficient of time integration scheme
        function c = MassCoefficient( obj )           
            c = (1 - obj.M_alpha_m) / ( obj.M_beta * obj.M_timestep^2) ;           
        end        
    end
end