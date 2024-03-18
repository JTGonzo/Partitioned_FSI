%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%                   2D Rigid Body Structural Solver                       %
%                                                                         % 
%             m.du/dt + c.u + k.(phi-X_0) = F in \Omega^s,                %                      
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   The equation is discretized by Generalized-alpha time integration     %
%   scheme and considers just 2-degrees of freedom translation of the     %
%   body. Here, u is the velocity, (phi-X_0) is the displacement, phi     %
%   is the position mapping and X_0 is the initial coordinates, F is      %
%   the force consisting of external and internal body forces.            %
%                                                                         %
%   Generalized-alpha scheme:                                             % 
%                                                                         %
%   u_alpha = alpha.u_(n+1) + (1-alpha).u_n                               %
%   dudt_alpham = alpham.dudt_(n+1) + (1-alpham).dudt_n                   %
%   phi_(n+1) = phi_n + dt.u_n + dt^2.((0.5-beta)dudt_n + beta.dudt_(n+1))%
%   u_(n+1) = u_n + dt.((1-gamma).dudt_n + gamma.dudt_(n+1))              %
%                                                                         %
%   Putting alpha = alpham = gamma = 0.5 and beta = 0.25, one can write   %
%   the whole equation in terms of u_(n+1) and the system can be formu-   %
%   lated as shown below.                                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sol] = rigidBody(Sol, solver, solid, Force)

% Sol.uSPrev is u_n and Sol.u is u_(n+1)
% Sol.dispS is the displacement at n and Sol.dispSPrev is at (n+1)

m = solid.mass ;
c = solid.damping ;
k = solid.stiffness ;
g = solid.gravity ;
dt = solver.dt ;

Force = Force.*solid.dom ;

% Right-hand side of the equation
RHS = Force + m.*g ...
    + (1/dt).*[m 0;0 m]*Sol.uSPrev ...
    - (1/2).*[c(1) c(3); c(3) c(2)]*Sol.uSPrev ...
    - [k(1) k(3); k(3) k(2)]*Sol.dispSPrev ...
    - (dt/4).*[k(1) k(3); k(3) k(2)]*Sol.uSPrev ;

% Left-hand side of the equation
LHS = (1/dt).*[m 0;0 m] ...
    + (1/2).*[c(1) c(3); c(3) c(2)] ...
    + (dt/4).*[k(1) k(3); k(3) k(2)] ;

% Solve the linear system
Sol.uS = LHS\RHS ;

Sol.uS = Sol.uS.*solid.dom ;

% Update the displacement
Sol.dispS = Sol.dispSPrev + (dt/2).*(Sol.uS + Sol.uSPrev) ;

end