%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Initialize Coupling solver setting, export files %%%%%%%%%%
%%%%%% convergence/extrapolation schemes and accleration method  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed-Point iteration parameters
Couple.rtol = 1e-4;                 % relative residual error
Couple.imin = 2;                    % min number of fix-point iterations
Couple.imax = 50;                   % max number of fix-point iterations
Couple.maxSteps = tf/dt;            % number of time steps
Couple.probeOutFreq = 1;            % physical variable output frequency
Couple.out = 1;                     % output frequency for coupling data

% Convergence acceleration parameters
algorithm = 'AndersonAcceleration';
Couple.omegaMax = 0.75; %0.1        % relaxation factor
Couple.reuse = 10;                   % time-steps to reuse
Couple.small = 1e-10;               % filtering thershold
Couple.filter = 2;                  % QR1 = 1; NM = 2; POD =3
Couple.count = 0;                   %| filtered columns

%% Initialize fixed-point iteration variables and data output files
r = zeros(length(uS),1);              % initial displacement residual 
Count  = zeros(tf/dt,2);              % iteration counter
omega  = Couple.omegaMax;             % initial relaxation factor
Wid = 2 ;  %5; - truek problem;                                 % watchpoint Id - linear elastic problem
WidF = 36;  %192; - turek problem;                              % watchpoint Id - turek problem

filename1 = sprintf('Results/FP_Iterations.oisd');              % iterations residual info data file
fileId1 = fopen(filename1,'w');
filename2 = sprintf('Results/%s_Disp.othd',problemString);      % variable info data file
fileId2 = fopen(filename2,'w');
filename3 = sprintf('Results/%s_Vel.othd',problemString);
fileId3 = fopen(filename3,'w');

%% Initalizing Vairables based on what coupling scheme is chosen
switch lower(algorithm)
    case 'explicit'
        Acc_flag = 0;
        model = Dummy();
    case 'gaussseidel'
        Acc_flag = 1;
        model = Dummy();
        omega  = 1;          
    case 'constant'  
        Acc_flag = 2;
        model = Dummy(); 
        omega = .5;    
    case 'aitken'  
        Acc_flag = 3;
        model = Aitken(); 
    case 'nifc'   
        Acc_flag = 4;
        p_it = ones(size(M.xr,1),1);
        q_it = zeros(size(M.xr,1),1);
        model = NIFC();   
%     case 'gmres'  
%         model=GMRES(small,reuse);
%         omega=0.01;
%         flag_t = 6;    
    case 'broydensecond'
        Acc_flag = 5;
        Couple.reuse = 0;                   
        model = Broyden(Couple.reuse);
    case 'genbroyden'
        Acc_flag = 6;
        Couple.reuse = 0;
        model = GB(Couple.small, Couple.reuse, Couple.count); %GB3
    case 'andersonacceleration'
        Acc_flag = 7;
        Couple.lindepen = 1e-7;     %| filtering threshold
        model = AA(Couple.lindepen,Couple.reuse, Couple.filter, Couple.count);
    case 'mvls'
        Acc_flag = 8;
        limit = 15;
        model = MVLS(Couple.small,limit, Couple.count);  
    case 'mvlss'
        Acc_flag = 9;
        limit = 15;
        model = MVLSS(Couple.small, limit, Couple.count);
    otherwise
        error('Unknown algorithm');
end

%% Additional functional tools (step's inital variable approx. / convergence criteria)
extrapolator = Extrapolator(zeros(length(MESH.Solid.Gamma_global),1));
convergence = Convergence(Couple.rtol,Couple.imin,Couple.imax,Couple.maxSteps,Couple.small,Couple.out, Acc_flag);
objects = {convergence,extrapolator,model};
