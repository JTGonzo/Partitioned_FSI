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
%algorithm = 'Constant';
algorithm = 'AndersonAcceleration';
Couple.omegaMax = 0.75; %0.1        % relaxation factor
Couple.reuse = 8;                   % time-steps to reuse
Couple.small = 1e-12;               % filtering thershold
Couple.filter = 1;                  % QR1 = 1; QR2 = 2; POD =3

%% Initialize fixed-point iteration variables and data output files
r = zeros(length(uS),1);              % initial displacement residual 
Count  = zeros(tf/dt,2);              % iteration counter
omega  = Couple.omegaMax;             % initial relaxation factor
Wid = 2;                              % watchpoint Id
WidF = 36;                              % watchpoint Id

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
    case 'gaussseidel'
        Acc_flag = 1;
        model = GSeidel();
        omega  = 1;          
    case 'constant'  
        Acc_flag = 2;
        model=Constant();   
        omega = .5;    
    case 'aitken'  
        Acc_flag = 3;
        model = Aitken(); 
    case 'nifc'   
        Acc_flag = 4;
        p_it = ones(size(M.xr,1),1);
        q_it = zeros(size(M.xr,1),1);
        model = NIFC();   
        %% WIP 
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
        model = GB(Couple.small, Couple.reuse);
    case 'andersonacceleration'
        Acc_flag = 7;
        model = AA(Couple.small,Couple.reuse, Couple.filter);
    case 'mvls'
        Acc_flag = 8;
        limit = 60;
        model = MVLS(Couple.small,Couple.limit, Couple.filter);  
    case 'mvlss'
        Acc_flag = 9;
        limit = 40;
        model = MVLSS(Couple.small,Couple.limit);
    otherwise
        error('Unknown algorithm');
end

%% Additional functional tools (step's inital variable approx. / convergence criteria)
extrapolator = Extrapolator(zeros(length(MESH.Solid.Gamma_global),1));
convergence = Convergence(Couple.rtol,Couple.imin,Couple.imax,Couple.maxSteps,Couple.small,Couple.out, Acc_flag);
objects = {convergence,extrapolator,model};
