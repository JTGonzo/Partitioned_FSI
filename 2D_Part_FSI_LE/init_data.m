%% Definte paths, folders and outputs locations
[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');

wrkDir = './' ;
problemString = 'Flap'; %'Flap' ; 'Mok'; 'Cavity'; 
addpath('C_Files/')
addpath('Mappers/')
addpath('Accelerators/')
addpath('Accelerators/Filters/')
addpath('import_export_time/')
vtk_filename = []; %'Figures/Mok_'; % 'Figures/PerpFlap_'; %[]; %'Figures/Turek_'; 'Figures/Mok_';

t = [];
param = [];
dim  =  2;  %3
MESH.dim  = dim; % problem dimensionality 

%% Load problem mesh and element properties
% Solid mesh
load('flap_S.mat') %'flap_S.mat' 'Mok_S.mat' 'Cavity_S.mat'
meshSolid.boundaries = boundaries;
meshSolid.elements = elements;
meshSolid.vertices = vertices;
meshSolid.rings = rings;

% Fluid mesh
load('flap_F.mat') ; %'flap_F.mat' 'Mok_S.mat' 'Cavity_S.mat'
meshFluid.boundaries = boundaries;
meshFluid.elements = elements;
meshFluid.vertices = vertices;
meshFluid.rings = rings;

use_SUPG = false;
quad_order   = 4; % for 3D problem quad_order = 5 
fem_F = {'P2', 'P1'}; 
fem_S = 'P2';
MESH.mapper = 1;                  %| NN = 1; Linear = 2; RBF =3

%% Load problem, solver and domain variable details
data_file_F = 'CFD_data';
data_file_S = 'CSM_data';

% read fluid domain parameters/settings
eval(data_file_F);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)    
     eval(['DATA.Fluid.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
end
DATA.Fluid.param = param;

clear data_fields_name data

% read solid domain parameters/settings
eval(data_file_S);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)    
    eval(['DATA.Solid.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
end
DATA.Solid.param = param;