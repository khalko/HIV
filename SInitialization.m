function SInitialization
% SInitialization Generates new SSFinder file for the model
%
% Latest revision 08.01.2021 
%
% Authors: M.Yu. Khristichenko (INM RAS)
%          Yu.M. Nechepurenko  (INM RAS)
%          E.V.  Sklyarova     (MIPT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING PARAMETERS:
%
% Folder name where the model files are located:
ModelFolderName     = 'LCMV'; 

% Name of the SSFinder file that will be generated in the given folder:
NewSSFinderFileName = 'SSFinder_LCMV.m';

% Name of the RHS file in the given folder:
ModelRHSFileName    = 'RHS_LCMV.m';

% Name of the SRHS file in the given folder (if you don't want use the SRHS 
% file to generate SSFinder write '' instead of the file name):
ModelSRHSFileName   = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Functions/Utilities');

[thisfilepath,~,~] = fileparts(mfilename('fullpath'));
modelFolderPath = [thisfilepath, '\', ModelFolderName, '\'];

paths = cellfun(@(a)[modelFolderPath, a], ...
    {NewSSFinderFileName ModelRHSFileName ModelSRHSFileName}, "uni", false);

if (strcmp(ModelSRHSFileName, ''))
    generateSSFinder(paths{1}, paths{2});
else 
    generateSSFinder(paths{1}, paths{2}, paths{3});
end 
     
end
