% initCobraToolbox(false)

clearvars
clc

changeCobraSolver('gurobi')
%%
% GAPSEQ models
Bacteroides_thetaiotaomicron_VPI_5482 = readCbModel('Bacteroides_thetaiotaomicron_VPI_5482.mat');
Bifidobacterium_animalis_lactis_PT33 = readCbModel('Bifidobacterium_animalis_lactis_PT33.mat');
Clostridium_innocuum_HFG2 = readCbModel('Clostridium_innocuum_HFG2.mat');
Clostridium_sp_7_2_43FAA = readCbModel('Clostridium_sp_7_2_43FAA.mat');
Clostridium_sp_M62 = readCbModel('Clostridium_sp_M62.mat');
Clostridium_symbiosum_WAL_14673 = readCbModel('Clostridium_symbiosum_WAL_14673.mat');
Lacticaseibacillus_paracasei_M38 = readCbModel('Lacticaseibacillus_paracasei_M38.mat');

% models = {Bacteroides_thetaiotaomicron_VPI_5482; Bifidobacterium_animalis_lactis_PT33; Clostridium_innocuum_HFG2; Clostridium_sp_7_2_43FAA; Clostridium_sp_M62; Clostridium_symbiosum_WAL_14673; Lacticaseibacillus_paracasei_M38};
Names = {'Bacteroides_thetaiotaomicron_VPI_5482'; 'Bifidobacterium_animalis_lactis_PT33'; 'Clostridium_innocuum_HFG2'; 'Clostridium_sp_7_2_43FAA'; 'Clostridium_sp_M62'; 'Clostridium_symbiosum_WAL_14673'; 'Lacticaseibacillus_paracasei_M38'};
%% Export to SBML format models
% Load COBRA Toolbox and SBML Toolbox

% Define the name of the new folder
folderName = 'SBML_Export_MEMOTE';

% Create the new folder
mkdir(folderName);

% Define the list of SBML file names (without extension)
sbmlFileNames = {'Bacteroides_thetaiotaomicron_VPI_5482'; 'Bifidobacterium_animalis_lactis_PT33'; 'Clostridium_innocuum_HFG2'; 'Clostridium_sp_7_2_43FAA'; 'Clostridium_sp_M62'; 'Clostridium_symbiosum_WAL_14673'; 'Lacticaseibacillus_paracasei_M38'};

% Export each GEM to SBML format in the new folder
for i = 1:numel(sbmlFileNames)
    % Load the model from *.mat file
    model = readCbModel([sbmlFileNames{i} '.mat']);
    
    % Export the model to SBML format in the new folder
    writeCbModel(model, 'format', 'sbml', 'fileName', fullfile(folderName, [sbmlFileNames{i} '.sbml']));
end

%%
% AGORA models
Bacteroides_thetaiotaomicron_VPI_5482 = readCbModel('Bacteroides_thetaiotaomicron_VPI_5482.mat');
Bifidobacterium_animalis_lactis_PT33 = readCbModel('Bifidobacterium_animalis_lactis_PT33.mat');
Clostridium_innocuum_HGF2 = readCbModel('Clostridium_innocuum_HGF2.mat');
Clostridium_sp_7_2_43FAA = readCbModel('Clostridium_sp_7_2_43FAA.mat');
Clostridium_sp_M62_1 = readCbModel('Clostridium_sp_M62_1.mat');
Clostridium_symbiosum_WAL_14673 = readCbModel('Clostridium_symbiosum_WAL_14673.mat');
Lacticaseibacillus_paracasei_M38 = readCbModel('Lacticaseibacillus_paracasei_M38.mat');

% models = {Bacteroides_thetaiotaomicron_VPI_5482; Bifidobacterium_animalis_lactis_PT33; Clostridium_innocuum_HFG2; Clostridium_sp_7_2_43FAA; Clostridium_sp_M62; Clostridium_symbiosum_WAL_14673; Lacticaseibacillus_paracasei_M38};
Names = {'Bacteroides_thetaiotaomicron_VPI_5482'; 'Bifidobacterium_animalis_lactis_PT33'; 'Clostridium_innocuum_HGF2'; 'Clostridium_sp_7_2_43FAA'; 'Clostridium_sp_M62_1'; 'Clostridium_symbiosum_WAL_14673'; 'Lacticaseibacillus_paracasei_M38'};

%% Export to SBML format models
% Load COBRA Toolbox and SBML Toolbox

% Define the name of the new folder
folderName = 'SBML_Export_MEMOTE';

% Create the new folder
mkdir(folderName);

% Define the list of SBML file names (without extension)
sbmlFileNames = {'Bacteroides_thetaiotaomicron_VPI_5482'; 'Bifidobacterium_animalis_lactis_PT33'; 'Clostridium_innocuum_HFG2'; 'Clostridium_sp_7_2_43FAA'; 'Clostridium_sp_M62'; 'Clostridium_symbiosum_WAL_14673'; 'Lacticaseibacillus_paracasei_M38'};

% Export each GEM to SBML format in the new folder
for i = 1:numel(sbmlFileNames)
    % Load the model from *.mat file
    model = readCbModel([sbmlFileNames{i} '.mat']);
    
    % Export the model to SBML format in the new folder
    writeCbModel(model, 'format', 'sbml', 'fileName', fullfile(folderName, [sbmlFileNames{i} '.sbml']));
end