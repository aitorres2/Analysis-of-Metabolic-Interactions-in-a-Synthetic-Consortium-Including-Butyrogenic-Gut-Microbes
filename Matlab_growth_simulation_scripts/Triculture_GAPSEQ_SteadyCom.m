%% SteadyCom pipeline
clearvars
clc

%% Models load gapseq or modelseed, kbase only
% initCobraToolbox(false)
changeCobraSolver('gurobi')

Clostridium_innocuum_HFG2 = readCbModel('Clostridium_innocuum_HFG2');
Bifidobacterium_animalis_lactis_PT33 =readCbModel('Bifidobacterium_animalis_lactis_PT33');
Bacteroides_thetaiotaomicron_VPI_5482 =readCbModel('Bacteroides_thetaiotaomicron_VPI_5482');

models = {Clostridium_innocuum_HFG2; Bifidobacterium_animalis_lactis_PT33; Bacteroides_thetaiotaomicron_VPI_5482};
Names = {'Clostridium_innocuum_HFG2'; 'Bifidobacterium_animalis_lactis_PT33'; 'Bacteroides_thetaiotaomicron_VPI_5482'};

%% Establish biomass a objective function per model
Clostridium_innocuum_HFG2 = changeObjective(Clostridium_innocuum_HFG2, 'bio1_biomass', 1);  % Set the biomass reaction coefficient to 1
Bifidobacterium_animalis_lactis_PT33 = changeObjective(Bifidobacterium_animalis_lactis_PT33, 'bio1_biomass', 1);  % Set the biomass reaction coefficient to 1
Bacteroides_thetaiotaomicron_VPI_5482 = changeObjective(Bacteroides_thetaiotaomicron_VPI_5482, 'bio1_biomass', 1);  % Set the biomass reaction coefficient to 1

%% Forcing that models are not missing any necessary reaction to grow
% for the gapseq models we use this function, as these are automatically
% generated models, and i am not very sure if the model will work or not on its own.

cultureMedium = 'cultureMediumZMBinulin_NOinulin_modelseed.xlsx';
carbonSource = 'EX_cpd11602_e0';

[Bifidobacterium_animalis_lactis_PT33,Bifidobacterium_animalis_lactis_PT33_ExcsAdded] = ...
    loadCultureMedium(Bifidobacterium_animalis_lactis_PT33, cultureMedium, carbonSource);

[Clostridium_innocuum_HFG2, Clostridium_innocuum_HFG2_ExcsAdded] = ...
    loadCultureMedium(Clostridium_innocuum_HFG2, cultureMedium, carbonSource);

[Bacteroides_thetaiotaomicron_VPI_5482,Bacteroides_thetaiotaomicron_VPI_5482_ExcsAdded] = ...
    loadCultureMedium(Bacteroides_thetaiotaomicron_VPI_5482, cultureMedium, carbonSource);

%% Condition of anaerobiosis
Clostridium_innocuum_HFG2 =changeRxnBounds(Clostridium_innocuum_HFG2,'EX_cpd00007_e0',0,'b');
Bifidobacterium_animalis_lactis_PT33 =changeRxnBounds(Bifidobacterium_animalis_lactis_PT33,'EX_cpd00007_e0',0,'b');
Bacteroides_thetaiotaomicron_VPI_5482 =changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482,'EX_cpd00007_e0',0,'b');

models = {Clostridium_innocuum_HFG2; Bifidobacterium_animalis_lactis_PT33; Bacteroides_thetaiotaomicron_VPI_5482};

%% Definition of culture medium
% maximum metabolite uptake rates

medioPUC = {
    'EX_cpd11602_e0';... % inulin
    'EX_cpd02298_e0';... FOS
    'EX_cpd00027_e0';... % D-glucose
    'EX_cpd00082_e0';... % D-fructose
    'EX_cpd00076_e0';... % Sucrose
    'EX_cpd00035_e0';...'EX_ala_L(e)'
    'EX_cpd00051_e0';...'EX_arg_L(e)'
    'EX_cpd00132_e0';...'EX_asn_L(e),
    'EX_cpd00041_e0';...'EX_asp_L(e)'
    'EX_cpd00084_e0';...'EX_cys_L(e)'
    'Ex_cpd00053_e0';...'EX_gln_L(e)'
    'EX_cpd00023_e0';...'EX_glu_L(e)' 
    'EX_cpd00033_e0';...'EX_gly(e)'
    'EX_cpd00119_e0';...'EX_his_L(e)'
    'EX_cpd00322_e0';...'EX_ile_L(e)'
    'EX_cpd00107_e0';...'EX_leu_L(e)'
    'EX_cpd00039_e0';...'EX_lys_L(e)'
    'EX_cpd00060_e0';...'EX_met_L(e)'
    'EX_cpd00066 _e0';...'EX_phe_L(e)'
    'EX_cpd00129_e0';...'EX_pro_L(e)'
    'EX_cpd00054_e0';...'EX_ser_L(e)'
    'EX_cpd00161_e0';...'EX_thr_L(e)'
    'EX_cpd00065_e0';...'EX_trp_L(e)'
    'EX_cpd00069_e0';...'EX_tyr_L(e)'
    'EX_cpd00156_e0';...'EX_val_L(e)'
    'EX_cpd00009_e0';...'EX_pi(e)',...
    'EX_cpd00205_e0';...'EX_k(e)',...
    'EX_cpd00254_e0';...'EX_mg2(e)',...
    'EX_cpd10516_e0';...'EX_fe3(e)',...
    'EX_cpd10515_e0';...'EX_fe2(e)',...
    'EX_cpd00971_e0';...'EX_na1(e)',...
    'EX_cpd00034_e0';...'EX_zn2(e)',...
    'EX_cpd00030_e0';... 'EX_mn2(e)',...
    'EX_cpd11574_e0';... 'EX_mobd(e)',...
    'EX_cpd19013_e0';... 'EX_nh4(e)',...
    'EX_cpd00063_e0';...'EX_ca2(e)'
    'EX_cpd00048_e0';... 'EX_so4(e)',...
    'EX_cpd00149_e0';... 'EX_cobalt2(e)',...
    'EX_cpd00058_e0';... 'EX_cu2(e)',...
    'EX_cpd00099_e0';... 'EX_cl(e)'
    'EX_cpd02201_e0';... 'EX_pnto_R(e)'%not registered in VMH. 
    'EX_cpd00393_e0';...'EX_fol(e)'
    'EX_cpd00218_e0';... %not registered in VMH 'EX_nac(e)'
    'EX_cpd00215_e0';...'EX_pydx(e)'
    'EX_cpd00443_e0';... %aminobenzoico
    'EX_cpd00309_e0';... 'EX_xan(e)';...
    'EX_cpd00121_e0';...'EX_inost(e)';...
    'EX_cpd00111_e0';... 'EX_gthox(e)';...
    'EX_cpd00042_e0';... 'EX_gthrd(e)';...
    'EX_cpd00104_e0';... 'EX_btn(e)',...
    'EX_cpd00220_e0';... 'EX_ribflv(e)',...
    'EX_cpd00305_e0';... 'EX_thm(e)',... not registered in VMH
    'EX_cpd00541_e0';... 'EX_lipoate(e)',... 
    'EX_cpd00028_e0';... 'EX_pheme(e)',... not registered in VMH
    'EX_cpd00557_e0';... 'EX_sheme(e)',...
    'EX_cpd00128_e0';... 'EX_ade(e)',...
    'EX_cpd00207_e0';... 'EX_gua(e)',...
    'EX_cpd00092_e0';... 'EX_ura(e)'
    'EX_cpd01741_e0';... 'EX_ddca(e)'
    'EX_cpd00239_e0'; ...'EX_h2s(e)',... 
    'EX_cpd00226_e0';... 'EX_hxan(e)'
    'EX_cpd11606_e0';... 'EX_mqn7(e)'
    'EX_cpd15500_e0';... 'EX_mqn8(e)'
    'EX_cpd00244_e0';... 'EX_ni2(e)'
    'EX_cpd01080_e0';... 'EX_ocdca(e)' 0.169280000000000
    'EX_cpd15560_e0';... 'EX_q8(e)'
    'EX_cpd00184_e0';... EX_thymd(e)'
    'EX_cpd00264_e0';... Spermidine, necessary for Bi. lactis PT33
    'EX_cpd00644_e0';... PANTOTHENATE, necessary for Clostridium sp M62
    };

ratesPUC = [1;0.1;0.05;0.05;0.05;0.251761678;0.290359287;0.225000000000000;0.571946391;...
    0.035088736;0.820199195;0.820199195;0.157022092;0.200883011;...
    0.392993839;0.669317632;0.571069172;0.206146322;0.358782322;...
    0.758793908;0.445626943;0.342992391;0.092107931;0.163162621;...
    0.483347333;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0.169280000000000;1;1;1;1];

iZMB1 = [medioPUC; {'EX_cpd00001_e0'}]; % H2O
ratesiZMB1 = [ratesPUC; 1000];

%% Growth test
% here the models are constrained with the medium, and then the FBAs are made
crecimiento = zeros(length(models),2);

% Change to add a different culture medium
media = {iZMB1};
    media_rates = {ratesiZMB1};
    sol_mono = {};

for i = 1:length(models)
    for j = 1:length(media)
        medium = media{j};
        rates = media_rates{j};

        model_test = changeRxnBounds(models{i}, models{i}.rxns(cellfun(@isempty, strfind(models{i}.rxns, 'EX_'))==0), 0, 'l');
        model_test = changeRxnBounds(model_test, medium, -rates, 'l');

        Clostridium_innocuum_HFG2 = changeRxnBounds(Clostridium_innocuum_HFG2, Clostridium_innocuum_HFG2.rxns(cellfun(@isempty, strfind(Clostridium_innocuum_HFG2.rxns, 'EX_'))==0), 0, 'l');
        Clostridium_innocuum_HFG2 = changeRxnBounds(Clostridium_innocuum_HFG2, medium, -rates, 'l');    
        
        fba_model_Clostridium_innocuum_HFG2= optimizeCbModel(Clostridium_innocuum_HFG2, 'max');
        exc_idx_Clostridium_innocuum_HFG2 = cellfun(@isempty, strfind(Clostridium_innocuum_HFG2.rxns, 'EX_'))==0;
        
        Bacteroides_thetaiotaomicron_VPI_5482 = changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482, Bacteroides_thetaiotaomicron_VPI_5482.rxns(cellfun(@isempty, strfind(Bacteroides_thetaiotaomicron_VPI_5482.rxns, 'EX_'))==0), 0, 'l');
        Bacteroides_thetaiotaomicron_VPI_5482 = changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482, medium, -rates, 'l');    
        
        fba_model_Bacteroides_thetaiotaomicron_VPI_5482= optimizeCbModel(Bacteroides_thetaiotaomicron_VPI_5482, 'max');
        exc_idx_Bacteroides_thetaiotaomicron_VPI_5482 = cellfun(@isempty, strfind(Bacteroides_thetaiotaomicron_VPI_5482.rxns, 'EX_'))==0;
        
        Bifidobacterium_animalis_lactis_PT33  = changeRxnBounds(Bifidobacterium_animalis_lactis_PT33 , Bifidobacterium_animalis_lactis_PT33 .rxns(cellfun(@isempty, strfind(Bifidobacterium_animalis_lactis_PT33 .rxns, 'EX_'))==0), 0, 'l');
        Bifidobacterium_animalis_lactis_PT33  = changeRxnBounds(Bifidobacterium_animalis_lactis_PT33 , medium, -rates, 'l');
        
        fba_model_Bifidobacterium_animalis_lactis_PT33 = optimizeCbModel(Bifidobacterium_animalis_lactis_PT33 , 'max');
        exc_idx_Bifidobacterium_animalis_lactis_PT33  = cellfun(@isempty, strfind(Bifidobacterium_animalis_lactis_PT33 .rxns, 'EX_'))==0; 
        
        fba_model= optimizeCbModel(model_test);
        exc_idx = cellfun(@isempty, strfind(model_test.rxns, 'EX_'))==0;
                
        Clostridium_innocuum_HFG2_mono_table = table(Clostridium_innocuum_HFG2.rxns(find(fba_model_Clostridium_innocuum_HFG2.v(:,1) ~= 0),1), ...
        fba_model_Clostridium_innocuum_HFG2.v(find(fba_model_Clostridium_innocuum_HFG2.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});
      
        Bifidobacterium_animalis_lactis_PT33_mono_table = table(Bifidobacterium_animalis_lactis_PT33 .rxns(find(fba_model_Bifidobacterium_animalis_lactis_PT33 .v(:,1) ~= 0),1), ...
        fba_model_Bifidobacterium_animalis_lactis_PT33.v(find(fba_model_Bifidobacterium_animalis_lactis_PT33.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});
    
        Bacteroides_thetaiotaomicron_VPI_5482_mono_table = table(Bacteroides_thetaiotaomicron_VPI_5482 .rxns(find(fba_model_Bacteroides_thetaiotaomicron_VPI_5482 .v(:,1) ~= 0),1), ...
        fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(find(fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});
        
        writetable(Clostridium_innocuum_HFG2_mono_table,'Clostridium_innocuum_HFG2_monocultivo_modelseed.xlsx');
        writetable(Bifidobacterium_animalis_lactis_PT33_mono_table,'Bifidobacterium_animalis_lactis_PT33_monocultivo_modelseed.xlsx');
        writetable(Bacteroides_thetaiotaomicron_VPI_5482_mono_table,'Bacteroides_thetaiotaomicron_VPI_5482_monocultivo_modelseed.xlsx');
        
        crecimiento(i, j) = fba_model.f;
        sol_mono{end+1} = fba_model.x;
    end   
end
%% Graphical fluxes of external reactions/metabolites of interest in monocultures
% cpd00023 = L-Glutamate
% cpd00281 = GABA
% cpd00211 = butyrate
% cpd00029 = acetate
% cpd11602 = inulin
% cpd00036 = succinate
% cpd00047 = formate
% cpd00159 = L-lactate
% cpd00221 = D-lactate
% cpd01022 = lactate
% cpd00141 = propionate

targets = {'bio1_biomass'; 'EX_cpd00211_e0'; 'EX_cpd00029_e0'; 'EX_cpd00141_e0'; 'EX_cpd00159_e0'; 'EX_cpd00221_e0'; 'EX_cpd00036_e0'; 'EX_cpd00047_e0'; 'EX_cpd11602_e0'};
ids_Clostridium_innocuum_HFG2 = findRxnIDs(Clostridium_innocuum_HFG2, targets);
targets_ids_Clostridium_innocuum_HFG2 = ids_Clostridium_innocuum_HFG2(ids_Clostridium_innocuum_HFG2 > 0);
targets_filt_Clostridium_innocuum_HFG2 = targets(ids_Clostridium_innocuum_HFG2 > 0);

ids_Bifidobacterium_animalis_lactis_PT33  = findRxnIDs(Bifidobacterium_animalis_lactis_PT33 , targets);
targets_ids_Bifidobacterium_animalis_lactis_PT33  = ids_Bifidobacterium_animalis_lactis_PT33 (ids_Bifidobacterium_animalis_lactis_PT33  > 0);
targets_filt_Bifidobacterium_animalis_lactis_PT33  = targets(ids_Bifidobacterium_animalis_lactis_PT33  > 0);

ids_Bacteroides_thetaiotaomicron_VPI_5482  = findRxnIDs(Bacteroides_thetaiotaomicron_VPI_5482 , targets);
targets_ids_Bacteroides_thetaiotaomicron_VPI_5482  = ids_Bacteroides_thetaiotaomicron_VPI_5482 (ids_Bacteroides_thetaiotaomicron_VPI_5482  > 0);
targets_filt_Bacteroides_thetaiotaomicron_VPI_5482  = targets(ids_Bacteroides_thetaiotaomicron_VPI_5482  > 0);

figure
hold on
bar(categorical(strrep(targets_filt_Clostridium_innocuum_HFG2, '_', '\_')), fba_model_Clostridium_innocuum_HFG2.x(targets_ids_Clostridium_innocuum_HFG2));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{1}, '_', '\_')})
title('mZMB + Inulin Mono-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Bifidobacterium_animalis_lactis_PT33 , '_', '\_')), fba_model_Bifidobacterium_animalis_lactis_PT33 .x(targets_ids_Bifidobacterium_animalis_lactis_PT33));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{2}, '_', '\_')})
title('mZMB + Inulin Mono-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Bacteroides_thetaiotaomicron_VPI_5482, '_', '\_')), fba_model_Bacteroides_thetaiotaomicron_VPI_5482.x(targets_ids_Bacteroides_thetaiotaomicron_VPI_5482));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{3}, '_', '\_')})
title('mZMB + Inulin Mono-culture')
hold off

metabolitos_mono_Clostridium_innocuum_HFG2= table(targets_filt_Clostridium_innocuum_HFG2, ...
    fba_model_Clostridium_innocuum_HFG2.x(targets_ids_Clostridium_innocuum_HFG2), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_mono_Bifidobacterium_animalis_lactis_PT33 = table(targets_filt_Bifidobacterium_animalis_lactis_PT33 , ...
    fba_model_Bifidobacterium_animalis_lactis_PT33 .x(targets_ids_Bifidobacterium_animalis_lactis_PT33 ), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_mono_Bacteroides_thetaiotaomicron_VPI_5482 = table(targets_filt_Bacteroides_thetaiotaomicron_VPI_5482 , ...
    fba_model_Bacteroides_thetaiotaomicron_VPI_5482 .x(targets_ids_Bacteroides_thetaiotaomicron_VPI_5482 ), ... 
    'VariableNames',{'Metabolite','v'});

%% Creation of the steadycom community model
% Names=Names';
EcCom = createMultipleSpeciesModel(models, Names);
% This is because the optimizer needs the csense field
EcCom.csense = char('E' * ones(1,numel(EcCom.mets))); 
% With this we will obtain the rates and reactions of each species, and those that are exclusive of the space between bacteria (lumen).
[EcCom.infoCom, EcCom.indCom] = getMultiSpeciesModelId(EcCom,Names);

EcCom.rxns=strrep(EcCom.rxns,'[u]','_e0');
EcCom.rxns=strrep(EcCom.rxns,'[u]tr','_e0');

%% Add biomass reaction information
rxnBiomass = EcCom.rxns(find(contains(EcCom.rxns, strcat(Names, 'bio1'))));
rxnBiomassId = findRxnIDs(EcCom, rxnBiomass);
EcCom.infoCom.spBm = rxnBiomass; 
EcCom.indCom.spBm = rxnBiomassId;

%% Definition of model additional constraint
% The steps to follow are:
% 1. Add the limits for the reactions of the bacteria separately
% 2. Add the limits for the reactions of the extracellular medium
% We will use iZMB1 to optimize

iZMB2=iZMB1;

% Reaction restrictions for each bacterium
for i=1:length(Names)
    
% We set the flows of each bacterium in an unrestricted way so as not to limit the exchange between the pairs
    EcCom = changeRxnBounds(EcCom, EcCom.infoCom.EXsp(:,i),-100, 'l');
    
end

% Here we set the bounds of our culture medium in the extracellular space
EcCom = changeRxnBounds(EcCom, EcCom.rxns(startsWith(EcCom.rxns, 'EX_')), 0, 'l');
EcCom = changeRxnBounds(EcCom, iZMB2, -ratesiZMB1, 'l');
EcCom = changeRxnBounds(EcCom, EcCom.infoCom.spBm{1}, 0.1*crecimiento(1,1), 'l')

%% Optimization
changeCobraSolver('gurobi')
options = struct();
options.GRguess = 0.1; % Initial value for growth rate
options.GRtol = 1e-4; % Tolerance for growth rate 
options.algorithm = 1; %  default
tic
[sol, result] = SteadyCom(EcCom, options);
toc

Clostridium_innocuum_HFG2_table = table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,1) ~= 0),1), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,1) ~= 0),1)), ... 
    'VariableNames',{'rxn','v'})

Bifidobacterium_animalis_lactis_PT33_table=table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,2) ~= 0),2), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,2) ~= 0),2)), ... 
    'VariableNames',{'rxn','v'})

Bacteroides_thetaiotaomicron_VPI_5482_table=table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,3) ~= 0),3), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,3) ~= 0),3)), ... 
    'VariableNames',{'rxn','v'})

com_table = table(EcCom.rxns(cellfun(@isempty, strfind(EcCom.rxns, 'EX_'))==0), ... 
    sol.full(cellfun(@isempty, strfind(EcCom.rxns, 'EX_'))==0), ... 
    'VariableNames',{'rxn','v'})

%% export excel tables with consumed/produced metabolites with steadycomFBA
Clostridium_innocuum_HFG2_table(~Clostridium_innocuum_HFG2_table.v,:) = []; % https://www.mathworks.com/matlabcentral/answers/407401-delete-zero-rows-from-a-table
Bifidobacterium_animalis_lactis_PT33_table(~Bifidobacterium_animalis_lactis_PT33_table.v,:) = [];
Bacteroides_thetaiotaomicron_VPI_5482_table(~Bacteroides_thetaiotaomicron_VPI_5482_table.v,:) = [];
com_table(~com_table.v,:) = [];

writetable(Clostridium_innocuum_HFG2_table,'Clostridium_innocuum_HFG2_consorcio_con_Bifidobacterium_animalis_lactis_PT33_Bacteroides_thetaiotaomicron_VPI_5482_modelseed.xlsx');
writetable(Bifidobacterium_animalis_lactis_PT33_table,'Bifidobacterium_animalis_lactis_PT33_consorcio_con_Clostridium_innocuum_HFG2_Bacteroides_thetaiotaomicron_VPI_5482_modelseed.xlsx');
writetable(Bacteroides_thetaiotaomicron_VPI_5482_table,'Bacteroides_thetaiotaomicron_VPI_5482_consorcio_con_Clostridium_innocuum_HFG2_Bifidobacterium_animalis_lactis_PT33_modelseed.xlsx');
writetable(com_table,'EX_comunidad_Clostridium_innocuum_HFG2_Bifidobacterium_animalis_lactis_PT33_Bacteroides_thetaiotaomicron_VPI_5482_modelseed.xlsx');

%% Graph for community external metabolites such as butyrate
targets_com_Clostridium_innocuum_HFG2 = {'Clostridium_innocuum_HFG2bio1_biomass';
    'Clostridium_innocuum_HFG2IEX_cpd00211_e0tr';
    'Clostridium_innocuum_HFG2IEX_cpd00029_e0tr';
    'Clostridium_innocuum_HFG2IEX_cpd00141_e0tr';
    'Clostridium_innocuum_HFG2IEX_cpd00159_e0tr'; 
    'Clostridium_innocuum_HFG2IEX_cpd00221_e0tr'; 
    'Clostridium_innocuum_HFG2IEX_cpd00036_e0tr'; 
    'Clostridium_innocuum_HFG2IEX_cpd00047_e0tr'; 
    'Clostridium_innocuum_HFG2IEX_cpd11602_e0tr'}; 
    
 targets_com_Bifidobacterium_animalis_lactis_PT33 = {'Bifidobacterium_animalis_lactis_PT33bio1_biomass';
    'Bifidobacterium_animalis_lactis_PT33IEX_cpd00211_e0tr';
    'Bifidobacterium_animalis_lactis_PT33IEX_cpd00029_e0tr';
    'Bifidobacterium_animalis_lactis_PT33IEX_cpd00141_e0tr';
    'Bifidobacterium_animalis_lactis_PT33IEX_cpd00159_e0tr'; 
    'Bifidobacterium_animalis_lactis_PT33IEX_cpd00221_e0tr'; 
    'Bifidobacterium_animalis_lactis_PT33IEX_cpd00036_e0tr'; 
    'Bifidobacterium_animalis_lactis_PT33IEX_cpd00047_e0tr'; 
    'Bifidobacterium_animalis_lactis_PT33IEX_cpd11602_e0tr'};

 targets_com_Bacteroides_thetaiotaomicron_VPI_5482 = {'Bacteroides_thetaiotaomicron_VPI_5482bio1_biomass';
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_cpd00211_e0tr';
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_cpd00029_e0tr';
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_cpd00141_e0tr';
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_cpd00159_e0tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_cpd00221_e0tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_cpd00036_e0tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_cpd00047_e0tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_cpd11602_e0tr'};
   
targets_com = {'EX_cpd00211_e0';
    'EX_cpd00029_e0';
    'EX_cpd00141_e0';
    'EX_cpd00159_e0'; 
    'EX_cpd00221_e0'; 
    'EX_cpd00036_e0'; 
    'EX_cpd00047_e0'; 
    'EX_cpd11602_e0'};
  
ids_Com_Clostridium_innocuum_HFG2 = findRxnIDs(EcCom, targets_com_Clostridium_innocuum_HFG2);
targets_ids_Com_Clostridium_innocuum_HFG2 = ids_Com_Clostridium_innocuum_HFG2(ids_Com_Clostridium_innocuum_HFG2 > 0);
targets_filt_Com_Clostridium_innocuum_HFG2 = targets_com_Clostridium_innocuum_HFG2(ids_Com_Clostridium_innocuum_HFG2 > 0);

ids_Com_Bifidobacterium_animalis_lactis_PT33 = findRxnIDs(EcCom, targets_com_Bifidobacterium_animalis_lactis_PT33);
targets_ids_Com_Bifidobacterium_animalis_lactis_PT33 = ids_Com_Bifidobacterium_animalis_lactis_PT33(ids_Com_Bifidobacterium_animalis_lactis_PT33 > 0);
targets_filt_Com_Bifidobacterium_animalis_lactis_PT33 = targets_com_Bifidobacterium_animalis_lactis_PT33(ids_Com_Bifidobacterium_animalis_lactis_PT33 > 0);

ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 = findRxnIDs(EcCom, targets_com_Bacteroides_thetaiotaomicron_VPI_5482);
targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 = ids_Com_Bacteroides_thetaiotaomicron_VPI_5482(ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 > 0);
targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482 = targets_com_Bacteroides_thetaiotaomicron_VPI_5482(ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 > 0);

ids_Com = findRxnIDs(EcCom, targets_com);
targets_ids_Com = ids_Com(ids_Com > 0);
targets_filt_Com = targets_com(ids_Com > 0);

figure
hold on
bar(categorical(strrep(targets_filt_Com_Clostridium_innocuum_HFG2, '_', '\_')), sol.full(targets_ids_Com_Clostridium_innocuum_HFG2));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{1}, '_', '\_')});
title('mZMB + Inulin Co-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com_Bifidobacterium_animalis_lactis_PT33, '_', '\_')), sol.full(targets_ids_Com_Bifidobacterium_animalis_lactis_PT33));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{2}, '_', '\_')});
title('mZMB + Inulin Co-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482, '_', '\_')), sol.full(targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{3}, '_', '\_')});
title('mZMB + Inulin Co-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com, '_', '\_')), sol.full(targets_ids_Com));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{1}, '_', '\_') strrep(Names{2}, '_', '\_') strrep(Names{3}, '_', '\_')});
title('mZMB + Inulin Co-culture')
hold off

metabolitos_com_Clostridium_innocuum_HFG2_table = table(targets_filt_Com_Clostridium_innocuum_HFG2, ...
    sol.full(targets_ids_Com_Clostridium_innocuum_HFG2), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_com_Bifidobacterium_animalis_lactis_PT33_table = table(targets_filt_Com_Bifidobacterium_animalis_lactis_PT33, ...
    sol.full(targets_ids_Com_Bifidobacterium_animalis_lactis_PT33), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_com_Bacteroides_thetaiotaomicron_VPI_5482_table = table(targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482, ...
    sol.full(targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_comunidad_totales_table = table(targets_filt_Com, ...
    sol.full(targets_ids_Com), ... 
    'VariableNames',{'Metabolite','v'});

%% SteadyCom-FVA
options.optBMpercent = 100;
options.rxnNameList = strcat('X_', Names);
options.optGRpercent = [0:1:100]; 
options.GRtol = 1e-6;
options.threads = 0;

[fvaComMin, fvaComMax, minFD, maxFD, GRvector, result, LP] = SteadyComFVA(EcCom, options);

%% Envelope plot for relative abundance by community growth rate
grComV = result.GRmax * options.optGRpercent / 100;
lgLabel = {strrep(Names{1}, '_', '\_');strrep(Names{2}, '_', '\_');;strrep(Names{3}, '_', '\_')};

x_lower = 10;
x_upper = 0;
close all
figure
hold on
grid on
x = [grComV(:); flipud(grComV(:))];
y_max = 0;
for j = 1:length(Names)
    y = [fvaComMin(j, :), fliplr(fvaComMax(j, :))];
    y = y(~isnan(y));
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

legend(lgLabel)
xlabel('Growth rate (h^{-1})');
ylabel('Relative abundance');
title('\underline{SteadyComFVA biomass}', 'Interpreter', 'latex');
hold off

%% FVA optimizing for metabolite of interest (butyrate)
options.optBMpercent = 100;
options.rxnNameList = EcCom.rxns(find(contains(EcCom.rxns, 'IEX_cpd00211')));
options.optGRpercent = [0:1:100];
options.GRtol = 1e-6;

[fvaComMin_B, fvaComMax_B, minFD_B, maxFD_B, GRvector_B, result_B, LP_B] = SteadyComFVA(EcCom, options);

%%
grComV = result.GRmax * options.optGRpercent / 100;
% lgLabel = {strrep(Names{2}, '_', '\_')};  % poner a la bacteria productora de butirato, nombre2
lgLabel = {strrep(Names{1}, '_', '\_');strrep(Names{3}, '_', '\_')}; % para 2 bacterias productoras del metabolito de interes, nombre
f = figure;
x_lower = 10;
x_upper = 0; 
hold on
grid on
x = [grComV(:); flipud(grComV(:))];
y_max = 0;
for j = 1:length(options.rxnNameList) % si son dos bacterias j = 1:length(options.rxnNameList); si es 1 bacteria: j = 1:length(options.rxnNameList(1)) 
    y = [fvaComMin_B(j, :), fliplr(fvaComMax_B(j, :))]; 
    y_max=max(y);
    y = y(~isnan(y));
    if crecimiento(j,1) == max(crecimiento(:,1))
        if min(y) < x_lower 
            x_lower = min(y);
        end 
        if max(y) >= x_upper
            x_upper = y_max;
        end 
        if max(y) >= y_max
            y_max = max(y)
        end 
    end 
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

tl(1) = title('\underline{SteadyComFVA butyrate}', 'Interpreter', 'latex');
lg = legend(lgLabel);
% lg.Box = 'off';
yl(1) = ylabel('Butyrate flux (mmol gDW-1 h-1)');
xl(1) = xlabel('Community growth rate (h^{-1})');
hold off


%% Combined figures
lgLabel1 = {strrep(Names{1}, '_', '\_'); strrep(Names{2}, '_', '\_'); strrep(Names{3}, '_', '\_')}; % poner todos los modelos mutantes

close all

% Create the main figure
figure('Position', [100, 100, 1200, 600])

% Define subplot positions
subplotPositions = [0.1 0.1 0.4 0.8; 0.55 0.1 0.4 0.8];

% Plot the first subplot
subplot('Position', subplotPositions(1, :))
hold on
grid on
x = [grComV(:); flipud(grComV(:))];
y_max = 0;
for j = 1:length(Names) % de 1 a n�mero de modelos mutantes
    y = [fvaComMin(j, :), fliplr(fvaComMax(j, :))];
    y = y(~isnan(y));
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

legend(lgLabel1)
xlabel('Growth rate (h^{-1})');
ylabel('Relative abundance'); % Relative abundance Bacteria growth rate (h^{-1})
title('\underline{SteadyComFVA biomass}', 'Interpreter', 'latex');

% Plot the second subplot
subplot('Position', subplotPositions(2, :))
hold on
grid on
x = [grComV(:); flipud(grComV(:))];
y_max = 0;
for j = 1:length(options.rxnNameList) % si son dos bacterias j = 1:length(options.rxnNameList); si es 1 bacteria: j = 1:length(options.rxnNameList(1))
    y = [fvaComMin_B(j, :), fliplr(fvaComMax_B(j, :))];
    y_max = max(y);
    y = y(~isnan(y));
    if crecimiento(j, 1) == max(crecimiento(:, 1))
        if min(y) < x_lower
            x_lower = min(y);
        end
        if max(y) >= x_upper
            x_upper = y_max;
        end
        if max(y) >= y_max
            y_max = max(y);
        end
    end
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

tl(1) = title('\underline{SteadyComFVA butyrate}', 'Interpreter', 'latex');
lgLabel2 = {strrep(Names{1}, '_', '\_')}; % poner a la bacteria productora de butirato, nombre2
lg = legend(lgLabel2);
yl(1) = ylabel('Butyrate flux (mmol gDW^{-1} h^{-1})');
xl(1) = xlabel('Community growth rate (h^{-1})');

% Save the combined figure
saveas(p(j, 1),'EX_combinado_Clostridium_innocuum_HFG2_Bifidobacterium_animalis_lactis_PT33_Bacteroides_thetaiotaomicron_VPI_5482_modelseed.png');
hold off
