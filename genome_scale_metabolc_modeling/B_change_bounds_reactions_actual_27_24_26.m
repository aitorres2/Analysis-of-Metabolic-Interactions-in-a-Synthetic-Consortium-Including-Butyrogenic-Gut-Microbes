clear all
clc

%% checkear raven, cobra, y gurobi:
%checkInstallation
% initCobraToolbox(false)
changeCobraSolver('gurobi')

% https://sysbiochalmers.github.io/RAVEN/doc/
% https://opencobra.github.io/cobratoolbox/stable/index.html

% agregar reacciones: 
% https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorial_modelManipulation.html

% conceptos basicos:
% https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorial_COBRAconcepts.html

% 1) FALTA UN ERROR MAS POR RESOLVER, ES QUE NINGUNA DE LAS BPBS PRODUCE
% BUTIRATO  Y EL BTHETHA SI PRODUCE BUTIRATO... SOLUCIONAR MANUALMENTE

% 2) AGREGAR LAS GPR RULES PARA GENES BUTIRATO, ACETATO, LACTATO, INULINA,
% FOS, TRASNPORTADORES, ETC, PARA INTEGRAR CORRECTAMENTE EL RNASEQ EN PT33
% Y HGF2

%% cargar modelo con cobra

cd /media/alexis/hdd2/postdoc/0_Simulaciones_MATLAB_paper_HGF2_PT33_correcciones/00_nuevos_modelos_metabolicos_carveme_gapseq/1_output_modelos_carveme
addpath(genpath(pwd));  % agregar subdirectorios al path para que readCbModel encuentre los .mat

% Cargar modelos limpios (output de fix_models_cobrapy.py)
% COBRApy exporta SBML con IDs limpios (M_xxx_c) y compartimentos correctos (C_c, C_p, C_e)
% writeCbModel de MATLAB manglea IDs y pierde compartimentos, por eso usamos COBRApy
% Todos los modelos: v3_clean_cobrapy.xml (output de fix_models_cobrapy.py)
% Bt v3_clean_cobrapy mantiene sinks de gapfill (necesarios para anaerobiosis)
Bacteroides_thetaiotaomicron_VPI_5482 = readCbModel('Bacteroides_thetaiotaomicron_VPI_5482_pangenome_metanex_bigg_polished_masschargecurated.mat');
Bifidobacterium_animalis_lactis_PT33 = readCbModel('Bifidobacterium_animalis_lactis_PT33_pangenome_metanex_bigg_polished_masschargecurated.mat');
Clostridium_sp_7_2_43FAA = readCbModel('Clostridium_sp_7_2_43FAA_pangenome_metanex_bigg_polished_masschargecurated.mat');
Clostridium_sp_HGF2 = readCbModel('Clostridium_sp_HGF2_v2_pangenome_metanex_bigg_polished_masschargecurated.mat');
Clostridium_sp_L2_50 = readCbModel('Clostridium_sp_L2_50_pangenome_metanex_bigg_polished_masschargecurated.mat');
Clostridium_sp_M62_1 = readCbModel('Clostridium_sp_M62_1_pangenome_metanex_bigg_polished_masschargecurated.mat');
Clostridium_symbiosum_WAL_14673 = readCbModel('Clostridium_symbiosum_WAL_14673_pangenome_metanex_bigg_polished_masschargecurated.mat');
Lacticaseibacillus_paracasei_M38 = readCbModel('Lacticaseibacillus_paracasei_M38_pangenome_metanex_bigg_polished_masschargecurated.mat');

modelos_raw = {Bacteroides_thetaiotaomicron_VPI_5482, ...
               Bifidobacterium_animalis_lactis_PT33, ...
               Clostridium_sp_7_2_43FAA, ...
               Clostridium_sp_HGF2, ...
               Clostridium_sp_L2_50, ...
               Clostridium_sp_M62_1, ...
               Clostridium_symbiosum_WAL_14673, ...
               Lacticaseibacillus_paracasei_M38};

nombres = {'Bacteroides_thetaiotaomicron_VPI_5482', ...
           'Bifidobacterium_animalis_lactis_PT33', ...
           'Clostridium_sp_7_2_43FAA', ...
           'Clostridium_sp_HGF2', ...
           'Clostridium_sp_L2_50', ...
           'Clostridium_sp_M62_1', ...
           'Clostridium_symbiosum_WAL_14673', ...
           'Lacticaseibacillus_paracasei_M38'};

%% Normalizar prefijos M_ y R_ en metabolitos y reacciones
% Los .mat exportados por COBRApy guardan IDs sin prefijo (accoa_c, EX_glc__D_e),
% pero todo el script B referencia M_accoa_c y R_EX_glc__D_e. Sin esta
% normalizacion, cada addReaction crea metabolitos M_xxx_c desconectados de los
% originales xxx_c, generando "islas" que no se conectan al metabolismo central
% y resultan en mu=0. Aqui forzamos prefijos M_ (mets) y R_ (rxns) en todos los
% modelos, idempotente: si ya tienen prefijo, no se duplica.
fprintf('\n=== NORMALIZANDO PREFIJOS M_ / R_ ===\n');
for i_norm = 1:length(modelos_raw)
    m_norm = modelos_raw{i_norm};
    n_mets_fix = 0;
    for j_norm = 1:length(m_norm.mets)
        if ~startsWith(m_norm.mets{j_norm}, 'M_')
            m_norm.mets{j_norm} = ['M_' m_norm.mets{j_norm}];
            n_mets_fix = n_mets_fix + 1;
        end
    end
    n_rxns_fix = 0;
    for j_norm = 1:length(m_norm.rxns)
        if ~startsWith(m_norm.rxns{j_norm}, 'R_')
            m_norm.rxns{j_norm} = ['R_' m_norm.rxns{j_norm}];
            n_rxns_fix = n_rxns_fix + 1;
        end
    end
    modelos_raw{i_norm} = m_norm;
    fprintf('  %s: %d mets + %d rxns con prefijo agregado\n', ...
            nombres{i_norm}, n_mets_fix, n_rxns_fix);
end
% Reasignar variables individuales (necesario porque modelos_raw es copia)
Bacteroides_thetaiotaomicron_VPI_5482 = modelos_raw{1};
Bifidobacterium_animalis_lactis_PT33  = modelos_raw{2};
Clostridium_sp_7_2_43FAA              = modelos_raw{3};
Clostridium_sp_HGF2                   = modelos_raw{4};
Clostridium_sp_L2_50                  = modelos_raw{5};
Clostridium_sp_M62_1                  = modelos_raw{6};
Clostridium_symbiosum_WAL_14673       = modelos_raw{7};
Lacticaseibacillus_paracasei_M38      = modelos_raw{8};
fprintf('=== PREFIJOS NORMALIZADOS ===\n\n');

%%
% surfNet(Bacteroides_thetaiotaomicron_VPI_5482, 'fru_c')
% surfNet(Bacteroides_thetaiotaomicron_VPI_5482, 'fru_e')
% surfNet(Bacteroides_thetaiotaomicron_VPI_5482, 'fru_p')

%% PARA AGREGAR INULINA  (METABOLITO) Y REACCIONES RELEVANTES UTILIZAR MISMAS REACCIONES Y METABOLITOS QUE EN VMH/AGORA
% Clostridium_innocuum_2959 = importModel(which('Clostridium_innocuum_2959.xml'));
% Bifidobacterium_animalis_subsp_lactis_BLC1 = importModel(which('Bifidobacterium_animalis_subsp_lactis_BLC1.xml'));

% ejemplo de gram negativa para fructosa:
% Escherichia_coli_str_K_12_substr_MG1655 = readCbModel(which('iML1515_Escherichia_coli_str_K_12_substr_MG1655.mat')); % BiGG model
% 
% surfNet(Escherichia_coli_str_K_12_substr_MG1655, 'fru_c')
% surfNet(Escherichia_coli_str_K_12_substr_MG1655, 'fru_e')
% surfNet(Escherichia_coli_str_K_12_substr_MG1655, 'fru_p')
% 
% %%
% % ejemplo de gram positiva para fructosa:
% Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 = readCbModel(which('iYS854_Staphylococcus_aureus_subsp_aureus_USA300_TCH1516.mat')); % BiGG model
% 
% surfNet(Staphylococcus_aureus_subsp_aureus_USA300_TCH1516, 'fru_c')
% surfNet(Staphylococcus_aureus_subsp_aureus_USA300_TCH1516, 'fru_e')
% surfNet(Staphylococcus_aureus_subsp_aureus_USA300_TCH1516, 'fru_p')

%%
% surfNet(Clostridium_innocuum_2959, 'EX_inulin(e)') % inulina en VMH de 29 fructosas
% surfNet(Clostridium_innocuum_2959, 'INULINabc')
% surfNet(Clostridium_innocuum_2959, 'INULINASE')
% surfNet(Clostridium_innocuum_2959, 'INULINASEe')
% 
% surfNet(Clostridium_innocuum_2959, 'EX_kestopt(e)') % FOS en VMH de 4 fructosas
% surfNet(Clostridium_innocuum_2959, 'KESTOPTASEe')
% 
% surfNet(Clostridium_innocuum_2959, 'EX_kesto(e)')
% 
% surfNet(Clostridium_innocuum_2959, 'EX_kestottr(e)')
% 
% % Con un metabolito
% surfNet(Clostridium_innocuum_2959, 'inulin[e]')
% surfNet(Clostridium_innocuum_2959, 'inulin[c]')
% 
% %%
% % Con una reacción
% surfNet(Bifidobacterium_animalis_subsp_lactis_BLC1, 'EX_inulin(e)') % inulina en VMH de 29 fructosas
% surfNet(Bifidobacterium_animalis_subsp_lactis_BLC1, 'INULINabc')
% surfNet(Bifidobacterium_animalis_subsp_lactis_BLC1, 'INULINASE')
% surfNet(Bifidobacterium_animalis_subsp_lactis_BLC1, 'INULINASEe')
% 
% surfNet(Bifidobacterium_animalis_subsp_lactis_BLC1, 'EX_kestopt(e)') % FOS en VMH de 4 fructosas
% surfNet(Bifidobacterium_animalis_subsp_lactis_BLC1, 'KESTOPTASEe')
% 
% surfNet(Bifidobacterium_animalis_subsp_lactis_BLC1, 'EX_kesto(e)')
% 
% surfNet(Bifidobacterium_animalis_subsp_lactis_BLC1, 'EX_kestottr(e)')
% 
% EX_kestottr(e)
% % Con un metabolito
% surfNet(Bifidobacterium_animalis_subsp_lactis_BLC1, 'inulin[e]')
% surfNet(Bifidobacterium_animalis_subsp_lactis_BLC1, 'inulin[c]')

%% Inulin reactions in VMH

% EX_inulin(e)
% Inulin exchange
% inulin[e] <=>
% Exchange/demand reaction

% INULINabc
% Inulin import through ABC transport system
% atp[c] + h2o[c] + inulin[e] -> adp[c] + h[c] + inulin[c] + pi[c]
% Transport, extracellular

% INULINASE
% Inulin degradation by beta-2->1-fructanase
% 29.0 h2o[c] + inulin[c] -> 29.0 fru[c] + glc_D[c]
% Plant polysaccharide degradation

% INULINASEe
% Inulin degradation, extracellular
% 25.0 h2o[e] + inulin[e] -> 25.0 fru[e] + kestopt[e]
% Plant polysaccharide degradation

% metabolites

% inulin
% Inulin
% C180H302O151

% kesto
% Kestose (2 fru, 1 glc inulin-type fructo-oligosaccharide)
% C18H32O16

% kestopt
% Kestopentaose (4 fru, 1 glc inulin-type fructo-oligosaccharide)
% C30H52O26

% kestottr
% Kestotetraose (3 fru, 1 glc inulin-type fructo-oligosaccharide)
% C24H42O21


%% FOS reactions in VMH

% EX_kestopt(e)
% Kestopentaose exchange
% kestopt[e] <=>
% Exchange/demand reaction
% 
% EX_kestottr(e)
% Kestotetraose exchange
% kestottr[e] <=>
% Exchange/demand reaction

% EX_kesto(e)
% Kestose exchange
% kesto[e] <=>
% Exchange/demand reaction

% KESTOASEe
% Beta-2-1 fructan degradation, extracellular
% 2.0 h2o[e] + kesto[e] -> 2.0 fru[e] + glc_D[e]
% Plant polysaccharide degradation

% KESTOTTRASEe
% Beta-2-1 fructan degradation, extracellular
% 3.0 h2o[e] + kestottr[e] -> 3.0 fru[e] + glc_D[e]
% Plant polysaccharide degradation

% KESTOPTASEe
% Beta-2-1 fructan degradation, extracellular
% 4.0 h2o[e] + kestopt[e] -> 4.0 fru[e] + glc_D[e]
% Plant polysaccharide degradation

% INULINASEe
% Inulin degradation, extracellular
% 25.0 h2o[e] + inulin[e] -> 25.0 fru[e] + kestopt[e]
% Plant polysaccharide degradation

%% AÑADIR METABOLITOS Y REACCIONES DE INULINA/FOS A MODELOS CARVEME
% Bt, PT33, C.symbiosum → kestose
% HGF2                  → inulina (sin FOS)
% M38                   → inulina + todos los FOS
%
% IMPORTANTE: Los modelos _v3_clean.xml (CarveMe + RAVEN import + COBRA export)
% usan formato BiGG con prefijo M_ y sufijo _c (citosol) / _p (periplamsa/extracelular).
% NO usan formato VMH con corchetes [c]/[e].
% Los metabolitos _e (boundaryCondition=true) fueron eliminados durante el import.
% El compartimento _p actua como extracelular en estos modelos.
%
% Mapeo VMH -> BiGG (CarveMe):
%   h2o[c]    -> M_h2o_c       h2o[e]    -> M_h2o_p
%   fru[c]    -> M_fru_c       fru[e]    -> M_fru_p
%   glc_D[c]  -> M_glc__D_c    glc_D[e]  -> M_glc__D_p
%   atp[c]    -> M_atp_c       adp[c]    -> M_adp_c
%   h[c]      -> M_h_c         pi[c]     -> M_pi_c
%   inulin[e] -> M_inulin_p    inulin[c] -> M_inulin_c
%   kesto[e]  -> M_kesto_p     kestopt[e]-> M_kestopt_p
%   kestottr[e] -> M_kestottr_p

% %% DIAGNOSTICO: verificar IDs de metabolitos existentes
% fprintf('=== DIAGNOSTICO DE IDs DE METABOLITOS ===\n');
% mets_check = {'h2o', 'fru', 'glc__D', 'atp', 'adp', 'pi', 'h'};
% modelo_test = Bifidobacterium_animalis_lactis_PT33;
% for j = 1:length(mets_check)
%     idx = find(contains(modelo_test.mets, mets_check{j}));
%     if ~isempty(idx)
%         fprintf('%s -> %s\n', mets_check{j}, strjoin(modelo_test.mets(idx), ', '));
%     else
%         fprintf('%s -> NO ENCONTRADO\n', mets_check{j});
%     end
% end
% fprintf('==========================================\n\n');

%% ADICION MANUAL DE REACCIONES DE TRANSPORTE DE FRUCTOSA
% =========================================================================
% Basado en comparacion con modelos BiGG de referencia:
%   iML1515 (E. coli, gram-negativa): e -> p -> c (3 compartimentos)
%   iYS854  (S. aureus, gram-positiva): e -> c (2 compartimentos)
%
% Todos los modelos CarveMe tienen 3 compartimentos (C_c, C_p, C_e),
% incluso las gram-positivas. Por lo tanto, todas siguen el patron:
%   R_EX_fru_e:  M_fru_e <=>              (exchange)
%   R_FRUtex:    M_fru_e <=> M_fru_p      (difusion outer membrane)
%   R_FRUptspp:  M_fru_p + M_pep_c -> M_f1p_c + M_pyr_c  (PTS periplasma)
%   R_FRUK:      M_atp_c + M_f1p_c -> M_adp_c + M_fdp_c + M_h_c  (kinasa)
%
% Estado de fructosa en _v3_clean_cobrapy.xml:
%                      Bt(G-)  PT33(G+) HGF2(G+) M62_1(G+) Csym(G+) M38(G+)
% M_fru_c               SI      SI       SI        SI        SI       SI
% M_fru_e              *NO*    *NO*      SI        SI        SI       SI
% M_fru_p              *NO*     SI       SI        SI        SI       SI
% R_EX_fru_e           *NO*    *NO*      SI        SI        SI       SI
% R_FRUtex             *NO*    *NO*      SI        SI        SI       SI
% R_FRUptspp           *NO*     SI       SI        --        SI       SI
% R_FRUK               *NO*     --       SI        --        SI       SI
%
% Correccion:
% - Bt: se agrega todo (EX, FRUtex, FRUptspp, FRUK)
% - PT33: se ELIMINA FRUptspp (Bifidobacterium no tiene PTS para fructosa)
%         se agrega EX, FRUtex, y FRUt2r (symport H+, basado en AGORA2 BB-12)
%         Ruta: fru[e] -> FRUt2r -> fru[c] -> HEX7 -> f6p[c] (bifid shunt)
% - HGF2, M62_1, Csym, M38: ya tienen pathway completo.
% =========================================================================

fprintf('\n=== ADICION DE TRANSPORTE DE FRUCTOSA ===\n');

% --- Bacteroides thetaiotaomicron (gram-negativa) ---
% Falta todo: M_fru_e, M_fru_p, exchange, transporte y kinasa
m = Bacteroides_thetaiotaomicron_VPI_5482;
fprintf('\nBacteroides_thetaiotaomicron_VPI_5482:\n');

% Metabolitos nuevos
m = addMetabolite(m, 'M_fru_e', 'metName', 'D-Fructose', ...
    'metFormula', 'C6H12O6');
m = addMetabolite(m, 'M_fru_p', 'metName', 'D-Fructose', ...
    'metFormula', 'C6H12O6');
fprintf('  Agregados: M_fru_e, M_fru_p\n');

% Exchange: M_fru_e <=>
m = addReaction(m, 'R_EX_fru_e', ...
    'reactionName', 'D-Fructose exchange', ...
    'metaboliteList', {'M_fru_e'}, ...
    'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

% Difusion outer membrane: M_fru_e <=> M_fru_p
m = addReaction(m, 'R_FRUtex', ...
    'reactionName', 'D-fructose transport via diffusion (extracellular to periplasm)', ...
    'metaboliteList', {'M_fru_e', 'M_fru_p'}, ...
    'stoichCoeffList', [-1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Transport, outer membrane');

% PTS periplasma: M_fru_p + M_pep_c -> M_f1p_c + M_pyr_c
m = addReaction(m, 'R_FRUptspp', ...
    'reactionName', 'D-fructose transport via PEP:Pyr PTS (periplasm)', ...
    'metaboliteList', {'M_fru_p', 'M_pep_c', 'M_f1p_c', 'M_pyr_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Transport, inner membrane');

% Fructose-1-phosphate kinase: M_atp_c + M_f1p_c -> M_adp_c + M_fdp_c + M_h_c
m = addReaction(m, 'R_FRUK', ...
    'reactionName', 'Fructose-1-phosphate kinase', ...
    'metaboliteList', {'M_atp_c', 'M_f1p_c', 'M_adp_c', 'M_fdp_c', 'M_h_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1, 1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Fructose and mannose metabolism');

fprintf('  Agregadas: R_EX_fru_e, R_FRUtex, R_FRUptspp, R_FRUK\n');
Bacteroides_thetaiotaomicron_VPI_5482 = m;

% --- Bifidobacterium animalis lactis PT33 (gram-positiva, SIN periplasma) ---
% Basado en AGORA2 BB-12: Bifidobacterium NO tiene PTS para fructosa.
% Ruta correcta: FRUt2r (symport H+, e->c directo) + HEX7 (fru->f6p).
% Gram-positiva: transporte directo e<->c, sin tex (no hay periplasma).
m = Bifidobacterium_animalis_lactis_PT33;
fprintf('\nBifidobacterium_animalis_lactis_PT33 (gram+, sin periplasma):\n');

% Metabolito nuevo: M_fru_e
if ~any(strcmp(m.mets, 'M_fru_e'))
    m = addMetabolite(m, 'M_fru_e', 'metName', 'D-Fructose', ...
        'metFormula', 'C6H12O6');
    fprintf('  Agregado: M_fru_e\n');
end

% Exchange: M_fru_e <=>
m = addReaction(m, 'R_EX_fru_e', ...
    'reactionName', 'D-Fructose exchange', ...
    'metaboliteList', {'M_fru_e'}, ...
    'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

% Symport H+ (AGORA2 BB-12): M_fru_e + M_h_e <=> M_fru_c + M_h_c
% Gram+: transporte directo e<->c (sin periplasma)
% Ruta: fru[e] -> FRUt2r -> fru[c] -> HEX7 -> f6p[c] -> bifid shunt
m = addReaction(m, 'R_FRUt2r', ...
    'reactionName', 'D-fructose reversible transport via proton symport', ...
    'metaboliteList', {'M_fru_e', 'M_h_e', 'M_fru_c', 'M_h_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Transport, inner membrane');

fprintf('  Agregadas: R_EX_fru_e, R_FRUt2r (symport AGORA2, e<->c directo)\n');
fprintf('  Ruta: fru[e] -> FRUt2r -> fru[c] -> HEX7 -> f6p[c]\n');
Bifidobacterium_animalis_lactis_PT33 = m;

% --- Verificacion de fructosa en todos los modelos ---
fprintf('\n--- Verificacion de fructosa en todos los modelos ---\n');
modelos_ver = {Bacteroides_thetaiotaomicron_VPI_5482, Bifidobacterium_animalis_lactis_PT33, ...
               Clostridium_sp_HGF2, Clostridium_sp_M62_1, ...
               Clostridium_symbiosum_WAL_14673, Lacticaseibacillus_paracasei_M38};
nombres_ver = {'Bt', 'PT33', 'HGF2', 'M62_1', 'Csym', 'M38'};
% Bt (gram-neg) usa FRUtex+FRUptspp; PT33 (gram+) usa FRUt2r directo
rxns_fru = {'R_EX_fru_e', 'R_FRUtex', 'R_FRUptspp', 'R_FRUK', 'R_FRUt2r'};
mets_fru = {'M_fru_c', 'M_fru_e'};

for i = 1:6
    mv = modelos_ver{i};
    rxns_ok = sum(ismember(rxns_fru, mv.rxns));
    mets_ok = sum(ismember(mets_fru, mv.mets));
    fprintf('  %6s: mets=%d/%d  rxns=%d/%d', nombres_ver{i}, mets_ok, length(mets_fru), rxns_ok, length(rxns_fru));
    % Listar lo que falta
    falta_rxns = rxns_fru(~ismember(rxns_fru, mv.rxns));
    falta_mets = mets_fru(~ismember(mets_fru, mv.mets));
    if isempty(falta_rxns) && isempty(falta_mets)
        fprintf('  OK\n');
    else
        if ~isempty(falta_mets)
            fprintf('  falta_mets: %s', strjoin(falta_mets, ', '));
        end
        if ~isempty(falta_rxns)
            fprintf('  falta_rxns: %s', strjoin(falta_rxns, ', '));
        end
        fprintf('\n');
    end
end
fprintf('==========================================\n');

%% ADICION MANUAL DE EXPORTACION DE BUTIRATO
% =========================================================================
% Los Clostridia (principales productores de butirato intestinal) tienen
% la maquinaria de sintesis (PBUTT + BUTKr -> M_but_c) pero NO pueden
% exportar butirato (faltan M_but_e, M_but_p, exchange y transporte).
% Bacteroides tiene R_BUTt4pp pero es import-only (irreversible).
%
% Referencia: iML1515 (E. coli, BiGG)
%   R_EX_but_e:  M_but_e <=>                               (exchange)
%   R_BUTtex:    M_but_e <=> M_but_p                       (difusion outer membrane)
%   R_BUTt2rpp:  M_but_p + M_h_p <=> M_but_c + M_h_c      (simport H+, reversible)
%
% Estado en _v3_clean_cobrapy.xml:
%                    Bt     HGF2   M62_1  C.sym  M38    PT33
% M_but_c            SI     SI     *NO*   SI     SI     NO
% R_PBUTT+R_BUTKr    SI     SI     *NO*   *NO*   SI     NO
% M_but_e/p          SI     *NO*   *NO*   *NO*   *NO*   NO
% R_EX_but_e         SI     *NO*   *NO*   *NO*   *NO*   NO
% R_BUTtex           SI     *NO*   *NO*   *NO*   *NO*   NO
% R_BUTt2rpp (rev)   *NO*   *NO*   *NO*   *NO*   *NO*   NO
%
% Notas biologicas:
%   - HGF2, M62_1, C.symbiosum: Clostridia, productores clasicos de butirato
%   - M38: Lacticaseibacillus, produce principalmente lactato
%           (CarveMe encontro genes de butirato, se agrega exportacion)
%   - PT33: Bifidobacterium, NO produce butirato -> sin cambios
% =========================================================================

fprintf('\n=== ADICION DE EXPORTACION DE BUTIRATO ===\n');

% --- Bacteroides: ya tiene M_but_e/p y R_EX_but_e, pero R_BUTt4pp es ---
% --- import-only. Agregar R_BUTt2rpp (reversible) para permitir export ---
m = Bacteroides_thetaiotaomicron_VPI_5482;
fprintf('\nBacteroides_thetaiotaomicron_VPI_5482:\n');

m = addReaction(m, 'R_BUTt2rpp', ...
    'reactionName', 'Butyrate transport via proton symport, reversible (periplasm)', ...
    'metaboliteList', {'M_but_p', 'M_h_p', 'M_but_c', 'M_h_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Transport, inner membrane');
fprintf('  Agregada: R_BUTt2rpp (transporte reversible c<=>p)\n');
Bacteroides_thetaiotaomicron_VPI_5482 = m;

% --- HGF2 (gram-positiva, SIN periplasma): falta exportacion butirato ---
m = Clostridium_sp_HGF2;
fprintf('\nClostridium_sp_HGF2 (gram+, sin periplasma):\n');

m = addMetabolite(m, 'M_but_e', 'metName', 'Butyrate (n-C4:0)', ...
    'metFormula', 'C4H7O2', 'Charge', -1);
fprintf('  Agregado: M_but_e\n');

m = addReaction(m, 'R_EX_but_e', ...
    'reactionName', 'Butyrate exchange', ...
    'metaboliteList', {'M_but_e'}, ...
    'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

% Gram+: transporte directo e<->c (sin tex ni periplasma)
m = addReaction(m, 'R_BUTt2r', ...
    'reactionName', 'Butyrate transport via proton symport, reversible', ...
    'metaboliteList', {'M_but_e', 'M_h_e', 'M_but_c', 'M_h_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Transport, inner membrane');

fprintf('  Agregadas: R_EX_but_e, R_BUTt2r (e<->c directo)\n');
Clostridium_sp_HGF2 = m;

% --- M62_1: NO tiene M_but_c ni R_PBUTT/R_BUTKr. Falta sintesis + export ---
m = Clostridium_sp_M62_1;
fprintf('\nClostridium_sp_M62_1:\n');

% Metabolitos de sintesis (tiene M_btcoa_c pero no M_but_c ni M_butpi_c)
m = addMetabolite(m, 'M_butpi_c', 'metName', 'Butanoyl phosphate', ...
    'metFormula', 'C4H7O5P');
m = addMetabolite(m, 'M_but_c', 'metName', 'Butyrate (n-C4:0)', ...
    'metFormula', 'C4H7O2');
fprintf('  Agregados: M_butpi_c, M_but_c\n');

% Sintesis: btcoa -> butpi -> but (ruta clasica Clostridial)
% R_PBUTT: btcoa_c + pi_c <=> butpi_c + coa_c (phosphotransbutyrylase)
m = addReaction(m, 'R_PBUTT', ...
    'reactionName', 'Phosphate butyryltransferase', ...
    'metaboliteList', {'M_btcoa_c', 'M_pi_c', 'M_butpi_c', 'M_coa_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Pyruvate metabolism');

% R_BUTKr: atp_c + but_c <=> adp_c + butpi_c (butyrate kinase)
m = addReaction(m, 'R_BUTKr', ...
    'reactionName', 'Butyrate kinase', ...
    'metaboliteList', {'M_atp_c', 'M_but_c', 'M_adp_c', 'M_butpi_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Pyruvate metabolism');

fprintf('  Agregadas: R_PBUTT, R_BUTKr (sintesis butirato)\n');

% Exportacion (gram+: sin periplasma, transporte directo e<->c)
m = addMetabolite(m, 'M_but_e', 'metName', 'Butyrate (n-C4:0)', ...
    'metFormula', 'C4H7O2', 'Charge', -1);

m = addReaction(m, 'R_EX_but_e', ...
    'reactionName', 'Butyrate exchange', ...
    'metaboliteList', {'M_but_e'}, ...
    'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

m = addReaction(m, 'R_BUTt2r', ...
    'reactionName', 'Butyrate transport via proton symport, reversible', ...
    'metaboliteList', {'M_but_e', 'M_h_e', 'M_but_c', 'M_h_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Transport, inner membrane');

fprintf('  Agregados: M_but_e + R_EX_but_e, R_BUTt2r (e<->c directo)\n');
Clostridium_sp_M62_1 = m;

% --- C.symbiosum (gram-positiva, SIN periplasma): falta export butirato ---
m = Clostridium_symbiosum_WAL_14673;
fprintf('\nClostridium_symbiosum_WAL_14673 (gram+, sin periplasma):\n');

m = addMetabolite(m, 'M_but_e', 'metName', 'Butyrate (n-C4:0)', ...
    'metFormula', 'C4H7O2', 'Charge', -1);
fprintf('  Agregado: M_but_e\n');

m = addReaction(m, 'R_EX_but_e', ...
    'reactionName', 'Butyrate exchange', ...
    'metaboliteList', {'M_but_e'}, ...
    'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

m = addReaction(m, 'R_BUTt2r', ...
    'reactionName', 'Butyrate transport via proton symport, reversible', ...
    'metaboliteList', {'M_but_e', 'M_h_e', 'M_but_c', 'M_h_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Transport, inner membrane');

fprintf('  Agregadas: R_EX_but_e, R_BUTt2r (e<->c directo)\n');
Clostridium_symbiosum_WAL_14673 = m;

% --- M38 (gram-positiva, SIN periplasma): falta exportacion butirato ---
m = Lacticaseibacillus_paracasei_M38;
fprintf('\nLacticaseibacillus_paracasei_M38 (gram+, sin periplasma):\n');

m = addMetabolite(m, 'M_but_e', 'metName', 'Butyrate (n-C4:0)', ...
    'metFormula', 'C4H7O2', 'Charge', -1);
fprintf('  Agregado: M_but_e\n');

m = addReaction(m, 'R_EX_but_e', ...
    'reactionName', 'Butyrate exchange', ...
    'metaboliteList', {'M_but_e'}, ...
    'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

m = addReaction(m, 'R_BUTt2r', ...
    'reactionName', 'Butyrate transport via proton symport, reversible', ...
    'metaboliteList', {'M_but_e', 'M_h_e', 'M_but_c', 'M_h_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Transport, inner membrane');

fprintf('  Agregadas: R_EX_but_e, R_BUTt2r (e<->c directo)\n');
Lacticaseibacillus_paracasei_M38 = m;

% PT33 (Bifidobacterium): NO produce butirato -> sin cambios
fprintf('\nBifidobacterium_animalis_lactis_PT33: sin cambios (no produce butirato)\n');

% --- Verificacion de butirato en todos los modelos ---
fprintf('\n--- Verificacion de butirato en todos los modelos ---\n');
modelos_but = {Bacteroides_thetaiotaomicron_VPI_5482, Bifidobacterium_animalis_lactis_PT33, ...
               Clostridium_sp_HGF2, Clostridium_sp_M62_1, ...
               Clostridium_symbiosum_WAL_14673, Lacticaseibacillus_paracasei_M38};
nombres_but = {'Bt', 'PT33', 'HGF2', 'M62_1', 'Csym', 'M38'};
% Bt (gram-neg) usa BUTtex+BUTt2rpp; gram+ usa BUTt2r (directo e<->c)
rxns_but = {'R_EX_but_e', 'R_BUTtex', 'R_BUTt2rpp', 'R_BUTt2r'};
mets_but = {'M_but_c', 'M_but_e'};

for i = 1:6
    mv = modelos_but{i};
    rxns_ok = sum(ismember(rxns_but, mv.rxns));
    mets_ok = sum(ismember(mets_but, mv.mets));
    fprintf('  %6s: mets=%d/%d  rxns=%d/%d', nombres_but{i}, mets_ok, length(mets_but), rxns_ok, length(rxns_but));
    if mets_ok == length(mets_but) && rxns_ok == length(rxns_but)
        fprintf('  OK (puede exportar butirato)\n');
    elseif ~any(ismember({'M_but_c'}, mv.mets))
        fprintf('  -- (no produce butirato)\n');
    else
        falta_r = rxns_but(~ismember(rxns_but, mv.rxns));
        falta_m = mets_but(~ismember(mets_but, mv.mets));
        if ~isempty(falta_m); fprintf('  falta: %s', strjoin(falta_m, ', ')); end
        if ~isempty(falta_r); fprintf('  falta: %s', strjoin(falta_r, ', ')); end
        fprintf('\n');
    end
end
fprintf('==========================================\n');

%% ========================================================================
%% REACCIONES DE INULINA Y FOS (basado en AGORA2)
%% ========================================================================
% Referencia: modelos AGORA2 de VMH (https://vmh.life)
%
% Resumen AGORA2:
%   Bt:          solo FOS (en periplasma [p]) - NO degrada inulina
%   L.paracasei: solo FOS (en [e], gram-positiva)
%   C.innocuum:  solo inulina intracelular
%   HGF2/Csym/M62/PT33: nada en AGORA2
%
% Modificaciones respecto a AGORA2:
%   - HGF2: inulina SOLO extracelular (INULINASEe) + EX_kestopt
%   - M38: inulina SOLO extracelular (INULINASEe) + todos los FOS
%   - Bt, PT33, Csym: solo degradacion de FOS (kesto, kestopt, kestottr)
%   - M62_1: sin reacciones de inulina/FOS
%   - NADIE tiene INULINabc ni INULINASE (inulina demasiado grande para transporte intracelular)
%
% Mecanismo de cross-feeding:
%   1. Degradador (HGF2/M38) degrada inulina extracelularmente (INULINASEe)
%      inulin_e + 25 h2o_e -> 25 fru_e + kestopt_e
%   2. fru_e disponible en medio para consumo directo por otras bacterias
%   3. kestopt_e disponible en medio via EX_kestopt(e)
%   4. Consumidores de FOS (PT33/Csym/M38/Bt) degradan kestopt/kesto/kestottr
%      a fru + glc en [e] (extracelular), luego entran por transporte normal
%   NOTA: Toda degradacion de FOS/inulina es EXTRACELULAR (enzimas secretadas),
%         por eso todos los metabolitos estan en _e (no en _p)
%
% Compartimentos en CarveMe:
%   _e = extracelular (usado por exchanges y reacciones extracelulares)
%   _p = periplasma (solo Bt, gram-neg: transporte e<->p<->c)
%   _c = citoplasma
%   Todas las reacciones de degradacion FOS/inulina usan _e (extracelular)
%% ========================================================================

%% Bt: solo FOS (NO degrada inulina, solo consume FOS liberados por degradadores)

fprintf('\n--- Bt: agregando solo FOS (NO inulina) ---\n');
m = Bacteroides_thetaiotaomicron_VPI_5482;

% Metabolitos nuevos en _e (solo FOS, NO inulina)
% La degradacion de FOS es EXTRACELULAR (enzimas secretadas), por eso _e
m = addMetabolite(m, 'M_kesto_e',    'metName', 'Kestose',       'metFormula', 'C18H32O16');
m = addMetabolite(m, 'M_kestopt_e',  'metName', 'Kestopentaose', 'metFormula', 'C30H52O26');
m = addMetabolite(m, 'M_kestottr_e', 'metName', 'Kestotetraose', 'metFormula', 'C24H42O21');

% FOS exchanges (metabolitos en _e, bidireccionales)
m = addReaction(m, 'EX_kesto(e)', 'reactionName', 'Kestose exchange', ...
    'metaboliteList', {'M_kesto_e'}, 'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');
m = addReaction(m, 'EX_kestopt(e)', 'reactionName', 'Kestopentaose exchange', ...
    'metaboliteList', {'M_kestopt_e'}, 'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');
m = addReaction(m, 'EX_kestottr(e)', 'reactionName', 'Kestotetraose exchange', ...
    'metaboliteList', {'M_kestottr_e'}, 'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

% FOS degradacion EXTRACELULAR (enzimas secretadas, todo en _e)
% Productos fru_e y glc__D_e entran a la celula por transportadores
% existentes: FRUtex (fru_e<->fru_p) + FRUptspp/FRUt3 (fru_p->c)
%             GLCtex (glc_e<->glc_p) + transportador glc (p->c)
m = addReaction(m, 'KESTOASEe', 'reactionName', 'Beta-2-1 fructan degradation (kestose)', ...
    'metaboliteList', {'M_h2o_e','M_kesto_e','M_fru_e','M_glc__D_e'}, ...
    'stoichCoeffList', [-2,-1,2,1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Plant polysaccharide degradation');
m = addReaction(m, 'KESTOPTASEe', 'reactionName', 'Beta-2-1 fructan degradation (kestopentaose)', ...
    'metaboliteList', {'M_h2o_e','M_kestopt_e','M_fru_e','M_glc__D_e'}, ...
    'stoichCoeffList', [-4,-1,4,1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Plant polysaccharide degradation');
m = addReaction(m, 'KESTOTTRASEe', 'reactionName', 'Beta-2-1 fructan degradation (kestotetraose)', ...
    'metaboliteList', {'M_h2o_e','M_kestottr_e','M_fru_e','M_glc__D_e'}, ...
    'stoichCoeffList', [-3,-1,3,1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Plant polysaccharide degradation');

Bacteroides_thetaiotaomicron_VPI_5482 = m;
fprintf('  Bt: 6 rxns FOS agregadas (sin inulina, todo en _e)\n');

%% HGF2 (gram+, sin periplasma): inulina SOLO extracelular
% INULINASEe: inulin_e -> 25 fru_e + kestopt_e (cross-feeding!)
% Gram+: todo en _e (sin periplasma)

fprintf('\n--- HGF2 (gram+): agregando inulina SOLO extracelular ---\n');
m = Clostridium_sp_HGF2;

% Metabolitos nuevos en _e (gram+ sin periplasma)
m = addMetabolite(m, 'M_inulin_e',  'metName', 'Inulin',        'metFormula', 'C180H302O151');
m = addMetabolite(m, 'M_kestopt_e', 'metName', 'Kestopentaose', 'metFormula', 'C30H52O26');

if ~any(strcmp(m.mets, 'M_fru_e'))
    m = addMetabolite(m, 'M_fru_e', 'metName', 'D-Fructose', 'metFormula', 'C6H12O6');
    fprintf('  NOTA: M_fru_e agregado (no existia)\n');
end

% Inulina exchange
m = addReaction(m, 'EX_inulin(e)', 'reactionName', 'Inulin exchange', ...
    'metaboliteList', {'M_inulin_e'}, 'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

% Inulina degradacion EXTRACELULAR (clave para cross-feeding!)
% inulin_e + 25 h2o_e -> 25 fru_e + kestopt_e
m = addReaction(m, 'INULINASEe', 'reactionName', 'Inulin degradation extracellular', ...
    'metaboliteList', {'M_h2o_e','M_inulin_e','M_fru_e','M_kestopt_e'}, ...
    'stoichCoeffList', [-25,-1,25,1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Plant polysaccharide degradation');

% Kestopentaose exchange (producto de INULINASEe, sale al medio)
m = addReaction(m, 'EX_kestopt(e)', 'reactionName', 'Kestopentaose exchange', ...
    'metaboliteList', {'M_kestopt_e'}, 'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

Clostridium_sp_HGF2 = m;
fprintf('  HGF2: %d rxns inulina agregadas (extracelular, _e)\n', 3);

%% M38 (gram+, sin periplasma): inulina SOLO extracelular + todos los FOS
% Gram+: todo en _e (sin periplasma)

fprintf('\n--- M38 (gram+): agregando inulina SOLO extracelular + FOS ---\n');
m = Lacticaseibacillus_paracasei_M38;

% Metabolitos nuevos en _e (gram+ sin periplasma)
m = addMetabolite(m, 'M_inulin_e',   'metName', 'Inulin',        'metFormula', 'C180H302O151');
m = addMetabolite(m, 'M_kesto_e',    'metName', 'Kestose',       'metFormula', 'C18H32O16');
m = addMetabolite(m, 'M_kestopt_e',  'metName', 'Kestopentaose', 'metFormula', 'C30H52O26');
m = addMetabolite(m, 'M_kestottr_e', 'metName', 'Kestotetraose', 'metFormula', 'C24H42O21');

if ~any(strcmp(m.mets, 'M_fru_e'))
    m = addMetabolite(m, 'M_fru_e', 'metName', 'D-Fructose', 'metFormula', 'C6H12O6');
    fprintf('  NOTA: M_fru_e agregado (no existia)\n');
end
if ~any(strcmp(m.mets, 'M_glc__D_e'))
    m = addMetabolite(m, 'M_glc__D_e', 'metName', 'D-Glucose', 'metFormula', 'C6H12O6');
    fprintf('  NOTA: M_glc__D_e agregado (no existia)\n');
end

% Inulina exchange
m = addReaction(m, 'EX_inulin(e)', 'reactionName', 'Inulin exchange', ...
    'metaboliteList', {'M_inulin_e'}, 'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

% Inulina degradacion EXTRACELULAR
% inulin_e + 25 h2o_e -> 25 fru_e + kestopt_e
m = addReaction(m, 'INULINASEe', 'reactionName', 'Inulin degradation extracellular', ...
    'metaboliteList', {'M_h2o_e','M_inulin_e','M_fru_e','M_kestopt_e'}, ...
    'stoichCoeffList', [-25,-1,25,1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Plant polysaccharide degradation');

% FOS exchanges
m = addReaction(m, 'EX_kesto(e)', 'reactionName', 'Kestose exchange', ...
    'metaboliteList', {'M_kesto_e'}, 'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');
m = addReaction(m, 'EX_kestopt(e)', 'reactionName', 'Kestopentaose exchange', ...
    'metaboliteList', {'M_kestopt_e'}, 'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');
m = addReaction(m, 'EX_kestottr(e)', 'reactionName', 'Kestotetraose exchange', ...
    'metaboliteList', {'M_kestottr_e'}, 'stoichCoeffList', [-1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Exchange/demand reaction');

% FOS degradacion extracelular (gram+: todo en _e)
m = addReaction(m, 'KESTOASEe', 'reactionName', 'Beta-2-1 fructan degradation (kestose)', ...
    'metaboliteList', {'M_h2o_e','M_kesto_e','M_fru_e','M_glc__D_e'}, ...
    'stoichCoeffList', [-2,-1,2,1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Plant polysaccharide degradation');
m = addReaction(m, 'KESTOPTASEe', 'reactionName', 'Beta-2-1 fructan degradation (kestopentaose)', ...
    'metaboliteList', {'M_h2o_e','M_kestopt_e','M_fru_e','M_glc__D_e'}, ...
    'stoichCoeffList', [-4,-1,4,1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Plant polysaccharide degradation');
m = addReaction(m, 'KESTOTTRASEe', 'reactionName', 'Beta-2-1 fructan degradation (kestotetraose)', ...
    'metaboliteList', {'M_h2o_e','M_kestottr_e','M_fru_e','M_glc__D_e'}, ...
    'stoichCoeffList', [-3,-1,3,1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Plant polysaccharide degradation');

% Gram+: NO necesita GLCtex (glc_D ya esta en _e, no en _p)

Lacticaseibacillus_paracasei_M38 = m;
fprintf('  M38: %d rxns inulina/FOS agregadas (extracelular, _e)\n', 8);

%% M38: reacciones biosintéticas faltantes (AGORA2 LC2W las tiene, CarveMe no)
% =========================================================================
% CarveMe no incluyo 4 reacciones clave que AGORA2 de L. paracasei LC2W si tiene.
% Sin estas reacciones, M38 no puede sintetizar isoleucina, leucina ni generar
% alfa-cetoglutarato via TCA parcial, y su crecimiento depende 100% de amino acids.
% Ref: AGORA2 Lactobacillus_paracasei_LC2W (VMH)
% =========================================================================
fprintf('\n--- M38: reacciones biosinteticas faltantes (AGORA2 LC2W) ---\n');
m = Lacticaseibacillus_paracasei_M38;

% R_ILETA: Isoleucine Transaminase (EC 2.6.1.42)
% akg_c + ile_L_c <=> 3mop_c + glu_L_c
if ~any(strcmp(m.rxns, 'R_ILETA'))
    m = addReaction(m, 'R_ILETA', ...
        'reactionName', 'Isoleucine Transaminase', ...
        'metaboliteList', {'M_akg_c', 'M_ile__L_c', 'M_3mop_c', 'M_glu__L_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Valine, leucine, and isoleucine metabolism');
    fprintf('  Agregada: R_ILETA (isoleucina transaminasa)\n');
else
    fprintf('  [OK] R_ILETA ya existe\n');
end

% R_LEUTA: Leucine Transaminase (EC 2.6.1.42, 2.6.1.6)
% akg_c + leu_L_c <=> 4mop_c + glu_L_c
if ~any(strcmp(m.rxns, 'R_LEUTA'))
    m = addReaction(m, 'R_LEUTA', ...
        'reactionName', 'Leucine Transaminase', ...
        'metaboliteList', {'M_akg_c', 'M_leu__L_c', 'M_4mop_c', 'M_glu__L_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Valine, leucine, and isoleucine metabolism');
    fprintf('  Agregada: R_LEUTA (leucina transaminasa)\n');
else
    fprintf('  [OK] R_LEUTA ya existe\n');
end

% R_VALTA: Valine Transaminase (EC 2.6.1.42)
% akg_c + val_L_c <=> 3mob_c + glu_L_c (BiGG: VALTA, MNXR96230)
if ~any(strcmp(m.rxns, 'R_VALTA'))
    m = addReaction(m, 'R_VALTA', ...
        'reactionName', 'Valine Transaminase', ...
        'metaboliteList', {'M_akg_c', 'M_val__L_c', 'M_3mob_c', 'M_glu__L_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Valine, leucine, and isoleucine metabolism');
    fprintf('  Agregada: R_VALTA (valina transaminasa)\n');
else
    fprintf('  [OK] R_VALTA ya existe\n');
end

% R_CS: Citrate Synthase (EC 2.3.3.1)
% accoa_c + h2o_c + oaa_c -> cit_c + coa_c + h_c (irreversible)
if ~any(strcmp(m.rxns, 'R_CS'))
    m = addReaction(m, 'R_CS', ...
        'reactionName', 'citrate synthase', ...
        'metaboliteList', {'M_accoa_c', 'M_h2o_c', 'M_oaa_c', 'M_cit_c', 'M_coa_c', 'M_h_c'}, ...
        'stoichCoeffList', [-1, -1, -1, 1, 1, 1], ...
        'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_CS (citrato sintasa)\n');
else
    fprintf('  [OK] R_CS ya existe\n');
end

% R_ACONT: Aconitase (EC 4.2.1.3)
% cit_c <=> icit_c (reversible)
if ~any(strcmp(m.rxns, 'R_ACONT'))
    m = addReaction(m, 'R_ACONT', ...
        'reactionName', 'Aconitase', ...
        'metaboliteList', {'M_cit_c', 'M_icit_c'}, ...
        'stoichCoeffList', [-1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_ACONT (aconitasa)\n');
else
    fprintf('  [OK] R_ACONT ya existe\n');
end

% R_PPCK: PEP carboxykinase (EC 4.1.1.32)
% atp_c + oaa_c <=> adp_c + co2_c + pep_c (BiGG: PPCK, MNXR103099)
if ~any(strcmp(m.rxns, 'R_PPCK'))
    m = addReaction(m, 'R_PPCK', ...
        'reactionName', 'Phosphoenolpyruvate carboxykinase (ATP)', ...
        'metaboliteList', {'M_atp_c', 'M_oaa_c', 'M_adp_c', 'M_co2_c', 'M_pep_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Pyruvate metabolism');
    fprintf('  Agregada: R_PPCK (PEP carboxikinasa)\n');
else
    fprintf('  [OK] R_PPCK ya existe\n');
end

% R_SUCD1: Succinate dehydrogenase (EC 1.3.5.2)
% fad_c + succ_c <=> fadh2_c + fum_c (BiGG: SUCD1, MNXR82874)
if ~any(strcmp(m.rxns, 'R_SUCD1'))
    m = addReaction(m, 'R_SUCD1', ...
        'reactionName', 'Succinate dehydrogenase', ...
        'metaboliteList', {'M_fad_c', 'M_succ_c', 'M_fadh2_c', 'M_fum_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_SUCD1 (succinato deshidrogenasa)\n');
else
    fprintf('  [OK] R_SUCD1 ya existe\n');
end

% R_ACS: Acetyl-CoA synthetase (EC 6.2.1.1)
% ac_c + atp_c + coa_c <=> accoa_c + amp_c + ppi_c (BiGG: ACS, MNXR190571)
if ~any(strcmp(m.rxns, 'R_ACS'))
    m = addReaction(m, 'R_ACS', ...
        'reactionName', 'Acetyl-CoA synthetase', ...
        'metaboliteList', {'M_ac_c', 'M_atp_c', 'M_coa_c', 'M_accoa_c', 'M_amp_c', 'M_ppi_c'}, ...
        'stoichCoeffList', [-1, -1, -1, 1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Pyruvate metabolism');
    fprintf('  Agregada: R_ACS (acetil-CoA sintetasa)\n');
else
    fprintf('  [OK] R_ACS ya existe\n');
end

Lacticaseibacillus_paracasei_M38 = m;
fprintf('  M38: BCAA (ILETA+LEUTA+VALTA) + TCA (CS+ACONT) + PPCK + SUCD1 + ACS completadas\n');

%% PT33 y Csym (gram+, sin periplasma): solo degradacion de FOS (cross-feeding)

fprintf('\n--- PT33, Csym (gram+): agregando degradacion FOS en _e ---\n');

modelos_fos = {Bifidobacterium_animalis_lactis_PT33, ...
               Clostridium_symbiosum_WAL_14673};
nombres_fos = {'Bifidobacterium_animalis_lactis_PT33', ...
               'Clostridium_symbiosum_WAL_14673'};
short_fos = {'PT33', 'Csym'};

for i = 1:2
    m = modelos_fos{i};

    % Metabolitos FOS en _e (gram+ sin periplasma)
    m = addMetabolite(m, 'M_kesto_e',    'metName', 'Kestose',       'metFormula', 'C18H32O16');
    m = addMetabolite(m, 'M_kestopt_e',  'metName', 'Kestopentaose', 'metFormula', 'C30H52O26');
    m = addMetabolite(m, 'M_kestottr_e', 'metName', 'Kestotetraose', 'metFormula', 'C24H42O21');

    if ~any(strcmp(m.mets, 'M_fru_e'))
        m = addMetabolite(m, 'M_fru_e', 'metName', 'D-Fructose', 'metFormula', 'C6H12O6');
        fprintf('  NOTA: M_fru_e agregado a %s\n', short_fos{i});
    end
    if ~any(strcmp(m.mets, 'M_glc__D_e'))
        m = addMetabolite(m, 'M_glc__D_e', 'metName', 'D-Glucose', 'metFormula', 'C6H12O6');
        fprintf('  NOTA: M_glc__D_e agregado a %s\n', short_fos{i});
    end

    % FOS exchanges (usan _e)
    m = addReaction(m, 'EX_kesto(e)', 'reactionName', 'Kestose exchange', ...
        'metaboliteList', {'M_kesto_e'}, 'stoichCoeffList', [-1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Exchange/demand reaction');
    m = addReaction(m, 'EX_kestopt(e)', 'reactionName', 'Kestopentaose exchange', ...
        'metaboliteList', {'M_kestopt_e'}, 'stoichCoeffList', [-1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Exchange/demand reaction');
    m = addReaction(m, 'EX_kestottr(e)', 'reactionName', 'Kestotetraose exchange', ...
        'metaboliteList', {'M_kestottr_e'}, 'stoichCoeffList', [-1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Exchange/demand reaction');

    % FOS degradacion extracelular (gram+: todo en _e)
    m = addReaction(m, 'KESTOASEe', 'reactionName', 'Beta-2-1 fructan degradation (kestose)', ...
        'metaboliteList', {'M_h2o_e','M_kesto_e','M_fru_e','M_glc__D_e'}, ...
        'stoichCoeffList', [-2,-1,2,1], ...
        'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
        'subSystem', 'Plant polysaccharide degradation');
    m = addReaction(m, 'KESTOPTASEe', 'reactionName', 'Beta-2-1 fructan degradation (kestopentaose)', ...
        'metaboliteList', {'M_h2o_e','M_kestopt_e','M_fru_e','M_glc__D_e'}, ...
        'stoichCoeffList', [-4,-1,4,1], ...
        'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
        'subSystem', 'Plant polysaccharide degradation');
    m = addReaction(m, 'KESTOTTRASEe', 'reactionName', 'Beta-2-1 fructan degradation (kestotetraose)', ...
        'metaboliteList', {'M_h2o_e','M_kestottr_e','M_fru_e','M_glc__D_e'}, ...
        'stoichCoeffList', [-3,-1,3,1], ...
        'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
        'subSystem', 'Plant polysaccharide degradation');

    % Gram+: NO necesita GLCtex (glc_D ya esta en _e, no en _p)

    modelos_fos{i} = m;
    fprintf('  %s: %d rxns FOS agregadas (extracelular, _e)\n', short_fos{i}, 6);
end

Bifidobacterium_animalis_lactis_PT33 = modelos_fos{1};
Clostridium_symbiosum_WAL_14673      = modelos_fos{2};

% M62_1: sin reacciones de inulina/FOS (tampoco en AGORA2)
fprintf('  M62_1: sin reacciones de inulina/FOS\n');

%% VERIFICACION: reacciones de inulina/FOS agregadas
fprintf('\n=== VERIFICACION DE REACCIONES INULINA/FOS ===\n');
check = {
    'Bt',   Bacteroides_thetaiotaomicron_VPI_5482, ...
        {'EX_kesto(e)','EX_kestopt(e)','EX_kestottr(e)', ...
         'KESTOASEe','KESTOPTASEe','KESTOTTRASEe'};
    'HGF2', Clostridium_sp_HGF2, ...
        {'EX_inulin(e)','INULINASEe','EX_kestopt(e)'};
    'M38',  Lacticaseibacillus_paracasei_M38, ...
        {'EX_inulin(e)','INULINASEe', ...
         'EX_kesto(e)','EX_kestopt(e)','EX_kestottr(e)', ...
         'KESTOASEe','KESTOPTASEe','KESTOTTRASEe'};
    'PT33', Bifidobacterium_animalis_lactis_PT33, ...
        {'EX_kesto(e)','EX_kestopt(e)','EX_kestottr(e)', ...
         'KESTOASEe','KESTOPTASEe','KESTOTTRASEe'};
    'Csym', Clostridium_symbiosum_WAL_14673, ...
        {'EX_kesto(e)','EX_kestopt(e)','EX_kestottr(e)', ...
         'KESTOASEe','KESTOPTASEe','KESTOTTRASEe'};
    'M62',  Clostridium_sp_M62_1, {};
};

for i = 1:size(check,1)
    if isempty(check{i,3})
        fprintf('%s: sin reacciones inulina/FOS (correcto)\n', check{i,1});
    else
        n_ok = sum(ismember(check{i,3}, check{i,2}.rxns));
        n_total = length(check{i,3});
        fprintf('%s: %d/%d rxns OK', check{i,1}, n_ok, n_total);
        if n_ok < n_total
            falta = check{i,3}(~ismember(check{i,3}, check{i,2}.rxns));
            fprintf(' --- FALTA: %s', strjoin(falta, ', '));
        end
        fprintf('\n');
    end
end

%% VERIFICACION: transporte de fructosa (necesario para cross-feeding)
fprintf('\n=== VERIFICACION: transporte de fructosa ===\n');
fprintf('Bt (gram-neg): FRUtex (e<->p) + FRUptspp/FRUt3 (p<->c)\n');
fprintf('Gram+: FRUt2r (e<->c directo) o fru_e ya en INULINASEe\n\n');

todos_modelos = {Bacteroides_thetaiotaomicron_VPI_5482, Bifidobacterium_animalis_lactis_PT33, ...
                 Clostridium_sp_HGF2, Clostridium_sp_M62_1, ...
                 Clostridium_symbiosum_WAL_14673, Lacticaseibacillus_paracasei_M38};
todos_nombres = {'Bt','PT33','HGF2','M62_1','Csym','M38'};

for i = 1:6
    mx = todos_modelos{i};
    tiene_FRUtex = any(strcmp(mx.rxns, 'R_FRUtex'));
    tiene_FRUt2r = any(strcmp(mx.rxns, 'R_FRUt2r'));
    tiene_EX_fru = any(strcmp(mx.rxns, 'R_EX_fru_e'));
    tiene_fru_e = any(strcmp(mx.mets, 'M_fru_e'));
    fprintf('%s: FRUtex=%d, FRUt2r=%d, EX_fru=%d, fru_e=%d\n', ...
        todos_nombres{i}, tiene_FRUtex, tiene_FRUt2r, tiene_EX_fru, tiene_fru_e);
end

%% Definicion medio de cultivo en formato BiGG
mZMB = {
'R_EX_glc__D_e';        % Glucosa (fuente principal de carbono)
'R_EX_ala__L_e';        % L-alanina
'R_EX_arg__L_e';        % L-arginina
'R_EX_asn__L_e';        % L-asparagina
'R_EX_asp__L_e';        % L-aspartato
'R_EX_cys__L_e';        % L-cisteína
'R_EX_gln__L_e';        % L-glutamina
'R_EX_glu__L_e';        % L-glutamato
'R_EX_gly_e';           % Glicina
'R_EX_his__L_e';        % L-histidina
'R_EX_ile__L_e';        % L-isoleucina
'R_EX_leu__L_e';        % L-leucina
'R_EX_lys__L_e';        % L-lisina
'R_EX_met__L_e';        % L-metionina
'R_EX_phe__L_e';        % L-fenilalanina
'R_EX_pro__L_e';        % L-prolina
'R_EX_ser__L_e';        % L-serina
'R_EX_thr__L_e';        % L-treonina
'R_EX_trp__L_e';        % L-triptófano
'R_EX_tyr__L_e';        % L-tirosina
'R_EX_val__L_e';        % L-valina
'R_EX_pi_e';            % Fosfato
'R_EX_k_e';             % Potasio
'R_EX_mg2_e';           % Magnesio
'R_EX_fe3_e';           % Hierro (III)
'R_EX_fe2_e';           % Hierro (II)
'R_EX_na1_e';           % Sodio
'R_EX_zn2_e';           % Zinc (cpd00034)
'R_EX_mn2_e';           % Manganeso
'R_EX_mobd_e';          % Molibdato
'R_EX_nh3_e';           % Amonio (nh3)
'R_EX_ca2_e';           % Calcio
'R_EX_so4_e';           % Sulfato
'R_EX_cobalt2_e';       % Cobalto
'R_EX_cu2_e';           % Cobre
'R_EX_cl_e';            % Cloruro
'R_EX_4ppan_e';         % Pantotenato (4ppan)
'R_EX_fol_e';           % Folato
'R_EX_nac_e';           % Niacina
'R_EX_pydx_e';          % Piridoxina
'R_EX_4abz_e';          % 4-aminobenzoato
'R_EX_xan_e';           % Xantina
'R_EX_inost_e';         % Inositol
'R_EX_gthox_e';         % Glutatión oxidado
'R_EX_gthrd_e';         % Glutatión reducido
'R_EX_btn_e';           % Biotina
'R_EX_ribflv_e';        % Riboflavina
'R_EX_thm_e';           % Tiamina
'R_EX_lipoate_e';       % Lipoato
'R_EX_pheme_e';         % Heme
'R_EX_sheme_e';         % Siroheme
'R_EX_ade_e';           % Adenina
'R_EX_gua_e';           % Guanina
'R_EX_ura_e';           % Uracilo
'R_EX_ddca_e';          % Laurato (cpd01741)
'R_EX_h2s_e';           % Sulfuro de hidrógeno (cpd00239 -> h2s)
'R_EX_hxan_e';          % Hipoxantina
'R_EX_mqn7_e';          % Menaquinona-7
'R_EX_mqn8_e';          % Menaquinona-8
'R_EX_ni2_e';           % Níquel
'R_EX_ocdca_e';         % Octadecanoato
'R_EX_q8_e';            % Ubiquinona-8
'R_EX_thymd_e';         % Timidina
'R_EX_spmd_e';          % Espermidina (cpd00264)
'R_EX_pnto__R_e';       % Pantotenato (cpd00644 -> pnto__R)
'R_EX_h2o__R_e';
'R_EX_nh4_e';           % Amonio (algunos modelos usan nh4 en vez de nh3)
};

% Rates individuales para cada componente del medio
% Glucosa: fuente principal de carbono, limitante a 10 mmol/gDW/h
% Aminoacidos: 0.5 mmol/gDW/h cada uno (~10 total, ~40 mmol C) — suficiente
%   para biosintesis y auxotrofias (LAB como M38 necesitan mas), pero fuerza
%   al modelo a usar glucosa (60 mmol C) como fuente principal de carbono.
%   Con rate=1000 o 1, Bt cataboliza aminoacidos en vez de glucosa (glucose=0).
%   Con rate=0.1, M38 queda infactible (auxotrofias). 0.5 es el balance.
% Iones/minerales/vitaminas/cofactores: 1000 (exceso, no limitantes)
rates_mZMB = zeros(length(mZMB), 1);
rates_mZMB(1) = 10;   % glucosa
% Aminoacidos (indices 2-21): 0.5 mmol/gDW/h cada uno
% Con rate=1000, los modelos catabolizan aminoacidos para ATP (ATPM=1000) en vez
% de usar glucosa como fuente principal de carbono.
% Con rate=0.1, M38 queda infactible (auxotrofias).
% 0.5 es el balance: suficiente para biosintesis, fuerza uso de glucosa para energia.
for idx_aa = 2:21
    rates_mZMB(idx_aa) = 1000;
end
% Iones, minerales, vitaminas, cofactores (indices 22-end): exceso
for idx_other = 22:length(mZMB)
    rates_mZMB(idx_other) = 1000;
end

%% 2. BIOMASA EXPERIMENTAL (ajustar a tus valores)
%% Define the old reaction name and the new reaction name
oldReactionName = 'R_Growth';
newReactionName = 'Growth_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Bacteroides_thetaiotaomicron_VPI_5482.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Bacteroides_thetaiotaomicron_VPI_5482.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end
%% Define the old reaction name and the new reaction name
oldReactionName = 'R_Growth';
newReactionName = 'Growth_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Bifidobacterium_animalis_lactis_PT33.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Bifidobacterium_animalis_lactis_PT33.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end
%% Define the old reaction name and the new reaction name
oldReactionName = 'R_Growth';
newReactionName = 'Growth_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Clostridium_sp_HGF2.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Clostridium_sp_HGF2.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end
%% Define the old reaction name and the new reaction name
oldReactionName = 'R_Growth';
newReactionName = 'Growth_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Clostridium_sp_M62_1.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Clostridium_sp_M62_1.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end
%% Define the old reaction name and the new reaction name
oldReactionName = 'R_Growth';
newReactionName = 'Growth_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Clostridium_symbiosum_WAL_14673.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Clostridium_symbiosum_WAL_14673.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end
%% Define the old reaction name and the new reaction name
oldReactionName = 'R_Growth';
newReactionName = 'Growth_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Lacticaseibacillus_paracasei_M38.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Lacticaseibacillus_paracasei_M38.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end
%%
Bacteroides_thetaiotaomicron_VPI_5482 = changeObjective(Bacteroides_thetaiotaomicron_VPI_5482, 'Growth_biomass', 1);
Bifidobacterium_animalis_lactis_PT33 = changeObjective(Bifidobacterium_animalis_lactis_PT33, 'Growth_biomass', 1);
Clostridium_sp_HGF2 = changeObjective(Clostridium_sp_HGF2, 'Growth_biomass', 1);
Clostridium_sp_M62_1 = changeObjective(Clostridium_sp_M62_1, 'Growth_biomass', 1);
Clostridium_symbiosum_WAL_14673 = changeObjective(Clostridium_symbiosum_WAL_14673, 'Growth_biomass', 1);
Lacticaseibacillus_paracasei_M38 = changeObjective(Lacticaseibacillus_paracasei_M38, 'Growth_biomass', 1);

biomasa_exp = [0.648; 0.357; 0.33; 0.141; 0.403; 0.21];  % h^-1

% Bacteroides_thetaiotaomicron_VPI_5482', '
% Bifidobacterium_animalis_lactis_PT33', 
% 'Clostridium_sp_HGF2', '
% Clostridium_sp_M62_1', 
% 'Clostridium_symbiosum_WAL_14673', 
% 'Lacticaseibacillus_paracasei_M38'};

%% CORRECCION FINAL BUTIRATO: BUTCT2 + eliminar BUTKr + Bacteroides
% =========================================================================
% Basado en analisis del modelo de referencia iCN900 (C. difficile, BiGG)
%
% Ruta de sintesis de butirato via acetil-CoA en los 3 Clostridium:
%   1. ACACT1r:    2 accoa_c -> aacoa_c + coa_c
%   2. HACD1/i:    aacoa_c + h_c + nadh_c -> 3hbcoa_c + nad_c
%   3. ECOAH1:     3hbcoa_c -> b2coa_c + h2o_c
%   4. ACOAD1/1fr: b2coa_c + cofactor_red -> btcoa_c + cofactor_ox
%   5. BUTCT2:     acac_c + btcoa_c -> aacoa_c + but_c  (PRODUCE BUTIRATO)
%
% Cambios:
%   - R_BUTKr: eliminar de todos los productores (no debe estar)
%   - R_BUTCT2: agregar como paso terminal (produce butirato)
%   - R_PBUTT: verificar presencia en los 3 (agregar si falta)
%   - Bacteroides: eliminar toda capacidad de butirato
% =========================================================================

fprintf('\n=== CORRECCION FINAL BUTIRATO ===\n');

% --- HGF2: eliminar BUTKr, agregar BUTCT2 ---
m = Clostridium_sp_HGF2;
fprintf('\nClostridium_sp_HGF2:\n');
if any(strcmp(m.rxns, 'R_BUTKr'))
    m = removeRxns(m, {'R_BUTKr'});
    fprintf('  Eliminada: R_BUTKr\n');
end
m = addReaction(m, 'R_BUTCT2', ...
    'reactionName', 'Acetoacetate:butyrate CoA-transferase', ...
    'metaboliteList', {'M_acac_c', 'M_btcoa_c', 'M_aacoa_c', 'M_but_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Pyruvate metabolism');
fprintf('  Agregada: R_BUTCT2 (acac + btcoa -> aacoa + but)\n');
Clostridium_sp_HGF2 = m;

% --- M62_1: eliminar BUTKr (agregada anteriormente en script), agregar BUTCT2 ---
m = Clostridium_sp_M62_1;
fprintf('\nClostridium_sp_M62_1:\n');
if any(strcmp(m.rxns, 'R_BUTKr'))
    m = removeRxns(m, {'R_BUTKr'});
    fprintf('  Eliminada: R_BUTKr\n');
end
m = addReaction(m, 'R_BUTCT2', ...
    'reactionName', 'Acetoacetate:butyrate CoA-transferase', ...
    'metaboliteList', {'M_acac_c', 'M_btcoa_c', 'M_aacoa_c', 'M_but_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Pyruvate metabolism');
fprintf('  Agregada: R_BUTCT2\n');
Clostridium_sp_M62_1 = m;

% --- C.symbiosum: agregar PBUTT + M_butpi_c (faltantes), agregar BUTCT2 ---
m = Clostridium_symbiosum_WAL_14673;
fprintf('\nClostridium_symbiosum_WAL_14673:\n');
if any(strcmp(m.rxns, 'R_BUTKr'))
    m = removeRxns(m, {'R_BUTKr'});
    fprintf('  Eliminada: R_BUTKr\n');
end
if ~any(strcmp(m.mets, 'M_butpi_c'))
    m = addMetabolite(m, 'M_butpi_c', 'metName', 'Butanoyl phosphate', ...
        'metFormula', 'C4H7O5P');
    fprintf('  Agregado: M_butpi_c\n');
end
if ~any(strcmp(m.rxns, 'R_PBUTT'))
    m = addReaction(m, 'R_PBUTT', ...
        'reactionName', 'Phosphate butyryltransferase', ...
        'metaboliteList', {'M_btcoa_c', 'M_pi_c', 'M_butpi_c', 'M_coa_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Pyruvate metabolism');
    fprintf('  Agregada: R_PBUTT\n');
end
m = addReaction(m, 'R_BUTCT2', ...
    'reactionName', 'Acetoacetate:butyrate CoA-transferase', ...
    'metaboliteList', {'M_acac_c', 'M_btcoa_c', 'M_aacoa_c', 'M_but_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
    'subSystem', 'Pyruvate metabolism');
fprintf('  Agregada: R_BUTCT2\n');
Clostridium_symbiosum_WAL_14673 = m;

% --- Bacteroides: eliminar reacciones de butirato ---
m = Bacteroides_thetaiotaomicron_VPI_5482;
fprintf('\nBacteroides_thetaiotaomicron_VPI_5482:\n');
rxns_but_remove = {'R_BUTKr', 'R_PBUTT', 'R_BUTt4pp', 'R_BUTtex', 'R_EX_but_e', 'R_BUTt2rpp'};
existing = rxns_but_remove(ismember(rxns_but_remove, m.rxns));
if ~isempty(existing)
    m = removeRxns(m, existing);
    fprintf('  Eliminadas: %s\n', strjoin(existing, ', '));
else
    fprintf('  Sin reacciones de butirato para eliminar\n');
end
Bacteroides_thetaiotaomicron_VPI_5482 = m;

% --- Verificacion final ---
fprintf('\n--- VERIFICACION FINAL BUTIRATO ---\n');
modelos_check = {Clostridium_sp_HGF2, Clostridium_sp_M62_1, ...
                 Clostridium_symbiosum_WAL_14673, Bacteroides_thetaiotaomicron_VPI_5482};
nombres_check = {'HGF2', 'M62_1', 'Csym', 'Bt'};
rxns_sintesis = {'R_ACACT1r', 'R_HACD1i', 'R_ACOAD1', 'R_ACOAD1fr', 'R_PBUTT', 'R_BUTCT2'};
% Bt (gram-neg, indice 4): BUTtex+BUTt2rpp; gram+ (indices 1-3): BUTt2r
for i = 1:4
    mv = modelos_check{i};
    fprintf('\n  %s:\n', nombres_check{i});
    for j = 1:length(rxns_sintesis)
        if any(strcmp(mv.rxns, rxns_sintesis{j}))
            fprintf('    [OK] %s\n', rxns_sintesis{j});
        else
            fprintf('    [--] %s (ausente)\n', rxns_sintesis{j});
        end
    end
    % Export: Bt usa BUTtex+BUTt2rpp, gram+ usa BUTt2r (directo e<->c)
    if i == 4  % Bt
        rxns_export = {'R_EX_but_e', 'R_BUTtex', 'R_BUTt2rpp'};
    else  % gram+
        rxns_export = {'R_EX_but_e', 'R_BUTt2r'};
    end
    for j = 1:length(rxns_export)
        if any(strcmp(mv.rxns, rxns_export{j}))
            fprintf('    [OK] %s (export)\n', rxns_export{j});
        else
            fprintf('    [--] %s (ausente)\n', rxns_export{j});
        end
    end
    if any(strcmp(mv.rxns, 'R_BUTKr'))
        fprintf('    [!!] R_BUTKr PRESENTE (deberia estar ausente!)\n');
    else
        fprintf('    [OK] R_BUTKr ausente (correcto)\n');
    end
end
fprintf('\n==========================================\n');

%% ADICION RUTA SINTESIS BUTIRATO: acetil-CoA -> butanoil-CoA -> butirato
% =========================================================================
% Ruta completa de butirato via acetil-CoA en Clostridium:
%   1. R_ACACT1r:   2 accoa_c <-> aacoa_c + coa_c            (tiolasa)
%   2. R_HACD1:     aacoa_c + h_c + nadh_c -> 3hbcoa_c + nad_c (3-HB-CoA DH, dir. sintesis)
%   3. R_ECOAH1:    3hbcoa_c -> b2coa_c + h2o_c               (crotonasa)
%   4. R_ACOAD1fr:  b2coa_c + fadh2_c <-> btcoa_c + fad_c     (butiril-CoA DH)
%   5. R_BUTCT2:    acac_c + btcoa_c -> aacoa_c + but_c       (CoA-transferasa, los 3)
%   6. R_BUTCT:     ac_c + btcoa_c -> accoa_c + but_c         (solo HGF2)
%
% R_BUTCT (BiGG) = BTCOAACCOAT = BUTCTr (VMH) = MNXR147175 (MetaNetX)
% EC 2.8.3.8 - Ruta principal en C. innocuum y Clostridium spp.
%
% ELIMINACIONES:
%   - R_PBUTT y R_BUTKr: eliminar de los 3 (ninguna tiene butirato kinasa)
%   - La ruta PBUTT+BUTKr (btcoa->butpi->but) no aplica a estas bacterias
%
% Produccion terminal de butirato:
%   - Los 3:  R_BUTCT2 (acac + btcoa -> aacoa + but)
%   - HGF2:   R_BUTCT2 + R_BUTCT (ac + btcoa -> accoa + but)
%
% Referencia: iCN900 (C. difficile, BiGG) + BiGG universal model + MetaNetX
% =========================================================================

fprintf('\n=== ADICION RUTA SINTESIS BUTIRATO (accoa -> btcoa -> but) ===\n');

productores_but = {Clostridium_sp_HGF2, Clostridium_sp_M62_1, Clostridium_symbiosum_WAL_14673};
nombres_prod = {'Clostridium_sp_HGF2', 'Clostridium_sp_M62_1', 'Clostridium_symbiosum_WAL_14673'};

for i = 1:3
    m = productores_but{i};
    fprintf('\n%s:\n', nombres_prod{i});

    % Asegurar metabolitos intermediarios (deberian existir de CarveMe)
    mets_check = {'M_3hbcoa_c', 'M_b2coa_c'};
    mets_names = {'(S)-3-Hydroxybutanoyl-CoA', 'Crotonoyl-CoA'};
    mets_formulas = {'C25H38N7O18P3S', 'C25H36N7O17P3S'};
    for j = 1:length(mets_check)
        if ~any(strcmp(m.mets, mets_check{j}))
            m = addMetabolite(m, mets_check{j}, 'metName', mets_names{j}, ...
                'metFormula', mets_formulas{j});
            fprintf('  Agregado metabolito: %s\n', mets_check{j});
        end
    end

    % --- ELIMINAR R_PBUTT y R_BUTKr (no aplican, sin butirato kinasa) ---
    rxns_eliminar = {'R_PBUTT', 'R_BUTKr'};
    for j = 1:length(rxns_eliminar)
        if any(strcmp(m.rxns, rxns_eliminar{j}))
            m = removeRxns(m, rxns_eliminar(j));
            fprintf('  Eliminada: %s\n', rxns_eliminar{j});
        end
    end

    % --- Paso 1: R_ACACT1r (tiolasa, reversible) ---
    % 2 accoa_c <-> aacoa_c + coa_c
    if ~any(strcmp(m.rxns, 'R_ACACT1r'))
        m = addReaction(m, 'R_ACACT1r', ...
            'reactionName', 'Acetyl-CoA C-acetyltransferase', ...
            'metaboliteList', {'M_accoa_c', 'M_aacoa_c', 'M_coa_c'}, ...
            'stoichCoeffList', [-2, 1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Fatty acid oxidation');
        fprintf('  Agregada: R_ACACT1r (2 accoa <-> aacoa + coa)\n');
    else
        fprintf('  [OK] R_ACACT1r ya existe\n');
    end

    % --- Paso 2: R_HACD1 (3-HB-CoA deshidrogenasa, dir. sintesis) ---
    % aacoa_c + h_c + nadh_c -> 3hbcoa_c + nad_c
    if ~any(strcmp(m.rxns, 'R_HACD1'))
        m = addReaction(m, 'R_HACD1', ...
            'reactionName', '3-hydroxyacyl-CoA dehydrogenase (acetoacetyl-CoA)', ...
            'metaboliteList', {'M_aacoa_c', 'M_h_c', 'M_nadh_c', 'M_3hbcoa_c', 'M_nad_c'}, ...
            'stoichCoeffList', [-1, -1, -1, 1, 1], ...
            'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
            'subSystem', 'Fatty acid oxidation');
        fprintf('  Agregada: R_HACD1 (aacoa + nadh -> 3hbcoa + nad)\n');
    else
        fprintf('  [OK] R_HACD1 ya existe\n');
    end

    % --- Paso 3: R_ECOAH1 (crotonasa / enoil-CoA hidratasa) ---
    % 3hbcoa_c -> b2coa_c + h2o_c
    if ~any(strcmp(m.rxns, 'R_ECOAH1'))
        m = addReaction(m, 'R_ECOAH1', ...
            'reactionName', '3-hydroxyacyl-CoA dehydratase (3-hydroxybutanoyl-CoA)', ...
            'metaboliteList', {'M_3hbcoa_c', 'M_b2coa_c', 'M_h2o_c'}, ...
            'stoichCoeffList', [-1, 1, 1], ...
            'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
            'subSystem', 'Fatty acid oxidation');
        fprintf('  Agregada: R_ECOAH1 (3hbcoa -> b2coa + h2o)\n');
    else
        fprintf('  [OK] R_ECOAH1 ya existe\n');
    end

    % --- Paso 4: R_BTCOADH (complejo Bcd-EtfAB, electron bifurcation) ---
    % VMH: BTCOADH = MetaNetX: buscar. Usa NADH + ferredoxina (no FADH2)
    % b2coa + fdxox + 2h + 2nadh -> btcoa + fdxrd + 2nad
    % Ref: AGORA2 HGF2, Csym, M62_1
    % Eliminar R_ACOAD1fr si existe (usaba FADH2, no disponible en anaerobiosis)
    if any(strcmp(m.rxns, 'R_ACOAD1fr'))
        m = removeRxns(m, 'R_ACOAD1fr');
        fprintf('  Eliminada: R_ACOAD1fr (usaba FADH2)\n');
    end
    if ~any(strcmp(m.rxns, 'R_BTCOADH'))
        % Verificar/agregar metabolito fdxrd (fdxo_2_2 ya existe en CarveMe como M_fdxo_2_2_c)
        if ~any(strcmp(m.mets, 'M_fdxrd_c'))
            m = addMetabolite(m, 'M_fdxrd_c', 'metName', 'Reduced ferredoxin', ...
                'metFormula', 'Fe2S2', 'Charge', 1);
            fprintf('  Agregado metabolito: M_fdxrd_c\n');
        end
        m = addReaction(m, 'R_BTCOADH', ...
            'reactionName', 'Butyryl-CoA dehydrogenase (Bcd-EtfAB electron bifurcation)', ...
            'metaboliteList', {'M_b2coa_c', 'M_fdxo_2_2_c', 'M_nadh_c', ...
                               'M_btcoa_c', 'M_fdxrd_c', 'M_nad_c'}, ...
            'stoichCoeffList', [-1, -2, -2, 1, 2, 2], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Pyruvate metabolism');
        fprintf('  Agregada: R_BTCOADH (b2coa + fdxox + 2h + 2nadh <-> btcoa + fdxrd + 2nad)\n');
    else
        fprintf('  [OK] R_BTCOADH ya existe\n');
    end

    % --- Paso 4b: R_FDNADOX_H (reciclaje de ferredoxina, tipo Rnf) ---
    % AGORA2: FDNADOX_H en los 3 Clostridia (HGF2, M62_1, Csym)
    % fdxrd + nad <-> fdxox + h[e] + nadh
    % SIN ESTA REACCION, fdxrd se acumula y fdxox se agota -> BTCOADH bloqueado
    if ~any(strcmp(m.rxns, 'R_FDNADOX_H'))
        m = addReaction(m, 'R_FDNADOX_H', ...
            'reactionName', 'Ferredoxin:NAD oxidoreductase (proton motive force)', ...
            'metaboliteList', {'M_fdxrd_c', 'M_nad_c', 'M_h_c', 'M_fdxo_2_2_c', 'M_h_e', 'M_nadh_c'}, ...
            'stoichCoeffList', [-2, -1, -2, 2, 1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Oxidative phosphorylation');
        fprintf('  Agregada: R_FDNADOX_H (fdxrd + nad <-> fdxox + h[e] + nadh) [Rnf complex]\n');
    else
        fprintf('  [OK] R_FDNADOX_H ya existe\n');
    end

    % --- Paso 4c: R_POR4 (piruvato:ferredoxina oxidorreductasa, PFOR) ---
    % AGORA2: POR4 en los 3 Clostridia
    % Escrito en BiGG como reversa: accoa + co2 + fdxrd <-> coa + fdxox + h + pyr
    % Fisiologicamente: pyr + coa + fdxox -> accoa + co2 + fdxrd
    % Reemplaza PDH en anaerobios (usa ferredoxina en vez de NAD)
    if ~any(strcmp(m.rxns, 'R_POR4'))
        m = addReaction(m, 'R_POR4', ...
            'reactionName', 'Pyruvate:ferredoxin oxidoreductase', ...
            'metaboliteList', {'M_accoa_c', 'M_co2_c', 'M_fdxrd_c', 'M_h_c', 'M_coa_c', 'M_fdxo_2_2_c', 'M_pyr_c'}, ...
            'stoichCoeffList', [-1, -1, -2, -1, 1, 2, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Pyruvate metabolism');
        fprintf('  Agregada: R_POR4 (accoa + co2 + fdxrd <-> coa + fdxox + h + pyr) [PFOR]\n');
    else
        fprintf('  [OK] R_POR4 ya existe\n');
    end

    % --- Paso 5: R_BUTCT2 (CoA-transferasa, los 3 productores) ---
    % acac_c + btcoa_c -> aacoa_c + but_c
    if ~any(strcmp(m.rxns, 'R_BUTCT2'))
        m = addReaction(m, 'R_BUTCT2', ...
            'reactionName', 'Acetoacetate:butyrate CoA-transferase', ...
            'metaboliteList', {'M_acac_c', 'M_btcoa_c', 'M_aacoa_c', 'M_but_c'}, ...
            'stoichCoeffList', [-1, -1, 1, 1], ...
            'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
            'subSystem', 'Pyruvate metabolism');
        fprintf('  Agregada: R_BUTCT2 (acac + btcoa -> aacoa + but)\n');
    else
        fprintf('  [OK] R_BUTCT2 ya existe\n');
    end

    % --- Paso 5b: R_ACACCT (produce acetoacetato libre para BUTCT2) ---
    % BiGG: ACACCT = VMH: ACACCTr = MetaNetX: MNXR95193
    % EC 2.8.3.8 / 2.8.3.9
    % acac_c + accoa_c <-> aacoa_c + ac_c
    % En reversa: aacoa -> acac (produce acetoacetato libre que BUTCT2 necesita)
    % Presente en AGORA2: HGF2 (ACACCTr), Csym (ACACCTr). M62_1 no la tiene.
    if ~any(strcmp(m.rxns, 'R_ACACCT'))
        m = addReaction(m, 'R_ACACCT', ...
            'reactionName', 'Acetyl-CoA:acetoacetyl-CoA transferase', ...
            'metaboliteList', {'M_acac_c', 'M_accoa_c', 'M_aacoa_c', 'M_ac_c'}, ...
            'stoichCoeffList', [-1, -1, 1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Pyruvate metabolism');
        fprintf('  Agregada: R_ACACCT (acac + accoa <-> aacoa + ac) [BiGG, EC 2.8.3.8/9]\n');
    else
        fprintf('  [OK] R_ACACCT ya existe\n');
    end

    % --- Paso 6: R_BUTCT (los 3 productores) ---
    % BiGG: BUTCT = VMH: BTCOAACCOAT = VMH: BUTCTr = MetaNetX: MNXR147175
    % EC 2.8.3.8 - ac_c + btcoa_c -> accoa_c + but_c
    % Ruta directa: no necesita acetoacetato libre
    if ~any(strcmp(m.rxns, 'R_BUTCT'))
        m = addReaction(m, 'R_BUTCT', ...
            'reactionName', 'Acetyl-CoA:butyrate-CoA transferase', ...
            'metaboliteList', {'M_ac_c', 'M_btcoa_c', 'M_accoa_c', 'M_but_c'}, ...
            'stoichCoeffList', [-1, -1, 1, 1], ...
            'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
            'subSystem', 'Pyruvate metabolism');
        fprintf('  Agregada: R_BUTCT (ac + btcoa -> accoa + but) [EC 2.8.3.8]\n');
    else
        fprintf('  [OK] R_BUTCT ya existe\n');
    end

    productores_but{i} = m;
end

Clostridium_sp_HGF2 = productores_but{1};
Clostridium_sp_M62_1 = productores_but{2};
Clostridium_symbiosum_WAL_14673 = productores_but{3};

% --- R_HYD4 solo para C. symbiosum (AGORA2) ---
% fdxrd -> fdxox + h2 (irreversible, hidrogenasa)
m = Clostridium_symbiosum_WAL_14673;
if ~any(strcmp(m.rxns, 'R_HYD4'))
    m = addReaction(m, 'R_HYD4', ...
        'reactionName', 'Hydrogen:ferredoxin oxidoreductase', ...
        'metaboliteList', {'M_fdxrd_c', 'M_h_c', 'M_fdxo_2_2_c', 'M_h2_c'}, ...
        'stoichCoeffList', [-2, -2, 2, 1], ...
        'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
        'subSystem', 'Pyruvate metabolism');
    fprintf('  Csym: Agregada R_HYD4 (fdxrd -> fdxox + h2) [hidrogenasa]\n');
end
Clostridium_symbiosum_WAL_14673 = m;

% --- Transporte directo de butirato e<->c para gram-positivas ---
% Las gram-positivas no usan ruta periplasmica eficientemente
% AGORA2 usa transporte directo. Agregar R_BUTt2r (BiGG) para los 3 productores
fprintf('\n--- TRANSPORTE DIRECTO BUTIRATO (gram-positivas) ---\n');
productores_but_t = {Clostridium_sp_HGF2, Clostridium_sp_M62_1, Clostridium_symbiosum_WAL_14673};
for i = 1:3
    m = productores_but_t{i};
    % R_BUTt2r: but_c + h_c <-> but_e + h_e (simporte H+, reversible, directo e<->c)
    if ~any(strcmp(m.rxns, 'R_BUTt2r'))
        if ~any(strcmp(m.mets, 'M_but_e'))
            m = addMetabolite(m, 'M_but_e', 'metName', 'Butyrate (n-C4:0)', 'metFormula', 'C4H7O2', 'Charge', -1);
        end
        m = addReaction(m, 'R_BUTt2r', ...
            'reactionName', 'Butyrate transport via proton symport (direct e<->c)', ...
            'metaboliteList', {'M_but_c', 'M_h_c', 'M_but_e', 'M_h_e'}, ...
            'stoichCoeffList', [-1, -1, 1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Transport');
        fprintf('  %s: Agregada R_BUTt2r (but_c + h_c <-> but_e + h_e)\n', nombres_prod{i});
    else
        fprintf('  %s: [OK] R_BUTt2r ya existe\n', nombres_prod{i});
    end
    % Asegurar EX_but_e sin R_ prefix tambien existe (para evitar duplicados en FVA)
    if ~any(strcmp(m.rxns, 'EX_but_e')) && ~any(strcmp(m.rxns, 'R_EX_but_e'))
        m = addReaction(m, 'R_EX_but_e', ...
            'reactionName', 'Butyrate exchange', ...
            'metaboliteList', {'M_but_e'}, ...
            'stoichCoeffList', [-1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Exchange/demand reaction');
    end
    productores_but_t{i} = m;
end
Clostridium_sp_HGF2 = productores_but_t{1};
Clostridium_sp_M62_1 = productores_but_t{2};
Clostridium_symbiosum_WAL_14673 = productores_but_t{3};

% --- Verificacion ruta completa ---
fprintf('\n--- VERIFICACION RUTA COMPLETA BUTIRATO ---\n');
rxns_sintesis = {'R_ACACT1r', 'R_HACD1', 'R_ECOAH1', 'R_BTCOADH'};
rxns_fdx = {'R_FDNADOX_H', 'R_POR4'};
rxns_terminal = {'R_ACACCT', 'R_BUTCT2', 'R_BUTCT'};
rxns_export = {'R_EX_but_e', 'R_BUTt2r'};
rxns_ausentes = {'R_BUTKr', 'R_PBUTT'};

modelos_verif_but = {Clostridium_sp_HGF2, Clostridium_sp_M62_1, Clostridium_symbiosum_WAL_14673};
for i = 1:3
    mv = modelos_verif_but{i};
    fprintf('\n  %s:\n', nombres_prod{i});
    fprintf('  Sintesis (accoa -> btcoa):\n');
    for j = 1:length(rxns_sintesis)
        if any(strcmp(mv.rxns, rxns_sintesis{j}))
            fprintf('    [OK] %s\n', rxns_sintesis{j});
        else
            fprintf('    [--] %s (FALTA!)\n', rxns_sintesis{j});
        end
    end
    fprintf('  Reciclaje ferredoxina:\n');
    for j = 1:length(rxns_fdx)
        if any(strcmp(mv.rxns, rxns_fdx{j}))
            fprintf('    [OK] %s\n', rxns_fdx{j});
        else
            fprintf('    [--] %s (FALTA!)\n', rxns_fdx{j});
        end
    end
    fprintf('  Produccion terminal (btcoa -> but):\n');
    for j = 1:length(rxns_terminal)
        if any(strcmp(mv.rxns, rxns_terminal{j}))
            fprintf('    [OK] %s\n', rxns_terminal{j});
        else
            fprintf('    [--] %s (FALTA!)\n', rxns_terminal{j});
        end
    end
    fprintf('  Exportacion:\n');
    for j = 1:length(rxns_export)
        if any(strcmp(mv.rxns, rxns_export{j}))
            fprintf('    [OK] %s\n', rxns_export{j});
        else
            fprintf('    [--] %s (FALTA!)\n', rxns_export{j});
        end
    end
    fprintf('  Deben estar ausentes:\n');
    for j = 1:length(rxns_ausentes)
        if any(strcmp(mv.rxns, rxns_ausentes{j}))
            fprintf('    [!!] %s PRESENTE (deberia estar ausente!)\n', rxns_ausentes{j});
        else
            fprintf('    [OK] %s ausente (correcto)\n', rxns_ausentes{j});
        end
    end
end
fprintf('\n==========================================\n');

%% CURACION EXPORTACION DE ACETATO
% =========================================================================
% Modelos que producen acetato pero NO pueden exportarlo:
%   - PT33 (Bifidobacterium): falta M_ac_e, M_ac_p, R_EX_ac_e, R_ACtex, R_ACt2rpp
%   - HGF2 (Clostridium): falta M_ac_e, M_ac_p, R_EX_ac_e, R_ACtex, R_ACt2rpp
%   - M62_1 (Clostridium): falta M_ac_e, M_ac_p, R_EX_ac_e, R_ACtex, R_ACt2rpp
%   - M38 (Lacticaseibacillus): falta M_ac_e, M_ac_p, R_EX_ac_e, R_ACtex, R_ACt2rpp
%
% Modelos OK (ya tienen acetato export): Bt, Csym
%
% Reacciones de referencia de iML1515 (E. coli, BiGG):
%   R_EX_ac_e:   M_ac_e <->                           (exchange, LB=-1000, UB=1000)
%   R_ACtex:     M_ac_e <-> M_ac_p                    (difusion membrana externa)
%   R_ACt2rpp:   M_ac_p + M_h_p <-> M_ac_c + M_h_c   (simporte proton, membrana interna)
% =========================================================================

fprintf('\n=== CURACION EXPORTACION DE ACETATO ===\n');

modelos_ac = {Bifidobacterium_animalis_lactis_PT33, Clostridium_sp_HGF2, ...
              Clostridium_sp_M62_1, Lacticaseibacillus_paracasei_M38};
nombres_ac = {'Bifidobacterium_animalis_lactis_PT33', 'Clostridium_sp_HGF2', ...
              'Clostridium_sp_M62_1', 'Lacticaseibacillus_paracasei_M38'};

for i = 1:4
    m = modelos_ac{i};
    fprintf('\n%s (gram+, sin periplasma):\n', nombres_ac{i});

    % Gram+: solo agregar M_ac_e y R_EX_ac_e (sin _p, sin tex/tpp)
    % El transporte directo e<->c (R_ACtr) se agrega en la seccion siguiente
    if ~any(strcmp(m.mets, 'M_ac_e'))
        m = addMetabolite(m, 'M_ac_e', 'metName', 'Acetate_e', ...
            'metFormula', 'C2H3O2', 'Charge', -1);
        fprintf('  Agregado metabolito: M_ac_e\n');
    end

    % R_EX_ac_e: ac_e <-> (exchange)
    if ~any(strcmp(m.rxns, 'R_EX_ac_e'))
        m = addReaction(m, 'R_EX_ac_e', ...
            'reactionName', 'Acetate exchange', ...
            'metaboliteList', {'M_ac_e'}, ...
            'stoichCoeffList', [-1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Exchange/demand reaction');
        fprintf('  Agregada: R_EX_ac_e\n');
    else
        fprintf('  [OK] R_EX_ac_e ya existe\n');
    end

    modelos_ac{i} = m;
end

Bifidobacterium_animalis_lactis_PT33 = modelos_ac{1};
Clostridium_sp_HGF2 = modelos_ac{2};
Clostridium_sp_M62_1 = modelos_ac{3};
Lacticaseibacillus_paracasei_M38 = modelos_ac{4};

%% CURACION TRANSPORTE DIRECTO e<->c PARA GRAM-POSITIVAS
% =========================================================================
% Las gram-positivas necesitan transporte directo e<->c (no via periplasma)
% AGORA2 usa R_ACtr para acetato y R_L_LACt2r / R_D_LACt2 para lactato
% Ref: AGORA2 de BLC1 (PT33), LC2W (M38), HGF2, M62_1, Csym
% =========================================================================
fprintf('\n=== TRANSPORTE DIRECTO e<->c PARA GRAM-POSITIVAS ===\n');

grampositivas = {Bifidobacterium_animalis_lactis_PT33, Clostridium_sp_HGF2, ...
                 Clostridium_sp_M62_1, Clostridium_symbiosum_WAL_14673, ...
                 Lacticaseibacillus_paracasei_M38};
nombres_gp = {'PT33', 'HGF2', 'M62_1', 'Csym', 'M38'};

for i = 1:5
    m = grampositivas{i};
    fprintf('\n%s:\n', nombres_gp{i});

    % R_ACtr: acetato transporte directo e<->c (BiGG valido, AGORA2: ACtr)
    if ~any(strcmp(m.rxns, 'R_ACtr'))
        % Asegurar que M_ac_e existe
        if ~any(strcmp(m.mets, 'M_ac_e'))
            m = addMetabolite(m, 'M_ac_e', 'metName', 'Acetate_e', 'metFormula', 'C2H3O2', 'Charge', -1);
        end
        m = addReaction(m, 'R_ACtr', ...
            'reactionName', 'Acetate transporter (direct e<->c)', ...
            'metaboliteList', {'M_ac_e', 'M_ac_c'}, ...
            'stoichCoeffList', [-1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Transport');
        fprintf('  Agregada: R_ACtr (acetato directo e<->c)\n');
    else
        fprintf('  [OK] R_ACtr ya existe\n');
    end

    % R_L_LACt2r: L-lactato transporte directo e<->c con H+ (BiGG valido)
    if ~any(strcmp(m.rxns, 'R_L_LACt2r'))
        if ~any(strcmp(m.mets, 'M_lac__L_e'))
            m = addMetabolite(m, 'M_lac__L_e', 'metName', 'L-Lactate_e', 'metFormula', 'C3H5O3', 'Charge', -1);
        end
        m = addReaction(m, 'R_L_LACt2r', ...
            'reactionName', 'L-Lactate reversible transport via proton symport (direct e<->c)', ...
            'metaboliteList', {'M_lac__L_e', 'M_h_e', 'M_lac__L_c', 'M_h_c'}, ...
            'stoichCoeffList', [-1, -1, 1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Transport');
        fprintf('  Agregada: R_L_LACt2r (L-lactato directo e<->c)\n');
    else
        fprintf('  [OK] R_L_LACt2r ya existe\n');
    end

    % R_D_LACt2: D-lactato transporte directo e<->c con H+ (BiGG valido)
    if ~any(strcmp(m.rxns, 'R_D_LACt2'))
        if ~any(strcmp(m.mets, 'M_lac__D_e'))
            m = addMetabolite(m, 'M_lac__D_e', 'metName', 'D-Lactate_e', 'metFormula', 'C3H5O3', 'Charge', -1);
        end
        m = addReaction(m, 'R_D_LACt2', ...
            'reactionName', 'D-Lactate transport via proton symport (direct e<->c)', ...
            'metaboliteList', {'M_lac__D_e', 'M_h_e', 'M_lac__D_c', 'M_h_c'}, ...
            'stoichCoeffList', [-1, -1, 1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Transport');
        fprintf('  Agregada: R_D_LACt2 (D-lactato directo e<->c)\n');
    else
        fprintf('  [OK] R_D_LACt2 ya existe\n');
    end

    % Asegurar que EX_ac_e, EX_lac__L_e, EX_lac__D_e existen
    if ~any(strcmp(m.rxns, 'R_EX_ac_e'))
        m = addReaction(m, 'R_EX_ac_e', 'reactionName', 'Acetate exchange', ...
            'metaboliteList', {'M_ac_e'}, 'stoichCoeffList', [-1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Exchange/demand reaction');
        fprintf('  Agregada: R_EX_ac_e\n');
    end
    if ~any(strcmp(m.rxns, 'R_EX_lac__L_e'))
        m = addReaction(m, 'R_EX_lac__L_e', 'reactionName', 'L-Lactate exchange', ...
            'metaboliteList', {'M_lac__L_e'}, 'stoichCoeffList', [-1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Exchange/demand reaction');
        fprintf('  Agregada: R_EX_lac__L_e\n');
    end
    if ~any(strcmp(m.rxns, 'R_EX_lac__D_e'))
        m = addReaction(m, 'R_EX_lac__D_e', 'reactionName', 'D-Lactate exchange', ...
            'metaboliteList', {'M_lac__D_e'}, 'stoichCoeffList', [-1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Exchange/demand reaction');
        fprintf('  Agregada: R_EX_lac__D_e\n');
    end

    grampositivas{i} = m;
end

Bifidobacterium_animalis_lactis_PT33 = grampositivas{1};
Clostridium_sp_HGF2 = grampositivas{2};
Clostridium_sp_M62_1 = grampositivas{3};
Clostridium_symbiosum_WAL_14673 = grampositivas{4};
Lacticaseibacillus_paracasei_M38 = grampositivas{5};

%% CURACION EXPORTACION DE L-LACTATO
% =========================================================================
% Modelos que necesitan exportacion de L-lactato:
%   - PT33 (Bifidobacterium): productor de lactato, falta export L y D
%   - M38 (Lacticaseibacillus): productor de lactato, falta export L y D
%   - Bt (Bacteroides): tiene D-lactato export pero falta L-lactato export
%
% Modelos OK: Csym (tiene todo), HGF2 y M62_1 (no producen L-lactato)
%
% Reacciones de referencia de iML1515 (E. coli, BiGG):
%   R_EX_lac__L_e:  M_lac__L_e <->                              (exchange)
%   R_L_LACtex:     M_lac__L_e <-> M_lac__L_p                   (difusion membrana ext)
%   R_L_LACt2rpp:   M_lac__L_p + M_h_p <-> M_lac__L_c + M_h_c  (simporte proton)
% =========================================================================

fprintf('\n=== CURACION EXPORTACION DE L-LACTATO ===\n');

% --- Gram-positivas (PT33, M38): solo M_lac__L_e + R_EX_lac__L_e ---
% El transporte directo e<->c (R_L_LACt2r) se agrega en seccion gram+ anterior
modelos_lacL_gp = {Bifidobacterium_animalis_lactis_PT33, Lacticaseibacillus_paracasei_M38};
nombres_lacL_gp = {'Bifidobacterium_animalis_lactis_PT33', 'Lacticaseibacillus_paracasei_M38'};

for i = 1:2
    m = modelos_lacL_gp{i};
    fprintf('\n%s (gram+, sin periplasma):\n', nombres_lacL_gp{i});

    if ~any(strcmp(m.mets, 'M_lac__L_e'))
        m = addMetabolite(m, 'M_lac__L_e', 'metName', 'L-Lactate_e', ...
            'metFormula', 'C3H5O3', 'Charge', -1);
        fprintf('  Agregado metabolito: M_lac__L_e\n');
    end

    if ~any(strcmp(m.rxns, 'R_EX_lac__L_e'))
        m = addReaction(m, 'R_EX_lac__L_e', ...
            'reactionName', 'L-Lactate exchange', ...
            'metaboliteList', {'M_lac__L_e'}, ...
            'stoichCoeffList', [-1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Exchange/demand reaction');
        fprintf('  Agregada: R_EX_lac__L_e\n');
    else
        fprintf('  [OK] R_EX_lac__L_e ya existe\n');
    end

    modelos_lacL_gp{i} = m;
end

Bifidobacterium_animalis_lactis_PT33 = modelos_lacL_gp{1};
Lacticaseibacillus_paracasei_M38 = modelos_lacL_gp{2};

% --- Bt (gram-negativa): transporte via periplasma (e<->p<->c) ---
m = Bacteroides_thetaiotaomicron_VPI_5482;
fprintf('\nBacteroides_thetaiotaomicron_VPI_5482 (gram-neg, con periplasma):\n');

if ~any(strcmp(m.mets, 'M_lac__L_p'))
    m = addMetabolite(m, 'M_lac__L_p', 'metName', 'L-Lactate_p', ...
        'metFormula', 'C3H5O3', 'Charge', -1);
    fprintf('  Agregado metabolito: M_lac__L_p\n');
end
if ~any(strcmp(m.mets, 'M_lac__L_e'))
    m = addMetabolite(m, 'M_lac__L_e', 'metName', 'L-Lactate_e', ...
        'metFormula', 'C3H5O3', 'Charge', -1);
    fprintf('  Agregado metabolito: M_lac__L_e\n');
end

if ~any(strcmp(m.rxns, 'R_L_LACt2rpp'))
    m = addReaction(m, 'R_L_LACt2rpp', ...
        'reactionName', 'L-lactate reversible transport via proton symport (periplasm)', ...
        'metaboliteList', {'M_lac__L_c', 'M_lac__L_p', 'M_h_c', 'M_h_p'}, ...
        'stoichCoeffList', [1, -1, 1, -1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Transport, Inner Membrane');
    fprintf('  Agregada: R_L_LACt2rpp\n');
else
    fprintf('  [OK] R_L_LACt2rpp ya existe\n');
end

if ~any(strcmp(m.rxns, 'R_L_LACtex'))
    m = addReaction(m, 'R_L_LACtex', ...
        'reactionName', 'L-Lactate transport via diffusion (extracellular to periplasm)', ...
        'metaboliteList', {'M_lac__L_e', 'M_lac__L_p'}, ...
        'stoichCoeffList', [-1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Transport, Outer Membrane');
    fprintf('  Agregada: R_L_LACtex\n');
else
    fprintf('  [OK] R_L_LACtex ya existe\n');
end

if ~any(strcmp(m.rxns, 'R_EX_lac__L_e'))
    m = addReaction(m, 'R_EX_lac__L_e', ...
        'reactionName', 'L-Lactate exchange', ...
        'metaboliteList', {'M_lac__L_e'}, ...
        'stoichCoeffList', [-1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Exchange/demand reaction');
    fprintf('  Agregada: R_EX_lac__L_e\n');
else
    fprintf('  [OK] R_EX_lac__L_e ya existe\n');
end

Bacteroides_thetaiotaomicron_VPI_5482 = m;

%% CURACION EXPORTACION DE D-LACTATO
% =========================================================================
% Modelos que necesitan exportacion de D-lactato:
%   - PT33 (Bifidobacterium): productor de lactato, falta export D
%   - HGF2 (Clostridium): falta export D
%   - M62_1 (Clostridium): falta export D
%   - M38 (Lacticaseibacillus): productor de lactato, falta export D
%
% Modelos OK: Bt (ya tiene D-lactato export), Csym (tiene todo)
%
% Reacciones de referencia de iML1515 (E. coli, BiGG):
%   R_EX_lac__D_e:  M_lac__D_e <->                              (exchange)
%   R_D_LACtex:     M_lac__D_e <-> M_lac__D_p                   (difusion membrana ext)
%   R_D_LACt2pp:    M_lac__D_p + M_h_p <-> M_lac__D_c + M_h_c  (simporte proton)
% =========================================================================

fprintf('\n=== CURACION EXPORTACION DE D-LACTATO ===\n');

modelos_lacD = {Bifidobacterium_animalis_lactis_PT33, Clostridium_sp_HGF2, ...
                Clostridium_sp_M62_1, Lacticaseibacillus_paracasei_M38};
nombres_lacD = {'Bifidobacterium_animalis_lactis_PT33', 'Clostridium_sp_HGF2', ...
                'Clostridium_sp_M62_1', 'Lacticaseibacillus_paracasei_M38'};

for i = 1:4
    m = modelos_lacD{i};
    fprintf('\n%s (gram+, sin periplasma):\n', nombres_lacD{i});

    % Gram+: solo agregar M_lac__D_e y R_EX_lac__D_e (sin _p, sin tex/tpp)
    % El transporte directo e<->c (R_D_LACt2) se agrega en seccion gram+ anterior
    if ~any(strcmp(m.mets, 'M_lac__D_e'))
        m = addMetabolite(m, 'M_lac__D_e', 'metName', 'D-Lactate_e', ...
            'metFormula', 'C3H5O3', 'Charge', -1);
        fprintf('  Agregado metabolito: M_lac__D_e\n');
    end

    % R_EX_lac__D_e: lac_D_e <-> (exchange)
    if ~any(strcmp(m.rxns, 'R_EX_lac__D_e'))
        m = addReaction(m, 'R_EX_lac__D_e', ...
            'reactionName', 'D-Lactate exchange', ...
            'metaboliteList', {'M_lac__D_e'}, ...
            'stoichCoeffList', [-1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Exchange/demand reaction');
        fprintf('  Agregada: R_EX_lac__D_e\n');
    else
        fprintf('  [OK] R_EX_lac__D_e ya existe\n');
    end

    modelos_lacD{i} = m;
end

Bifidobacterium_animalis_lactis_PT33 = modelos_lacD{1};
Clostridium_sp_HGF2 = modelos_lacD{2};
Clostridium_sp_M62_1 = modelos_lacD{3};
Lacticaseibacillus_paracasei_M38 = modelos_lacD{4};

%% VERIFICACION EXPORTACION ACETATO Y LACTATO
% =========================================================================
fprintf('\n--- VERIFICACION EXPORTACION ACETATO Y LACTATO ---\n');

todos_modelos = {Bacteroides_thetaiotaomicron_VPI_5482, ...
                 Bifidobacterium_animalis_lactis_PT33, ...
                 Clostridium_sp_HGF2, ...
                 Clostridium_sp_M62_1, ...
                 Clostridium_symbiosum_WAL_14673, ...
                 Lacticaseibacillus_paracasei_M38};
todos_nombres = {'Bt', 'PT33', 'HGF2', 'M62_1', 'Csym', 'M38'};

% Acetato: Bt usa tex+tpp (via periplasma); gram+ usa ACtr (directo e<->c)
% Bt = indice 1, gram+ = indices 2-6
fprintf('\n  ACETATO (export):\n');
for i = 1:6
    mv = todos_modelos{i};
    if i == 1  % Bt (gram-neg)
        rxns_check = {'R_ACt2rpp', 'R_ACtex', 'R_EX_ac_e', 'R_ACt2ipp'};
    else  % gram+
        rxns_check = {'R_ACtr', 'R_EX_ac_e'};
    end
    n_ok = sum(ismember(rxns_check, mv.rxns));
    if n_ok == length(rxns_check)
        fprintf('    [OK] %6s: %d/%d rxns acetato\n', todos_nombres{i}, n_ok, length(rxns_check));
    else
        falta = rxns_check(~ismember(rxns_check, mv.rxns));
        fprintf('    [!!] %6s: %d/%d rxns acetato FALTA: %s\n', todos_nombres{i}, n_ok, length(rxns_check), strjoin(falta, ', '));
    end
end

% L-Lactato: Bt usa tex+tpp; gram+ usa L_LACt2r (directo e<->c)
fprintf('\n  L-LACTATO (export):\n');
for i = 1:6
    mv = todos_modelos{i};
    if i == 1  % Bt
        rxns_check = {'R_L_LACt2rpp', 'R_L_LACtex', 'R_EX_lac__L_e'};
    else  % gram+
        rxns_check = {'R_L_LACt2r', 'R_EX_lac__L_e'};
    end
    n_ok = sum(ismember(rxns_check, mv.rxns));
    if n_ok == length(rxns_check)
        fprintf('    [OK] %6s: %d/%d rxns L-lactato\n', todos_nombres{i}, n_ok, length(rxns_check));
    else
        falta = rxns_check(~ismember(rxns_check, mv.rxns));
        fprintf('    [--] %6s: %d/%d rxns L-lactato FALTA: %s\n', todos_nombres{i}, n_ok, length(rxns_check), strjoin(falta, ', '));
    end
end

% D-Lactato: Bt usa tex+tpp; gram+ usa D_LACt2 (directo e<->c)
fprintf('\n  D-LACTATO (export):\n');
for i = 1:6
    mv = todos_modelos{i};
    if i == 1  % Bt
        rxns_check = {'R_D_LACt2pp', 'R_D_LACtex', 'R_EX_lac__D_e'};
    else  % gram+
        rxns_check = {'R_D_LACt2', 'R_EX_lac__D_e'};
    end
    n_ok = sum(ismember(rxns_check, mv.rxns));
    if n_ok == length(rxns_check)
        fprintf('    [OK] %6s: %d/%d rxns D-lactato\n', todos_nombres{i}, n_ok, length(rxns_check));
    else
        falta = rxns_check(~ismember(rxns_check, mv.rxns));
        fprintf('    [--] %6s: %d/%d rxns D-lactato FALTA: %s\n', todos_nombres{i}, n_ok, length(rxns_check), strjoin(falta, ', '));
    end
end

fprintf('\n==========================================\n');

%% SUCCINATO EXPORT: Bacteroides thetaiotaomicron
% =========================================================================
% El modelo CarveMe puede importar succinato (R_HMR_9610 Na+/succ e->c,
% R_SUCCt2_3pp 3H+/succ p->c) pero NO puede exportarlo (ambas irreversibles).
% Se agregan las reacciones equivalentes al modelo AGORA2, usando IDs BiGG:
%   R_SUCCtex: succ[e] <=> succ[p]  (difusion membrana externa)
%   R_SUCCtpp: h[p] + succ[p] <=> h[c] + succ[c]  (simporte 1H+, reversible)
% AGORA2 usa SUCCt2rpp (VMH) = SUCCtpp (BiGG) — MetaNetX: MNXR104622
% Ref: AGORA2 Bacteroides_thetaiotaomicron_VPI_5482, BiGG iJN678/iSynCJ816
% =========================================================================
fprintf('\n--- SUCCINATO EXPORT: Bacteroides_thetaiotaomicron ---\n');

m = Bacteroides_thetaiotaomicron_VPI_5482;
fprintf('\nBacteroides_thetaiotaomicron_VPI_5482:\n');

% R_SUCCtex: difusion membrana externa (e <=> p)
m = addReaction(m, 'R_SUCCtex', ...
    'reactionName', 'Succinate transport via diffusion (extracellular to periplasm)', ...
    'metaboliteList', {'M_succ_e', 'M_succ_p'}, ...
    'stoichCoeffList', [-1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Transport, outer membrane');

% R_SUCCtpp: simporte reversible con 1H+ (p <=> c) — BiGG ID (equiv. AGORA2 SUCCt2rpp)
m = addReaction(m, 'R_SUCCtpp', ...
    'reactionName', 'Succinate transport via proton symport, reversible (periplasm)', ...
    'metaboliteList', {'M_succ_p', 'M_h_p', 'M_succ_c', 'M_h_c'}, ...
    'stoichCoeffList', [-1, -1, 1, 1], ...
    'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
    'subSystem', 'Transport, inner membrane');

fprintf('  Agregadas: R_SUCCtex, R_SUCCtpp (transporte succinato reversible e<=>p<=>c)\n');
Bacteroides_thetaiotaomicron_VPI_5482 = m;

% Verificacion succinato (busca con y sin prefijo R_)
fprintf('\n  Verificacion succinato Bt:\n');
rxns_succ = {'R_SUCCtex', 'R_SUCCtpp', 'R_EX_succ_e'};
for j = 1:length(rxns_succ)
    rxn_id = rxns_succ{j};
    % Buscar con prefijo R_ y sin el
    rxn_id_alt = regexprep(rxn_id, '^R_', '');
    if any(strcmp(m.rxns, rxn_id)) || any(strcmp(m.rxns, rxn_id_alt))
        fprintf('    [OK] %s presente\n', rxn_id);
    else
        fprintf('    [!!] %s FALTA\n', rxn_id);
    end
end
fprintf('\n==========================================\n');

%% SUCCINATO EXPORT: Clostridium sp HGF2 (gram-positiva)
% =========================================================================
% HGF2 no tiene M_succ_e, ni exchange, ni transporte e<->c para succinato.
% Se usa la ruta directa e<->c de AGORA2 (C. innocuum 2959), sin periplasma:
%   R_EX_succ_e: succ[e] <=>                    (exchange)
%   R_SUCCt2r:   h[e] + succ[e] <=> h[c] + succ[c]  (simporte 1H+, reversible)
% Ambos son IDs BiGG validos (iJR904, iNF517, iCN718, etc.)
% Ref: AGORA2 Clostridium_innocuum_2959 (gram-positiva, 2 compartimentos)
% =========================================================================
fprintf('\n--- SUCCINATO EXPORT: Clostridium_sp_HGF2 (gram-positiva, ruta directa e<->c) ---\n');

m = Clostridium_sp_HGF2;
fprintf('\nClostridium_sp_HGF2:\n');

% Agregar metabolito M_succ_e si no existe
if ~any(strcmp(m.mets, 'M_succ_e'))
    m = addMetabolite(m, 'M_succ_e', ...
        'metName', 'Succinate', ...
        'metFormula', 'C4H4O4', ...
        'charge', -2, ...
        'csense', 'E');
    fprintf('  Agregado metabolito: M_succ_e\n');
else
    fprintf('  [OK] M_succ_e ya existe\n');
end

% R_EX_succ_e: exchange
if ~any(strcmp(m.rxns, 'R_EX_succ_e')) && ~any(strcmp(m.rxns, 'EX_succ_e'))
    m = addReaction(m, 'R_EX_succ_e', ...
        'reactionName', 'Succinate exchange', ...
        'metaboliteList', {'M_succ_e'}, ...
        'stoichCoeffList', [-1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Exchange/demand reaction');
    fprintf('  Agregada: R_EX_succ_e\n');
else
    fprintf('  [OK] R_EX_succ_e ya existe\n');
end

% R_SUCCt2r: transporte directo e <=> c con 1H+ (ruta gram-positiva AGORA2)
if ~any(strcmp(m.rxns, 'R_SUCCt2r'))
    m = addReaction(m, 'R_SUCCt2r', ...
        'reactionName', 'Succinate transport via proton symport', ...
        'metaboliteList', {'M_succ_e', 'M_h_e', 'M_succ_c', 'M_h_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Transport, inner membrane');
    fprintf('  Agregada: R_SUCCt2r (transporte directo e<=>c, 1H+ simporte)\n');
else
    fprintf('  [OK] R_SUCCt2r ya existe\n');
end

Clostridium_sp_HGF2 = m;

% Verificacion succinato HGF2
fprintf('\n  Verificacion succinato HGF2:\n');
rxns_succ_hgf2 = {'R_SUCCt2r', 'R_EX_succ_e'};
for j = 1:length(rxns_succ_hgf2)
    rxn_id = rxns_succ_hgf2{j};
    rxn_id_alt = regexprep(rxn_id, '^R_', '');
    if any(strcmp(m.rxns, rxn_id)) || any(strcmp(m.rxns, rxn_id_alt))
        fprintf('    [OK] %s presente\n', rxn_id);
    else
        fprintf('    [!!] %s FALTA\n', rxn_id);
    end
end
fprintf('\n==========================================\n');

%% CURACION EXPORTACION: Clostridium symbiosum (gram-positiva, ruta directa e<->c)
% =========================================================================
% En AGORA2, C. symbiosum tiene transporte directo e<->c reversible para:
%   Acetato:   R_ACtr     ac[e] <=> ac[c]              (difusion, reversible)
%   L-Lactato: R_L_LACt   lac_L[e] <=> lac_L[c]        (difusion, reversible)
%   D-Lactato: R_D_LACt2  h[e]+lac_D[e] <=> h[c]+lac_D[c] (1H+, reversible)
%   Succinato: R_SUCCt2r  h[e]+succ[e] <=> h[c]+succ[c]   (1H+, reversible)
% En CarveMe, las reacciones de transporte existen pero son import-only (LB=0).
% Se agregan transportadores reversibles y se corrigen bounds.
% Ref: AGORA2 Clostridium_symbiosum_WAL_14673
% =========================================================================
fprintf('\n--- CURACION EXPORTACION: Clostridium_symbiosum (gram-positiva) ---\n');

m = Clostridium_symbiosum_WAL_14673;
fprintf('\nClostridium_symbiosum_WAL_14673:\n');

% 1) ACETATO: agregar R_ACtr (difusion directa e<=>c, sin H+) — BiGG valido
%    CarveMe tiene R_ACt2r (import-only) y R_ACt2r_1 (UB=0), no sirven para exportar
if ~any(strcmp(m.rxns, 'R_ACtr'))
    m = addReaction(m, 'R_ACtr', ...
        'reactionName', 'Acetate transporter', ...
        'metaboliteList', {'M_ac_e', 'M_ac_c'}, ...
        'stoichCoeffList', [-1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Transport, inner membrane');
    fprintf('  Agregada: R_ACtr (acetato, difusion directa e<=>c)\n');
else
    fprintf('  [OK] R_ACtr ya existe\n');
end

% 2) L-LACTATO: agregar R_L_LACt2r (simporte 1H+, reversible e<=>c) — BiGG valido
%    AGORA2 tiene R_L_LACt (difusion sin H+), pero ese ID no existe en BiGG.
%    Usamos R_L_LACt2r que es BiGG valido (iNF517, iYS854, iCN900, etc.)
%    CarveMe tiene R_L_LACt2r pero import-only (LB=0). Corregimos bounds.
rxn_idx = find(strcmp(m.rxns, 'R_L_LACt2r'));
if ~isempty(rxn_idx)
    m = changeRxnBounds(m, 'R_L_LACt2r', -1000, 'l');
    fprintf('  Corregida: R_L_LACt2r bounds -> [-1000, 1000] (L-lactato, ahora reversible)\n');
else
    m = addReaction(m, 'R_L_LACt2r', ...
        'reactionName', 'L lactate reversible transport via proton symport', ...
        'metaboliteList', {'M_lac__L_e', 'M_h_e', 'M_lac__L_c', 'M_h_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Transport, inner membrane');
    fprintf('  Agregada: R_L_LACt2r (L-lactato, 1H+ simporte reversible e<=>c)\n');
end

% 3) D-LACTATO: R_D_LACt2 existe pero import-only (LB=0). Hacer reversible.
%    En AGORA2 es reversible (-1000 a 1000). BiGG: D_LACt2 valido.
rxn_idx = find(strcmp(m.rxns, 'R_D_LACt2'));
if ~isempty(rxn_idx)
    m = changeRxnBounds(m, 'R_D_LACt2', -1000, 'l');
    fprintf('  Corregida: R_D_LACt2 bounds -> [-1000, 1000] (D-lactato, ahora reversible)\n');
else
    fprintf('  [!!] R_D_LACt2 no encontrada\n');
end

% 4) SUCCINATO: R_SUCCt2r existe pero import-only (LB=0). Hacer reversible.
%    En AGORA2 es reversible (-1000 a 1000). BiGG: SUCCt2r valido.
rxn_idx = find(strcmp(m.rxns, 'R_SUCCt2r'));
if ~isempty(rxn_idx)
    m = changeRxnBounds(m, 'R_SUCCt2r', -1000, 'l');
    fprintf('  Corregida: R_SUCCt2r bounds -> [-1000, 1000] (succinato, ahora reversible)\n');
else
    fprintf('  [!!] R_SUCCt2r no encontrada\n');
end

Clostridium_symbiosum_WAL_14673 = m;

% Verificacion Csym
fprintf('\n  Verificacion Csym exportacion:\n');
rxns_check_csym = {'R_ACtr', 'R_EX_ac_e', 'R_L_LACt2r', 'R_EX_lac__L_e', ...
                    'R_D_LACt2', 'R_EX_lac__D_e', 'R_SUCCt2r', 'R_EX_succ_e'};
for j = 1:length(rxns_check_csym)
    rxn_id = rxns_check_csym{j};
    rxn_id_alt = regexprep(rxn_id, '^R_', '');
    if any(strcmp(m.rxns, rxn_id)) || any(strcmp(m.rxns, rxn_id_alt))
        fprintf('    [OK] %s\n', rxn_id);
    else
        fprintf('    [!!] %s FALTA\n', rxn_id);
    end
end

fprintf('\n==========================================\n');

%% CURACION EXPORTACION: Bacteroides thetaiotaomicron (gram-negativa, ruta e<->p<->c)
% =========================================================================
% Bt es gram-negativa, usa 3 compartimentos. Problemas detectados:
% 1) Acetato: tiene ACtex (e<->p, reversible) y ACt2ipp (c->p, UB=0 bloqueado).
%    ACt2ipp: ac_c + h_c -> ac_p + h_p con UB=0 (no puede exportar c->p).
%    Solucion: cambiar UB de ACt2ipp a 1000 para permitir exportacion.
%    ACt2ipp es BiGG valido (iAF987). Con ACtex reversible, la ruta queda completa.
% 2) D-Lactato: tiene D_LACt2 (e->c directo, import-only LB=0) y EX_lac__D_e.
%    AGORA2 usa ruta periplasmica: D_LACtex (e<->p) + D_LACt2pp (p<->c).
%    Solucion: agregar D_LACtex + D_LACt2pp como en AGORA2/iML1515.
% Ref: AGORA2 Bacteroides_thetaiotaomicron_VPI_5482, BiGG iML1515
% =========================================================================
fprintf('\n--- CURACION EXPORTACION: Bacteroides_thetaiotaomicron (gram-negativa) ---\n');

m = Bacteroides_thetaiotaomicron_VPI_5482;
fprintf('\nBacteroides_thetaiotaomicron_VPI_5482:\n');

% 1) ACETATO: desbloquear R_ACt2ipp (cambiar UB de 0 a 1000)
%    ACt2ipp: ac_c + h_c -> ac_p + h_p. Con UB=1000, permite exportar c->p.
%    Combinado con ACtex (e<->p reversible) y EX_ac_e, la ruta queda completa.
rxn_idx = find(strcmp(m.rxns, 'R_ACt2ipp'));
if isempty(rxn_idx)
    rxn_idx = find(strcmp(m.rxns, 'ACt2ipp'));
end
if ~isempty(rxn_idx)
    m = changeRxnBounds(m, m.rxns{rxn_idx}, 1000, 'u');
    fprintf('  Corregida: %s UB -> 1000 (acetato, export c->p desbloqueado)\n', m.rxns{rxn_idx});
else
    fprintf('  [!!] ACt2ipp no encontrada\n');
end

% 2) D-LACTATO: agregar ruta periplasmica AGORA2
%    Agregar M_lac__D_p si no existe
if ~any(strcmp(m.mets, 'M_lac__D_p'))
    m = addMetabolite(m, 'M_lac__D_p', ...
        'metName', 'D-Lactate', ...
        'metFormula', 'C3H5O3', ...
        'charge', -1, ...
        'csense', 'E');
    fprintf('  Agregado metabolito: M_lac__D_p\n');
else
    fprintf('  [OK] M_lac__D_p ya existe\n');
end

%    R_D_LACtex: difusion membrana externa (e <=> p) — BiGG valido (iML1515)
if ~any(strcmp(m.rxns, 'R_D_LACtex'))
    m = addReaction(m, 'R_D_LACtex', ...
        'reactionName', 'D-lactate transport via diffusion (extracellular to periplasm)', ...
        'metaboliteList', {'M_lac__D_e', 'M_lac__D_p'}, ...
        'stoichCoeffList', [-1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Transport, outer membrane');
    fprintf('  Agregada: R_D_LACtex (D-lactato, difusion e<=>p)\n');
else
    fprintf('  [OK] R_D_LACtex ya existe\n');
end

%    R_D_LACt2pp: simporte 1H+ membrana interna (p <=> c) — BiGG valido (iML1515)
if ~any(strcmp(m.rxns, 'R_D_LACt2pp'))
    m = addReaction(m, 'R_D_LACt2pp', ...
        'reactionName', 'D-lactate transport via proton symport (periplasm)', ...
        'metaboliteList', {'M_lac__D_p', 'M_h_p', 'M_lac__D_c', 'M_h_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Transport, inner membrane');
    fprintf('  Agregada: R_D_LACt2pp (D-lactato, 1H+ simporte reversible p<=>c)\n');
else
    fprintf('  [OK] R_D_LACt2pp ya existe\n');
end

Bacteroides_thetaiotaomicron_VPI_5482 = m;

% Verificacion Bt
fprintf('\n  Verificacion Bt exportacion:\n');
rxns_check_bt = {'R_ACtex', 'R_EX_ac_e', 'R_D_LACtex', 'R_D_LACt2pp', 'R_EX_lac__D_e'};
for j = 1:length(rxns_check_bt)
    rxn_id = rxns_check_bt{j};
    rxn_id_alt = regexprep(rxn_id, '^R_', '');
    if any(strcmp(m.rxns, rxn_id)) || any(strcmp(m.rxns, rxn_id_alt))
        fprintf('    [OK] %s\n', rxn_id);
    else
        fprintf('    [!!] %s FALTA\n', rxn_id);
    end
end

fprintf('\n==========================================\n');

%% CURACION SUCCINATO: Clostridium sp M62_1 (gram-positiva, ruta directa e<->c)
% =========================================================================
% M62_1 tiene SUCCtex (e<->p) y EX_succ_e, pero NO tiene transporte p<->c
% ni transporte directo e<->c para succinato.
% AGORA2 M62_1 usa ruta directa: R_SUCCt2r h[e]+succ[e] <=> h[c]+succ[c]
% Solucion: agregar R_SUCCt2r (gram-positiva, ruta directa e<->c)
% Ref: AGORA2 Clostridium_sp_M62_1
% =========================================================================
fprintf('\n--- CURACION SUCCINATO: Clostridium_sp_M62_1 (gram-positiva) ---\n');

m = Clostridium_sp_M62_1;
fprintf('\nClostridium_sp_M62_1:\n');

if ~any(strcmp(m.rxns, 'R_SUCCt2r'))
    m = addReaction(m, 'R_SUCCt2r', ...
        'reactionName', 'Succinate transport via proton symport', ...
        'metaboliteList', {'M_succ_e', 'M_h_e', 'M_succ_c', 'M_h_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Transport, inner membrane');
    fprintf('  Agregada: R_SUCCt2r (succinato, transporte directo e<=>c, 1H+ simporte)\n');
else
    fprintf('  [OK] R_SUCCt2r ya existe\n');
end

Clostridium_sp_M62_1 = m;

% Verificacion M62_1 succinato
fprintf('\n  Verificacion M62_1 succinato:\n');
rxns_check_m62 = {'R_SUCCt2r', 'R_EX_succ_e'};
for j = 1:length(rxns_check_m62)
    rxn_id = rxns_check_m62{j};
    rxn_id_alt = regexprep(rxn_id, '^R_', '');
    if any(strcmp(m.rxns, rxn_id)) || any(strcmp(m.rxns, rxn_id_alt))
        fprintf('    [OK] %s\n', rxn_id);
    else
        fprintf('    [!!] %s FALTA\n', rxn_id);
    end
end

fprintf('\n==========================================\n');

%% CURACION EXPORTACION ETANOL: todas las bacterias
% =========================================================================
% Todos los modelos tienen ALCD2x (produccion de etanol en citoplasma):
%   acald[c] + h[c] + nadh[c] <=> etoh[c] + nad[c]
% Pero la mayoria NO puede exportar etanol al medio extracelular.
%
% Estado actual:
%   Bt:    OK (ETOHtex e<->p, ETOHtrpp p<->c, EX_etoh_e)
%   HGF2:  tiene ETOHtrpp (p<->c) pero NO ETOHtex ni EX_etoh_e
%   M62_1: NO tiene transporte ni exchange de etanol
%   Csym:  NO tiene transporte ni exchange de etanol
%   PT33:  NO tiene transporte ni exchange de etanol
%   M38:   NO tiene transporte ni exchange de etanol
%
% Solucion gram-positivas: ruta directa e<->c (AGORA2 usa ETOHt2r)
%   R_ETOHt2r: etoh[e] + h[e] <=> etoh[c] + h[c] (symport H+, reversible)
%   R_EX_etoh_e: etoh[e] <=> (exchange)
% Ref: AGORA2 HGF2, M62_1, Csym, L.paracasei LC2W, BB-12
% =========================================================================
fprintf('\n=== CURACION EXPORTACION ETANOL ===\n');

% --- Gram-positivas: transporte directo e<->c ---
grampositivas_etoh = {Bifidobacterium_animalis_lactis_PT33, Clostridium_sp_HGF2, ...
                      Clostridium_sp_M62_1, Clostridium_symbiosum_WAL_14673, ...
                      Lacticaseibacillus_paracasei_M38};
nombres_gp_etoh = {'PT33', 'HGF2', 'M62_1', 'Csym', 'M38'};

for i = 1:5
    m = grampositivas_etoh{i};
    fprintf('\n%s:\n', nombres_gp_etoh{i});

    % Asegurar metabolito M_etoh_e
    if ~any(strcmp(m.mets, 'M_etoh_e'))
        m = addMetabolite(m, 'M_etoh_e', 'metName', 'Ethanol', ...
            'metFormula', 'C2H6O', 'Charge', 0);
        fprintf('  Agregado metabolito: M_etoh_e\n');
    end

    % R_ETOHt2r: transporte directo e<->c con H+ (AGORA2)
    if ~any(strcmp(m.rxns, 'R_ETOHt2r'))
        m = addReaction(m, 'R_ETOHt2r', ...
            'reactionName', 'Ethanol reversible transport via proton symport (direct e<->c)', ...
            'metaboliteList', {'M_etoh_e', 'M_h_e', 'M_etoh_c', 'M_h_c'}, ...
            'stoichCoeffList', [-1, -1, 1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Transport');
        fprintf('  Agregada: R_ETOHt2r (etanol directo e<->c, symport H+)\n');
    else
        fprintf('  [OK] R_ETOHt2r ya existe\n');
    end

    % R_EX_etoh_e: exchange
    if ~any(strcmp(m.rxns, 'R_EX_etoh_e')) && ~any(strcmp(m.rxns, 'EX_etoh_e'))
        m = addReaction(m, 'R_EX_etoh_e', ...
            'reactionName', 'Ethanol exchange', ...
            'metaboliteList', {'M_etoh_e'}, ...
            'stoichCoeffList', [-1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Exchange/demand reaction');
        fprintf('  Agregada: R_EX_etoh_e\n');
    else
        fprintf('  [OK] R_EX_etoh_e ya existe\n');
    end

    grampositivas_etoh{i} = m;
end

Bifidobacterium_animalis_lactis_PT33 = grampositivas_etoh{1};
Clostridium_sp_HGF2 = grampositivas_etoh{2};
Clostridium_sp_M62_1 = grampositivas_etoh{3};
Clostridium_symbiosum_WAL_14673 = grampositivas_etoh{4};
Lacticaseibacillus_paracasei_M38 = grampositivas_etoh{5};

% --- Bt (gram-negativa): ya tiene ETOHtex + ETOHtrpp + EX_etoh_e ---
fprintf('\nBt: [OK] ya tiene ETOHtex + ETOHtrpp + EX_etoh_e\n');

fprintf('\n==========================================\n');

%% CURACION EXPORTACION FORMIATO: todas las bacterias
% =========================================================================
% Todos los modelos tienen PFL (pyruvate formate lyase):
%   coa[c] + pyr[c] -> accoa[c] + for[c]
% Pero la mayoria NO puede exportar formiato al medio extracelular.
%
% Estado actual:
%   Bt:    FORtex (e<->p) + EX_for_e, pero FALTA FORtpp (p<->c)
%   HGF2:  NO tiene transporte ni exchange de formiato (solo for_c y for_p)
%   M62_1: FORt2 (import-only e->c) + FORtppi (c->p irreversible) + EX_for_e
%   Csym:  FORt2 (import-only e->c) + EX_for_e
%   PT33:  NO tiene transporte ni exchange de formiato
%   M38:   NO tiene transporte ni exchange de formiato
%
% Solucion gram-positivas: ruta directa e<->c (AGORA2 usa FORt)
%   R_FORt: for[e] <=> for[c] (difusion simple, reversible)
%   R_EX_for_e: for[e] <=> (exchange, si no existe)
% Solucion Bt (gram-neg): agregar FORtppi (p<->c) — ID correcto BiGG
%   R_FORtppi: for[p] <=> for[c] (difusion periplasma, reversible)
%   NOTA: VMH usa "FORtpp", BiGG usa "FORtppi" (ambos MNXR99620)
% Ref: AGORA2 todos los modelos
% =========================================================================
fprintf('\n=== CURACION EXPORTACION FORMIATO ===\n');

% --- Gram-positivas: transporte directo e<->c ---
grampositivas_for = {Bifidobacterium_animalis_lactis_PT33, Clostridium_sp_HGF2, ...
                     Clostridium_sp_M62_1, Clostridium_symbiosum_WAL_14673, ...
                     Lacticaseibacillus_paracasei_M38};
nombres_gp_for = {'PT33', 'HGF2', 'M62_1', 'Csym', 'M38'};

for i = 1:5
    m = grampositivas_for{i};
    fprintf('\n%s:\n', nombres_gp_for{i});

    % Asegurar metabolito M_for_e
    if ~any(strcmp(m.mets, 'M_for_e'))
        m = addMetabolite(m, 'M_for_e', 'metName', 'Formate', ...
            'metFormula', 'CHO2', 'Charge', -1);
        fprintf('  Agregado metabolito: M_for_e\n');
    end

    % R_FORt: transporte directo e<->c (difusion simple, AGORA2)
    if ~any(strcmp(m.rxns, 'R_FORt'))
        m = addReaction(m, 'R_FORt', ...
            'reactionName', 'Formate transport via diffusion (direct e<->c)', ...
            'metaboliteList', {'M_for_e', 'M_for_c'}, ...
            'stoichCoeffList', [-1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Transport');
        fprintf('  Agregada: R_FORt (formiato directo e<->c, difusion)\n');
    else
        fprintf('  [OK] R_FORt ya existe\n');
    end

    % R_EX_for_e: exchange
    if ~any(strcmp(m.rxns, 'R_EX_for_e')) && ~any(strcmp(m.rxns, 'EX_for_e'))
        m = addReaction(m, 'R_EX_for_e', ...
            'reactionName', 'Formate exchange', ...
            'metaboliteList', {'M_for_e'}, ...
            'stoichCoeffList', [-1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Exchange/demand reaction');
        fprintf('  Agregada: R_EX_for_e\n');
    else
        fprintf('  [OK] R_EX_for_e ya existe\n');
    end

    grampositivas_for{i} = m;
end

Bifidobacterium_animalis_lactis_PT33 = grampositivas_for{1};
Clostridium_sp_HGF2 = grampositivas_for{2};
Clostridium_sp_M62_1 = grampositivas_for{3};
Clostridium_symbiosum_WAL_14673 = grampositivas_for{4};
Lacticaseibacillus_paracasei_M38 = grampositivas_for{5};

% --- Bt (gram-negativa): agregar FORtppi (p<->c), ya tiene FORtex + EX_for_e ---
% NOTA: VMH usa "FORtpp", BiGG usa "FORtppi" (MNXR99620)
fprintf('\nBt:\n');
m = Bacteroides_thetaiotaomicron_VPI_5482;

% Asegurar metabolito M_for_p
if ~any(strcmp(m.mets, 'M_for_p'))
    m = addMetabolite(m, 'M_for_p', 'metName', 'Formate', ...
        'metFormula', 'CHO2', 'Charge', -1);
    fprintf('  Agregado metabolito: M_for_p\n');
end

% R_FORtppi: difusion periplasma p<->c (BiGG ID correcto, MNXR99620)
if ~any(strcmp(m.rxns, 'R_FORtppi')) && ~any(strcmp(m.rxns, 'R_FORtpp'))
    m = addReaction(m, 'R_FORtppi', ...
        'reactionName', 'Formate transport via diffusion (periplasm to cytoplasm)', ...
        'metaboliteList', {'M_for_p', 'M_for_c'}, ...
        'stoichCoeffList', [-1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Transport, Inner Membrane');
    fprintf('  Agregada: R_FORtppi (formiato p<->c, difusion periplasma)\n');
else
    fprintf('  [OK] R_FORtppi ya existe\n');
end
fprintf('  [OK] FORtex (e<->p) y EX_for_e ya existen\n');

Bacteroides_thetaiotaomicron_VPI_5482 = m;

% --- Verificacion etanol y formiato ---
fprintf('\n--- VERIFICACION EXPORTACION ETANOL Y FORMIATO ---\n');
todos_modelos_ef = {Bacteroides_thetaiotaomicron_VPI_5482, ...
                    Bifidobacterium_animalis_lactis_PT33, ...
                    Clostridium_sp_HGF2, ...
                    Clostridium_sp_M62_1, ...
                    Clostridium_symbiosum_WAL_14673, ...
                    Lacticaseibacillus_paracasei_M38};
todos_nombres_ef = {'Bt', 'PT33', 'HGF2', 'M62_1', 'Csym', 'M38'};

% Etanol
fprintf('\n  ETANOL (export):\n');
for i = 1:6
    mv = todos_modelos_ef{i};
    % Buscar produccion
    tiene_alcd = any(strcmp(mv.rxns, 'R_ALCD2x')) || any(strcmp(mv.rxns, 'ALCD2x'));
    % Buscar transporte (directo o via periplasma)
    tiene_transp = any(strcmp(mv.rxns, 'R_ETOHt2r')) || ...
                   (any(strcmp(mv.rxns, 'R_ETOHtex')) && any(strcmp(mv.rxns, 'R_ETOHtrpp')));
    % Buscar exchange
    tiene_ex = any(strcmp(mv.rxns, 'R_EX_etoh_e')) || any(strcmp(mv.rxns, 'EX_etoh_e'));
    if tiene_alcd && tiene_transp && tiene_ex
        fprintf('    [OK] %6s: produccion + transporte + exchange\n', todos_nombres_ef{i});
    else
        fprintf('    [!!] %6s: ALCD=%d transp=%d EX=%d\n', todos_nombres_ef{i}, tiene_alcd, tiene_transp, tiene_ex);
    end
end

% Formiato
fprintf('\n  FORMIATO (export):\n');
for i = 1:6
    mv = todos_modelos_ef{i};
    % Buscar produccion (PFL)
    tiene_pfl = any(strcmp(mv.rxns, 'R_PFL')) || any(strcmp(mv.rxns, 'PFL'));
    % Buscar transporte (directo o via periplasma)
    tiene_transp = any(strcmp(mv.rxns, 'R_FORt')) || ...
                   (any(strcmp(mv.rxns, 'R_FORtex')) && (any(strcmp(mv.rxns, 'R_FORtppi')) || any(strcmp(mv.rxns, 'R_FORtpp'))));
    % Buscar exchange
    tiene_ex = any(strcmp(mv.rxns, 'R_EX_for_e')) || any(strcmp(mv.rxns, 'EX_for_e'));
    if tiene_pfl && tiene_transp && tiene_ex
        fprintf('    [OK] %6s: produccion + transporte + exchange\n', todos_nombres_ef{i});
    else
        fprintf('    [!!] %6s: PFL=%d transp=%d EX=%d\n', todos_nombres_ef{i}, tiene_pfl, tiene_transp, tiene_ex);
    end
end

fprintf('\n==========================================\n');

%% CURACION METABOLISMO CENTRAL: reacciones faltantes vs AGORA2
% =========================================================================
% Comparacion sistematica CarveMe vs AGORA2 (metabolismo central).
% Solo se agregan reacciones que existen en el homologo AGORA2.
% Todos los IDs verificados contra BiGG database y MetaNetX.
% =========================================================================
fprintf('\n=== CURACION METABOLISMO CENTRAL (vs AGORA2) ===\n');

% --- Bt: Bacteroides thetaiotaomicron ---
% Faltantes vs AGORA2: SUCDi, FRD2, FRD7, FRUt3 (p<->c fructosa)
% NOTA: ACONTa/ACONTb no se agregan porque Bt ya tiene CS+ACONT+ICDHyr+AKGDH
% NOTA: ALCD19 es gliceraldehido<->glicerol (NO etanol), no relevante
fprintf('\n--- Bt: reacciones faltantes metabolismo central ---\n');
m = Bacteroides_thetaiotaomicron_VPI_5482;

% R_SUCDi: Succinate dehydrogenase irreversible (EC 1.3.5.1)
% q8_c + succ_c -> fum_c + q8h2_c (BiGG: SUCDi, MNXR99641)
if ~any(strcmp(m.rxns, 'R_SUCDi'))
    m = addReaction(m, 'R_SUCDi', ...
        'reactionName', 'Succinate dehydrogenase (ubiquinone-8)', ...
        'metaboliteList', {'M_q8_c', 'M_succ_c', 'M_fum_c', 'M_q8h2_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', false, 'lowerBound', 0, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_SUCDi (succinato DH, ubiquinona)\n');
else
    fprintf('  [OK] R_SUCDi ya existe\n');
end

% R_FRD2: Fumarate reductase (menaquinone) (EC 1.3.5.4)
% fum_c + mql8_c <=> mqn8_c + succ_c (BiGG: FRD2, MNXR198970)
if ~any(strcmp(m.rxns, 'R_FRD2'))
    m = addReaction(m, 'R_FRD2', ...
        'reactionName', 'Fumarate reductase (menaquinol-8)', ...
        'metaboliteList', {'M_fum_c', 'M_mql8_c', 'M_mqn8_c', 'M_succ_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_FRD2 (fumarato reductasa, menaquinona)\n');
else
    fprintf('  [OK] R_FRD2 ya existe\n');
end

% R_FRD7: Fumarate reductase (ubiquinone) (EC 1.3.5.4)
% fum_c + q8h2_c <=> q8_c + succ_c (BiGG: FRD7, MNXR99641)
if ~any(strcmp(m.rxns, 'R_FRD7'))
    m = addReaction(m, 'R_FRD7', ...
        'reactionName', 'Fumarate reductase (ubiquinol-8)', ...
        'metaboliteList', {'M_fum_c', 'M_q8h2_c', 'M_q8_c', 'M_succ_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_FRD7 (fumarato reductasa, ubiquinona)\n');
else
    fprintf('  [OK] R_FRD7 ya existe\n');
end

% R_FRUt3: Fructose transport p<->c con H+ (BiGG: FRUt3)
% h_p + fru_p <=> h_c + fru_c
% NOTA: VMH usa FRUt2rpp, BiGG usa FRUt3 (MNXR190257)
if ~any(strcmp(m.rxns, 'R_FRUt3')) && ~any(strcmp(m.rxns, 'R_FRUt2rpp'))
    m = addReaction(m, 'R_FRUt3', ...
        'reactionName', 'Fructose transport via proton symport (periplasm)', ...
        'metaboliteList', {'M_h_p', 'M_fru_p', 'M_h_c', 'M_fru_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Transport, Inner Membrane');
    fprintf('  Agregada: R_FRUt3 (fructosa p<->c, symport H+)\n');
else
    fprintf('  [OK] R_FRUt3/FRUt2rpp ya existe\n');
end

Bacteroides_thetaiotaomicron_VPI_5482 = m;

% --- PT33: Bifidobacterium animalis ---
% Faltantes vs AGORA2 BB-12: PPCK, FUM, FRD7
% NOTA: PFK REMOVIDA - Bifidobacterium NO tiene PFK en genoma, usa bifid shunt
%       (PKETF/PKETX, gen RVCYZ000630) como ruta principal
% NOTA: ACONTa/ACONTb no se agregan, PT33 ya tiene CS+ACONT+ICDHyr
% NOTA: SUCDi no se agrega, PT33 ya tiene SUCD1
fprintf('\n--- PT33: reacciones faltantes metabolismo central ---\n');
m = Bifidobacterium_animalis_lactis_PT33;

% R_PPCK: PEP carboxykinase (EC 4.1.1.32)
% atp_c + oaa_c <=> adp_c + co2_c + pep_c (BiGG: PPCK, MNXR103099)
if ~any(strcmp(m.rxns, 'R_PPCK'))
    m = addReaction(m, 'R_PPCK', ...
        'reactionName', 'Phosphoenolpyruvate carboxykinase (ATP)', ...
        'metaboliteList', {'M_atp_c', 'M_oaa_c', 'M_adp_c', 'M_co2_c', 'M_pep_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Pyruvate metabolism');
    fprintf('  Agregada: R_PPCK (PEP carboxikinasa)\n');
else
    fprintf('  [OK] R_PPCK ya existe\n');
end

% R_FUM: Fumarase (EC 4.2.1.2)
% fum_c + h2o_c <=> mal_L_c (BiGG: FUM, MNXR99705)
if ~any(strcmp(m.rxns, 'R_FUM'))
    m = addReaction(m, 'R_FUM', ...
        'reactionName', 'Fumarase', ...
        'metaboliteList', {'M_fum_c', 'M_h2o_c', 'M_mal__L_c'}, ...
        'stoichCoeffList', [-1, -1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_FUM (fumarasa)\n');
else
    fprintf('  [OK] R_FUM ya existe\n');
end

% R_FRD7: Fumarate reductase (ubiquinone)
% fum_c + q8h2_c <=> q8_c + succ_c (BiGG: FRD7)
if ~any(strcmp(m.rxns, 'R_FRD7'))
    m = addReaction(m, 'R_FRD7', ...
        'reactionName', 'Fumarate reductase (ubiquinol-8)', ...
        'metaboliteList', {'M_fum_c', 'M_q8h2_c', 'M_q8_c', 'M_succ_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_FRD7 (fumarato reductasa, ubiquinona)\n');
else
    fprintf('  [OK] R_FRD7 ya existe\n');
end

Bifidobacterium_animalis_lactis_PT33 = m;

% --- HGF2: Clostridium sp. HGF2 ---
% Faltantes vs AGORA2: FRD2, FRD7
% NOTA: PPPNt2r es 3-fenilpropanoato (NO propionato), no relevante
fprintf('\n--- HGF2: reacciones faltantes metabolismo central ---\n');
m = Clostridium_sp_HGF2;

% R_FRD2: Fumarate reductase (menaquinone)
if ~any(strcmp(m.rxns, 'R_FRD2'))
    m = addReaction(m, 'R_FRD2', ...
        'reactionName', 'Fumarate reductase (menaquinol-8)', ...
        'metaboliteList', {'M_fum_c', 'M_mql8_c', 'M_mqn8_c', 'M_succ_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_FRD2 (fumarato reductasa, menaquinona)\n');
else
    fprintf('  [OK] R_FRD2 ya existe\n');
end

% R_FRD7: Fumarate reductase (ubiquinone)
if ~any(strcmp(m.rxns, 'R_FRD7'))
    m = addReaction(m, 'R_FRD7', ...
        'reactionName', 'Fumarate reductase (ubiquinol-8)', ...
        'metaboliteList', {'M_fum_c', 'M_q8h2_c', 'M_q8_c', 'M_succ_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_FRD7 (fumarato reductasa, ubiquinona)\n');
else
    fprintf('  [OK] R_FRD7 ya existe\n');
end

Clostridium_sp_HGF2 = m;

% --- M62_1: Clostridium sp. M62_1 ---
% Faltantes vs AGORA2: SUCD1, FRD7, ACALD
fprintf('\n--- M62_1: reacciones faltantes metabolismo central ---\n');
m = Clostridium_sp_M62_1;

% R_SUCD1: Succinate dehydrogenase (FAD)
% fad_c + succ_c <=> fadh2_c + fum_c (BiGG: SUCD1, MNXR82874)
if ~any(strcmp(m.rxns, 'R_SUCD1'))
    m = addReaction(m, 'R_SUCD1', ...
        'reactionName', 'Succinate dehydrogenase', ...
        'metaboliteList', {'M_fad_c', 'M_succ_c', 'M_fadh2_c', 'M_fum_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_SUCD1 (succinato DH, FAD)\n');
else
    fprintf('  [OK] R_SUCD1 ya existe\n');
end

% R_FRD7: Fumarate reductase (ubiquinone)
if ~any(strcmp(m.rxns, 'R_FRD7'))
    m = addReaction(m, 'R_FRD7', ...
        'reactionName', 'Fumarate reductase (ubiquinol-8)', ...
        'metaboliteList', {'M_fum_c', 'M_q8h2_c', 'M_q8_c', 'M_succ_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_FRD7 (fumarato reductasa, ubiquinona)\n');
else
    fprintf('  [OK] R_FRD7 ya existe\n');
end

% R_ACALD: Acetaldehyde dehydrogenase, acetylating (EC 1.2.1.10)
% acald_c + coa_c + nad_c <=> accoa_c + h_c + nadh_c (BiGG: ACALD, MNXR190509)
if ~any(strcmp(m.rxns, 'R_ACALD'))
    m = addReaction(m, 'R_ACALD', ...
        'reactionName', 'Acetaldehyde dehydrogenase (acetylating)', ...
        'metaboliteList', {'M_acald_c', 'M_coa_c', 'M_nad_c', 'M_accoa_c', 'M_h_c', 'M_nadh_c'}, ...
        'stoichCoeffList', [-1, -1, -1, 1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Pyruvate metabolism');
    fprintf('  Agregada: R_ACALD (acetaldehido DH, acetilante)\n');
else
    fprintf('  [OK] R_ACALD ya existe\n');
end

Clostridium_sp_M62_1 = m;

% --- Csym: Clostridium symbiosum ---
% Faltantes vs AGORA2: PPCK, SUCD1, FRD7
% NOTA: PPPNt2r es 3-fenilpropanoato, no relevante
fprintf('\n--- Csym: reacciones faltantes metabolismo central ---\n');
m = Clostridium_symbiosum_WAL_14673;

% R_PPCK: PEP carboxykinase
if ~any(strcmp(m.rxns, 'R_PPCK'))
    m = addReaction(m, 'R_PPCK', ...
        'reactionName', 'Phosphoenolpyruvate carboxykinase (ATP)', ...
        'metaboliteList', {'M_atp_c', 'M_oaa_c', 'M_adp_c', 'M_co2_c', 'M_pep_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Pyruvate metabolism');
    fprintf('  Agregada: R_PPCK (PEP carboxikinasa)\n');
else
    fprintf('  [OK] R_PPCK ya existe\n');
end

% R_SUCD1: Succinate dehydrogenase (FAD)
if ~any(strcmp(m.rxns, 'R_SUCD1'))
    m = addReaction(m, 'R_SUCD1', ...
        'reactionName', 'Succinate dehydrogenase', ...
        'metaboliteList', {'M_fad_c', 'M_succ_c', 'M_fadh2_c', 'M_fum_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_SUCD1 (succinato DH, FAD)\n');
else
    fprintf('  [OK] R_SUCD1 ya existe\n');
end

% R_FRD7: Fumarate reductase (ubiquinone)
if ~any(strcmp(m.rxns, 'R_FRD7'))
    m = addReaction(m, 'R_FRD7', ...
        'reactionName', 'Fumarate reductase (ubiquinol-8)', ...
        'metaboliteList', {'M_fum_c', 'M_q8h2_c', 'M_q8_c', 'M_succ_c'}, ...
        'stoichCoeffList', [-1, -1, 1, 1], ...
        'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
        'subSystem', 'Citric acid cycle');
    fprintf('  Agregada: R_FRD7 (fumarato reductasa, ubiquinona)\n');
else
    fprintf('  [OK] R_FRD7 ya existe\n');
end

Clostridium_symbiosum_WAL_14673 = m;

% --- Verificacion metabolismo central ---
fprintf('\n--- VERIFICACION METABOLISMO CENTRAL ---\n');
todos_mc = {Bacteroides_thetaiotaomicron_VPI_5482, ...
            Bifidobacterium_animalis_lactis_PT33, ...
            Clostridium_sp_HGF2, ...
            Clostridium_sp_M62_1, ...
            Clostridium_symbiosum_WAL_14673, ...
            Lacticaseibacillus_paracasei_M38};
nombres_mc = {'Bt', 'PT33', 'HGF2', 'M62_1', 'Csym', 'M38'};

rxns_check = {'R_CS', 'R_ACONT', 'R_ICDHyr', 'R_FUM', 'R_MDH', ...
              'R_SUCD1', 'R_SUCDi', 'R_FRD2', 'R_FRD7', ...
              'R_PPCK', ...
              'R_ACALD', 'R_ACS', ...
              'R_ILETA', 'R_LEUTA', 'R_VALTA'};

fprintf('  %8s', '');
for j = 1:length(rxns_check)
    fprintf(' %7s', strrep(rxns_check{j}, 'R_', ''));
end
fprintf('\n');
for i = 1:6
    mv = todos_mc{i};
    fprintf('  %8s', nombres_mc{i});
    for j = 1:length(rxns_check)
        if any(strcmp(mv.rxns, rxns_check{j}))
            fprintf('  %6s', 'OK');
        else
            fprintf('  %6s', '-');
        end
    end
    fprintf('\n');
end

fprintf('\n==========================================\n');

%% VERIFICAR CONSISTENCIA ESTEQUIOMÉTRICA
% modelo = Bifidobacterium_animalis_lactis_PT33;
% 
% % Encontrar reacciones inconsistentes
% [SConsistentMetBool, SConsistentRxnBool, SInConsistentMetBool, SInConsistentRxnBool] = ...
%     findStoichConsistentSubset(modelo);
% 
% fprintf('Metabolitos consistentes: %d/%d\n', sum(SConsistentMetBool), length(modelo.mets));
% fprintf('Reacciones consistentes: %d/%d\n', sum(SConsistentRxnBool), length(modelo.rxns));
% fprintf('Reacciones INCONSISTENTES: %d\n', sum(SInConsistentRxnBool));
% 
% % Ver las inconsistentes
% idx_incons = find(SInConsistentRxnBool);
% fprintf('\nPrimeras 10 reacciones inconsistentes:\n');
% for i = 1:min(10, length(idx_incons))
%     fprintf('%s - %s\n', modelo.rxns{idx_incons(i)}, modelo.rxnNames{idx_incons(i)});
% end

% % Guardar
% T = table(modelo.rxns(idx_incons), modelo.rxnNames(idx_incons), ...
%     'VariableNames', {'Reaccion', 'Nombre'});
% writetable(T, 'reacciones_inconsistentes.xlsx');

%% BUSCAR REACCIONES ASOCIADAS A UN GEN (utilidad comentada)
% modelo = Clostridium_sp_HGF2;
% gen_id = 'RCHGF000837';
% idx_gen = find(strcmp(modelo.genes, gen_id));
% if ~isempty(idx_gen)
%     rxns_asociadas = find(modelo.rxnGeneMat(:, idx_gen));
%     for i = 1:length(rxns_asociadas)
%         fprintf('%s - %s\n', modelo.rxns{rxns_asociadas(i)}, modelo.rxnNames{rxns_asociadas(i)});
%     end
% end

%% =========================================================================
%% TRANSPORTE VIA PERIPLASMA (e<->p<->c) PARA GRAM+ MODELOS
%% =========================================================================
% CarveMe asigna periplasma a TODOS los modelos uniformemente (convencion BiGG/iML1515).
% Esta seccion agrega los pares tex (e<->p) + pp (p<->c) para fru, but, lac_L, lac_D,
% succ y glc en los gram+, MANTENIENDO los transportadores directos e<->c que CarveMe
% ya pone (no se elimina nada, son rutas paralelas).
%
% IDs BiGG validados (todos en bigg_models_reactions.tsv):
%   FRUtex      fru_e <-> fru_p                        (difusion membrana exterior)
%   FRUt3       fru_p + h_p <-> fru_c + h_c            (symport 1H+, p<->c)
%   FRUptspp    pep_c + fru_p <-> f1p_c + pyr_c        (PTS, p<->c) -- ya en HGF2/M62/Csym/M38 via CarveMe
%   BUTtex      but_e <-> but_p                        (difusion)
%   BUTt2rpp    but_p + h_p <-> but_c + h_c            (symport 1H+ reversible)
%   L_LACtex    lac_L_e <-> lac_L_p                    (difusion)
%   L_LACt2rpp  lac_L_p + h_p <-> lac_L_c + h_c        (symport 1H+ reversible)
%   D_LACtex    lac_D_e <-> lac_D_p                    (difusion)
%   D_LACt2pp   lac_D_p + h_p <-> lac_D_c + h_c        (symport 1H+)
%   SUCCtex     succ_e <-> succ_p                      (difusion)
%   SUCCtpp     succ_p + h_p <-> succ_c + h_c          (symport 1H+ reversible)
%   GLCtex      glc_D_e <-> glc_D_p                    (difusion)
%   GLCt2pp     glc_D_p + h_p <-> glc_D_c + h_c        (symport 1H+, gram+ sin PTS)
%
% Justificacion paper: CarveMe (Machado 2018) usa periplasma como compartimento
% uniforme de modelado para gram+ y gram-. Las reacciones tex/pp son IDs BiGG
% canonicos presentes en iML1515, iJO1366, iAF692.
%
% Asignacion por bacteria:
%   PT33  (Bifidobacterium, sin PTS para fru/glc): FRUtex+FRUt3 + GLCtex+GLCt2pp
%                                                  + BUTtex+BUTt2rpp + LACtex+LACt2rpp (L/D)
%   HGF2  (gram+ con PTS): BUTtex+BUTt2rpp + LACtex+LACt2rpp (L/D) + SUCCtex+SUCCtpp
%   M62_1 (gram+, productor butirato): BUTtex+BUTt2rpp + LACtex+LACt2rpp (L/D)
%   Csym  (gram+ con PTS): BUTtex+BUTt2rpp + LACtex+LACt2rpp (L/D) + SUCCtex+SUCCtpp
%   M38   (LAB homofermentativa, solo lactato): LACtex+LACt2rpp (L/D)
%
% Kesto/inulin/kestopt/kestottr: NO se agrega tex/pp (BiGG no los tiene).
% Se mantienen en _e — la hidrolisis es extracelular por enzimas secretadas
% (KESTOASEe, INULINASEe), convencion AGORA2/VMH ya implementada arriba.
%
% Idempotencia: si una reaccion o metabolito _p ya existe, NO se duplica.
% =========================================================================

fprintf('\n\n=========================================================\n');
fprintf('  TRANSPORTE VIA PERIPLASMA PARA GRAM+ (paralelo a directo)\n');
fprintf('=========================================================\n');

% --- PT33 (Bifidobacterium, sin PTS para fru ni glc) ---
fprintf('\n--- PT33 ---\n');
m = Bifidobacterium_animalis_lactis_PT33;
m = addPeripPair(m, 'fru',    'FRUtex',   'FRUt3',      'D-fructose transport via diffusion (e<->p)',                'D-fructose transport via H+ symport (p<->c)');
m = addPeripPair(m, 'glc__D', 'GLCtex',   'GLCt2pp',    'D-glucose transport via diffusion (e<->p)',                 'D-glucose transport via H+ symport (p<->c)');
m = addPeripPair(m, 'but',    'BUTtex',   'BUTt2rpp',   'Butyrate transport via diffusion (e<->p)',                  'Butyrate transport via H+ symport, reversible (p<->c)');
m = addPeripPair(m, 'lac__L', 'L_LACtex', 'L_LACt2rpp', 'L-lactate transport via diffusion (e<->p)',                 'L-lactate transport via H+ symport, reversible (p<->c)');
m = addPeripPair(m, 'lac__D', 'D_LACtex', 'D_LACt2pp',  'D-lactate transport via diffusion (e<->p)',                 'D-lactate transport via H+ symport (p<->c)');
Bifidobacterium_animalis_lactis_PT33 = m;

% --- HGF2 (Clostridium, gram+ con PTS — fru ya tiene tex+ptspp en MCC) ---
fprintf('\n--- HGF2 ---\n');
m = Clostridium_sp_HGF2;
m = addPeripPair(m, 'but',    'BUTtex',   'BUTt2rpp',   'Butyrate transport via diffusion (e<->p)',                  'Butyrate transport via H+ symport, reversible (p<->c)');
m = addPeripPair(m, 'lac__L', 'L_LACtex', 'L_LACt2rpp', 'L-lactate transport via diffusion (e<->p)',                 'L-lactate transport via H+ symport, reversible (p<->c)');
m = addPeripPair(m, 'lac__D', 'D_LACtex', 'D_LACt2pp',  'D-lactate transport via diffusion (e<->p)',                 'D-lactate transport via H+ symport (p<->c)');
m = addPeripPair(m, 'succ',   'SUCCtex',  'SUCCtpp',    'Succinate transport via diffusion (e<->p)',                 'Succinate transport via H+ symport, reversible (p<->c)');
Clostridium_sp_HGF2 = m;

% --- M62_1 (Clostridium, gram+, productor butirato; sin evidencia para succ) ---
fprintf('\n--- M62_1 ---\n');
m = Clostridium_sp_M62_1;
m = addPeripPair(m, 'but',    'BUTtex',   'BUTt2rpp',   'Butyrate transport via diffusion (e<->p)',                  'Butyrate transport via H+ symport, reversible (p<->c)');
m = addPeripPair(m, 'lac__L', 'L_LACtex', 'L_LACt2rpp', 'L-lactate transport via diffusion (e<->p)',                 'L-lactate transport via H+ symport, reversible (p<->c)');
m = addPeripPair(m, 'lac__D', 'D_LACtex', 'D_LACt2pp',  'D-lactate transport via diffusion (e<->p)',                 'D-lactate transport via H+ symport (p<->c)');
Clostridium_sp_M62_1 = m;

% --- Csym (Clostridium symbiosum, gram+ con PTS) ---
fprintf('\n--- Csym ---\n');
m = Clostridium_symbiosum_WAL_14673;
m = addPeripPair(m, 'but',    'BUTtex',   'BUTt2rpp',   'Butyrate transport via diffusion (e<->p)',                  'Butyrate transport via H+ symport, reversible (p<->c)');
m = addPeripPair(m, 'lac__L', 'L_LACtex', 'L_LACt2rpp', 'L-lactate transport via diffusion (e<->p)',                 'L-lactate transport via H+ symport, reversible (p<->c)');
m = addPeripPair(m, 'lac__D', 'D_LACtex', 'D_LACt2pp',  'D-lactate transport via diffusion (e<->p)',                 'D-lactate transport via H+ symport (p<->c)');
m = addPeripPair(m, 'succ',   'SUCCtex',  'SUCCtpp',    'Succinate transport via diffusion (e<->p)',                 'Succinate transport via H+ symport, reversible (p<->c)');
Clostridium_symbiosum_WAL_14673 = m;

% --- M38 (Lacticaseibacillus paracasei, LAB homofermentativa) ---
fprintf('\n--- M38 ---\n');
m = Lacticaseibacillus_paracasei_M38;
m = addPeripPair(m, 'lac__L', 'L_LACtex', 'L_LACt2rpp', 'L-lactate transport via diffusion (e<->p)',                 'L-lactate transport via H+ symport, reversible (p<->c)');
m = addPeripPair(m, 'lac__D', 'D_LACtex', 'D_LACt2pp',  'D-lactate transport via diffusion (e<->p)',                 'D-lactate transport via H+ symport (p<->c)');
Lacticaseibacillus_paracasei_M38 = m;

fprintf('\n=========================================================\n');
fprintf('  TRANSPORTE VIA PERIPLASMA PARA GRAM+ AGREGADO\n');
fprintf('=========================================================\n');

%% =========================================================================
%% SETEAR MEDIO ZMB + ANAEROBIOSIS EN TODOS LOS MODELOS
%% =========================================================================
fprintf('\n\n=======================================================\n');
fprintf('  SETEANDO MEDIO ZMB + GLUCOSA + ANAEROBIOSIS\n');
fprintf('=======================================================\n');

Bacteroides_thetaiotaomicron_VPI_5482 = setearMedioZMB(Bacteroides_thetaiotaomicron_VPI_5482, mZMB, rates_mZMB, 'glucosa');
Bifidobacterium_animalis_lactis_PT33  = setearMedioZMB(Bifidobacterium_animalis_lactis_PT33, mZMB, rates_mZMB, 'glucosa');
Clostridium_sp_HGF2                  = setearMedioZMB(Clostridium_sp_HGF2, mZMB, rates_mZMB, 'glucosa');
Clostridium_sp_M62_1                 = setearMedioZMB(Clostridium_sp_M62_1, mZMB, rates_mZMB, 'glucosa');
Clostridium_symbiosum_WAL_14673      = setearMedioZMB(Clostridium_symbiosum_WAL_14673, mZMB, rates_mZMB, 'glucosa');
Lacticaseibacillus_paracasei_M38     = setearMedioZMB(Lacticaseibacillus_paracasei_M38, mZMB, rates_mZMB, 'glucosa');

fprintf('  [OK] Medio ZMB + glucosa seteado en los 6 modelos\n');

% %% =========================================================================
% %% RESTRINGIR ATPM: valor de E. coli (iML1515, BiGG Models)
% %% =========================================================================
% % NOTA: Se probo restringir ATPM a 8.39 (E. coli iML1515, BiGG) pero con
% % ATPM capado, el solver no produce butirato porque la demanda de ATP es
% % demasiado baja y se satisface solo con acetato/lactato (rutas cortas).
% % Se deja ATPM en su valor original (LB=0, UB=1000) y se maximiza como
% % objetivo. ATPM=1000 es un artefacto de ciclos internos, pero los flujos
% % de fermentacion (butirato, acetato, lactato) son cualitativamente correctos.
% % Referencia: AGORA2 tambien usa ATPM LB=0, UB=1000 sin restriccion.

%% =========================================================================
%% SIMULACIONES: fijar mu experimental, maximizar ATPM + pFBA L2
%% =========================================================================

fprintf('\n\n=======================================================\n');
fprintf('  SIMULACIONES: mu fijo + loopless max ATPM + pFBA L2\n');
fprintf('=======================================================\n');

nombres_sim = {'Bt', 'PT33', 'HGF2', 'M62_1', 'Csym', 'M38'};
modelos_sim = {Bacteroides_thetaiotaomicron_VPI_5482, ...
               Bifidobacterium_animalis_lactis_PT33, ...
               Clostridium_sp_HGF2, ...
               Clostridium_sp_M62_1, ...
               Clostridium_symbiosum_WAL_14673, ...
               Lacticaseibacillus_paracasei_M38};
mu_exp_sim = [0.648; 0.357; 0.33; 0.141; 0.403; 0.21];

condiciones = {'glucosa', 'fructosa'};

for c = 1:length(condiciones)
    cond = condiciones{c};
    fprintf('\n=== %s ===\n', upper(cond));

    % ---- PASO 1: FBA normal (max ATPM, con ciclos) ----
    fprintf('\n--- PASO 1: FBA normal (max ATPM, con ciclos) ---\n');
    fprintf('%-8s | %8s | %8s | %10s\n', 'Modelo', 'mu_max', 'mu_real', 'ATPM_raw');
    fprintf('%s\n', repmat('-', 1, 45));

    sol_raw_all = cell(6,1);
    m_all = cell(6,1);

    for i = 1:6
        m = modelos_sim{i};
        m = setearMedioZMB(m, mZMB, rates_mZMB, cond);
        % Biomasa relajada: crece hasta un maximo de mu_exp (no forzada exacta)
        m = changeRxnBounds(m, 'Growth_biomass', 0, 'l');
        m = changeRxnBounds(m, 'Growth_biomass', mu_exp_sim(i), 'u');

        atpm_idx = find(strcmp(m.rxns, 'R_ATPM') | strcmp(m.rxns, 'ATPM'));
        if ~isempty(atpm_idx)
            m = changeObjective(m, m.rxns{atpm_idx(1)}, 1);
        end

        sol_raw = optimizeCbModel(m, 'max');

        % Obtener mu real alcanzado
        bio_idx = find(strcmp(m.rxns, 'Growth_biomass'));
        if sol_raw.stat == 1 && ~isempty(atpm_idx)
            mu_real = sol_raw.v(bio_idx(1));
            fprintf('%-8s | %8.4f | %8.4f | %10.4f\n', nombres_sim{i}, mu_exp_sim(i), mu_real, sol_raw.v(atpm_idx(1)));
        elseif sol_raw.stat == 1
            mu_real = sol_raw.v(bio_idx(1));
            fprintf('%-8s | %8.4f | %8.4f | %10s\n', nombres_sim{i}, mu_exp_sim(i), mu_real, 'NO ATPM');
        else
            fprintf('%-8s | %8.4f | %8s | %10s\n', nombres_sim{i}, mu_exp_sim(i), 'INFACT', 'INFACT');
        end

        sol_raw_all{i} = sol_raw;
        m_all{i} = m;
    end

    % ---- PASO 2: cycleFreeFlux (eliminar ciclos del vector de flujos) ----
    % Primero eliminar metabolitos huerfanos (sin reacciones) que causan
    % filas de ceros en la matriz S y hacen fallar findStoichConsistentSubset
    for i = 1:6
        m = m_all{i};
        orphans = find(~any(m.S ~= 0, 2));  % filas con todos ceros
        if ~isempty(orphans)
            fprintf('  %s: eliminando %d metabolitos huerfanos\n', nombres_sim{i}, length(orphans));
            m = removeMetabolites(m, m.mets(orphans), false);
            m_all{i} = m;
            % Actualizar vector de flujos (solo se eliminaron mets, no rxns)
        end
    end

    fprintf('\n--- PASO 2: cycleFreeFlux (elimina ciclos internos) ---\n');
    fprintf('%-8s | %10s | %10s\n', 'Modelo', 'ATPM_raw', 'ATPM_CF');
    fprintf('%s\n', repmat('-', 1, 35));

    sol_cf_all = cell(6,1);

    for i = 1:6
        m = m_all{i};
        atpm_idx = find(strcmp(m.rxns, 'R_ATPM') | strcmp(m.rxns, 'ATPM'));

        if sol_raw_all{i}.stat == 1
            % cycleFreeFlux: toma el vector de flujos y elimina loops internos
            v_cf = cycleFreeFlux(sol_raw_all{i}.v, m.c, m);
            sol_cf = sol_raw_all{i};
            sol_cf.v = v_cf;

            if ~isempty(atpm_idx)
                fprintf('%-8s | %10.4f | %10.4f\n', nombres_sim{i}, ...
                    sol_raw_all{i}.v(atpm_idx(1)), v_cf(atpm_idx(1)));
            else
                fprintf('%-8s | %10s | %10s\n', nombres_sim{i}, 'NO ATPM', 'NO ATPM');
            end
        else
            sol_cf = sol_raw_all{i};
            fprintf('%-8s | %10s | %10s\n', nombres_sim{i}, 'INFACT', 'INFACT');
        end

        sol_cf_all{i} = sol_cf;
    end

    % ---- PASO 3: pFBA L2 con ATPM fijado al valor cycle-free ----
    fprintf('\n--- PASO 3: pFBA L2 (ATPM fijado al valor cycle-free) ---\n');
    fprintf('%-8s | %10s | %10s\n', 'Modelo', 'ATPM_CF', 'stat_L2');
    fprintf('%s\n', repmat('-', 1, 35));

    for i = 1:6
        m = m_all{i};
        atpm_idx = find(strcmp(m.rxns, 'R_ATPM') | strcmp(m.rxns, 'ATPM'));

        if ~isempty(atpm_idx) && sol_cf_all{i}.stat == 1
            atpm_cf = sol_cf_all{i}.v(atpm_idx(1));
            m = changeRxnBounds(m, m.rxns{atpm_idx(1)}, atpm_cf, 'b');
            sol = optimizeCbModel(m, 'max', 1e-6);
            if sol.stat == 1
                fprintf('%-8s | %10.4f | %10s\n', nombres_sim{i}, atpm_cf, 'OK');
            else
                fprintf('%-8s | %10.4f | %10s\n', nombres_sim{i}, atpm_cf, 'INFACT');
            end
        else
            sol = sol_cf_all{i};
            fprintf('%-8s | %10s | %10s\n', nombres_sim{i}, '---', 'SKIP');
        end
    end

    % ---- Flujos clave: RAW vs CYCLE-FREE vs L2 ----
    fprintf('\n--- Flujos clave: RAW vs CF vs L2 ---\n');
    fprintf('%-12s', 'Metabolito');
    for i = 1:6
        fprintf(' | %8s', nombres_sim{i});
    end
    fprintf('\n%s\n', repmat('-', 1, 12 + 6*11));

    mets_ll = {'R_EX_ac_e','EX_ac_e','Acetato'; ...
               'R_EX_lac__L_e','EX_lac__L_e','L-Lactato'; ...
               'R_EX_but_e','EX_but_e','Butirato'; ...
               'R_EX_succ_e','EX_succ_e','Succinato'; ...
               'R_EX_glc__D_e','EX_glc__D_e','Glucosa'};
    % Solo mostrar cycle-free (el mas relevante)
    for j = 1:size(mets_ll,1)
        fprintf('%-12s', mets_ll{j,3});
        for i = 1:6
            sol_i = sol_cf_all{i};
            m_i = m_all{i};
            idx = find(strcmp(m_i.rxns, mets_ll{j,1}) | strcmp(m_i.rxns, mets_ll{j,2}));
            if ~isfield(sol_i, 'v') || isempty(sol_i.v)
                fprintf(' | %8s', 'INF');
            elseif ~isempty(idx)
                v = sol_i.v(idx(1));
                if abs(v) < 1e-6, fprintf(' | %8s', '0');
                else, fprintf(' | %8.2f', v); end
            else
                fprintf(' | %8s', '---');
            end
        end
        fprintf('\n');
    end
    fprintf('(flujos de cycleFreeFlux, sin ciclos internos)\n');
end

%% =========================================================================
%% FBA PARA EXPLORAR: guardar soluciones en workspace
%% =========================================================================
% Guarda en workspace: fba_Bt, fba_PT33, etc. (soluciones FBA)
% y tabla_Bt, tabla_PT33, etc. (tablas con flujos != 0)
% Condicion: ZMB + glucosa, mu fijo al experimental, max ATPM + pFBA L2

fprintf('\n\n=======================================================\n');
fprintf('  FBA PARA EXPLORAR EN WORKSPACE (cycleFreeFlux + L2)\n');
fprintf('=======================================================\n');

for i = 1:6
    m = modelos_sim{i};
    m = setearMedioZMB(m, mZMB, rates_mZMB, 'glucosa');
    % Biomasa relajada: crece hasta un maximo de mu_exp (no forzada exacta)
    m = changeRxnBounds(m, 'Growth_biomass', 0, 'l');
    m = changeRxnBounds(m, 'Growth_biomass', mu_exp_sim(i), 'u');

    % Paso 1: FBA normal max ATPM
    atpm_idx = find(strcmp(m.rxns, 'R_ATPM') | strcmp(m.rxns, 'ATPM'));
    if ~isempty(atpm_idx)
        m = changeObjective(m, m.rxns{atpm_idx(1)}, 1);
    end
    sol_raw = optimizeCbModel(m, 'max');

    % Eliminar metabolitos huerfanos antes de cycleFreeFlux
    orphans = find(~any(m.S ~= 0, 2));
    if ~isempty(orphans)
        m = removeMetabolites(m, m.mets(orphans), false);
    end

    % Paso 2: cycleFreeFlux — eliminar ciclos del vector de flujos
    if sol_raw.stat == 1
        v_cf = cycleFreeFlux(sol_raw.v, m.c, m);
        sol_cf = sol_raw;
        sol_cf.v = v_cf;
    else
        sol_cf = sol_raw;
    end

    % Guardar solucion cycle-free
    eval(['fba_cf_' nombres_sim{i} ' = sol_cf;']);

    if ~isempty(atpm_idx) && sol_raw.stat == 1
        fprintf('  [CF] %s: ATPM raw=%.1f -> CF=%.4f\n', ...
            nombres_sim{i}, sol_raw.v(atpm_idx(1)), v_cf(atpm_idx(1)));
    end

    % Paso 3: Fijar ATPM al valor cycle-free + pFBA L2
    if ~isempty(atpm_idx) && sol_cf.stat == 1
        m = changeRxnBounds(m, m.rxns{atpm_idx(1)}, sol_cf.v(atpm_idx(1)), 'b');
    end
    sol = optimizeCbModel(m, 'max', 1e-6);  % L2 con ATPM fijo

    % Tabla con flujos != 0
    idx_nz = find(abs(sol.v) > 1e-6);
    t = table(m.rxns(idx_nz), sol.v(idx_nz), 'VariableNames', {'rxn', 'flux'});
    t = sortrows(t, 'flux');

    % Guardar en workspace con nombre corto
    eval(['fba_' nombres_sim{i} ' = sol;']);
    eval(['tabla_' nombres_sim{i} ' = t;']);
    eval(['modelo_' nombres_sim{i} ' = m;']);

    fprintf('  [OK] fba_%s, tabla_%s, modelo_%s (%d rxns con flujo)\n', ...
        nombres_sim{i}, nombres_sim{i}, nombres_sim{i}, height(t));
end

%% =========================================================================
%% PRODUCCION DE METABOLITOS CLAVE (butirato, acetato, lactato, succinato, etc)
%% =========================================================================
% Flujo positivo en EX_ = exportacion (produccion)
% Flujo negativo en EX_ = importacion (consumo)

fprintf('\n\n=======================================================\n');
fprintf('  PRODUCCION DE METABOLITOS CLAVE (glucosa, ZMB)\n');
fprintf('=======================================================\n');

% Metabolitos de interes: {ID exchange con R_, ID sin R_, nombre}
mets_clave = {
    'R_EX_but_e',      'EX_but_e',      'Butirato';
    'R_EX_ac_e',       'EX_ac_e',       'Acetato';
    'R_EX_lac__L_e',   'EX_lac__L_e',   'L-Lactato';
    'R_EX_lac__D_e',   'EX_lac__D_e',   'D-Lactato';
    'R_EX_succ_e',     'EX_succ_e',     'Succinato';
    'R_EX_for_e',      'EX_for_e',      'Formiato';
    'R_EX_etoh_e',     'EX_etoh_e',     'Etanol';
    'R_EX_co2_e',      'EX_co2_e',      'CO2';
    'R_EX_h2_e',       'EX_h2_e',       'H2';
    'R_EX_ppa_e',      'EX_ppa_e',      'Propionato';
    'R_EX_glc__D_e',   'EX_glc__D_e',   'Glucosa';
};

% Encabezado
fprintf('\n%-12s', 'Metabolito');
for i = 1:6
    fprintf(' | %8s', nombres_sim{i});
end
fprintf('\n%s\n', repmat('-', 1, 12 + 6*11));

% Para cada metabolito, buscar flujo en cada modelo
for j = 1:size(mets_clave, 1)
    fprintf('%-12s', mets_clave{j,3});
    for i = 1:6
        sol = eval(['fba_' nombres_sim{i}]);
        m = eval(['modelo_' nombres_sim{i}]);
        idx = find(strcmp(m.rxns, mets_clave{j,1}) | strcmp(m.rxns, mets_clave{j,2}));
        if isempty(sol.v)
            fprintf(' | %8s', 'INF');
        elseif ~isempty(idx)
            flux = sol.v(idx(1));
            if abs(flux) < 1e-6
                fprintf(' | %8s', '0');
            else
                fprintf(' | %8.4f', flux);
            end
        else
            fprintf(' | %8s', '---');
        end
    end
    fprintf('\n');
end

fprintf('\n(+) = produce/exporta, (-) = consume/importa, (---) = no tiene la reaccion\n');

%% =========================================================================
%% COMPARACION: RAW vs CYCLE-FREE vs L2 (flujos clave y ATPM)
%% =========================================================================
fprintf('\n\n=======================================================\n');
fprintf('  COMPARACION: RAW vs CF vs L2\n');
fprintf('=======================================================\n');
fprintf('  RAW = FBA con ciclos, CF = cycleFreeFlux, L2 = pFBA L2 (ATPM=CF)\n\n');

% ATPM
fprintf('%-12s', 'ATPM');
for i = 1:6
    fprintf(' | %8s', nombres_sim{i});
end
fprintf('\n');
fprintf('%-12s', '  RAW');
for i = 1:6
    sol_cf_i = eval(['fba_cf_' nombres_sim{i}]);
    m_i = eval(['modelo_' nombres_sim{i}]);
    atpm_idx = find(strcmp(m_i.rxns, 'R_ATPM') | strcmp(m_i.rxns, 'ATPM'));
    % Buscar raw: fba_cf tiene .v ya con CF, necesitamos el raw
    % Lo reportamos del paso 1 de la seccion anterior (sol_raw guardado)
    % Pero aqui solo tenemos fba_cf y fba. Usemos la info del header.
    fprintf(' | %8s', '(ver P1)');
end
fprintf('\n');
fprintf('%-12s', '  CF');
for i = 1:6
    sol_cf_i = eval(['fba_cf_' nombres_sim{i}]);
    m_i = eval(['modelo_' nombres_sim{i}]);
    atpm_idx = find(strcmp(m_i.rxns, 'R_ATPM') | strcmp(m_i.rxns, 'ATPM'));
    if ~isempty(atpm_idx) && ~isempty(sol_cf_i.v)
        fprintf(' | %8.1f', sol_cf_i.v(atpm_idx(1)));
    else
        fprintf(' | %8s', '---');
    end
end
fprintf('\n');
fprintf('%-12s', '  L2');
for i = 1:6
    sol_l2_i = eval(['fba_' nombres_sim{i}]);
    m_i = eval(['modelo_' nombres_sim{i}]);
    atpm_idx = find(strcmp(m_i.rxns, 'R_ATPM') | strcmp(m_i.rxns, 'ATPM'));
    if ~isempty(atpm_idx) && ~isempty(sol_l2_i.v)
        fprintf(' | %8.1f', sol_l2_i.v(atpm_idx(1)));
    else
        fprintf(' | %8s', '---');
    end
end
fprintf('\n\n');

% Metabolitos clave: CF vs L2
fprintf('%-12s', 'Metabolito');
for i = 1:6
    fprintf(' | %17s', nombres_sim{i});
end
fprintf('\n%-12s', '');
for i = 1:6
    fprintf(' | %8s %8s', 'CF', 'L2');
end
fprintf('\n%s\n', repmat('-', 1, 12 + 6*20));

mets_comp = {'R_EX_ac_e','EX_ac_e','Acetato'; ...
             'R_EX_lac__L_e','EX_lac__L_e','L-Lactato'; ...
             'R_EX_but_e','EX_but_e','Butirato'; ...
             'R_EX_succ_e','EX_succ_e','Succinato'; ...
             'R_EX_glc__D_e','EX_glc__D_e','Glucosa'};
for j = 1:size(mets_comp,1)
    fprintf('%-12s', mets_comp{j,3});
    for i = 1:6
        sol_cf_i = eval(['fba_cf_' nombres_sim{i}]);
        sol_l2_i = eval(['fba_' nombres_sim{i}]);
        m_i = eval(['modelo_' nombres_sim{i}]);
        idx = find(strcmp(m_i.rxns, mets_comp{j,1}) | strcmp(m_i.rxns, mets_comp{j,2}));
        if isempty(idx)
            fprintf(' | %8s %8s', '---', '---');
        elseif isempty(sol_cf_i.v) || isempty(sol_l2_i.v)
            fprintf(' | %8s %8s', 'INF', 'INF');
        else
            v_cf = sol_cf_i.v(idx(1)); v_l2 = sol_l2_i.v(idx(1));
            if abs(v_cf)<1e-6, s_cf='0'; else, s_cf=sprintf('%.2f',v_cf); end
            if abs(v_l2)<1e-6, s_l2='0'; else, s_l2=sprintf('%.2f',v_l2); end
            fprintf(' | %8s %8s', s_cf, s_l2);
        end
    end
    fprintf('\n');
end

fprintf('\nCF = cycleFreeFlux (sin ciclos), L2 = pFBA L2 (ATPM fijado al CF)\n');
fprintf('Si ATPM_CF << ATPM_RAW, cycleFreeFlux elimino ciclos exitosamente.\n');

%% =========================================================================
%% DIAGNOSTICO Bt: verificar que exchanges se cerraron bien
%% =========================================================================
fprintf('\n--- DIAGNOSTICO Bt: exchanges con LB < 0 que NO estan en ZMB ---\n');
m_bt = modelo_Bt;
ex_bt = unique([m_bt.rxns(contains(m_bt.rxns, 'EX_')); m_bt.rxns(findExcRxns(m_bt))]);
zmb_alt = cellfun(@(x) regexprep(x, '^R_', ''), mZMB, 'UniformOutput', false);
zmb_all = [mZMB; zmb_alt; {'R_EX_h2o_e';'EX_h2o_e';'R_EX_co2_e';'EX_co2_e';'R_EX_h_e';'EX_h_e'}];
for j = 1:length(ex_bt)
    idx = find(strcmp(m_bt.rxns, ex_bt{j}));
    if m_bt.lb(idx) < -1e-6 && ~any(strcmp(zmb_all, ex_bt{j}))
        fprintf('  [ABIERTO] %s  LB=%.1f\n', ex_bt{j}, m_bt.lb(idx));
    end
end

%% =========================================================================
%% FVA: rango posible de produccion de metabolitos clave
%% =========================================================================
% FVA muestra el rango [min, max] de flujo posible manteniendo mu fijo
% IMPORTANTE: usamos Growth_biomass como objetivo (ya fijado por LB=UB=mu_exp)
% Asi el FVA explora todo el rango factible sin forzar ATPM al maximo

fprintf('\n\n=======================================================\n');
fprintf('  FVA: rango posible de produccion (mu fijo)\n');
fprintf('=======================================================\n');

fprintf('\n%-12s', 'Metabolito');
for i = 1:6
    fprintf(' | %14s', nombres_sim{i});
end
fprintf('\n%s\n', repmat('-', 1, 12 + 6*17));

mets_fva = {
    'R_EX_but_e',    'EX_but_e',    'Butirato';
    'R_EX_ac_e',     'EX_ac_e',     'Acetato';
    'R_EX_lac__L_e', 'EX_lac__L_e', 'L-Lactato';
    'R_EX_lac__D_e', 'EX_lac__D_e', 'D-Lactato';
    'R_EX_succ_e',   'EX_succ_e',   'Succinato';
    'R_EX_for_e',    'EX_for_e',    'Formiato';
};

for j = 1:size(mets_fva, 1)
    fprintf('%-12s', mets_fva{j,3});
    for i = 1:6
        m = eval(['modelo_' nombres_sim{i}]);

        % Cambiar objetivo a Growth_biomass para el FVA
        m = changeObjective(m, 'Growth_biomass', 1);

        % Buscar la reaccion (con o sin R_), usar solo la primera coincidencia
        idx = find(strcmp(m.rxns, mets_fva{j,1}) | strcmp(m.rxns, mets_fva{j,2}));
        if ~isempty(idx)
            idx = idx(1);  % evitar duplicados
            try
                [minF, maxF] = fluxVariability(m, 100, 'max', m.rxns(idx));
                fprintf(' | %6.1f/%6.1f', minF, maxF);
            catch
                fprintf(' | %14s', 'INFACT');
            end
        else
            fprintf(' | %14s', '---');
        end
    end
    fprintf('\n');
end
fprintf('\nFormato: min/max (positivo=exporta, si max>0 PUEDE producir)\n');

%% =========================================================================
%% DIAGNOSTICO PREFIJO: verificar si los modelos usan M_ o no
%% =========================================================================
fprintf('\n\n=======================================================\n');
fprintf('  DIAGNOSTICO PREFIJO METABOLITOS\n');
fprintf('=======================================================\n');
for i = 1:6
    m = eval(['modelo_' nombres_sim{i}]);
    % Mostrar los primeros 5 metabolitos
    fprintf('\n  %s (primeros 5 mets):\n', nombres_sim{i});
    for j = 1:min(5, length(m.mets))
        fprintf('    %s\n', m.mets{j});
    end
    % Contar cuantos empiezan con M_
    n_M = sum(startsWith(m.mets, 'M_'));
    fprintf('  Total mets: %d, con M_: %d, sin M_: %d\n', ...
        length(m.mets), n_M, length(m.mets) - n_M);
    % Verificar accoa_c
    has_M = any(strcmp(m.mets, 'M_accoa_c'));
    has_no = any(strcmp(m.mets, 'accoa_c'));
    fprintf('  M_accoa_c: %s (%d rxns) | accoa_c: %s (%d rxns)\n', ...
        string(has_M), nnz(m.S(strcmp(m.mets, 'M_accoa_c'), :)), ...
        string(has_no), nnz(m.S(strcmp(m.mets, 'accoa_c'), :)));
end

%% =========================================================================
%% DIAGNOSTICO BUTIRATO: conectividad de metabolitos en la ruta
%% =========================================================================
% Si un metabolito esta en solo 1 reaccion, es un dead-end (bloqueado)
fprintf('\n\n=======================================================\n');
fprintf('  DIAGNOSTICO BUTIRATO: conectividad metabolitos\n');
fprintf('=======================================================\n');

mets_ruta_but = {'M_accoa_c', 'M_aacoa_c', 'M_3hbcoa_c', 'M_b2coa_c', ...
                 'M_btcoa_c', 'M_but_c', 'M_but_e', ...
                 'M_fdxo_2_2_c', 'M_fdxrd_c', 'M_ac_c', 'M_acac_c', 'M_h2_c'};

for i = [3, 4, 5]  % HGF2, M62_1, Csym
    m = eval(['modelo_' nombres_sim{i}]);
    fprintf('\n  %s:\n', nombres_sim{i});
    for j = 1:length(mets_ruta_but)
        idx_met = find(strcmp(m.mets, mets_ruta_but{j}));
        if ~isempty(idx_met)
            nrxns = nnz(m.S(idx_met, :));
            if nrxns <= 1
                fprintf('    [DEAD-END] %s (%d rxn)\n', mets_ruta_but{j}, nrxns);
            else
                fprintf('    [OK] %s (%d rxns)\n', mets_ruta_but{j}, nrxns);
            end
        else
            % Buscar sin prefijo M_
            met_sin_M = regexprep(mets_ruta_but{j}, '^M_', '');
            idx_met2 = find(strcmp(m.mets, met_sin_M));
            if ~isempty(idx_met2)
                nrxns = nnz(m.S(idx_met2, :));
                fprintf('    [PREFIJO!] %s NO existe, pero %s SI (%d rxns)\n', ...
                    mets_ruta_but{j}, met_sin_M, nrxns);
            else
                fprintf('    [MISSING] %s\n', mets_ruta_but{j});
            end
        end
    end

    % Tambien verificar si las reacciones de la ruta pueden llevar flujo
    fprintf('  Flujo maximo posible por cada reaccion (obj=biomasa, mu libre):\n');
    m_test = m;
    m_test = changeObjective(m_test, 'Growth_biomass', 0);
    m_test.lb(strcmp(m_test.rxns, 'Growth_biomass')) = 0;
    m_test.ub(strcmp(m_test.rxns, 'Growth_biomass')) = 1000;
    rxns_check = {'R_ACACT1r', 'R_HACD1', 'R_ECOAH1', 'R_BTCOADH', ...
                  'R_FDNADOX_H', 'R_POR4', 'R_BUTCT', 'R_BUTt2r', 'R_EX_but_e'};
    for j = 1:length(rxns_check)
        ridx = find(strcmp(m_test.rxns, rxns_check{j}));
        if ~isempty(ridx)
            m_tmp = changeObjective(m_test, rxns_check{j}, 1);
            sol = optimizeCbModel(m_tmp, 'max');
            fprintf('    max(%s) = %.4f\n', rxns_check{j}, sol.f);
        end
    end
end

%% =========================================================================
%% DIAGNOSTICO PT33/M38: conectividad acetato/lactato
%% =========================================================================
fprintf('\n\n=======================================================\n');
fprintf('  DIAGNOSTICO PT33/M38: conectividad acetato/lactato\n');
fprintf('=======================================================\n');

mets_ac_lac = {'M_ac_c', 'M_ac_e', 'M_lac__L_c', 'M_lac__L_e', ...
               'M_lac__D_c', 'M_lac__D_e'};

for i = [2, 6]  % PT33, M38
    m = eval(['modelo_' nombres_sim{i}]);
    fprintf('\n  %s:\n', nombres_sim{i});
    for j = 1:length(mets_ac_lac)
        idx_met = find(strcmp(m.mets, mets_ac_lac{j}));
        if ~isempty(idx_met)
            nrxns = nnz(m.S(idx_met, :));
            if nrxns <= 1
                fprintf('    [DEAD-END] %s (%d rxn)\n', mets_ac_lac{j}, nrxns);
            else
                fprintf('    [OK] %s (%d rxns)\n', mets_ac_lac{j}, nrxns);
            end
        else
            met_sin_M = regexprep(mets_ac_lac{j}, '^M_', '');
            idx_met2 = find(strcmp(m.mets, met_sin_M));
            if ~isempty(idx_met2)
                nrxns = nnz(m.S(idx_met2, :));
                fprintf('    [PREFIJO!] %s NO existe, pero %s SI (%d rxns)\n', ...
                    mets_ac_lac{j}, met_sin_M, nrxns);
            else
                fprintf('    [MISSING] %s\n', mets_ac_lac{j});
            end
        end
    end

    % Test flujo maximo acetato y lactato
    fprintf('  Flujo maximo posible (obj=biomasa libre):\n');
    m_test = m;
    m_test = changeObjective(m_test, 'Growth_biomass', 0);
    m_test.lb(strcmp(m_test.rxns, 'Growth_biomass')) = 0;
    m_test.ub(strcmp(m_test.rxns, 'Growth_biomass')) = 1000;
    rxns_check = {'R_ACtr', 'R_EX_ac_e', 'R_L_LACt2r', 'R_EX_lac__L_e', ...
                  'R_D_LACt2', 'R_EX_lac__D_e'};
    for j = 1:length(rxns_check)
        ridx = find(strcmp(m_test.rxns, rxns_check{j}));
        if ~isempty(ridx)
            m_tmp = changeObjective(m_test, rxns_check{j}, 1);
            sol = optimizeCbModel(m_tmp, 'max');
            fprintf('    max(%s) = %.4f\n', rxns_check{j}, sol.f);
        end
    end
end

%% =========================================================================
%% CORRECCION DE REACCIONES CarveMe CON BUGS DE MASA/CARGA
%% =========================================================================
% Estas reacciones vienen desbalanceadas de CarveMe.
% MEMOTE las detecta y rompe stoichiometric consistency.
% Fixes:
%   CMCBTFL       -> falta 1 h_c como reactante (H:+1, charge:-3)
%   SALCHS4FEabcpp-> salchs4fe_c formula H47 deberia ser H46
%   DTBTt         -> h_c espurio como producto
%   UACCpts       -> 3 h_c espurios como productos
%   TAGabc        -> 1 h_c espurio como producto (2->1)

fprintf('\n=== Corrigiendo reacciones CarveMe con bugs de masa ===\n');

% Lista de modelos y sus nombres para iterar
modelos_fix = {Bifidobacterium_animalis_lactis_PT33, ...
               Clostridium_sp_HGF2, ...
               Clostridium_sp_M62_1, ...
               Clostridium_symbiosum_WAL_14673, ...
               Lacticaseibacillus_paracasei_M38};
nombres_fix = {'PT33', 'HGF2', 'M62_1', 'Csym', 'M38'};

for i = 1:length(modelos_fix)
    m = modelos_fix{i};
    nom = nombres_fix{i};

    % --- CMCBTFL: agregar 1 h_c como reactante ---
    % Reaccion: cmcbtt_c + fe3_e -> fcmcbtt_c (falta h_c)
    % Presente en: PT33, M62_1, Csym, M38
    ridx = find(strcmp(m.rxns, 'R_CMCBTFL') | strcmp(m.rxns, 'CMCBTFL'));
    if ~isempty(ridx)
        midx = find(strcmp(m.mets, 'M_h_c'));
        if ~isempty(midx)
            m.S(midx, ridx) = m.S(midx, ridx) - 1;  % agregar 1 h_c como reactante
            fprintf('  [%s] CMCBTFL: agregado h_c como reactante\n', nom);
        end
    end

    % --- SALCHS4FEabcpp: corregir formula de salchs4fe_c ---
    % salchs4fe_c tiene H47 pero deberia ser H46 (como salchs4fe_p)
    % Presente en: PT33, M38
    midx = find(strcmp(m.mets, 'M_salchs4fe_c'));
    if ~isempty(midx)
        old_formula = m.metFormulas{midx};
        if contains(old_formula, 'H47')
            m.metFormulas{midx} = strrep(old_formula, 'H47', 'H46');
            fprintf('  [%s] SALCHS4FEabcpp: salchs4fe_c formula H47->H46\n', nom);
        end
    end

    % --- DTBTt: quitar h_c espurio de productos ---
    % Reaccion tiene 1 h_c de mas como producto
    % Presente en: M62_1
    ridx = find(strcmp(m.rxns, 'R_DTBTt') | strcmp(m.rxns, 'DTBTt'));
    if ~isempty(ridx)
        midx = find(strcmp(m.mets, 'M_h_c'));
        if ~isempty(midx) && m.S(midx, ridx) > 0
            m.S(midx, ridx) = m.S(midx, ridx) - 1;  % quitar 1 h_c de productos
            fprintf('  [%s] DTBTt: quitado h_c espurio de productos\n', nom);
        end
    end

    % --- UACCpts: quitar 3 h_c espurios de productos ---
    % Reaccion tiene 3 h_c de mas como productos
    % Presente en: Csym
    ridx = find(strcmp(m.rxns, 'R_UACCpts') | strcmp(m.rxns, 'UACCpts'));
    if ~isempty(ridx)
        midx = find(strcmp(m.mets, 'M_h_c'));
        if ~isempty(midx) && m.S(midx, ridx) > 0
            m.S(midx, ridx) = m.S(midx, ridx) - 3;  % quitar 3 h_c de productos
            fprintf('  [%s] UACCpts: quitados 3 h_c espurios de productos\n', nom);
        end
    end

    % --- TAGabc: quitar 1 h_c espurio de productos ---
    % Reaccion tiene 2 h_c pero deberia tener 1
    % Presente en: HGF2 (baja prioridad, ya pasa MEMOTE, pero corregimos igual)
    ridx = find(strcmp(m.rxns, 'R_TAGabc') | strcmp(m.rxns, 'TAGabc'));
    if ~isempty(ridx)
        midx = find(strcmp(m.mets, 'M_h_c'));
        if ~isempty(midx) && m.S(midx, ridx) > 1
            m.S(midx, ridx) = m.S(midx, ridx) - 1;  % quitar 1 h_c de productos
            fprintf('  [%s] TAGabc: quitado 1 h_c espurio (2->1)\n', nom);
        end
    end

    % Guardar modelo corregido de vuelta
    modelos_fix{i} = m;
end

% Copiar modelos corregidos de vuelta a las variables originales
Bifidobacterium_animalis_lactis_PT33 = modelos_fix{1};
Clostridium_sp_HGF2 = modelos_fix{2};
Clostridium_sp_M62_1 = modelos_fix{3};
Clostridium_symbiosum_WAL_14673 = modelos_fix{4};
Lacticaseibacillus_paracasei_M38 = modelos_fix{5};

fprintf('=== Correcciones CarveMe completadas ===\n');

%% =========================================================================
%% EXPORTACION DE MODELOS FINALES (sin constraints de simulacion)
%% =========================================================================
% Se exportan los modelos curados (con reacciones agregadas, bounds corregidos)
% pero SIN las restricciones de simulacion (mu libre, ATPM libre, medio no forzado).
% Los modelos exportados son los originales tras curacion (pasos anteriores),
% NO los modelo_Bt/modelo_PT33 del workspace (esos tienen mu fijo + medio ZMB).

fprintf('\n\n=======================================================\n');
fprintf('  EXPORTACION DE MODELOS FINALES (SBML para MEMOTE)\n');
fprintf('=======================================================\n');

nombres_export = {'Bacteroides_thetaiotaomicron_VPI_5482', ...
                  'Bifidobacterium_animalis_lactis_PT33', ...
                  'Clostridium_sp_HGF2', ...
                  'Clostridium_sp_M62_1', ...
                  'Clostridium_symbiosum_WAL_14673', ...
                  'Lacticaseibacillus_paracasei_M38'};
modelos_export = {Bacteroides_thetaiotaomicron_VPI_5482, ...
                  Bifidobacterium_animalis_lactis_PT33, ...
                  Clostridium_sp_HGF2, ...
                  Clostridium_sp_M62_1, ...
                  Clostridium_symbiosum_WAL_14673, ...
                  Lacticaseibacillus_paracasei_M38};

output_dir = pwd;

for i = 1:6
    m = modelos_export{i};
    nombre = nombres_export{i};

    % Objetivo = Growth_biomass (sin fijar mu)
    m = changeObjective(m, 'Growth_biomass', 1);

    % Exportar .mat
    mat_file = fullfile(output_dir, [nombre '.mat']);
    eval([nombre ' = m;']);
    save(mat_file, nombre);
    fprintf('  [OK] %s.mat\n', nombre);

    % Exportar .xml (SBML, para MEMOTE)
    xml_file = fullfile(output_dir, [nombre '_curado.xml']);
    m.modelAnnotation = {};  % fix compatibilidad writeSBML
    writeCbModel(m, 'format', 'sbml', 'fileName', xml_file);
    fprintf('  [OK] %s_curado.xml\n', nombre);
end

fprintf('\n=== EXPORTACION COMPLETA ===\n');
fprintf('Archivos exportados en: %s\n', output_dir);
fprintf('Sufijo _curado.xml: modelos con reacciones agregadas y bounds corregidos\n');
fprintf('Para MEMOTE: memote run --filename report.html modelo_curado.xml\n');

%% =========================================================================
%% FUNCION LOCAL: setearMedioZMB
%% =========================================================================
% Funcion auxiliar para setear medio ZMB en un modelo
% Retorna el modelo con bounds configurados para anaerobiosis + ZMB
% carbon_source: 'glucosa', 'fructosa', o 'ambas'

function m = setearMedioZMB(m, mZMB, rates_mZMB, carbon_source)
    % 1) Cerrar TODAS las reacciones de exchange (LB=0, UB=1000)
    % Buscar por nombre (EX_) y por estructura (findExcRxns)
    ex_by_name = m.rxns(contains(m.rxns, 'EX_'));
    ex_by_struct = m.rxns(findExcRxns(m));
    ex_rxns = unique([ex_by_name; ex_by_struct]);
    for j = 1:length(ex_rxns)
        m = changeRxnBounds(m, ex_rxns{j}, 0, 'l');
        m = changeRxnBounds(m, ex_rxns{j}, 1000, 'u');
    end

    % 2) Abrir medio ZMB (importacion de nutrientes)
    for j = 1:length(mZMB)
        rxn_id = mZMB{j};
        rxn_id_alt = regexprep(rxn_id, '^R_', '');
        idx = find(strcmp(m.rxns, rxn_id) | strcmp(m.rxns, rxn_id_alt));
        if ~isempty(idx)
            m = changeRxnBounds(m, m.rxns{idx(1)}, -rates_mZMB(j), 'l');
        end
    end

    % 3) Fuente de carbono segun condicion
    rxn_glc = find(strcmp(m.rxns, 'R_EX_glc__D_e') | strcmp(m.rxns, 'EX_glc__D_e'));
    rxn_fru = find(strcmp(m.rxns, 'R_EX_fru_e') | strcmp(m.rxns, 'EX_fru_e'));

    % NOTA: Se fuerza consumo exacto de la fuente de carbono (LB=UB=-10)
    % porque con LB=-10,UB=1000 el modelo puede elegir no consumir glucosa
    % (Bt cataboliza aminoacidos en vez de glucosa si no se fuerza).
    % NOTA: La fuente de carbono activa se fuerza con LB=UB=-10 ('b').
    % La fuente INACTIVA solo se cierra para importacion (LB=0, 'l'),
    % manteniendo UB=1000 (del paso 1). Asi el exchange queda (0,1000)
    % en el .mat, permitiendo que SteadyCom cree IEX funcionales.
    % Antes se usaba 'b' (LB=UB=0) que bloqueaba completamente el exchange.
    switch carbon_source
        case 'glucosa'
            if ~isempty(rxn_glc)
                m = changeRxnBounds(m, m.rxns{rxn_glc(1)}, -10, 'b');
            end
            % fru: solo cerrar importacion (LB=0), UB=1000 ya esta del paso 1
            % (no usar 'b' que pone UB=0 y rompe IEX en SteadyCom)
        case 'fructosa'
            % glc: solo cerrar importacion (LB=0), UB=1000 ya esta del paso 1
            if ~isempty(rxn_fru)
                m = changeRxnBounds(m, m.rxns{rxn_fru(1)}, -10, 'b');
            end
        case 'ambas'
            if ~isempty(rxn_glc)
                m = changeRxnBounds(m, m.rxns{rxn_glc(1)}, -10, 'b');
            end
            if ~isempty(rxn_fru)
                m = changeRxnBounds(m, m.rxns{rxn_fru(1)}, -10, 'b');
            end
    end

    % 4) Anaerobiosis: cerrar O2
    rxn_o2 = find(strcmp(m.rxns, 'R_EX_o2_e') | strcmp(m.rxns, 'EX_o2_e'));
    if ~isempty(rxn_o2)
        m = changeRxnBounds(m, m.rxns{rxn_o2(1)}, 0, 'l');
    end

    % 5) H2O, CO2, H+: exchange abiertos bidireccional
    rxns_libres = {'R_EX_h2o_e', 'R_EX_co2_e', 'R_EX_h_e'};
    for j = 1:length(rxns_libres)
        rxn_id = rxns_libres{j};
        rxn_id_alt = regexprep(rxn_id, '^R_', '');
        idx = find(strcmp(m.rxns, rxn_id) | strcmp(m.rxns, rxn_id_alt));
        if ~isempty(idx)
            m = changeRxnBounds(m, m.rxns{idx(1)}, -1000, 'l');
        end
    end
end

%% =========================================================================
%% FUNCION LOCAL: addPeripPair
%% =========================================================================
% Agrega un par de transporte tex (e<->p, difusion) + pp (p<->c, symport 1H+)
% para un metabolito dado. Idempotente: si la reaccion o metabolito _p ya
% existe, no se duplica. Si M_<metStem>_p no existe, se crea copiando
% formula/carga/nombre desde M_<metStem>_c.
%
% Inputs:
%   m         - modelo COBRA
%   metStem   - raiz del metabolito sin prefijo M_ ni sufijo (ej. 'fru', 'lac__L')
%   texId     - ID BiGG de la reaccion tex (ej. 'FRUtex')
%   ppId      - ID BiGG de la reaccion pp (ej. 'FRUt3', 'BUTt2rpp')
%   texName   - nombre legible para tex
%   ppName    - nombre legible para pp
% Output:
%   m         - modelo modificado
function m = addPeripPair(m, metStem, texId, ppId, texName, ppName)
    met_e = ['M_' metStem '_e'];
    met_p = ['M_' metStem '_p'];
    met_c = ['M_' metStem '_c'];

    % Pre-check: necesitamos met_c (si no existe, no podemos avanzar)
    idx_c = find(strcmp(m.mets, met_c), 1);
    if isempty(idx_c)
        fprintf('  [SKIP %s] %s no existe en citosol\n', metStem, met_c);
        return;
    end

    % Crear M_<met>_e si no existe (heredando propiedades de _c)
    if ~any(strcmp(m.mets, met_e))
        m = addMetabolite(m, met_e, ...
            'metName',    [m.metNames{idx_c} ' (extracellular)'], ...
            'metFormula', m.metFormulas{idx_c}, ...
            'charge',     m.metCharges(idx_c), ...
            'csense',     'E');
        fprintf('  Agregado metabolito: %s\n', met_e);
    end

    % Crear M_<met>_p si no existe (heredando propiedades de _c)
    if ~any(strcmp(m.mets, met_p))
        m = addMetabolite(m, met_p, ...
            'metName',    [m.metNames{idx_c} ' (periplasm)'], ...
            'metFormula', m.metFormulas{idx_c}, ...
            'charge',     m.metCharges(idx_c), ...
            'csense',     'E');
        fprintf('  Agregado metabolito: %s\n', met_p);
    end

    % Agregar tex (e<->p, difusion)
    if ~any(strcmp(m.rxns, ['R_' texId])) && ~any(strcmp(m.rxns, texId))
        m = addReaction(m, ['R_' texId], ...
            'reactionName',  texName, ...
            'metaboliteList', {met_e, met_p}, ...
            'stoichCoeffList', [-1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Transport, outer membrane');
        fprintf('  Agregada: R_%s\n', texId);
    else
        fprintf('  [OK] %s ya existe\n', texId);
    end

    % Agregar pp (p<->c, symport 1H+ reversible)
    if ~any(strcmp(m.rxns, ['R_' ppId])) && ~any(strcmp(m.rxns, ppId))
        m = addReaction(m, ['R_' ppId], ...
            'reactionName',  ppName, ...
            'metaboliteList', {met_p, 'M_h_p', met_c, 'M_h_c'}, ...
            'stoichCoeffList', [-1, -1, 1, 1], ...
            'reversible', true, 'lowerBound', -1000, 'upperBound', 1000, ...
            'subSystem', 'Transport, inner membrane');
        fprintf('  Agregada: R_%s\n', ppId);
    else
        fprintf('  [OK] %s ya existe\n', ppId);
    end
end

% /media/alexis/hdd2/postdoc/0_Simulaciones_MATLAB_paper_HGF2_PT33_correcciones/00_nuevos_modelos_metabolicos_carveme_gapseq/1_output_modelos_carveme/
% 
% nombres de gems exportados
% Bacteroides_thetaiotaomicron_VPI_5482.mat
% Bacteroides_thetaiotaomicron_VPI_5482_curado.xml
% 
% Bifidobacterium_animalis_lactis_PT33.mat
% Bifidobacterium_animalis_lactis_PT33_curado.xml
% 
% Clostridium_sp_HGF2.mat
% Clostridium_sp_HGF2_curado.xml
% 
% Clostridium_sp_M62_1.mat
% Clostridium_sp_M62_1_curado.xml
% 
% Clostridium_symbiosum_WAL_14673.mat
% Clostridium_symbiosum_WAL_14673_curado.xml
% 
% Lacticaseibacillus_paracasei_M38.mat
% Lacticaseibacillus_paracasei_M38_curado.xml

