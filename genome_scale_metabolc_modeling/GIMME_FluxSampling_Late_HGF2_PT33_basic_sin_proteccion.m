%% =====================================================================
%  GIMME_FluxSampling_Late_HGF2_PT33_basic_sin_proteccion.m
%  Version simplificada SIN proteccion de reacciones GIMME
%  =====================================================================
%
%  Misma logica que GIMME_FluxSampling_Late_HGF2_PT33_basic.m pero:
%   - SIN lista de reacciones protegidas
%   - SIN recovery reverse-greedy
%   - Solo chequeo de growth con error claro si GIMME deja biomasa = 0
%
%  Pasos:
%   1.  GIMME de HGF2 con counts de Consortium_Late
%   2.  GIMME de PT33 con counts de Consortium_Late
%   3.  Polish HGF2 para SteadyCom (formato BiGG)
%   4.  Polish PT33 para SteadyCom
%   5.  Crear modelo comunitario con createMultipleSpeciesModel
%   6.  Configurar medio comunitario (ZMB + inulina + anaerobiosis)
%   7.  SteadyCom (verificacion libre)
%   8.  Fijar biomasas individuales (mu_exp + fracciones RNA-seq Late)
%   9.  Restriccion HGF2 acetato/lactato UB=0
%  10.  Flux sampling CHRR
%  11.  Guardar resultados
%
%  Requiere: COBRA Toolbox inicializado, Gurobi como solver
% =====================================================================

clear; clc;

% initCobraToolbox(false)
changeCobraSolver('gurobi');

cd /media/alexis/hdd2/postdoc/0_Simulaciones_MATLAB_paper_HGF2_PT33_correcciones/00_nuevos_modelos_metabolicos_carveme_gapseq/1_output_modelos_carveme/1_contextualizacion_trasncriptomica

baseDir = fileparts(mfilename('fullpath'));
if isempty(baseDir)
    baseDir = pwd;
end


%% =====================================================================
%  CONFIGURACION GENERAL
%  =====================================================================

% --- Archivos de entrada ---
file_modelo_HGF2  = fullfile(baseDir, 'Clostridium_sp_HGF2.mat');
file_modelo_PT33  = fullfile(baseDir, 'Bifidobacterium_animalis_lactis_PT33.mat');
file_counts_HGF2  = fullfile(baseDir, 'Raw_Counts_HGF2_only_kallisto.csv');
file_counts_PT33  = fullfile(baseDir, 'Raw_Counts_PT33_only_kallisto.csv');
file_annot_HGF2   = fullfile(baseDir, 'RNA_seq_HGF2', 'Genome_annotation', 'Clostridium_HGF2_annotations.csv');
file_annot_PT33   = fullfile(baseDir, 'RNA_seq_PT33', 'Genome_annotation', 'Bi_animalis_lactis_PT33_annotations.csv');

% --- Columnas de Consortium_Late en los CSVs de counts ---
cols_late_HGF2 = {'Consortium_Late_R1', 'Consortium_Late_R2', 'Consortium_Late_R3'};
cols_late_PT33 = {'Consortium_Late_R1', 'Consortium_Late_R2', 'Consortium_Late_R3'};

% --- Parametros experimentales ---
mu_exp_HGF2 = 0.073;   % h^-1 (HGF2 Late)
mu_exp_PT33 = 0.357;   % h^-1 (PT33 Late)

% Fracciones de RNA-seq Late: estimadas como la proporcion de CDS-mapped reads
% atribuidas a cada organismo en el cocultivo Late (promedio por replica:
% HGF2 = 24114 reads, PT33 = 4490600 reads -> 0.53% / 99.47%).
% Estos valores son los reportados en la seccion de resultados del paper.
frac_HGF2_Late = 0.0053;
frac_PT33_Late = 0.9947;

% --- Parametros GIMME ---
percentil_TPM = 40;    % percentil del TPM positivo

% --- Parametros flux sampling ---
nStepsPerPoint  = 500;       % thinning alto para reducir autocorrelacion
nPointsReturned = 200000;
nombre_escenario = 'Consortium_Late_Inulina';


%% =====================================================================
%  MEDIO ZMB (sin fuente de carbono — la inulina se agrega aparte)
%  =====================================================================
mZMB = {
    'R_EX_ala__L_e';     % L-alanina
    'R_EX_arg__L_e';     % L-arginina
    'R_EX_asn__L_e';     % L-asparagina
    'R_EX_asp__L_e';     % L-aspartato
    'R_EX_cys__L_e';     % L-cisteina
    'R_EX_gln__L_e';     % L-glutamina
    'R_EX_glu__L_e';     % L-glutamato
    'R_EX_gly_e';        % Glicina
    'R_EX_his__L_e';     % L-histidina
    'R_EX_ile__L_e';     % L-isoleucina
    'R_EX_leu__L_e';     % L-leucina
    'R_EX_lys__L_e';     % L-lisina
    'R_EX_met__L_e';     % L-metionina
    'R_EX_phe__L_e';     % L-fenilalanina
    'R_EX_pro__L_e';     % L-prolina
    'R_EX_ser__L_e';     % L-serina
    'R_EX_thr__L_e';     % L-treonina
    'R_EX_trp__L_e';     % L-triptofano
    'R_EX_tyr__L_e';     % L-tirosina
    'R_EX_val__L_e';     % L-valina
    'R_EX_pi_e';         % Fosfato
    'R_EX_k_e';          % Potasio
    'R_EX_mg2_e';        % Magnesio
    'R_EX_fe3_e';        % Hierro (III)
    'R_EX_fe2_e';        % Hierro (II)
    'R_EX_na1_e';        % Sodio
    'R_EX_zn2_e';        % Zinc
    'R_EX_mn2_e';        % Manganeso
    'R_EX_mobd_e';       % Molibdato
    'R_EX_nh3_e';        % Amonio
    'R_EX_ca2_e';        % Calcio
    'R_EX_so4_e';        % Sulfato
    'R_EX_cobalt2_e';    % Cobalto
    'R_EX_cu2_e';        % Cobre
    'R_EX_cl_e';         % Cloruro
    'R_EX_4ppan_e';      % Pantotenato
    'R_EX_fol_e';        % Folato
    'R_EX_nac_e';        % Niacina
    'R_EX_pydx_e';       % Piridoxina
    'R_EX_4abz_e';       % 4-aminobenzoato
    'R_EX_xan_e';        % Xantina
    'R_EX_inost_e';      % Inositol
    'R_EX_gthox_e';      % Glutation oxidado
    'R_EX_gthrd_e';      % Glutation reducido
    'R_EX_btn_e';        % Biotina
    'R_EX_ribflv_e';     % Riboflavina
    'R_EX_thm_e';        % Tiamina
    'R_EX_lipoate_e';    % Lipoato
    'R_EX_pheme_e';      % Heme
    'R_EX_sheme_e';      % Siroheme
    'R_EX_ade_e';        % Adenina
    'R_EX_gua_e';        % Guanina
    'R_EX_ura_e';        % Uracilo
    'R_EX_ddca_e';       % Laurato
    'R_EX_h2s_e';        % H2S
    'R_EX_hxan_e';       % Hipoxantina
    'R_EX_mqn7_e';       % Menaquinona-7
    'R_EX_mqn8_e';       % Menaquinona-8
    'R_EX_ni2_e';        % Niquel
    'R_EX_ocdca_e';      % Octadecanoato
    'R_EX_q8_e';         % Ubiquinona-8
    'R_EX_thymd_e';      % Timidina
    'R_EX_spmd_e';       % Espermidina
    'R_EX_pnto__R_e';    % Pantotenato R
    'R_EX_h2o__R_e';
    'R_EX_nh4_e';        % Amonio (nh4)
};

% Rates: AAs=1, resto=1000. Inulina (fuente de C) se fija a -10 mas abajo.
n_aa = 20;
rates_mZMB = repmat(1000, length(mZMB), 1);
rates_mZMB(1:n_aa) = 1;

%% #####################################################################
%  PARTE 1: GIMME PARA HGF2 (counts de Consortium_Late)
%  #####################################################################
fprintf('\n#######################################################\n');
fprintf('  GIMME para HGF2 — Consortium_Late + inulina\n');
fprintf('#######################################################\n');


% --- 1.1) Cargar modelo HGF2 ---
fprintf('\n--- Cargando modelo HGF2 ---\n');
tmp_HGF2 = load(file_modelo_HGF2);
fnames_HGF2 = fieldnames(tmp_HGF2);
model_HGF2 = tmp_HGF2.(fnames_HGF2{1});
fprintf('  Rxns: %d, Mets: %d, Genes: %d\n', ...
    length(model_HGF2.rxns), length(model_HGF2.mets), length(model_HGF2.genes));

if ~isfield(model_HGF2, 'rules') || isempty(model_HGF2.rules)
    model_HGF2 = generateRules(model_HGF2);
    fprintf('  Campo "rules" generado desde grRules\n');
end


% --- 1.2) Configurar medio HGF2 (ZMB + inulina + anaerobiosis) ---
fprintf('\n--- Configurando medio HGF2 (ZMB + inulina) ---\n');

% (a) Cerrar TODOS los exchanges (LB=0, UB=1000)
ex_by_name_HGF2   = model_HGF2.rxns(contains(model_HGF2.rxns, 'EX_'));
ex_by_struct_HGF2 = model_HGF2.rxns(findExcRxns(model_HGF2));
ex_rxns_HGF2 = unique([ex_by_name_HGF2; ex_by_struct_HGF2]);
for j = 1:length(ex_rxns_HGF2)
    model_HGF2 = changeRxnBounds(model_HGF2, ex_rxns_HGF2{j}, 0, 'l');
    model_HGF2 = changeRxnBounds(model_HGF2, ex_rxns_HGF2{j}, 1000, 'u');
end
fprintf('  Exchanges cerrados: %d\n', length(ex_rxns_HGF2));

% (b) Abrir medio ZMB
nAbiertas_HGF2 = 0;
for j = 1:length(mZMB)
    rxn_id     = mZMB{j};
    rxn_id_alt = regexprep(rxn_id, '^R_', '');
    idx = find(strcmp(model_HGF2.rxns, rxn_id) | strcmp(model_HGF2.rxns, rxn_id_alt));
    if ~isempty(idx)
        model_HGF2 = changeRxnBounds(model_HGF2, model_HGF2.rxns{idx(1)}, -rates_mZMB(j), 'l');
        nAbiertas_HGF2 = nAbiertas_HGF2 + 1;
    end
end
fprintf('  Nutrientes ZMB abiertos: %d / %d\n', nAbiertas_HGF2, length(mZMB));

% (c) Cerrar glucosa y fructosa libres (forzar uso de inulina)
idx_glc_HGF2 = find(strcmp(model_HGF2.rxns, 'R_EX_glc__D_e') | strcmp(model_HGF2.rxns, 'EX_glc__D_e'));
idx_fru_HGF2 = find(strcmp(model_HGF2.rxns, 'R_EX_fru_e')    | strcmp(model_HGF2.rxns, 'EX_fru_e'));
if ~isempty(idx_glc_HGF2)
    model_HGF2 = changeRxnBounds(model_HGF2, model_HGF2.rxns{idx_glc_HGF2(1)}, 0, 'l');
end
if ~isempty(idx_fru_HGF2)
    model_HGF2 = changeRxnBounds(model_HGF2, model_HGF2.rxns{idx_fru_HGF2(1)}, 0, 'l');
end

% (d) Abrir inulina (LB=UB=-10)
idx_inu_HGF2 = find(strcmp(model_HGF2.rxns, 'EX_inulin(e)') | ...
                    strcmp(model_HGF2.rxns, 'R_EX_inulin_e') | ...
                    strcmp(model_HGF2.rxns, 'EX_inulin_e'));
if ~isempty(idx_inu_HGF2)
    model_HGF2 = changeRxnBounds(model_HGF2, model_HGF2.rxns{idx_inu_HGF2(1)}, -10, 'b');
    fprintf('  Inulina abierta: %s LB=UB=-10\n', model_HGF2.rxns{idx_inu_HGF2(1)});
else
    fprintf('  ADVERTENCIA: reaccion de inulina NO encontrada en HGF2!\n');
end

% (e) Anaerobiosis: cerrar O2
idx_o2_HGF2 = find(strcmp(model_HGF2.rxns, 'R_EX_o2_e') | strcmp(model_HGF2.rxns, 'EX_o2_e'));
if ~isempty(idx_o2_HGF2)
    model_HGF2 = changeRxnBounds(model_HGF2, model_HGF2.rxns{idx_o2_HGF2(1)}, 0, 'l');
end

% (f) H2O, CO2, H+ libres (LB=-1000)
rxns_libres_HGF2 = {'R_EX_h2o_e', 'R_EX_co2_e', 'R_EX_h_e'};
for j = 1:length(rxns_libres_HGF2)
    rxn_id     = rxns_libres_HGF2{j};
    rxn_id_alt = regexprep(rxn_id, '^R_', '');
    idx = find(strcmp(model_HGF2.rxns, rxn_id) | strcmp(model_HGF2.rxns, rxn_id_alt));
    if ~isempty(idx)
        model_HGF2 = changeRxnBounds(model_HGF2, model_HGF2.rxns{idx(1)}, -1000, 'l');
    end
end


% --- 1.3) Leer raw counts HGF2 ---
fprintf('\n--- Leyendo raw counts HGF2 ---\n');
T_counts_HGF2 = readtable(file_counts_HGF2, 'VariableNamingRule', 'preserve');
fprintf('  Genes en counts: %d\n', height(T_counts_HGF2));

geneIDs_counts_HGF2 = T_counts_HGF2{:, 1};
if isstring(geneIDs_counts_HGF2)
    geneIDs_counts_HGF2 = cellstr(geneIDs_counts_HGF2);
end

% Extraer columnas Consortium_Late
countMatrix_HGF2 = zeros(height(T_counts_HGF2), length(cols_late_HGF2));
for c = 1:length(cols_late_HGF2)
    colIdx = find(strcmp(T_counts_HGF2.Properties.VariableNames, cols_late_HGF2{c}));
    if isempty(colIdx)
        error('Columna %s no encontrada en %s', cols_late_HGF2{c}, file_counts_HGF2);
    end
    countMatrix_HGF2(:, c) = T_counts_HGF2{:, colIdx};
end
fprintf('  Columnas: %s\n', strjoin(cols_late_HGF2, ', '));


% --- 1.4) Leer gene lengths HGF2 ---
fprintf('\n--- Leyendo gene lengths HGF2 ---\n');
T_annot_HGF2 = readtable(file_annot_HGF2, 'VariableNamingRule', 'modify');
fprintf('  Genes en anotaciones: %d\n', height(T_annot_HGF2));

geneIDs_annot_HGF2 = T_annot_HGF2{:, 1};
if isstring(geneIDs_annot_HGF2)
    geneIDs_annot_HGF2 = cellstr(geneIDs_annot_HGF2);
end

% Buscar columna de length (puede llamarse length_nt_, length_nt, etc.)
lenColIdx_HGF2 = find(strcmp(T_annot_HGF2.Properties.VariableNames, 'length_nt_'));
if isempty(lenColIdx_HGF2)
    allCols = T_annot_HGF2.Properties.VariableNames;
    lenColIdx_HGF2 = find(contains(allCols, 'length') & contains(allCols, 'nt'));
    lenColIdx_HGF2 = lenColIdx_HGF2(1);
end

rawLen_HGF2 = T_annot_HGF2{:, lenColIdx_HGF2};
if iscell(rawLen_HGF2) || isstring(rawLen_HGF2)
    geneLengths_annot_HGF2 = str2double(rawLen_HGF2);
else
    geneLengths_annot_HGF2 = double(rawLen_HGF2);
end

% Mapear lengths a la tabla de counts (mismo orden que counts)
[~, idxAnnot_HGF2, idxCounts_HGF2] = intersect(geneIDs_annot_HGF2, geneIDs_counts_HGF2, 'stable');
geneLengths_HGF2 = nan(height(T_counts_HGF2), 1);
geneLengths_HGF2(idxCounts_HGF2) = geneLengths_annot_HGF2(idxAnnot_HGF2);
fprintf('  Genes con length mapeado: %d / %d\n', sum(~isnan(geneLengths_HGF2)), height(T_counts_HGF2));


% --- 1.5) Calcular TPM HGF2 (vectorizado, sin loop interno por gen) ---
fprintf('\n--- Calculando TPM HGF2 ---\n');
length_kb_HGF2 = geneLengths_HGF2 / 1000;
length_kb_HGF2(length_kb_HGF2 <= 0 | isnan(length_kb_HGF2)) = NaN;

% RPK = counts / (length en kb)  — broadcast por columna
rpk_HGF2 = countMatrix_HGF2 ./ length_kb_HGF2;
rpk_HGF2(isnan(rpk_HGF2)) = 0;

% TPM = (RPK / sumRPK) * 1e6 por columna
sum_rpk_HGF2 = sum(rpk_HGF2, 1);
tpmMatrix_HGF2 = (rpk_HGF2 ./ sum_rpk_HGF2) * 1e6;
tpm_avg_HGF2 = mean(tpmMatrix_HGF2, 2);

fprintf('  TPM — min: %.2f, median: %.2f, max: %.2f\n', ...
    min(tpm_avg_HGF2), median(tpm_avg_HGF2), max(tpm_avg_HGF2));
fprintf('  Genes con TPM > 0: %d / %d\n', sum(tpm_avg_HGF2 > 0), length(tpm_avg_HGF2));


% --- 1.6) Mapear gene IDs al modelo HGF2 (CHGFRNA -> RCHGF) ---
fprintf('\n--- Mapeando gene IDs al modelo HGF2 ---\n');
geneIDs_mapped_HGF2 = regexprep(geneIDs_counts_HGF2, '^CHGFRNA', 'RCHGF');
fprintf('  Mapeo HGF2: CHGFRNA -> RCHGF\n');

[genesEnModelo_HGF2, idxMapped_HGF2, idxModel_HGF2] = ...
    intersect(geneIDs_mapped_HGF2, model_HGF2.genes, 'stable');
fprintf('  Genes mapeados al modelo: %d / %d\n', ...
    length(genesEnModelo_HGF2), length(model_HGF2.genes));


% --- 1.7) Construir expressionData HGF2 ---
expressionData_HGF2.gene  = model_HGF2.genes;
expressionData_HGF2.value = nan(length(model_HGF2.genes), 1);
for g = 1:length(genesEnModelo_HGF2)
    expressionData_HGF2.value(idxModel_HGF2(g)) = tpm_avg_HGF2(idxMapped_HGF2(g));
end


% --- 1.8) Mapear expresion a reacciones HGF2 ---
fprintf('\n--- Mapeando expresion a reacciones HGF2 ---\n');
[expressionRxns_HGF2, ~] = mapExpressionToReactions(model_HGF2, expressionData_HGF2);
fprintf('  Rxns con expresion: %d, sin expresion: %d\n', ...
    sum(expressionRxns_HGF2 >= 0), sum(expressionRxns_HGF2 < 0));


% --- 1.9) Threshold y GIMME HGF2 ---
tpm_positivos_HGF2 = tpm_avg_HGF2(tpm_avg_HGF2 > 0);
threshold_HGF2 = prctile(tpm_positivos_HGF2, percentil_TPM);

fprintf('\n--- Corriendo GIMME HGF2 (threshold=%.2f, percentil %d) ---\n', ...
    threshold_HGF2, percentil_TPM);

options_gimme_HGF2 = struct();
options_gimme_HGF2.solver         = 'GIMME';
options_gimme_HGF2.expressionRxns = expressionRxns_HGF2;
options_gimme_HGF2.threshold      = threshold_HGF2;
options_gimme_HGF2.obj_frac       = 0.9;

tissueModel_HGF2 = createTissueSpecificModel(model_HGF2, options_gimme_HGF2);
fprintf('  Rxns: %d -> %d (eliminadas: %d)\n', ...
    length(model_HGF2.rxns), length(tissueModel_HGF2.rxns), ...
    length(model_HGF2.rxns) - length(tissueModel_HGF2.rxns));


% --- 1.10) Verificar growth HGF2 (sin recovery, solo error si falla) ---
solGIMME_HGF2 = optimizeCbModel(tissueModel_HGF2, 'max');
solOrig_HGF2  = optimizeCbModel(model_HGF2, 'max');
fprintf('  Growth: GIMME=%.4f, original=%.4f (%.1f%%)\n', ...
    solGIMME_HGF2.f, solOrig_HGF2.f, (solGIMME_HGF2.f / solOrig_HGF2.f) * 100);

if solGIMME_HGF2.f < 1e-6
    error(['HGF2 GIMME dejo growth = 0. ' ...
           'Considera bajar percentil_TPM o revisar reacciones eliminadas. ' ...
           'Si necesitas el recovery automatico, usa el script con proteccion.']);
end


% --- 1.11) Guardar modelo GIMME HGF2 Late ---
HGF2_GIMME_Late = tissueModel_HGF2;
save(fullfile(baseDir, 'HGF2_GIMME_Late_sin_proteccion.mat'), 'HGF2_GIMME_Late');
fprintf('\n  Guardado: HGF2_GIMME_Late_sin_proteccion.mat\n');




%% #####################################################################
%  PARTE 2: GIMME PARA PT33 (counts de Consortium_Late)
%  #####################################################################
fprintf('\n#######################################################\n');
fprintf('  GIMME para PT33 — Consortium_Late + inulina\n');
fprintf('#######################################################\n');


% --- 2.1) Cargar modelo PT33 ---
fprintf('\n--- Cargando modelo PT33 ---\n');
tmp_PT33 = load(file_modelo_PT33);
fnames_PT33 = fieldnames(tmp_PT33);
model_PT33 = tmp_PT33.(fnames_PT33{1});
fprintf('  Rxns: %d, Mets: %d, Genes: %d\n', ...
    length(model_PT33.rxns), length(model_PT33.mets), length(model_PT33.genes));

if ~isfield(model_PT33, 'rules') || isempty(model_PT33.rules)
    model_PT33 = generateRules(model_PT33);
    fprintf('  Campo "rules" generado desde grRules\n');
end


% --- 2.2) Configurar medio PT33 (ZMB + inulina + cross-feeding fru/glc/kestopt) ---
fprintf('\n--- Configurando medio PT33 (ZMB + inulina + cross-feeding) ---\n');

% (a) Cerrar TODOS los exchanges
ex_by_name_PT33   = model_PT33.rxns(contains(model_PT33.rxns, 'EX_'));
ex_by_struct_PT33 = model_PT33.rxns(findExcRxns(model_PT33));
ex_rxns_PT33 = unique([ex_by_name_PT33; ex_by_struct_PT33]);
for j = 1:length(ex_rxns_PT33)
    model_PT33 = changeRxnBounds(model_PT33, ex_rxns_PT33{j}, 0, 'l');
    model_PT33 = changeRxnBounds(model_PT33, ex_rxns_PT33{j}, 1000, 'u');
end
fprintf('  Exchanges cerrados: %d\n', length(ex_rxns_PT33));

% (b) Abrir medio ZMB
nAbiertas_PT33 = 0;
for j = 1:length(mZMB)
    rxn_id     = mZMB{j};
    rxn_id_alt = regexprep(rxn_id, '^R_', '');
    idx = find(strcmp(model_PT33.rxns, rxn_id) | strcmp(model_PT33.rxns, rxn_id_alt));
    if ~isempty(idx)
        model_PT33 = changeRxnBounds(model_PT33, model_PT33.rxns{idx(1)}, -rates_mZMB(j), 'l');
        nAbiertas_PT33 = nAbiertas_PT33 + 1;
    end
end
fprintf('  Nutrientes ZMB abiertos: %d / %d\n', nAbiertas_PT33, length(mZMB));

% (c) Cerrar glucosa y fructosa libres por defecto (las abrimos como cross-feeding abajo)
idx_glc_PT33 = find(strcmp(model_PT33.rxns, 'R_EX_glc__D_e') | strcmp(model_PT33.rxns, 'EX_glc__D_e'));
idx_fru_PT33 = find(strcmp(model_PT33.rxns, 'R_EX_fru_e')    | strcmp(model_PT33.rxns, 'EX_fru_e'));
if ~isempty(idx_glc_PT33)
    model_PT33 = changeRxnBounds(model_PT33, model_PT33.rxns{idx_glc_PT33(1)}, 0, 'l');
end
if ~isempty(idx_fru_PT33)
    model_PT33 = changeRxnBounds(model_PT33, model_PT33.rxns{idx_fru_PT33(1)}, 0, 'l');
end

% (d) Abrir inulina (PT33 no tiene inulinasa, pero la dejamos por consistencia)
idx_inu_PT33 = find(strcmp(model_PT33.rxns, 'EX_inulin(e)') | ...
                    strcmp(model_PT33.rxns, 'R_EX_inulin_e') | ...
                    strcmp(model_PT33.rxns, 'EX_inulin_e'));
if ~isempty(idx_inu_PT33)
    model_PT33 = changeRxnBounds(model_PT33, model_PT33.rxns{idx_inu_PT33(1)}, -10, 'b');
    fprintf('  Inulina abierta: %s LB=UB=-10\n', model_PT33.rxns{idx_inu_PT33(1)});
end

% (e) CROSS-FEEDING: PT33 recibe fru, glc y kestopt de HGF2 (que degrada inulina)
%     Sin esto, GIMME elimina los transportadores de azucar simples de PT33.
idx_fru_cf = find(strcmp(model_PT33.rxns, 'R_EX_fru_e') | strcmp(model_PT33.rxns, 'EX_fru_e') | strcmp(model_PT33.rxns, 'EX_fru(e)'));
idx_glc_cf = find(strcmp(model_PT33.rxns, 'R_EX_glc__D_e') | strcmp(model_PT33.rxns, 'EX_glc__D_e') | strcmp(model_PT33.rxns, 'EX_glc__D(e)'));
idx_kpt_cf = find(strcmp(model_PT33.rxns, 'EX_kestopt(e)') | strcmp(model_PT33.rxns, 'R_EX_kestopt_e') | strcmp(model_PT33.rxns, 'EX_kestopt_e'));

if ~isempty(idx_fru_cf)
    model_PT33 = changeRxnBounds(model_PT33, model_PT33.rxns{idx_fru_cf(1)}, -10, 'l');
    fprintf('  Cross-feeding: %s LB=-10 (fructosa de HGF2)\n', model_PT33.rxns{idx_fru_cf(1)});
end
if ~isempty(idx_glc_cf)
    model_PT33 = changeRxnBounds(model_PT33, model_PT33.rxns{idx_glc_cf(1)}, -10, 'l');
    fprintf('  Cross-feeding: %s LB=-10 (glucosa de HGF2)\n', model_PT33.rxns{idx_glc_cf(1)});
end
if ~isempty(idx_kpt_cf)
    model_PT33 = changeRxnBounds(model_PT33, model_PT33.rxns{idx_kpt_cf(1)}, -10, 'l');
    fprintf('  Cross-feeding: %s LB=-10 (kestopentaosa de HGF2)\n', model_PT33.rxns{idx_kpt_cf(1)});
end

% (f) Anaerobiosis
idx_o2_PT33 = find(strcmp(model_PT33.rxns, 'R_EX_o2_e') | strcmp(model_PT33.rxns, 'EX_o2_e'));
if ~isempty(idx_o2_PT33)
    model_PT33 = changeRxnBounds(model_PT33, model_PT33.rxns{idx_o2_PT33(1)}, 0, 'l');
end

% (g) H2O, CO2, H+ libres
rxns_libres_PT33 = {'R_EX_h2o_e', 'R_EX_co2_e', 'R_EX_h_e'};
for j = 1:length(rxns_libres_PT33)
    rxn_id     = rxns_libres_PT33{j};
    rxn_id_alt = regexprep(rxn_id, '^R_', '');
    idx = find(strcmp(model_PT33.rxns, rxn_id) | strcmp(model_PT33.rxns, rxn_id_alt));
    if ~isempty(idx)
        model_PT33 = changeRxnBounds(model_PT33, model_PT33.rxns{idx(1)}, -1000, 'l');
    end
end


% --- 2.3) Leer raw counts PT33 ---
fprintf('\n--- Leyendo raw counts PT33 ---\n');
T_counts_PT33 = readtable(file_counts_PT33, 'VariableNamingRule', 'preserve');
fprintf('  Genes en counts: %d\n', height(T_counts_PT33));

geneIDs_counts_PT33 = T_counts_PT33{:, 1};
if isstring(geneIDs_counts_PT33)
    geneIDs_counts_PT33 = cellstr(geneIDs_counts_PT33);
end

countMatrix_PT33 = zeros(height(T_counts_PT33), length(cols_late_PT33));
for c = 1:length(cols_late_PT33)
    colIdx = find(strcmp(T_counts_PT33.Properties.VariableNames, cols_late_PT33{c}));
    if isempty(colIdx)
        error('Columna %s no encontrada en %s', cols_late_PT33{c}, file_counts_PT33);
    end
    countMatrix_PT33(:, c) = T_counts_PT33{:, colIdx};
end
fprintf('  Columnas: %s\n', strjoin(cols_late_PT33, ', '));


% --- 2.4) Leer gene lengths PT33 ---
fprintf('\n--- Leyendo gene lengths PT33 ---\n');
T_annot_PT33 = readtable(file_annot_PT33, 'VariableNamingRule', 'modify');
fprintf('  Genes en anotaciones: %d\n', height(T_annot_PT33));

geneIDs_annot_PT33 = T_annot_PT33{:, 1};
if isstring(geneIDs_annot_PT33)
    geneIDs_annot_PT33 = cellstr(geneIDs_annot_PT33);
end

lenColIdx_PT33 = find(strcmp(T_annot_PT33.Properties.VariableNames, 'length_nt'));
if isempty(lenColIdx_PT33)
    allCols = T_annot_PT33.Properties.VariableNames;
    lenColIdx_PT33 = find(contains(allCols, 'length') & contains(allCols, 'nt'));
    lenColIdx_PT33 = lenColIdx_PT33(1);
end

rawLen_PT33 = T_annot_PT33{:, lenColIdx_PT33};
if iscell(rawLen_PT33) || isstring(rawLen_PT33)
    geneLengths_annot_PT33 = str2double(rawLen_PT33);
else
    geneLengths_annot_PT33 = double(rawLen_PT33);
end

[~, idxAnnot_PT33, idxCounts_PT33] = intersect(geneIDs_annot_PT33, geneIDs_counts_PT33, 'stable');
geneLengths_PT33 = nan(height(T_counts_PT33), 1);
geneLengths_PT33(idxCounts_PT33) = geneLengths_annot_PT33(idxAnnot_PT33);
fprintf('  Genes con length mapeado: %d / %d\n', sum(~isnan(geneLengths_PT33)), height(T_counts_PT33));


% --- 2.5) Calcular TPM PT33 (vectorizado) ---
fprintf('\n--- Calculando TPM PT33 ---\n');
length_kb_PT33 = geneLengths_PT33 / 1000;
length_kb_PT33(length_kb_PT33 <= 0 | isnan(length_kb_PT33)) = NaN;

rpk_PT33 = countMatrix_PT33 ./ length_kb_PT33;
rpk_PT33(isnan(rpk_PT33)) = 0;

sum_rpk_PT33 = sum(rpk_PT33, 1);
tpmMatrix_PT33 = (rpk_PT33 ./ sum_rpk_PT33) * 1e6;
tpm_avg_PT33 = mean(tpmMatrix_PT33, 2);

fprintf('  TPM — min: %.2f, median: %.2f, max: %.2f\n', ...
    min(tpm_avg_PT33), median(tpm_avg_PT33), max(tpm_avg_PT33));
fprintf('  Genes con TPM > 0: %d / %d\n', sum(tpm_avg_PT33 > 0), length(tpm_avg_PT33));


% --- 2.6) Mapear gene IDs PT33 (directos, sin renombrar) ---
fprintf('\n--- Mapeando gene IDs al modelo PT33 ---\n');
geneIDs_mapped_PT33 = geneIDs_counts_PT33;

[genesEnModelo_PT33, idxMapped_PT33, idxModel_PT33] = ...
    intersect(geneIDs_mapped_PT33, model_PT33.genes, 'stable');
fprintf('  Genes mapeados al modelo: %d / %d\n', ...
    length(genesEnModelo_PT33), length(model_PT33.genes));


% --- 2.7) Construir expressionData PT33 ---
expressionData_PT33.gene  = model_PT33.genes;
expressionData_PT33.value = nan(length(model_PT33.genes), 1);
for g = 1:length(genesEnModelo_PT33)
    expressionData_PT33.value(idxModel_PT33(g)) = tpm_avg_PT33(idxMapped_PT33(g));
end


% --- 2.8) Mapear expresion a reacciones PT33 ---
fprintf('\n--- Mapeando expresion a reacciones PT33 ---\n');
[expressionRxns_PT33, ~] = mapExpressionToReactions(model_PT33, expressionData_PT33);
fprintf('  Rxns con expresion: %d, sin expresion: %d\n', ...
    sum(expressionRxns_PT33 >= 0), sum(expressionRxns_PT33 < 0));


% --- 2.9) Threshold y GIMME PT33 ---
tpm_positivos_PT33 = tpm_avg_PT33(tpm_avg_PT33 > 0);
threshold_PT33 = prctile(tpm_positivos_PT33, percentil_TPM);

fprintf('\n--- Corriendo GIMME PT33 (threshold=%.2f, percentil %d) ---\n', ...
    threshold_PT33, percentil_TPM);

options_gimme_PT33 = struct();
options_gimme_PT33.solver         = 'GIMME';
options_gimme_PT33.expressionRxns = expressionRxns_PT33;
options_gimme_PT33.threshold      = threshold_PT33;
options_gimme_PT33.obj_frac       = 0.9;

tissueModel_PT33 = createTissueSpecificModel(model_PT33, options_gimme_PT33);
fprintf('  Rxns: %d -> %d (eliminadas: %d)\n', ...
    length(model_PT33.rxns), length(tissueModel_PT33.rxns), ...
    length(model_PT33.rxns) - length(tissueModel_PT33.rxns));


% --- 2.10) Verificar growth PT33 (sin recovery, solo error si falla) ---
solGIMME_PT33 = optimizeCbModel(tissueModel_PT33, 'max');
solOrig_PT33  = optimizeCbModel(model_PT33, 'max');
fprintf('  Growth: GIMME=%.4f, original=%.4f (%.1f%%)\n', ...
    solGIMME_PT33.f, solOrig_PT33.f, (solGIMME_PT33.f / solOrig_PT33.f) * 100);

if solGIMME_PT33.f < 1e-6
    error(['PT33 GIMME dejo growth = 0. ' ...
           'Considera bajar percentil_TPM o revisar reacciones eliminadas. ' ...
           'Si necesitas el recovery automatico, usa el script con proteccion.']);
end


% --- 2.11) Guardar modelo GIMME PT33 Late ---
PT33_GIMME_Late = tissueModel_PT33;
save(fullfile(baseDir, 'PT33_GIMME_Late_sin_proteccion.mat'), 'PT33_GIMME_Late');
fprintf('\n  Guardado: PT33_GIMME_Late_sin_proteccion.mat\n');




%% #####################################################################
%  PARTE 3: POLISH HGF2 PARA STEADYCOM
%  Convierte IDs al formato BiGG (xxx[c], xxx(e), sin prefijo R_)
%  #####################################################################
fprintf('\n#######################################################\n');
fprintf('  Polish HGF2 para SteadyCom\n');
fprintf('#######################################################\n');

model_HGF2_polish = tissueModel_HGF2;

% --- 3.1) Metabolitos: M_xxx_c -> xxx[c] ---
model_HGF2_polish.mets = regexprep(model_HGF2_polish.mets, '^M_(.+)_([cep])$', '$1[$2]');

% Mets sin bracket: convertir formato xxx_c -> xxx[c]
sin_bracket_HGF2 = cellfun(@(x) ~contains(x, '['), model_HGF2_polish.mets);
if any(sin_bracket_HGF2)
    model_HGF2_polish.mets(sin_bracket_HGF2) = ...
        regexprep(model_HGF2_polish.mets(sin_bracket_HGF2), '^(.+)_([cep])$', '$1[$2]');
end

% --- 3.2) Reacciones: quitar prefijo R_ ---
model_HGF2_polish.rxns = regexprep(model_HGF2_polish.rxns, '^R_', '');

% --- 3.3) Eliminar reacciones duplicadas ---
[~, idx_uniq_HGF2] = unique(model_HGF2_polish.rxns, 'last');
if length(idx_uniq_HGF2) < length(model_HGF2_polish.rxns)
    n_dups_HGF2 = length(model_HGF2_polish.rxns) - length(idx_uniq_HGF2);
    fprintf('  Eliminando %d reacciones duplicadas en HGF2\n', n_dups_HGF2);
    idx_uniq_HGF2 = sort(idx_uniq_HGF2);
    idx_dups_HGF2 = setdiff(1:length(model_HGF2_polish.rxns), idx_uniq_HGF2);
    for dd = 1:length(idx_dups_HGF2)
        model_HGF2_polish.rxns{idx_dups_HGF2(dd)} = ['DUP_' num2str(dd) '_' model_HGF2_polish.rxns{idx_dups_HGF2(dd)}];
    end
    model_HGF2_polish = removeRxns(model_HGF2_polish, model_HGF2_polish.rxns(idx_dups_HGF2), 'metFlag', false);
end

% --- 3.4) Sufijos de compartimento en reacciones (xxx_e -> xxx(e)) ---
model_HGF2_polish.rxns = regexprep(model_HGF2_polish.rxns, '_e$', '(e)');
model_HGF2_polish.rxns = regexprep(model_HGF2_polish.rxns, '_c$', '(c)');
model_HGF2_polish.rxns = regexprep(model_HGF2_polish.rxns, '_p$', '(p)');

% --- 3.5) Renombrar biomasa: Growth_biomass -> bio ---
model_HGF2_polish.rxns = strrep(model_HGF2_polish.rxns, 'Growth_biomass', 'bio');

% --- 3.6) Llenar celdas vacias en campos opcionales ---
campos_opt = {'metFormulas'; 'genes'; 'metNames'; 'rxnNames'; 'subSystems'};
for ff = 1:length(campos_opt)
    if isfield(model_HGF2_polish, campos_opt{ff})
        vacias = cellfun(@isempty, model_HGF2_polish.(campos_opt{ff}));
        model_HGF2_polish.(campos_opt{ff})(vacias) = {''};
    end
end

% --- 3.7) Objetivo + anaerobiosis ---
model_HGF2_polish = changeObjective(model_HGF2_polish, 'bio', 1);
if ismember('EX_o2(e)', model_HGF2_polish.rxns)
    model_HGF2_polish = changeRxnBounds(model_HGF2_polish, 'EX_o2(e)', 0, 'l');
end

fprintf('  HGF2 polished: %d rxns, %d mets\n', ...
    length(model_HGF2_polish.rxns), length(model_HGF2_polish.mets));




%% #####################################################################
%  PARTE 4: POLISH PT33 PARA STEADYCOM
%  #####################################################################
fprintf('\n#######################################################\n');
fprintf('  Polish PT33 para SteadyCom\n');
fprintf('#######################################################\n');

model_PT33_polish = tissueModel_PT33;

% --- 4.1) Metabolitos ---
model_PT33_polish.mets = regexprep(model_PT33_polish.mets, '^M_(.+)_([cep])$', '$1[$2]');
sin_bracket_PT33 = cellfun(@(x) ~contains(x, '['), model_PT33_polish.mets);
if any(sin_bracket_PT33)
    model_PT33_polish.mets(sin_bracket_PT33) = ...
        regexprep(model_PT33_polish.mets(sin_bracket_PT33), '^(.+)_([cep])$', '$1[$2]');
end

% --- 4.2) Reacciones: quitar R_ ---
model_PT33_polish.rxns = regexprep(model_PT33_polish.rxns, '^R_', '');

% --- 4.3) Duplicadas ---
[~, idx_uniq_PT33] = unique(model_PT33_polish.rxns, 'last');
if length(idx_uniq_PT33) < length(model_PT33_polish.rxns)
    n_dups_PT33 = length(model_PT33_polish.rxns) - length(idx_uniq_PT33);
    fprintf('  Eliminando %d reacciones duplicadas en PT33\n', n_dups_PT33);
    idx_uniq_PT33 = sort(idx_uniq_PT33);
    idx_dups_PT33 = setdiff(1:length(model_PT33_polish.rxns), idx_uniq_PT33);
    for dd = 1:length(idx_dups_PT33)
        model_PT33_polish.rxns{idx_dups_PT33(dd)} = ['DUP_' num2str(dd) '_' model_PT33_polish.rxns{idx_dups_PT33(dd)}];
    end
    model_PT33_polish = removeRxns(model_PT33_polish, model_PT33_polish.rxns(idx_dups_PT33), 'metFlag', false);
end

% --- 4.4) Sufijos de compartimento ---
model_PT33_polish.rxns = regexprep(model_PT33_polish.rxns, '_e$', '(e)');
model_PT33_polish.rxns = regexprep(model_PT33_polish.rxns, '_c$', '(c)');
model_PT33_polish.rxns = regexprep(model_PT33_polish.rxns, '_p$', '(p)');

% --- 4.5) Biomasa ---
model_PT33_polish.rxns = strrep(model_PT33_polish.rxns, 'Growth_biomass', 'bio');

% --- 4.6) Celdas vacias ---
for ff = 1:length(campos_opt)
    if isfield(model_PT33_polish, campos_opt{ff})
        vacias = cellfun(@isempty, model_PT33_polish.(campos_opt{ff}));
        model_PT33_polish.(campos_opt{ff})(vacias) = {''};
    end
end

% --- 4.7) Objetivo + anaerobiosis ---
model_PT33_polish = changeObjective(model_PT33_polish, 'bio', 1);
if ismember('EX_o2(e)', model_PT33_polish.rxns)
    model_PT33_polish = changeRxnBounds(model_PT33_polish, 'EX_o2(e)', 0, 'l');
end

fprintf('  PT33 polished: %d rxns, %d mets\n', ...
    length(model_PT33_polish.rxns), length(model_PT33_polish.mets));




%% #####################################################################
%  PARTE 5: CREAR MODELO COMUNITARIO Y CONFIGURAR MEDIO
%  #####################################################################
fprintf('\n#######################################################\n');
fprintf('  Creando modelo comunitario HGF2 + PT33\n');
fprintf('#######################################################\n');

% Medio ZMB en formato BiGG (compartimento [u])
mZMB_cobra = regexprep(mZMB, '^R_EX_(.+)_e$', 'EX_$1(e)');
iZMB_cobra = [mZMB_cobra; {'EX_h2o(e)'}];
rates_iZMB = [rates_mZMB; 1000];
iZMB_u     = strrep(iZMB_cobra, '(e)', '[u]');


% --- 5.1) createMultipleSpeciesModel ---
Names  = {'HGF2'; 'PT33'};
models = {model_HGF2_polish; model_PT33_polish};

EcCom = createMultipleSpeciesModel(models, Names);
EcCom.csense = repmat('E', numel(EcCom.mets), 1);
[EcCom.infoCom, EcCom.indCom] = getMultiSpeciesModelId(EcCom, Names);

% Asignar reacciones de biomasa
rxnBio   = strcat(Names, 'bio');
rxnBioId = findRxnIDs(EcCom, rxnBio);
EcCom.infoCom.spBm = rxnBio;
EcCom.indCom.spBm  = rxnBioId;

fprintf('  Modelo comunitario: %d rxns, %d mets\n', length(EcCom.rxns), length(EcCom.mets));


% --- 5.2) Configurar medio comunitario ---
fprintf('\n--- Configurando medio comunitario (ZMB + inulina) ---\n');

% (a) Abrir IEX (transferencias entre especies y pool)
exsp1 = EcCom.infoCom.EXsp(:,1);  exsp1 = exsp1(~cellfun(@isempty, exsp1));
exsp2 = EcCom.infoCom.EXsp(:,2);  exsp2 = exsp2(~cellfun(@isempty, exsp2));
EcCom = changeRxnBounds(EcCom, exsp1, -100, 'l');
EcCom = changeRxnBounds(EcCom, exsp2, -100, 'l');
EcCom = changeRxnBounds(EcCom, exsp1,  100, 'u');
EcCom = changeRxnBounds(EcCom, exsp2,  100, 'u');

% (b) Cerrar todos los exchanges comunitarios [u]
rxns_EX_u = EcCom.rxns(startsWith(EcCom.rxns, 'EX_') & contains(EcCom.rxns, '[u]'));
EcCom = changeRxnBounds(EcCom, rxns_EX_u, 0, 'l');

% (c) Abrir medio ZMB en [u]
[rxns_validas, ~, ~] = intersect(iZMB_u, EcCom.rxns);
for jj = 1:length(rxns_validas)
    pos_original = find(strcmp(iZMB_u, rxns_validas{jj}));
    EcCom = changeRxnBounds(EcCom, rxns_validas{jj}, -rates_iZMB(pos_original(1)), 'l');
end
fprintf('  Medio ZMB abierto en [u]: %d nutrientes\n', length(rxns_validas));

% (d) Inulina como fuente de carbono
if ismember('EX_inulin[u]', EcCom.rxns)
    EcCom = changeRxnBounds(EcCom, 'EX_inulin[u]', -10, 'l');
    fprintf('  Inulina: EX_inulin[u] LB=-10\n');
else
    fprintf('  ADVERTENCIA: EX_inulin[u] no encontrada!\n');
end

% (e) Anaerobiosis: cerrar O2 comunitario
if ismember('EX_o2[u]', EcCom.rxns)
    EcCom = changeRxnBounds(EcCom, 'EX_o2[u]', 0, 'b');
end

% (f) H2O, CO2, H+ libres
rxns_free_u = {'EX_h2o[u]'; 'EX_co2[u]'; 'EX_h[u]'};
for jj = 1:length(rxns_free_u)
    if ismember(rxns_free_u{jj}, EcCom.rxns)
        EcCom = changeRxnBounds(EcCom, rxns_free_u{jj}, -1000, 'l');
    end
end




%% #####################################################################
%  PARTE 6: STEADYCOM (verificacion sin restriccion de abundancias)
%  #####################################################################
fprintf('\n#######################################################\n');
fprintf('  SteadyCom (libre, verificacion de comunidad)\n');
fprintf('#######################################################\n');

options_sc = struct();
options_sc.GRguess   = 0.1;
options_sc.GRtol     = 1e-3;
options_sc.algorithm = 1;

[sol_sc, result_sc] = SteadyCom(EcCom, options_sc);
fprintf('  Status: %s, GRmax: %.6f\n', result_sc.stat, result_sc.GRmax);

if isnan(result_sc.GRmax) || result_sc.GRmax < 1e-6
    error('SteadyCom sin crecimiento — revisar configuracion del medio comunitario');
end

totalBM_free = sum(result_sc.BM);
fprintf('  Biomasa total: %.6f\n', totalBM_free);
fprintf('  HGF2: BM=%.6f (%.1f%%)\n', result_sc.BM(1), result_sc.BM(1)/totalBM_free*100);
fprintf('  PT33: BM=%.6f (%.1f%%)\n', result_sc.BM(2), result_sc.BM(2)/totalBM_free*100);




%% #####################################################################
%  PARTE 7: FIJAR BIOMASAS SEGUN RNA-seq + mu_exp
%  #####################################################################
fprintf('\n#######################################################\n');
fprintf('  Fijando biomasas individuales (RNA-seq Late + mu_exp)\n');
fprintf('#######################################################\n');

% mu de la comunidad ponderado por las fracciones RNA-seq
mu_com_Late = mu_exp_HGF2 * frac_HGF2_Late + mu_exp_PT33 * frac_PT33_Late;
vBM_HGF2    = mu_com_Late * frac_HGF2_Late;
vBM_PT33    = mu_com_Late * frac_PT33_Late;

fprintf('  mu_exp HGF2 = %.4f, PT33 = %.4f\n', mu_exp_HGF2, mu_exp_PT33);
fprintf('  Fracciones RNA-seq Late: HGF2=%.3f%%, PT33=%.3f%%\n', ...
    frac_HGF2_Late*100, frac_PT33_Late*100);
fprintf('  mu comunidad ponderado = %.6f h^-1\n', mu_com_Late);
fprintf('  vBM_HGF2 = %.6f, vBM_PT33 = %.6f\n', vBM_HGF2, vBM_PT33);

% Tolerancia minima para evitar el error "x0 not interior" en CHRR.
% Fijar LB == UB exacto crea un politopo de volumen 0 que rompe el rounding
% del mve_solver_cobra. Usar un rango +/- tol da volumen biologicamente
% indistinguible pero matematicamente bien definido.
tol_bm = 1e-4;   % relajado de 1e-6 a 1e-4 para que ACHR genere warmup points sin LPs infactibles

% Fijar reacciones de biomasa de cada especie (con tolerancia)
idx_bio_HGF2 = findRxnIDs(EcCom, 'HGF2bio');
idx_bio_PT33 = findRxnIDs(EcCom, 'PT33bio');

if idx_bio_HGF2 > 0
    EcCom = changeRxnBounds(EcCom, 'HGF2bio', vBM_HGF2 - tol_bm, 'l');
    EcCom = changeRxnBounds(EcCom, 'HGF2bio', vBM_HGF2 + tol_bm, 'u');
    fprintf('  Fijado HGF2bio = %.6f (+/- %.0e)\n', vBM_HGF2, tol_bm);
end
if idx_bio_PT33 > 0
    EcCom = changeRxnBounds(EcCom, 'PT33bio', vBM_PT33 - tol_bm, 'l');
    EcCom = changeRxnBounds(EcCom, 'PT33bio', vBM_PT33 + tol_bm, 'u');
    fprintf('  Fijado PT33bio = %.6f (+/- %.0e)\n', vBM_PT33, tol_bm);
end

% Fijar IEX de biomasa con la misma tolerancia
iex_bio_rxns = EcCom.rxns(contains(EcCom.rxns, 'IEX_biomass'));
for jj = 1:length(iex_bio_rxns)
    rxn_iex_b = iex_bio_rxns{jj};
    if startsWith(rxn_iex_b, 'HGF2')
        EcCom = changeRxnBounds(EcCom, rxn_iex_b, vBM_HGF2 - tol_bm, 'l');
        EcCom = changeRxnBounds(EcCom, rxn_iex_b, vBM_HGF2 + tol_bm, 'u');
        fprintf('  Fijado %s = %.6f (+/- %.0e)\n', rxn_iex_b, vBM_HGF2, tol_bm);
    elseif startsWith(rxn_iex_b, 'PT33')
        EcCom = changeRxnBounds(EcCom, rxn_iex_b, vBM_PT33 - tol_bm, 'l');
        EcCom = changeRxnBounds(EcCom, rxn_iex_b, vBM_PT33 + tol_bm, 'u');
        fprintf('  Fijado %s = %.6f (+/- %.0e)\n', rxn_iex_b, vBM_PT33, tol_bm);
    end
end




%% #####################################################################
%  PARTE 8: RESTRICCION HGF2 ACETATO + LACTATO (UB=0, solo consumo)
%  COMENTADA — corrida sin restriccion para comparar con el caso AcConstrained
%  #####################################################################
fprintf('\n#######################################################\n');
fprintf('  PARTE 8 desactivada: SIN restriccion HGF2 ac/lac\n');
fprintf('#######################################################\n');

% % Acetato HGF2: solo consumo
% idx_ac_HGF2 = find(strcmp(EcCom.rxns, 'HGF2IEX_ac[u]tr'));
% if ~isempty(idx_ac_HGF2)
%     old_ub = EcCom.ub(idx_ac_HGF2);
%     EcCom.ub(idx_ac_HGF2) = 0;
%     fprintf('  HGF2IEX_ac[u]tr: UB cambiado de %.1f a 0\n', old_ub);
% else
%     fprintf('  WARNING: HGF2IEX_ac[u]tr no encontrado!\n');
% end
%
% % L-lactato HGF2: solo consumo
% idx_lacL_HGF2 = find(strcmp(EcCom.rxns, 'HGF2IEX_lac__L[u]tr'));
% if ~isempty(idx_lacL_HGF2)
%     old_ub = EcCom.ub(idx_lacL_HGF2);
%     EcCom.ub(idx_lacL_HGF2) = 0;
%     fprintf('  HGF2IEX_lac__L[u]tr: UB cambiado de %.1f a 0\n', old_ub);
% else
%     fprintf('  WARNING: HGF2IEX_lac__L[u]tr no encontrado!\n');
% end
%
% % D-lactato HGF2: solo consumo
% idx_lacD_HGF2 = find(strcmp(EcCom.rxns, 'HGF2IEX_lac__D[u]tr'));
% if ~isempty(idx_lacD_HGF2)
%     old_ub = EcCom.ub(idx_lacD_HGF2);
%     EcCom.ub(idx_lacD_HGF2) = 0;
%     fprintf('  HGF2IEX_lac__D[u]tr: UB cambiado de %.1f a 0\n', old_ub);
% else
%     fprintf('  WARNING: HGF2IEX_lac__D[u]tr no encontrado!\n');
% end


% --- Verificar factibilidad antes del sampling ---
sol_check = optimizeCbModel(EcCom, 'max');
if sol_check.stat ~= 1
    fprintf('\n  Modelo infactible. Relajando feasTol a 1e-4...\n');
    changeCobraSolverParams('LP', 'feasTol', 1e-4);
    sol_check = optimizeCbModel(EcCom, 'max');
    if sol_check.stat ~= 1
        error('Modelo sigue infactible despues de relajar feasTol');
    end
end
fprintf('\n  Modelo factible (obj=%.6f)\n', sol_check.f);




%% #####################################################################
%  PARTE 9: FLUX SAMPLING CHRR
%  #####################################################################
fprintf('\n#######################################################\n');
fprintf('  Flux Sampling CHRR (nPoints=%d, nSteps=%d)\n', nPointsReturned, nStepsPerPoint);
fprintf('#######################################################\n');

changeCobraSolverParams('LP', 'feasTol', 1e-4);
fprintf('  feasTol relajado a 1e-4 para CHRR\n');

opts = struct();
opts.nStepsPerPoint  = nStepsPerPoint;
opts.nPointsReturned = nPointsReturned;

% tic;
% [modelSampling, samples] = sampleCbModel(EcCom, ...
%     sprintf('Sampling_%s_noAcConstraint_sin_proteccion', nombre_escenario), ...
%     'CHRR', opts);
% t_sampling = toc;

tic;
[modelSampling, samples] = sampleCbModel(EcCom, ...
    sprintf('Sampling_%s_noAcConstraint_sin_proteccion', nombre_escenario), ...
    'CHRR', opts);
t_sampling = toc;
                    

fprintf('  Sampling completado en %.1f seg\n', t_sampling);
fprintf('  Dimensiones: %d rxns x %d puntos\n', size(samples, 1), size(samples, 2));




%% #####################################################################
%  PARTE 10: GUARDAR RESULTADOS
%  #####################################################################
fprintf('\n#######################################################\n');
fprintf('  Guardando resultados\n');
fprintf('#######################################################\n');

outDir_sampling = fullfile(baseDir, 'flux_sampling');
if ~exist(outDir_sampling, 'dir'), mkdir(outDir_sampling); end

nombre_out = [nombre_escenario '_noAcConstraint_sin_proteccion'];

% .mat completo
sampling_result = struct();
sampling_result.samples         = samples;
sampling_result.rxns            = EcCom.rxns;
sampling_result.model           = EcCom;
sampling_result.result_sc       = result_sc;
sampling_result.nPointsReturned = nPointsReturned;
sampling_result.nStepsPerPoint  = nStepsPerPoint;
sampling_result.mu_com          = mu_com_Late;
sampling_result.vBM_HGF2        = vBM_HGF2;
sampling_result.vBM_PT33        = vBM_PT33;
sampling_result.frac_HGF2       = frac_HGF2_Late;
sampling_result.frac_PT33       = frac_PT33_Late;

save(fullfile(outDir_sampling, sprintf('sampling_%s.mat', nombre_out)), ...
    'sampling_result', '-v7.3');
fprintf('  Guardado: sampling_%s.mat\n', nombre_out);

% CSVs
writematrix(samples, fullfile(outDir_sampling, ...
    sprintf('samples_%s_%dk.csv', nombre_out, round(nPointsReturned/1000))));
writecell(EcCom.rxns, fullfile(outDir_sampling, ...
    sprintf('rxn_names_%s.csv', nombre_out)));
fprintf('  CSVs guardados\n');




%% #####################################################################
%  ESTADISTICAS RAPIDAS DE FLUJOS IEX
%  #####################################################################
fprintf('\n--- Estadisticas de flujos IEX (sampling) ---\n');
mets_stats     = {'ac', 'but', 'lac__L', 'lac__D', 'for', 'etoh', ...
                  'succ', 'co2', 'fru', 'glc__D', 'inulin', 'kestopt', 'glu__L'};
nombres_stats  = {'Acetato', 'Butirato', 'L-Lactato', 'D-Lactato', ...
                  'Formato', 'Etanol', 'Succinato', 'CO2', ...
                  'Fructosa', 'Glucosa', 'Inulina', 'Kestopentaosa', 'Glutamato'};

fprintf('  %-15s  %-8s  %10s  %10s  %10s  %s\n', ...
    'Metabolito', 'Bacteria', 'Media', 'Mediana', 'SD', '[Min, Max]');
fprintf('  %s\n', repmat('-', 1, 80));

for mm = 1:length(mets_stats)
    % HGF2
    rxn_name_H = ['HGF2IEX_' mets_stats{mm} '[u]tr'];
    idx_H = find(strcmp(EcCom.rxns, rxn_name_H));
    if ~isempty(idx_H)
        vals = samples(idx_H, :);
        if any(abs(vals) > 1e-9)
            fprintf('  %-15s  %-8s  %10.4f  %10.4f  %10.4f  [%.4f, %.4f]\n', ...
                nombres_stats{mm}, 'HGF2', mean(vals), median(vals), std(vals), min(vals), max(vals));
        end
    end
    % PT33
    rxn_name_P = ['PT33IEX_' mets_stats{mm} '[u]tr'];
    idx_P = find(strcmp(EcCom.rxns, rxn_name_P));
    if ~isempty(idx_P)
        vals = samples(idx_P, :);
        if any(abs(vals) > 1e-9)
            fprintf('  %-15s  %-8s  %10.4f  %10.4f  %10.4f  [%.4f, %.4f]\n', ...
                nombres_stats{mm}, 'PT33', mean(vals), median(vals), std(vals), min(vals), max(vals));
        end
    end
end

fprintf('\n#######################################################\n');
fprintf('  GIMME_FluxSampling_Late_HGF2_PT33_basic_sin_proteccion.m FINALIZADO\n');
fprintf('#######################################################\n');


%% Outputs quedan en la siguiente carpeta:
% /media/alexis/hdd2/postdoc/0_Simulaciones_MATLAB_paper_HGF2_PT33_correcciones/00_nuevos_modelos_metabolicos_carveme_gapseq/1_output_modelos_carveme/1_contextualizacion_trasncriptomica/flux_sampling

% archcivos generados:

% sampling_Consortium_Late_Inulina_noAcConstraint_sin_proteccion.mat
% samples_Consortium_Late_Inulina_noAcConstraint_sin_proteccion_200k.csv
% rxn_names_Consortium_Late_Inulina_noAcConstraint_sin_proteccion.csv

% HGF2_GIMME_Late_sin_proteccion.mat (HGF2_GIMME_Late); 1786 → 1432 (354 eliminadas, 19.8%) 
% PT33_GIMME_Late_sin_proteccion.mat (PT33_GIMME_Late); 1638 → 1406 (232 eliminadas, 14.2%)

load('HGF2_GIMME_Late_sin_proteccion.mat')   
load('PT33_GIMME_Late_sin_proteccion.mat')   

