% =========================================================================
% D_curar_GPR.m
% =========================================================================
% Curacion manual de GPR (Gene-Protein-Reaction rules) para HGF2 y PT33.
%
% Sin GPR, GIMME no puede usar datos RNA-seq para contextualizar el modelo.
% Este script agrega GPR a reacciones clave del metabolismo de butirato (HGF2),
% complejo Rnf, degradacion de FOS (PT33) y otras reacciones por cruce de
% numeros EC con eggnog-mapper.
%
% Pipeline:
%   1. B_change_bounds_reactions_actual_27_24_26.m  -> .mat finales
%   2. >>> D_curar_GPR.m (este script) <<<           -> .mat con GPR
%   3. C_monocultivos_cocultivos_inulina.m           -> simulaciones
%
% Input:  Clostridium_sp_HGF2.mat, Bifidobacterium_animalis_lactis_PT33.mat
% Output: los mismos archivos, sobrescritos con GPR agregadas
%
% Fuentes de anotacion GPR:
%   1. eggnog-mapper: cruce de EC numbers entre modelo y anotacion del genoma
%   2. Revision manual de literatura: rutas de butirato (Rnf, BTCOADH, BUTCT)
%      y degradacion FOS (KESTOASE, beta-fructofuranosidasa)
%   3. BLASTp contra genomas de referencia para EC sin cruce directo
% =========================================================================
clear all; clc;
changeCobraSolver('gurobi', 'LP');

cd /media/alexis/hdd2/postdoc/0_Simulaciones_MATLAB_paper_HGF2_PT33_correcciones/00_nuevos_modelos_metabolicos_carveme_gapseq/1_output_modelos_carveme

%% =========================================================================
%% HGF2: Clostridium sp. HGF2
%% =========================================================================
% GPR a agregar:
%   - Ruta del butirato (4): HACD1, ECOAH1, BTCOADH (complejo 3 sub.), BUTCT
%   - Inulinasa extracelular: INULINASEe
%   - Piruvato:ferredoxina oxidorreductasa: POR4
%   - Complejo Rnf (6 subunidades): FDNADOX_H
%   - Cruces por EC (14): reacciones donde el numero EC del modelo coincide
%     con la anotacion de eggnog-mapper
fprintf('\n========================================\n');
fprintf('  HGF2: Clostridium sp. HGF2\n');
fprintf('========================================\n');

s = load('Clostridium_sp_HGF2.mat');
hgf2_var = fieldnames(s); hgf2_var = hgf2_var{1};
m_hgf2 = s.(hgf2_var);
fprintf('Cargado: %s (%d rxns, %d mets, %d genes)\n', ...
        hgf2_var, length(m_hgf2.rxns), length(m_hgf2.mets), length(m_hgf2.genes));
fprintf('GPR pre-existentes: %d/%d (%.1f%%)\n', ...
        sum(~cellfun(@isempty, m_hgf2.grRules)), length(m_hgf2.rxns), ...
        100*sum(~cellfun(@isempty, m_hgf2.grRules))/length(m_hgf2.rxns));

% Cell array {rxn_id, gpr_string} - rxn_id sin prefijo R_, el helper lo busca con/sin
hgf2_gpr = {
    % --- Ruta del butirato ---
    'HACD1',        'RCHGF003502';
    'ECOAH1',       'RCHGF003501';
    'BTCOADH',      'RCHGF003503 and RCHGF003504 and RCHGF003505';
    'BUTCT',        'RCHGF003500';
    % --- Inulina ---
    'INULINASEe',   'RCHGF003366';
    % --- Piruvato y ferredoxina ---
    'POR4',         'RCHGF003035';
    % --- Complejo Rnf (6 subunidades) ---
    'FDNADOX_H',    'RCHGF003725 and RCHGF003727 and RCHGF003728 and RCHGF003729 and RCHGF003730 and RCHGF003731';
    % --- Cruces por EC (eggnog-mapper) ---
    '4ABUTD',       'RCHGF002900';
    'ACACT6r',      'RCHGF003506 or RCHGF003507';
    'ACOAD1f',      'RCHGF002689 or RCHGF003505';
    'G3PAT120',     'RCHGF003970';
    'G3PAT161',     'RCHGF003970';
    'G3PAT180',     'RCHGF003970';
    'G3PAT181',     'RCHGF003970';
    'GLDBRAN2',     'RCHGF000050';
    'GTHRDHpp',     'RCHGF001629';
    'HMBS',         'RCHGF000195';
    'LGTHL',        'RCHGF000811 or RCHGF001113';
    'NP1',          'RCHGF002144 or RCHGF002686 or RCHGF003340';
    'FEENTERabcpp', 'RCHGF000332';
    'UAG4Ei',       'RCHGF000111';
};

m_hgf2 = aplicarGPR(m_hgf2, hgf2_gpr, 'HGF2');
fprintf('GPR despues: %d/%d (%.1f%%)\n', ...
        sum(~cellfun(@isempty, m_hgf2.grRules)), length(m_hgf2.rxns), ...
        100*sum(~cellfun(@isempty, m_hgf2.grRules))/length(m_hgf2.rxns));

% Verificar FBA
m_hgf2 = changeObjective(m_hgf2, 'Growth_biomass', 1);
sol = optimizeCbModel(m_hgf2);
fprintf('FBA mu = %.4f (status=%d)\n', sol.f, sol.stat);

% Guardar conservando el nombre de variable original
eval([hgf2_var ' = m_hgf2;']);
save('Clostridium_sp_HGF2.mat', hgf2_var);
fprintf('[OK] Guardado Clostridium_sp_HGF2.mat (variable: %s)\n', hgf2_var);

%% =========================================================================
%% PT33: Bifidobacterium animalis subsp. lactis PT33
%% =========================================================================
% GPR a agregar:
%   - Degradacion de FOS (3): KESTOASEe, KESTOTTRASEe, KESTOPTASEe
%     (todas asignadas al mismo gen de beta-fructofuranosidasa)
%   - Cruces por EC (32): reacciones donde el numero EC del modelo coincide
%     con la anotacion de eggnog-mapper
%
% Ademas se ELIMINA PFK (fosfofructoquinasa). Bifidobacterium NO tiene este
% gen en su genoma; usa el bifid shunt (via fosfocetoalasa, F6PPK), que es
% la ruta caracteristica del genero.
fprintf('\n========================================\n');
fprintf('  PT33: Bifidobacterium animalis subsp. lactis PT33\n');
fprintf('========================================\n');

s = load('Bifidobacterium_animalis_lactis_PT33.mat');
pt33_var = fieldnames(s); pt33_var = pt33_var{1};
m_pt33 = s.(pt33_var);
fprintf('Cargado: %s (%d rxns, %d mets, %d genes)\n', ...
        pt33_var, length(m_pt33.rxns), length(m_pt33.mets), length(m_pt33.genes));
fprintf('GPR pre-existentes: %d/%d (%.1f%%)\n', ...
        sum(~cellfun(@isempty, m_pt33.grRules)), length(m_pt33.rxns), ...
        100*sum(~cellfun(@isempty, m_pt33.grRules))/length(m_pt33.rxns));

% Eliminar PFK (Bifidobacterium no tiene gen, usa bifid shunt)
idx_pfk = find(strcmp(m_pt33.rxns, 'PFK') | strcmp(m_pt33.rxns, 'R_PFK'), 1);
if ~isempty(idx_pfk)
    m_pt33 = removeRxns(m_pt33, m_pt33.rxns(idx_pfk), 'metFlag', false);
    fprintf('PFK eliminada (Bifidobacterium no tiene gen, usa bifid shunt)\n');
else
    fprintf('PFK no encontrada (ya fue eliminada previamente)\n');
end

pt33_gpr = {
    % --- Degradacion de FOS (beta-fructofuranosidasa) ---
    'KESTOASEe',    'RVCYZ001167';
    'KESTOTTRASEe', 'RVCYZ001167';
    'KESTOPTASEe',  'RVCYZ001167';
    % --- Cruces por EC (eggnog-mapper) ---
    '3SALATAi',     'RVCYZ000360 or RVCYZ000409 or RVCYZ000751 or RVCYZ000753 or RVCYZ000963 or RVCYZ001195';
    'ACGS',         'RVCYZ000109';
    'AGPAT120',     'RVCYZ000351 or RVCYZ001226';
    'AGPAT141',     'RVCYZ000351 or RVCYZ001226';
    'AGPAT161',     'RVCYZ000351 or RVCYZ001226';
    'AGPAT180',     'RVCYZ000351 or RVCYZ001226';
    'AGPAT181',     'RVCYZ000351 or RVCYZ001226';
    'AIRC2',        'RVCYZ000238';
    'DHFR2i',       'RVCYZ000554';
    'THFOR2',       'RVCYZ000554';
    'DRBK',         'RVCYZ001386';
    'FACOAL80',     'RVCYZ000063 or RVCYZ000154 or RVCYZ000695 or RVCYZ001047';
    'FACOAL160',    'RVCYZ000063 or RVCYZ000154 or RVCYZ000695 or RVCYZ001047';
    'FACOAL180_2',  'RVCYZ000063 or RVCYZ000154 or RVCYZ000695 or RVCYZ001047';
    'FEENTERabcpp', 'RVCYZ000711';
    'G1PACT',       'RVCYZ000086';
    'GDPTPDP',      'RVCYZ001243';
    'PPGPPDP',      'RVCYZ001243';
    'GLBRAN2',      'RVCYZ000608';
    'GPDDA3pp',     'RVCYZ000797';
    'GPDDA4',       'RVCYZ000797';
    'GRTT',         'RVCYZ001228';
    'MDDCP1pp',     'RVCYZ000330 or RVCYZ001129';
    'MDDCP4pp',     'RVCYZ000330 or RVCYZ001129';
    'MI3PP',        'RVCYZ000193';
    'PGSA120',      'RVCYZ000009 or RVCYZ000040 or RVCYZ001215';
    'PGSA161',      'RVCYZ000009 or RVCYZ000040 or RVCYZ001215';
    'PGSA180',      'RVCYZ000009 or RVCYZ000040 or RVCYZ001215';
    'PGSA181',      'RVCYZ000009 or RVCYZ000040 or RVCYZ001215';
    'PYK2',         'RVCYZ000619';
    'PYK4',         'RVCYZ000619';
    'UGLT',         'RVCYZ000993';
};

m_pt33 = aplicarGPR(m_pt33, pt33_gpr, 'PT33');
fprintf('GPR despues: %d/%d (%.1f%%)\n', ...
        sum(~cellfun(@isempty, m_pt33.grRules)), length(m_pt33.rxns), ...
        100*sum(~cellfun(@isempty, m_pt33.grRules))/length(m_pt33.rxns));

% Verificar FBA
m_pt33 = changeObjective(m_pt33, 'Growth_biomass', 1);
sol = optimizeCbModel(m_pt33);
fprintf('FBA mu = %.4f (status=%d)\n', sol.f, sol.stat);

% Guardar conservando el nombre de variable original
eval([pt33_var ' = m_pt33;']);
save('Bifidobacterium_animalis_lactis_PT33.mat', pt33_var);
fprintf('[OK] Guardado Bifidobacterium_animalis_lactis_PT33.mat (variable: %s)\n', pt33_var);

fprintf('\n========================================\n');
fprintf('  CURACION GPR COMPLETA\n');
fprintf('========================================\n');
fprintf('Listo para correr C_monocultivos_cocultivos_inulina.m\n');

%% =========================================================================
%% FUNCION LOCAL: aplicarGPR
%% =========================================================================
% Aplica una lista de reglas GPR a un modelo. Para cada par {rxn_id, gpr},
% busca la rxn equivalente en m.rxns probando con/sin prefijo R_ y con/sin
% sufijo 'i' (ej: 'HACD1', 'R_HACD1', 'HACD1i', 'R_HACD1i'). Si encuentra
% la rxn y NO tiene GPR pre-existente, asigna la GPR. Si ya tiene GPR la
% deja como esta (idempotente). Reporta no encontradas al final.
function m = aplicarGPR(m, gpr_pairs, etiqueta)
    n_added = 0; n_skipped = 0; not_found = {};
    for i = 1:size(gpr_pairs, 1)
        rid = gpr_pairs{i, 1};
        gpr = gpr_pairs{i, 2};
        % Buscar con/sin R_ prefix y con/sin sufijo 'i'
        candidatos = {rid, ['R_' rid], [rid 'i'], ['R_' rid 'i']};
        idx = [];
        for c = 1:length(candidatos)
            idx = find(strcmp(m.rxns, candidatos{c}), 1);
            if ~isempty(idx), break; end
        end
        if isempty(idx)
            not_found{end+1} = rid; %#ok<AGROW>
            continue;
        end
        if isempty(m.grRules{idx})
            m.grRules{idx} = gpr;
            n_added = n_added + 1;
            gpr_short = gpr;
            if length(gpr_short) > 50
                gpr_short = [gpr_short(1:50) '...'];
            end
            fprintf('  + %-20s -> %s\n', m.rxns{idx}, gpr_short);
        else
            n_skipped = n_skipped + 1;
        end
    end
    fprintf('\n[%s] agregadas=%d, ya tenian GPR=%d, no_encontradas=%d\n', ...
            etiqueta, n_added, n_skipped, length(not_found));
    if ~isempty(not_found)
        fprintf('[!!] No encontradas: %s\n', strjoin(not_found, ', '));
    end
end
