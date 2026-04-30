%% =====================================================================
%  G_fix_tex_t2r_bounds.m
%
%  Aplica fix de bounds a transportadores reversibles + fixes puntuales
%  metabolicos justificados, en los .mat de los 6 modelos del paper.
%
%  Justificacion (BiGG: bigg_models_reactions.tsv):
%    - sufijo 'r' (t2r, t4rpp): reversible explicito en BiGG
%    - tpp, t3pp: difusion / antiport (BiGG <-> en todos los casos)
%    - tex: outer-membrane / porinas (ya suelen estar OK con LB=-1000)
%    - abc: NO se tocan (consumen ATP, son irreversibles biologicamente)
%
%  CarveMe a veces deja estos transportadores con LB=0 cuando deberian
%  ser LB=-1000 segun BiGG. Este fix abre la reversibilidad solo para
%  rxns con LB=0 cuya entrada en BiGG diga "<->".
%
%  IMPORTANTE: NO se hace fix general de TODAS las rxns metabolicas. Eso
%  causa loops termodinamicos disipativos en pFBA L2 (aunque BiGG diga
%  <->, muchas rxns son fisiologicamente irreversibles).
%
%  Fixes metabolicos puntuales (lista blanca, justificados caso por caso):
%    - Csym R_BUTCT: escrito al reves (accoa+but -> ac+btcoa, LB=0).
%      Solo CONSUME butirato. Fix: LB=-1000 (BiGG dice <->) para que
%      pueda producir butirato corriendo en reverso.
%
%  IMPORTANTE: este script SOLO escribe .mat (writeCbModel rompe SBML).
%  Para los .xml usar el notebook G_fix_tex_t2r_bounds.ipynb.
%
%  Output: 2_modelos_t2r_tpp_fixed/<nombre>_fixed.mat + fix_log_mat.txt
% =====================================================================

clear; clc;

% --- Inicializar COBRA si hace falta ---
% initCobraToolbox(false)
try
    changeCobraSolver('gurobi');
catch
    fprintf('  (gurobi no disponible — se intenta solver default)\n');
end

baseDir = '/media/alexis/hdd2/postdoc/0_Simulaciones_MATLAB_paper_HGF2_PT33_correcciones/00_nuevos_modelos_metabolicos_carveme_gapseq/1_output_modelos_carveme';
outDir  = fullfile(baseDir, '2_modelos_t2r_tpp_fixed');
if ~exist(outDir, 'dir'), mkdir(outDir); end

biggTSV = '/media/alexis/hdd2/BIGG_database/bigg_models_reactions.tsv';

modelFiles = {
    'Bacteroides_thetaiotaomicron_VPI_5482.mat';
    'Bifidobacterium_animalis_lactis_PT33.mat';
    'Clostridium_sp_HGF2.mat';
    'Clostridium_sp_M62_1.mat';
    'Clostridium_symbiosum_WAL_14673.mat';
    'Lacticaseibacillus_paracasei_M38.mat';
};
shortNames = {'Bt'; 'PT33'; 'HGF2'; 'M62_1'; 'Csym'; 'M38'};

% Sufijos a parchar (longest first para que endsWith no confunda 'tpp' con 't4rpp')
suffixes = {'t4rpp', 't3pp', 't2r', 'tpp'};

% Lista blanca de fixes metabolicos puntuales por modelo
% Formato: container con clave = shortName, valor = cell array de rxn IDs (sin R_)
metabolic_fixes = containers.Map();
metabolic_fixes('Csym') = {'BUTCT'};   % escrito al reves -> abrir reversibilidad


%% --- Cargar tabla BiGG (id -> reaction_string) ---
fprintf('Cargando tabla BiGG (%s)...\n', biggTSV);
fid = fopen(biggTSV, 'r');
if fid == -1
    error('No se pudo abrir BiGG TSV en %s', biggTSV);
end
header = fgetl(fid); %#ok<NASGU>  % skip header
bigg_ids = {};
bigg_rstr = {};
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    parts = strsplit(line, sprintf('\t'), 'CollapseDelimiters', false);
    if length(parts) >= 3
        bigg_ids{end+1, 1} = parts{1};
        bigg_rstr{end+1, 1} = parts{3};
    end
end
fclose(fid);
fprintf('  %d reacciones BiGG cargadas\n', length(bigg_ids));
biggMap = containers.Map(bigg_ids, bigg_rstr);


%% --- Inicializar log ---
logPath = fullfile(outDir, 'fix_log_mat.txt');
flog = fopen(logPath, 'w');
fprintf(flog, 'G_fix_tex_t2r_bounds.m — fix de bounds .mat\n');
fprintf(flog, 'Fecha: %s\n', datestr(now));
fprintf(flog, 'Sufijos transporte: %s\n', strjoin(suffixes, ', '));
fprintf(flog, 'Regla: LB=0 + UB>0 + BiGG dice "<->"  ->  LB=-1000\n');
fprintf(flog, '+ fixes metabolicos puntuales por modelo (lista blanca)\n\n');


%% --- Procesar cada modelo ---
total_changes = 0;
for k = 1:length(modelFiles)
    sn = shortNames{k};
    matIn = fullfile(baseDir, modelFiles{k});
    fprintf('\n=========================================\n');
    fprintf('  %s  (%s)\n', sn, modelFiles{k});
    fprintf('=========================================\n');
    fprintf(flog, '\n=========================================\n  %s\n=========================================\n', sn);

    if ~isfile(matIn)
        fprintf('  ! No existe: %s — saltando\n', matIn);
        fprintf(flog, '  ! No existe — saltado\n');
        continue;
    end

    tmp = load(matIn);
    fnames = fieldnames(tmp);
    model = tmp.(fnames{1});
    fprintf('  Variable: %s, rxns=%d, mets=%d\n', fnames{1}, length(model.rxns), length(model.mets));

    % FBA antes
    fba_before = NaN;
    try
        sol = optimizeCbModel(model, 'max');
        fba_before = sol.f;
        fprintf('  FBA antes: %.6f\n', fba_before);
    catch ME
        fprintf('  (no se pudo correr FBA antes: %s)\n', ME.message);
    end

    n_fix = struct('t2r',0,'t4rpp',0,'t3pp',0,'tpp',0,'metab_puntual',0, ...
                   'total',0,'skipped_no_bigg',0,'skipped_no_rev',0, ...
                   'already_open',0);

    % --- Fase 1: fix de transportadores por sufijo ---
    for j = 1:length(model.rxns)
        rid = model.rxns{j};
        if startsWith(rid, 'R_')
            base = rid(3:end);
        else
            base = rid;
        end
        % Buscar sufijo
        suf = '';
        for s = 1:length(suffixes)
            if endsWith(base, suffixes{s})
                suf = suffixes{s};
                break;
            end
        end
        if isempty(suf), continue; end

        % Verificar contra BiGG
        if ~isKey(biggMap, base)
            n_fix.skipped_no_bigg = n_fix.skipped_no_bigg + 1;
            continue;
        end
        rstr = biggMap(base);
        if ~contains(rstr, '<->') && ~contains(rstr, '<=>')
            n_fix.skipped_no_rev = n_fix.skipped_no_rev + 1;
            continue;
        end

        % Aplicar fix solo si LB=0 y UB>0
        if model.lb(j) == 0 && model.ub(j) > 0
            old_lb = model.lb(j);
            model.lb(j) = -1000;
            n_fix.(suf) = n_fix.(suf) + 1;
            n_fix.total  = n_fix.total + 1;
            fprintf('  FIX %-6s %-30s LB %g -> -1000\n', suf, rid, old_lb);
            fprintf(flog, 'FIX %-6s %-30s LB %g -> -1000\n', suf, rid, old_lb);
        elseif model.lb(j) < 0
            n_fix.already_open = n_fix.already_open + 1;
        end
    end

    % --- Fase 2: fixes metabolicos puntuales (lista blanca por modelo) ---
    if isKey(metabolic_fixes, sn)
        target_rxns = metabolic_fixes(sn);
        for tt = 1:length(target_rxns)
            target = target_rxns{tt};
            % Buscar la rxn (con o sin prefijo R_)
            idx = find(strcmp(model.rxns, target) | strcmp(model.rxns, ['R_' target]));
            if isempty(idx)
                fprintf('  WARN: %s no encontrada en %s\n', target, sn);
                fprintf(flog, 'WARN: %s no encontrada en %s\n', target, sn);
                continue;
            end
            j = idx(1);
            % Verificar BiGG <->
            if ~isKey(biggMap, target)
                fprintf('  WARN: %s no esta en BiGG, no se parcha\n', target);
                continue;
            end
            rstr = biggMap(target);
            if ~contains(rstr, '<->') && ~contains(rstr, '<=>')
                fprintf('  WARN: %s no es <-> en BiGG, no se parcha\n', target);
                continue;
            end
            if model.lb(j) == 0 && model.ub(j) > 0
                old_lb = model.lb(j);
                model.lb(j) = -1000;
                n_fix.metab_puntual = n_fix.metab_puntual + 1;
                n_fix.total = n_fix.total + 1;
                fprintf('  FIX METAB %-30s LB %g -> -1000  [BiGG: %s]\n', ...
                    model.rxns{j}, old_lb, rstr(1:min(60,length(rstr))));
                fprintf(flog, 'FIX METAB %-30s LB %g -> -1000  [BiGG: %s]\n', ...
                    model.rxns{j}, old_lb, rstr);
            else
                fprintf('  %s ya tiene LB=%g (no necesita fix)\n', target, model.lb(j));
            end
        end
    end

    fprintf('  Total fix: %d  (t2r=%d, t4rpp=%d, t3pp=%d, tpp=%d, metab=%d)\n', ...
        n_fix.total, n_fix.t2r, n_fix.t4rpp, n_fix.t3pp, n_fix.tpp, n_fix.metab_puntual);
    fprintf('  Ya con LB<0: %d  |  No en BiGG: %d  |  No <-> en BiGG: %d\n', ...
        n_fix.already_open, n_fix.skipped_no_bigg, n_fix.skipped_no_rev);
    fprintf(flog, 'Total fix: %d  (t2r=%d t4rpp=%d t3pp=%d tpp=%d metab=%d)\n', ...
        n_fix.total, n_fix.t2r, n_fix.t4rpp, n_fix.t3pp, n_fix.tpp, n_fix.metab_puntual);
    fprintf(flog, 'Ya abiertos: %d, no BiGG: %d, no <->: %d\n', ...
        n_fix.already_open, n_fix.skipped_no_bigg, n_fix.skipped_no_rev);
    total_changes = total_changes + n_fix.total;

    % FBA despues
    try
        sol = optimizeCbModel(model, 'max');
        fba_after = sol.f;
        fprintf('  FBA despues: %.6f\n', fba_after);
        fprintf(flog, 'FBA antes=%.6f, despues=%.6f\n', fba_before, fba_after);
        if ~isnan(fba_before) && fba_before > 1e-6
            rel = fba_after / fba_before;
            fprintf('  Razon despues/antes: %.4f\n', rel);
            if abs(rel - 1) > 0.05
                fprintf('  *** ALERTA: cambio de biomasa > 5%% ***\n');
                fprintf(flog, 'ALERTA: cambio biomasa > 5%% (rel=%.4f)\n', rel);
            end
        end
    catch ME
        fprintf('  (no se pudo correr FBA despues: %s)\n', ME.message);
    end

    % Guardar .mat conservando el nombre original de la variable
    [~, basemat, ~] = fileparts(modelFiles{k});
    matOut = fullfile(outDir, [basemat '_fixed.mat']);
    saveStruct = struct();
    saveStruct.(fnames{1}) = model;
    save(matOut, '-struct', 'saveStruct');
    fprintf('  Guardado: %s\n', matOut);
    fprintf(flog, 'Guardado: %s\n', matOut);
end

fprintf('\n=========================================\n');
fprintf('  Listo. Total cambios: %d\n', total_changes);
fprintf('  Output .mat: %s\n', outDir);
fprintf('  Para los .xml usar G_fix_tex_t2r_bounds.ipynb\n');
fprintf('=========================================\n');
fprintf(flog, '\nTOTAL CAMBIOS: %d\n', total_changes);
fclose(flog);
