%% SCRIPT C: MONOCULTIVOS FOS/INULINA (modelos CarveMe/BiGG)
% C_monocultivos_cocultivos_inulina.m
%
% debes correr D2_reexportar_mat_desde_xml.m previamente en los gems que vienen de cobraprython
%
% 6 monocultivos con fuente de carbono especifica:
%   Bt, PT33, M62_1, Csym -> FOS (kesto) LB=UB=-10
%   HGF2, M38             -> inulina LB=UB=-10
%
% Estrategia: max biomasa (UB=mu_exp) -> pFBA L2
% M62_1 no tiene reacciones de FOS ni inulina, pero se le pone FOS igualmente

clearvars
clc

% initCobraToolbox(false)
changeCobraSolver('gurobi');

% debes correr D2_reexportar_mat_desde_xml.m previamente en los gems que vienen de cobraprython

%% ========================================================================
%% CONFIGURACION GENERAL
%% ========================================================================

cd /media/alexis/hdd2/postdoc/0_Simulaciones_MATLAB_paper_HGF2_PT33_correcciones/00_nuevos_modelos_metabolicos_carveme_gapseq/1_output_modelos_carveme/3_modelos_curados_fix_v1

modelFiles = {
    'Bacteroides_thetaiotaomicron_VPI_5482.mat';
    'Bifidobacterium_animalis_lactis_PT33.mat';
    'Clostridium_sp_HGF2.mat';
    'Clostridium_sp_M62_1.mat';
    'Clostridium_symbiosum_WAL_14673.mat';
    'Lacticaseibacillus_paracasei_M38.mat'
};

M38_carveme = readCbModel('Lacticaseibacillus_paracasei_M38.mat')

shortNames = {'Bt'; 'PT33'; 'HGF2'; 'M62_1'; 'Csym'; 'M38'};

% mu experimentales
% HGF2 = 0.073 (fase tardia, consumo inulina)
mu_exp = [0.648; 0.357; 0.073; 0.141; 0.403; 0.21];

% Fuente de carbono por bacteria:
%   'kesto'   = FOS (Bt, PT33, M62_1, Csym)
%   'inulina' = inulina (HGF2, M38)
fuenteC = {'kesto'; 'kesto'; 'inulina'; 'kesto'; 'kesto'; 'inulina'};
% Nombres para graficar (ingles, nombre real del metabolito)
fuenteC_label = {'Kestose'; 'Kestose'; 'Inulin'; 'Kestose'; 'Kestose'; 'Inulin'};

% Quien tiene la capacidad de degradar cada fuente
% (M62_1 no degrada FOS, pero se le pone igualmente)
tiene_inulina = [false; false; true;  false; false; true];
tiene_kesto   = [true;  true;  false; false; true;  true];

%% MEDIO ZMB (sin fuente de carbono — se agrega aparte)
% aa=1 para TODAS las condiciones
mZMB = {
'R_EX_ala__L_e';        % L-alanina
'R_EX_arg__L_e';        % L-arginina
'R_EX_asn__L_e';        % L-asparagina
'R_EX_asp__L_e';        % L-aspartato
'R_EX_cys__L_e';        % L-cisteina
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
'R_EX_trp__L_e';        % L-triptofano
'R_EX_tyr__L_e';        % L-tirosina
'R_EX_val__L_e';        % L-valina
'R_EX_ca2_e';           % Calcio
'R_EX_cl_e';            % Cloruro
'R_EX_cobalt2_e';       % Cobalto
'R_EX_cu2_e';           % Cobre
'R_EX_fe2_e';           % Hierro(II)
'R_EX_fe3_e';           % Hierro(III)
'R_EX_k_e';             % Potasio
'R_EX_mg2_e';           % Magnesio
'R_EX_mn2_e';           % Manganeso
'R_EX_mobd_e';          % Molibdato
'R_EX_na1_e';           % Sodio
'R_EX_zn2_e';           % Zinc
'R_EX_so4_e';           % Sulfato
'R_EX_pi_e';            % Fosfato inorganico
'R_EX_h2o_e';           % Agua
'R_EX_h_e';             % Proton
'R_EX_btn_e';           % Biotina
'R_EX_fol_e';           % Folato
'R_EX_ncam_e';          % Nicotinamida
'R_EX_pydxn_e';         % Piridoxina
'R_EX_ribflv_e';        % Riboflavina
'R_EX_thm_e';           % Tiamina
'R_EX_4abz_e';          % 4-aminobenzoato
'R_EX_nac_e';           % Nicotinato
'R_EX_cbl1_e';          % Cianocobalamina
'R_EX_adocbl_e';        % Adenosilcobalamina
'R_EX_pheme_e';         % Protohemo
'R_EX_sheme_e';         % Sirohemo
'R_EX_2dmmq8_e';        % 2-demetilmenaquinona-8
'R_EX_hxan_e';          % Hipoxantina
'R_EX_mqn7_e';          % Menaquinona-7
'R_EX_mqn8_e';          % Menaquinona-8
'R_EX_ni2_e';           % Niquel
'R_EX_ocdca_e';         % Octadecanoato
'R_EX_q8_e';            % Ubiquinona-8
'R_EX_thymd_e';         % Timidina
'R_EX_spmd_e';          % Espermidina
'R_EX_pnto__R_e';       % Pantotenato
'R_EX_h2o__R_e';
'R_EX_nh4_e';           % Amonio (alternativo)
};

% Rates: aa=1, minerales/vitaminas=1000
n_aa = 20;  % primeros 20 entries de mZMB son aminoacidos
rates_mZMB = repmat(1000, length(mZMB), 1);
rates_mZMB(1:n_aa) = 1;

% Metabolitos de interes para reportar (SCFAs + gases) — labels en ingles
mets_clave = {
    'R_EX_but_e',      'EX_but_e',      'Butyrate';
    'R_EX_ac_e',       'EX_ac_e',       'Acetate';
    'R_EX_lac__D_e',   'EX_lac__D_e',   'D-Lactate';
    'R_EX_lac__L_e',   'EX_lac__L_e',   'L-Lactate';
    'R_EX_succ_e',     'EX_succ_e',     'Succinate';
    'R_EX_ppa_e',      'EX_ppa_e',      'Propionate';
    'R_EX_for_e',      'EX_for_e',      'Formate';
    'R_EX_etoh_e',     'EX_etoh_e',     'Ethanol';
    'R_EX_co2_e',      'EX_co2_e',      'CO2';
    'R_EX_h2_e',       'EX_h2_e',       'H2';
};

% Aminoacidos (L-aminoacidos, los biologicos)
aa_clave = {
    'R_EX_ala__L_e',   'EX_ala__L_e',   'Ala';
    'R_EX_arg__L_e',   'EX_arg__L_e',   'Arg';
    'R_EX_asn__L_e',   'EX_asn__L_e',   'Asn';
    'R_EX_asp__L_e',   'EX_asp__L_e',   'Asp';
    'R_EX_cys__L_e',   'EX_cys__L_e',   'Cys';
    'R_EX_gln__L_e',   'EX_gln__L_e',   'Gln';
    'R_EX_glu__L_e',   'EX_glu__L_e',   'Glu';
    'R_EX_gly_e',      'EX_gly_e',      'Gly';
    'R_EX_his__L_e',   'EX_his__L_e',   'His';
    'R_EX_ile__L_e',   'EX_ile__L_e',   'Ile';
    'R_EX_leu__L_e',   'EX_leu__L_e',   'Leu';
    'R_EX_lys__L_e',   'EX_lys__L_e',   'Lys';
    'R_EX_met__L_e',   'EX_met__L_e',   'Met';
    'R_EX_phe__L_e',   'EX_phe__L_e',   'Phe';
    'R_EX_pro__L_e',   'EX_pro__L_e',   'Pro';
    'R_EX_ser__L_e',   'EX_ser__L_e',   'Ser';
    'R_EX_thr__L_e',   'EX_thr__L_e',   'Thr';
    'R_EX_trp__L_e',   'EX_trp__L_e',   'Trp';
    'R_EX_tyr__L_e',   'EX_tyr__L_e',   'Tyr';
    'R_EX_val__L_e',   'EX_val__L_e',   'Val';
};

% Combinar todo para la simulacion
todos_clave = [mets_clave; aa_clave];

mets_nombres = todos_clave(:,3);
n_mets = size(todos_clave, 1);
n_scfa = size(mets_clave, 1);   % 10 primeros = SCFAs/gases
n_aa_plot = size(aa_clave, 1);  % 20 siguientes = aminoacidos
aa_nombres = aa_clave(:,3);

%% ========================================================================
%% MONOCULTIVOS FOS/INULINA
%% ========================================================================

fprintf('\n========================================\n');
fprintf('  MONOCULTIVOS FOS/INULINA\n');
fprintf('========================================\n');
fprintf('Estrategia: max biomasa (UB=mu_exp) -> pFBA L2\n');
fprintf('Bt,PT33,M62_1,Csym: FOS (kesto) LB=UB=-10\n');
fprintf('HGF2,M38: inulina LB=UB=-10\n');
fprintf('HGF2 mu=0.073 (fase tardia)\n');
fprintf('aa=1\n\n');

mono_tabla = table();

for i = 1:6
    [fila, sol] = simularMonocultivo(modelFiles{i}, shortNames{i}, mu_exp(i), ...
        mZMB, rates_mZMB, fuenteC{i}, tiene_inulina(i), tiene_kesto(i), ...
        todos_clave, mets_nombres, n_mets);
    mono_tabla = [mono_tabla; fila];
end

fprintf('\n=== TABLA MONOCULTIVOS ===\n');
disp(mono_tabla);

%% ========================================================================
%% EXPORTAR RESULTADOS
%% ========================================================================

writetable(mono_tabla, 'resultados_monocultivos_FOS_inulina_C.xlsx');
fprintf('Guardado: resultados_monocultivos_FOS_inulina_C.xlsx\n');

save('resultados_monocultivos_FOS_inulina_C.mat', 'mono_tabla');

%% ========================================================================
%% GRAFICOS DE BARRAS
%% ========================================================================

% SCFAs a graficar (sin formato, con propionate) — labels en ingles
mets_plot = {'Butyrate'; 'Acetate'; 'D-Lactate'; 'L-Lactate'; 'Succinate'; 'Propionate'};
n_plot = length(mets_plot);

% Colores: Butyrate=rojo, Acetate=azul, D-lac=verde oscuro, L-lac=verde claro,
%          Succinate=naranjo, Propionate=amarillo
colores_scfa = [
    0.8  0.0  0.0;   % Butyrate - rojo
    0.0  0.3  0.8;   % Acetate - azul
    0.0  0.5  0.0;   % D-Lactate - verde oscuro
    0.4  0.8  0.2;   % L-Lactate - verde claro
    1.0  0.5  0.0;   % Succinate - naranjo
    0.9  0.8  0.0;   % Propionate - amarillo
];

% Extraer datos
bacterias = mono_tabla.Bacteria;
n_bac = length(bacterias);
datos_flujos = zeros(n_bac, n_plot);
for j = 1:n_plot
    datos_flujos(:, j) = mono_tabla.(mets_plot{j});
end
mu_real = mono_tabla.mu_real;

% Aminoacidos: extraer flujos
datos_aa = zeros(n_bac, n_aa_plot);
for j = 1:n_aa_plot
    datos_aa(:, j) = mono_tabla.(aa_nombres{j});
end

% Etiquetas con fuente de carbono
etiquetas_bac = cell(n_bac, 1);
for i = 1:n_bac
    etiquetas_bac{i} = sprintf('%s (%s)', bacterias{i}, fuenteC_label{i});
end

% --- Filtrar aminoacidos con flujo positivo ---
aa_pos = datos_aa;
aa_pos(aa_pos < 0) = 0;
idx_aa_pos = find(any(aa_pos > 1e-6, 1));
n_aa_pos = length(idx_aa_pos);

% Colores aminoacidos (tonos violeta/rosa)
colores_aa_base = [
    0.6  0.2  0.8;   % violeta
    0.8  0.4  0.9;   % lila
    0.5  0.0  0.5;   % purpura oscuro
    0.9  0.2  0.6;   % rosa fuerte
    0.7  0.5  1.0;   % lavanda
    0.4  0.1  0.6;   % indigo
    0.9  0.6  0.8;   % rosa claro
    0.6  0.0  0.4;   % magenta oscuro
    0.8  0.3  0.7;   % fucsia
    0.5  0.3  0.9;   % azul-violeta
];
colores_aa = colores_aa_base(mod((1:n_aa_pos)-1, size(colores_aa_base,1))+1, :);

% --- Figura 1: dos paneles lado a lado ---
%   Izquierda : SCFAs / fermentaciones
%   Derecha   : aminoacidos secretados (flujo positivo)
fig1 = figure('Position', [50 50 1500 550], 'Color', 'w');

% Panel izquierdo: SCFAs
ax1 = subplot(1, 2, 1);
b1 = bar(ax1, datos_flujos);
for j = 1:n_plot
    b1(j).FaceColor = colores_scfa(j, :);
end
set(ax1, 'XTickLabel', etiquetas_bac, 'XTickLabelRotation', 35, ...
    'FontSize', 12, 'FontName', 'Arial', 'LineWidth', 1, 'Box', 'on', ...
    'TickLabelInterpreter', 'none');
ylabel(ax1, 'pFBA flux (mmol/gDW/h)', 'FontSize', 13, 'FontName', 'Arial');
title(ax1, 'Fermentation products', 'FontSize', 14, 'FontName', 'Arial');
legend(ax1, mets_plot, 'Location', 'best', 'FontSize', 10, 'FontName', 'Arial');
grid(ax1, 'on');

% Panel derecho: aminoacidos secretados
ax2 = subplot(1, 2, 2);
if n_aa_pos > 0
    datos_aa_filt   = max(datos_aa(:, idx_aa_pos), 0);  % solo flujos positivos (secrecion)
    nombres_aa_filt = aa_nombres(idx_aa_pos);
    b2 = bar(ax2, datos_aa_filt);
    for j = 1:n_aa_pos
        b2(j).FaceColor = colores_aa(j, :);
    end
    set(ax2, 'XTickLabel', etiquetas_bac, 'XTickLabelRotation', 35, ...
        'FontSize', 12, 'FontName', 'Arial', 'LineWidth', 1, 'Box', 'on', ...
        'TickLabelInterpreter', 'none');
    ylabel(ax2, 'pFBA flux (mmol/gDW/h)', 'FontSize', 13, 'FontName', 'Arial');
    title(ax2, 'Amino acid secretion', 'FontSize', 14, 'FontName', 'Arial');
    legend(ax2, nombres_aa_filt, 'Location', 'best', 'FontSize', 10, ...
        'FontName', 'Arial');
    grid(ax2, 'on');
    fprintf('  Aminoacidos secretados (%d): %s\n', n_aa_pos, ...
        strjoin(nombres_aa_filt, ', '));
else
    text(ax2, 0.5, 0.5, 'No amino acids secreted', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized', ...
        'FontSize', 14, 'FontName', 'Arial');
    title(ax2, 'Amino acid secretion', 'FontSize', 14, 'FontName', 'Arial');
    set(ax2, 'XTick', [], 'YTick', []);
    fprintf('  Ningun aminoacido se secreta (todos LB=0 / consumo)\n');
end

exportgraphics(fig1, 'fig_flujos_monocultivos_C.pdf', 'ContentType', 'vector');
exportgraphics(fig1, 'fig_flujos_monocultivos_C.png', 'Resolution', 600);
fprintf('Figura 1: SCFAs + AA secretados exportada (PDF + PNG 600dpi)\n');

% --- Figura 2: Tasa de crecimiento ---
fig2 = figure('Position', [100 100 600 400], 'Color', 'w');
b3 = bar(mu_real, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'k', 'LineWidth', 1);
set(gca, 'XTickLabel', etiquetas_bac, 'XTickLabelRotation', 30, ...
    'FontSize', 12, 'FontName', 'Arial', 'LineWidth', 1, 'Box', 'on', ...
    'TickLabelInterpreter', 'none');
ylabel('Growth rate (h^{-1})', 'FontSize', 13, 'FontName', 'Arial');
grid on;
hold on;
for i = 1:n_bac
    plot([i-0.4 i+0.4], [mu_exp(i) mu_exp(i)], 'r--', 'LineWidth', 2);
end
legend({'\mu_{predicted}', '\mu_{experimental}'}, 'Location', 'northeast', ...
    'FontSize', 12, 'FontName', 'Arial');
hold off;

exportgraphics(fig2, 'fig_crecimiento_monocultivos_C.pdf', 'ContentType', 'vector');
exportgraphics(fig2, 'fig_crecimiento_monocultivos_C.png', 'Resolution', 600);
fprintf('Figura 2: Crecimiento (ver GUI)\n');

% --- Tabla resumen completa para graficar ---
% SCFAs + aminoacidos en una sola tabla
varnames = ['Bacteria', 'FuenteC', 'mu_exp', 'mu_real', mets_plot', aa_nombres'];
datos_todos = [datos_flujos, datos_aa];
tabla_plot = table(bacterias, mono_tabla.FuenteC, mu_exp, mu_real, ...
    'VariableNames', {'Bacteria', 'FuenteC', 'mu_exp', 'mu_real'});
for j = 1:n_plot
    tabla_plot.(mets_plot{j}) = datos_flujos(:, j);
end
for j = 1:n_aa_plot
    tabla_plot.(aa_nombres{j}) = datos_aa(:, j);
end
fprintf('\n=== TABLA RESUMEN PARA GRAFICO ===\n');
disp(tabla_plot);
writetable(tabla_plot, 'tabla_resumen_grafico_C.xlsx');
fprintf('Guardado: tabla_resumen_grafico_C.xlsx\n');

%% RESUMEN
fprintf('\n========================================\n');
fprintf('RESUMEN\n');
fprintf('========================================\n');
fprintf('Monocultivos: %d\n', height(mono_tabla));
for i = 1:height(mono_tabla)
    fprintf('  %s: fuente=%s, mu_exp=%.3f, mu_real=%.4f\n', ...
        mono_tabla.Bacteria{i}, mono_tabla.FuenteC{i}, ...
        mono_tabla.mu_exp(i), mono_tabla.mu_real(i));
end

%% ========================================================================
%% FUNCIONES LOCALES
%% ========================================================================

function m = configurarMedio(m, mZMB, rates_mZMB, fuenteC, tieneInulina, tieneKesto)
    % Configura medio ZMB con fuente de carbono especificada
    % fuenteC: 'inulina' o 'kesto'

    % 1) Cerrar TODAS las reacciones de exchange (LB=0, UB=1000)
    ex_by_name = m.rxns(contains(m.rxns, 'EX_'));
    ex_by_struct = m.rxns(findExcRxns(m));
    ex_rxns = unique([ex_by_name; ex_by_struct]);
    for j = 1:length(ex_rxns)
        m = changeRxnBounds(m, ex_rxns{j}, 0, 'l');
        m = changeRxnBounds(m, ex_rxns{j}, 1000, 'u');
    end

    % 2) Abrir medio ZMB (aminoacidos aa=1, iones/vitaminas=1000)
    for j = 1:length(mZMB)
        rxn_id = mZMB{j};
        rxn_id_alt = regexprep(rxn_id, '^R_', '');
        idx = find(strcmp(m.rxns, rxn_id) | strcmp(m.rxns, rxn_id_alt));
        if ~isempty(idx)
            m = changeRxnBounds(m, m.rxns{idx(1)}, -rates_mZMB(j), 'l');
        end
    end

    % 3) Fuente de carbono
    if strcmp(fuenteC, 'inulina')
        if tieneInulina
            rxn_inu = find(strcmp(m.rxns, 'EX_inulin(e)') | ...
                           strcmp(m.rxns, 'R_EX_inulin_e') | ...
                           strcmp(m.rxns, 'EX_inulin_e'));
            if ~isempty(rxn_inu)
                m = changeRxnBounds(m, m.rxns{rxn_inu(1)}, -10, 'b');
                fprintf('  Inulina: %s LB=UB=-10\n', m.rxns{rxn_inu(1)});
            else
                fprintf('  ADVERTENCIA: reaccion de inulina no encontrada\n');
            end
        else
            fprintf('  No degrada inulina (fuente presente pero no consumible)\n');
        end

    elseif strcmp(fuenteC, 'kesto')
        if tieneKesto
            rxn_kesto = find(strcmp(m.rxns, 'EX_kesto(e)') | ...
                             strcmp(m.rxns, 'R_EX_kesto_e') | ...
                             strcmp(m.rxns, 'EX_kesto_e'));
            if ~isempty(rxn_kesto)
                m = changeRxnBounds(m, m.rxns{rxn_kesto(1)}, -10, 'b');
                fprintf('  Kesto/FOS: %s LB=UB=-10\n', m.rxns{rxn_kesto(1)});
            else
                fprintf('  ADVERTENCIA: reaccion de kesto no encontrada\n');
            end
        else
            fprintf('  No degrada kesto/FOS (fuente presente pero no consumible)\n');
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

function [fila, sol_out] = simularMonocultivo(modelFile, nombre, mu, ...
    mZMB, rates_mZMB, fuenteC, tieneInulina, tieneKesto, mets_clave, mets_nombres, n_mets)
    % Simula monocultivo: max biomasa (UB=mu_exp) -> pFBA L2

    fprintf('--- %s (fuente: %s, mu_exp=%.3f) ---\n', nombre, fuenteC, mu);

    % Cargar modelo
    data = load(modelFile);
    fn = fieldnames(data);
    m = data.(fn{1});

    % Configurar medio
    m = configurarMedio(m, mZMB, rates_mZMB, fuenteC, tieneInulina, tieneKesto);

    % Biomasa: LB=0, UB=mu_exp
    m = changeRxnBounds(m, 'Growth_biomass', 0, 'l');
    m = changeRxnBounds(m, 'Growth_biomass', mu, 'u');

    % Objetivo: max biomasa
    m = changeObjective(m, 'Growth_biomass', 1);

    % PASO 1: FBA max biomasa
    sol_raw = optimizeCbModel(m, 'max');

    if sol_raw.stat ~= 1
        fprintf('  FBA infactible\n\n');
        fila = crearFilaMono(nombre, fuenteC, mu, NaN, NaN(n_mets,1), mets_nombres);
        sol_out = sol_raw;
        return;
    end

    mu_opt = sol_raw.v(findRxnIDs(m, 'Growth_biomass'));
    fprintf('  PASO 1 (max bio): mu=%.4f\n', mu_opt);

    % PASO 2: Fijar biomasa al optimo -> pFBA L2
    m = changeRxnBounds(m, 'Growth_biomass', mu_opt, 'b');
    sol = optimizeCbModel(m, 'max', 1e-6);  % L2

    if sol.stat ~= 1
        fprintf('  pFBA L2 infactible\n\n');
        fila = crearFilaMono(nombre, fuenteC, mu, mu_opt, NaN(n_mets,1), mets_nombres);
        sol_out = sol;
        return;
    end

    mu_L2 = sol.v(findRxnIDs(m, 'Growth_biomass'));
    fprintf('  PASO 2 (pFBA L2): mu=%.4f\n', mu_L2);

    % Extraer flujos de metabolitos clave
    flujos_mets = zeros(n_mets, 1);
    for j = 1:n_mets
        idx_r = find(strcmp(m.rxns, mets_clave{j,1}) | strcmp(m.rxns, mets_clave{j,2}));
        if ~isempty(idx_r)
            flujos_mets(j) = sol.v(idx_r(1));
        end
    end

    % Mostrar metabolitos con flujo
    for j = 1:n_mets
        if abs(flujos_mets(j)) > 1e-6
            fprintf('    %-12s: %8.4f\n', mets_nombres{j}, flujos_mets(j));
        end
    end

    % ATPM
    idx_atpm = find(strcmp(m.rxns, 'R_ATPM') | strcmp(m.rxns, 'ATPM'));
    if ~isempty(idx_atpm)
        fprintf('    ATPM:         %8.4f (LB=%.2f, UB=%.2f)\n', ...
            sol.v(idx_atpm(1)), m.lb(idx_atpm(1)), m.ub(idx_atpm(1)));
    end

    fprintf('\n');

    fila = crearFilaMono(nombre, fuenteC, mu, mu_L2, flujos_mets, mets_nombres);
    sol_out = sol;
end

function fila = crearFilaMono(nombre, fuenteC, mu_exp, mu_real, flujos, nombres_mets)
    fila = table();
    fila.Bacteria = {nombre};
    fila.FuenteC = {fuenteC};
    fila.mu_exp = mu_exp;
    fila.mu_real = mu_real;
    for j = 1:length(nombres_mets)
        fila.(nombres_mets{j}) = flujos(j);
    end
end
