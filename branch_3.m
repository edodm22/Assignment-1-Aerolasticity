%% ============================================================
%  ANALISI PRESSIONE DI DIVERGENZA E LUNGHEZZA (Ritz-Galerkin N-DoF)
%% ============================================================
% xB fissato, punto A (attacco a terra) fisso.
% Analisi del trade-off tra prestazioni aeroelastiche e massa/ingombro.
clear; clc; close all;

% Parametri strutturali e aerodinamici
param.EJ = 1e7;
param.GJ = 1e6;
param.EA = 1e8;
param.CL_a = 2*pi;
param.b = 14;
param.c = 1;
param.e = param.c/4;
param.gamma = 15; % Baseline
xB = -0.1;
N = 5; % Livello di approssimazione Ritz-Galerkin (4 garantisce convergenza)

% Definizione vincoli geometrici (A fisso)
gamma_ref = param.gamma;
yB_ref = param.b/2;
zA = yB_ref * tan(deg2rad(gamma_ref)); % Quota fissa al muro

% Vettori di iterazione (MODIFICATO QUI)
gamma_vec_deg = 15:1:40; % Discretizzato da 15 a 40 ogni 1 grado
q_div_vec = zeros(size(gamma_vec_deg));
lunghezza_barra = zeros(size(gamma_vec_deg));

for step = 1:length(gamma_vec_deg)
    g = deg2rad(gamma_vec_deg(step));
    
    % 1. Geometria e Vincolo yB
    yB = zA / tan(g);
    lunghezza_barra(step) = zA / sin(g); 
    
    % Se il puntone finisce fuori dall'ala, scarta il dato (con tolleranza numerica)
    if yB > (param.b/2 + 1e-6) || yB <= 0
        q_div_vec(step) = NaN;
        lunghezza_barra(step) = NaN;
        continue;
    end
    
    etaB = yB / param.b;
    keq_z = (param.EA * sin(g)^3) / zA;
    
    % 2. Inizializzazione Matrici N x N
    K_ww = zeros(N, N); K_tt = zeros(N, N); K_wt = zeros(N, N);
    A_wt = zeros(N, N); A_tt = zeros(N, N);
    
    % 3. Riempimento Matrici Galerkin
    for i = 1:N
        for j = 1:N
            pot_ww = (i-1) + (j-1); 
            K_ww(i,j) = param.EJ * ((i+1)*i/param.b^2) * ((j+1)*j/param.b^2) * (param.b / (pot_ww + 1));
            
            pot_tt = (i-1) + (j-1);
            K_tt(i,j) = param.GJ * (i/param.b) * (j/param.b) * (param.b / (pot_tt + 1));
            
            pot_wt = (i+1) + j; 
            A_wt(i,j) = param.c * param.CL_a * (param.b / (pot_wt + 1));
            
            pot_tt_aero = i + j;
            A_tt(i,j) = param.c * param.e * param.CL_a * (param.b / (pot_tt_aero + 1));
            
            phi_wi_B = etaB^(i+1); phi_wj_B = etaB^(j+1);
            phi_ti_B = etaB^i;     phi_tj_B = etaB^j;
            
            K_ww(i,j) = K_ww(i,j) + keq_z * phi_wi_B * phi_wj_B;
            K_tt(i,j) = K_tt(i,j) + keq_z * (xB^2) * phi_ti_B * phi_tj_B;
            K_wt(i,j) = -keq_z * xB * phi_wi_B * phi_tj_B;
        end
    end
    
    % 4. Assemblaggio Globale e Risoluzione
    K_global = [K_ww, K_wt; K_wt', K_tt];
    A_global = [zeros(N,N), A_wt; zeros(N,N), A_tt];
    
    [~, D] = eig(K_global, A_global);
    autovalori = diag(D);
    
    % Estrazione della radice critica positiva
    valori_fisici = autovalori(imag(autovalori) == 0 & real(autovalori) > 0 & ~isinf(autovalori));
    
    if isempty(valori_fisici)
        q_div_vec(step) = NaN;
    else
        q_div_vec(step) = min(valori_fisici);
    end
end

%% Plotting con Doppio Asse Y
figure('Color','w','Name','Analisi Multidisciplinare Strallo');
grid on; hold on;

% --- ASSE SINISTRO: Pressione di Divergenza ---
yyaxis left
plot(gamma_vec_deg, q_div_vec, '-ob', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'DisplayName', 'Pressione di Divergenza q_D');
ylabel('Divergence Pressure q_D [Pa]', 'FontWeight', 'bold');
set(gca, 'YColor', 'b'); 

% --- ASSE DESTRO: Lunghezza della barra ---
yyaxis right
plot(gamma_vec_deg, lunghezza_barra, '--or', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', 'Lunghezza Strallo L_s');
ylabel('Strut Length L_s [m]', 'FontWeight', 'bold');
set(gca, 'YColor', 'r'); 

% --- Formattazione Generale ---
xlabel('Strut Angle \gamma [deg]', 'FontWeight', 'bold');
title(['Divergenza e Geometria vs \gamma (x_B = ', num2str(xB), ', N = ', num2str(N), ')'], 'FontSize', 12);
legend('Location', 'north');

% Limite fisico apertura alare (corrisponde esattamente a 15 gradi)
xline(rad2deg(atan(zA/(param.b/2))), ':k', 'y_B = b/2 Limit', 'LineWidth', 1.5);