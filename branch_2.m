% =========================================================================
% Assignment 1 - Punto 1.e: Studio di Convergenza Ritz-Galerkin
% =========================================================================
clear; clc; close all;

%% 1. Parametri di input (Dati dal problema)
EJ = 1e7;              % [Nm^2] Rigidezza flessionale
GJ = 1e6;              % [Nm^2] Rigidezza torsionale
EA = 1e8;              % [N] Rigidezza assiale puntone
C_Lalpha = 2*pi;       % [-] Pendenza retta portanza
b = 14;                % [m] Apertura alare
c = 1;                 % [m] Corda
e = c/4;               % [m] Distanza AC-EA (Positivo se AC avanti a EA)
gamma = 15 * pi/180;   % [rad] Angolo del puntone
xB = -0.33;            % [m] Posizione della cerniera B

% Calcolo rigidezza equivalente del puntone lungo l'asse z
L_strut = (b/2) / cos(gamma);            % Lunghezza fisica del puntone
k_eq_z = (EA / L_strut) * sin(gamma)^2;  % Rigidezza trasversale equivalente

%% 2. Impostazioni dello studio di convergenza
N_max = 13; % Numero massimo di funzioni di forma da testare (per w e per theta)
qD_history = zeros(1, N_max); % Vettore per salvare le pressioni di divergenza
N_values = 1:N_max;

%% 3. Ciclo Iterativo sulle funzioni di forma
for N = N_values
    
    % Inizializzazione delle sottomatrici (N x N)
    K_ww = zeros(N, N); % Rigidezza flessionale
    K_tt = zeros(N, N); % Rigidezza torsionale
    K_wt = zeros(N, N); % Accoppiamento strut
    
    A_wt = zeros(N, N); % Lavoro portanza su flessione
    A_tt = zeros(N, N); % Lavoro momento su torsione
    
    % Riempimento delle sottomatrici
    for i = 1:N
        for j = 1:N
            % Nota sugli integrali: usiamo la coordinata adimensionale eta = y/b
            % L'integrale da 0 a b di (y/b)^A dy equivale a b / (A + 1)
            
            % -- MATRICI STRUTTURALI (Ala Pulita) --
            pot_ww = (i-1) + (j-1); 
            K_ww(i,j) = EJ * ((i+1)*i/b^2) * ((j+1)*j/b^2) * (b / (pot_ww + 1));
            
            pot_tt_strut = (i-1) + (j-1);
            K_tt(i,j) = GJ * (i/b) * (j/b) * (b / (pot_tt_strut + 1));
            
            % -- MATRICI AERODINAMICHE --
            pot_wt = (i+1) + j; 
            A_wt(i,j) = c * C_Lalpha * (b / (pot_wt + 1));
            
            pot_tt_aero = i + j;
            A_tt(i,j) = c * e * C_Lalpha * (b / (pot_tt_aero + 1));
            
            % -- CONTRIBUTO DEL PUNTONE (Valutato in y = b/2 -> eta = 0.5) --
            phi_wi_b2 = (0.5)^(i+1);
            phi_wj_b2 = (0.5)^(j+1);
            phi_ti_b2 = (0.5)^i;
            phi_tj_b2 = (0.5)^j;
            
            K_ww(i,j) = K_ww(i,j) + k_eq_z * phi_wi_b2 * phi_wj_b2;
            K_tt(i,j) = K_tt(i,j) + k_eq_z * (xB^2) * phi_ti_b2 * phi_tj_b2;
            K_wt(i,j) = -k_eq_z * xB * phi_wi_b2 * phi_tj_b2;
        end
    end
    
    % 4. Assemblaggio del Sistema Globale (Matrici 2N x 2N)
    K_global = [K_ww,       K_wt; 
                K_wt',      K_tt];  % La trasposta assicura la simmetria
            
    A_global = [zeros(N,N), A_wt; 
                zeros(N,N), A_tt];
            
    % 5. Risoluzione del Problema agli Autovalori (Divergenza)
    [V, D] = eig(K_global, A_global);
    autovalori = diag(D);
    
    % Filtra solo gli autovalori REALI, NEGATIVI e scarta gli INFINITI (~isinf)
    valori_negativi = autovalori(imag(autovalori) == 0 & real(autovalori) < 0 & ~isinf(autovalori));
    
    if isempty(valori_negativi)
        qD_history(N) = NaN; % Nessuna radice negativa trovata per questo N
    else
        % Tra i valori negativi, cerchiamo quello con valore assoluto più piccolo,
        % ovvero quello più vicino a zero. Per numeri negativi usiamo max()
        % (es: tra -5000 e -150, max() restituisce -150).
        qD_history(N) = max(valori_negativi); 
    end
end

%% 6. Plot dei Risultati
figure('Color', 'w');
plot(N_values, qD_history, '-o', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'Color', 'r');
grid on;
set(gca, 'FontSize', 12);
xlabel('Numero di funzioni di forma (N) per DoF', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Pressione di Divergenza q_D (Ramo Negativo) [Pa]', 'FontSize', 14, 'FontWeight', 'bold');
title('Studio di Convergenza per q_D negativa (Ritz-Galerkin)', 'FontSize', 16);

% Stampa in console i valori
disp('--- Risultati dello Studio di Convergenza (Radici Negative) ---');
for i = 1:N_max
    if isnan(qD_history(i))
        fprintf('N = %d -> q_D = Nessuna radice negativa trovata\n', i);
    else
        fprintf('N = %d -> q_D = %.2f Pa\n', i, qD_history(i));
    end
end