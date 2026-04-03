clear all; 
close all; 
clc;

%% ============================================================
%                    PARAMETERS STRUCT
%% ============================================================
% All physical and geometrical parameters are grouped here
% to keep the code clean and easily modifiable.

param.EJ = 1e7;          % Bending stiffness of the wing [Nm^2]
param.GJ = 1e6;          % Torsional stiffness of the wing [Nm^2]
param.EA = 1e8;          % Axial stiffness of the strut [N]

param.CL_a = 2*pi;       % Lift curve slope [1/rad]
param.CL_beta = pi;      % Lift variation due to control deflection
param.CM_beta = -0.1;    % Moment coefficient due to control

param.b = 14;            % Wing span [m]
param.c = 1;             % Chord length [m]
param.e = param.c/4;     % Distance between Aerodynamic Center and Elastic Axis

param.gamma = deg2rad(15); % Reference strut angle [rad]
param.ba = 2.5;           % Aileron span [m]
beta0 = deg2rad(10);     % Input for Point 3


%% ============================================================
%                      POINT 1
%             AEROELASTIC DIVERGENCE ANALYSIS
%% ============================================================
% Compute the divergence dynamic pressure q_D as a function
% of the chordwise position of the strut attachment (x_B).

% ------------------------------------------------------------
% Equivalent stiffness of the strut in vertical direction
% ------------------------------------------------------------
% The strut provides an additional elastic contribution.
% Only the vertical component contributes to bending-torsion coupling.

k_eqz = (param.EA * cos(param.gamma) / (param.b/2)) * (sin(param.gamma)^2);

% ------------------------------------------------------------
% Sweep of strut position along the chord
% ------------------------------------------------------------
% xB = -0.5 → Leading Edge
% xB =  0   → Elastic Axis
% xB =  0.5 → Trailing Edge

xB = linspace(-0.5, 0.5, 1800);

% ------------------------------------------------------------
% Structural stiffness terms (independent of q)
% ------------------------------------------------------------
K11 = 4*param.EJ/param.b^3 + k_eqz/16;   % bending stiffness term
K22_base = param.GJ/param.b;             % torsional stiffness (base)
K12_coeff = -k_eqz/8;                    % coupling coefficient

% ------------------------------------------------------------
% Aerodynamic influence coefficients
% ------------------------------------------------------------
% These terms come from the aerodynamic matrix
% and represent how lift and moment depend on deformation.

A22 = param.e * param.b * param.c * param.CL_a / 3;
A12 = param.b * param.c * param.CL_a / 4;

% ------------------------------------------------------------
% Divergence condition
% ------------------------------------------------------------
% Divergence occurs when det(K - q*A) = 0
% Here we compute q explicitly as a ratio:
%   q_D = (structural terms) / (aerodynamic terms)

Num = K11*(K22_base + k_eqz*xB.^2/4) - (K12_coeff*xB).^2;
D   = K11*(A22) - (K12_coeff*xB).*(A12);

q_div = Num ./ D;

% ------------------------------------------------------------
% Plot results
% ------------------------------------------------------------
figure('Color','w','Name','Divergence Pressure');
plot(xB, q_div, 'b','LineWidth',2);
grid on; hold on;

yline(0,'k','HandleVisibility','off');
xline(0,'--r','Elastic Axis');

xlabel('x_B [-]');
ylabel('q_D [Pa]');
title('Divergence Pressure vs Strut Position');

text(-0.48, 4e5, 'Leading Edge','FontWeight','bold');
text(0.3, 4e5, 'Trailing Edge','FontWeight','bold');

% ------------------------------------------------------------
% Critical position
% ------------------------------------------------------------
% This corresponds to a singularity (denominator = 0),
% i.e. where divergence pressure tends to infinity.

x_critic = (K11 * A22) / (K12_coeff * A12);

%% ============================================================
%                      POINT 2
%        CONTROL EFFECTIVENESS (FIXED POINT A)
%% ============================================================
% Study how control effectiveness E_C varies with dynamic pressure
% and with the strut angle gamma.

% ------------------------------------------------------------
% Fixed point A condition
% ------------------------------------------------------------
% The geometry is constrained so that point A is fixed.
% This introduces a relation between zA and gamma.

gamma_ref = param.gamma;
yB_ref = param.b/2;
zA = yB_ref * tan(gamma_ref);

% ------------------------------------------------------------
% Choose a specific chordwise position for the strut
% ------------------------------------------------------------
xB = -0.33; 

% ------------------------------------------------------------
% Define parameter ranges
% ------------------------------------------------------------
gamma_deg = [10, 15, 20];       % different strut inclinations
gamma_deg = linspace(10,20,5); 
q_vec = linspace(0, 3e5, 1000); % dynamic pressure range

figure('Color','w'); hold on; grid on;
colors = {'r','b','g', 'c', 'm'};

for j = 1:length(gamma_deg)

    g = deg2rad(gamma_deg(j));

    % --------------------------------------------------------
    % Geometry update due to gamma variation
    % --------------------------------------------------------
    % Since A is fixed, changing gamma changes yB

    yB = zA / tan(g);
    etaB = yB / param.b;   % nondimensional position

    % --------------------------------------------------------
    % Equivalent stiffness (nonlinear with gamma)
    % --------------------------------------------------------
    keq_z = (param.EA * sin(g)^3) / zA;

    % --------------------------------------------------------
    % Structural stiffness matrix K
    % --------------------------------------------------------
    % Includes coupling between bending and torsion

    K11 = 4*param.EJ/param.b^3 + keq_z * etaB^4;
    K22 = param.GJ/param.b + keq_z * xB^2 * etaB^2;
    K12 = -keq_z * xB * etaB^3;

    K = [K11, K12;
         K12, K22];

    % --------------------------------------------------------
    % Aerodynamic matrix
    % --------------------------------------------------------
    A = [0, param.b*param.c*param.CL_a/4;
         0, param.e*param.b*param.c*param.CL_a/3];

    % --------------------------------------------------------
    % Aileron forcing vector
    % --------------------------------------------------------
    % Represents aerodynamic loads due to control deflection

    F1b = (param.c * param.CL_beta / (3*param.b^2)) * ...
          (param.b^3 - (param.b-param.ba)^3);

    F2b = (param.c^2 * param.CM_beta / (2*param.b) + ...
           param.e*param.c*param.CL_beta / (2*param.b)) * ...
          (param.b^2 - (param.b-param.ba)^2);

    F_beta = [F1b; F2b];

    % --------------------------------------------------------
    % Solve system for each dynamic pressure
    % --------------------------------------------------------
    % Equation:
    %   (K - q*A) * q_vec = q * F_beta

    Ec = zeros(size(q_vec));

    for i = 1:length(q_vec)

        q = q_vec(i);

        % Check stability using eigenvalues of the aeroelastic system
        eig_vals = eig(K - q*A);


        % Only compute solution if system is still stable
        % metto questa condizoni in quanto voglio analizzare quando per
        % ciascun caso di angoo gamma ottengo la divergenza (quindi quando
        % il determiante della matrice [K - q*A] diventa zero => ovvero
        % quando almeno uno degli autovalori smette di essere positivo) 
        

        % !! decommenta se vuoigrafici completi senza intrruzione (e
        % commenta sotto) !!
        % if all(eig_vals > 0) 
        %     sol = (K - q*A) \ (q * F_beta);
        %     theta = sol(2); % torsional response
        % 
        %     % Control effectiveness
        %     Ec(i) = 1 + (param.b * param.CL_a * theta) / ...
        %                 (2 * param.CL_beta * param.ba);
        % 
        % else
        %     % Outside physical region (post-divergence)
        %     Ec(i) = NaN;
        % end


       
        % !! decommenta se vuoigrafici completi senza intrruzione (e
        % commenta sopra) !! 
        sol = (K - q*A) \ (q * F_beta);
            theta = sol(2); % torsional response

            
            % Control Effectiveness expression
            Ec(i) = 1 + (param.b * param.CL_a * theta) / ...
            (2 * param.CL_beta * beta0 * param.ba);
    end

    % Plot curve
    plot(q_vec, Ec, 'Color', colors{j}, 'LineWidth',2, ...
        'DisplayName',['\gamma = ', num2str(gamma_deg(j)), '°']);
end

% ------------------------------------------------------------
% Final plot formatting
% ------------------------------------------------------------
yline(0,'k','HandleVisibility','off');
yline(1,'k--','Rigid Wing','HandleVisibility','off');

xlabel('Dynamic Pressure q [Pa]');
ylabel('Control Effectiveness E_C');
title(['Control Effectiveness (x_B = ', num2str(xB), ')']);

legend('Location','best');


%% ============================================================
%                      POINT 3
%             DEFORMED SHAPE COMPARISON
%% ============================================================
% 1. Reference: No Strut Divergence
qD_no_strut = (3 * param.GJ) / (param.e * param.b^2 * param.c * param.CL_a); %formula ricavata a mano 
q_test = 0.4 * qD_no_strut; 

% 2. System WITHOUT Strut
K_no = [4*param.EJ/param.b^3, 0; 0, param.GJ/param.b];
sol_no = (K_no - q_test*A) \ (q_test * F_beta); % A rimane invariata

% 3. System WITH Strut (using gamma = 15 deg)
g15 = deg2rad(15); % per gambiare gamma cambia qua! 
yB15 = zA / tan(g15); 
etaB15 = yB15 / param.b;
keq15 = (param.EA * sin(g15)^3) / zA;
K_strut = [4*param.EJ/param.b^3 + keq15*etaB15^4, -keq15*xB*etaB15^3;
           -keq15*xB*etaB15^3, param.GJ/param.b + keq15*xB^2*etaB15^2];
sol_strut = (K_strut - q_test*A) \ (q_test * F_beta);

% 4. Plotting Deformed Shapes
y_span = linspace(0, param.b, 100);
w_no = (y_span/param.b).^2 * sol_no(1) * beta0;
t_no = (y_span/param.b) * sol_no(2) * beta0;
w_strut = (y_span/param.b).^2 * sol_strut(1) * beta0;
t_strut = (y_span/param.b) * sol_strut(2) * beta0;

figure('Color','w','Name','Deformed Shape Comparison');
subplot(2,1,1);
plot(y_span, w_no, 'r--', 'LineWidth', 2); hold on;
plot(y_span, w_strut, 'b', 'LineWidth', 2);
ylabel('Bending w(y) [m]'); grid on; legend('No Strut','With Strut');
title(['Wing Deformation at q = ', num2str(round(q_test)), ' Pa']);

subplot(2,1,2);
plot(y_span, rad2deg(t_no), 'r--', 'LineWidth', 2); hold on;
plot(y_span, rad2deg(t_strut), 'b', 'LineWidth', 2);
ylabel('Twist \theta(y) [deg]'); xlabel('y [m]'); grid on;


%% ============================================================
%                      POINT 4
%       CONTROL REVERSAL vs STRUT POSITION (N=1)
%% ============================================================ 
% Find the dynamic pressure q_R at which control effectiveness 
% becomes zero, as a function of the strut position xB.

% 1. Setup Base Matrices (N=1)
Kw_clean = [4*param.EJ/param.b^3, 0; 0, param.GJ/param.b];
Ka = [0, param.b*param.c*param.CL_a/4; 
      0, param.e*param.b*param.c*param.CL_a/3];

% Aileron forcing vector (Normalized per unit beta)
F1b_unit = (param.c * param.CL_beta / (3*param.b^2)) * ...
           (param.b^3 - (param.b-param.ba)^3);
F2b_unit = (param.c^2 * param.CM_beta / (2*param.b) + ...
           param.e*param.c*param.CL_beta / (2*param.b)) * ...
           (param.b^2 - (param.b-param.ba)^2);
Fa_unit = [F1b_unit; F2b_unit];

% 2. Reversal Aerodynamic Matrix (Ka_rev)
v_Lalpha = [0; param.c * param.CL_a * param.b / 2];
L_beta_unit = param.c * param.CL_beta * param.ba;
Ka_rev = Ka - (Fa_unit * v_Lalpha') / L_beta_unit;

% 3. Geometry (Fixed Point A logic)
zA = (param.b/2) * tan(param.gamma);
yB_base = param.b/2;
etaB = yB_base / param.b;
keq_z = (param.EA * sin(param.gamma)^3) / zA;

% 4. Analytic Asymptote (N=1 logic)
Ks_11 = Kw_clean(1,1) + keq_z * etaB^4;
K_slope = -keq_z * etaB^3; 
xB_crit_R = (Ks_11 * Ka_rev(2,2)) / (K_slope * Ka_rev(1,2));

% 5. Baseline: Unbraced Wing Reversal
[~, D_unbraced] = eig(Kw_clean, Ka_rev);
q_eig_un = diag(D_unbraced);
valid_q_un = real(q_eig_un(abs(imag(q_eig_un)) < 1e-5 & q_eig_un > 0));
qR_unbraced = min(valid_q_un);

% 6. Sweep over Strut Position xB
xB_vec_R = linspace(-param.c/2, param.c/2, 1000);
qR_vec = zeros(size(xB_vec_R));

for i = 1:length(xB_vec_R)
    xB_curr = xB_vec_R(i);
    K_strut = [keq_z * etaB^4,            -keq_z * xB_curr * etaB^3;
              -keq_z * xB_curr * etaB^3,   keq_z * xB_curr^2 * etaB^2];
    Ks = Kw_clean + K_strut; 
    [~, D_braced] = eig(Ks, Ka_rev);
    q_eig_br = diag(D_braced);
    valid_q_br = real(q_eig_br(abs(imag(q_eig_br)) < 1e-5));
    if isempty(valid_q_br)
        qR_vec(i) = NaN; 
    else
        [~, idx_br] = min(abs(valid_q_br));
        qR_vec(i) = valid_q_br(idx_br);
    end
end

% 7. Plotting (Light Mode - White Background)
figure('Name', 'Point 4 - Control Reversal (N=1)', 'Position', [200, 200, 800, 550], 'Color', 'w');
ax = gca; hold on; grid on;
% Impostiamo lo sfondo bianco e gli assi neri
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
    'GridColor', [0.8 0.8 0.8], 'TickLabelInterpreter', 'latex');

% Baseline (Unbraced)
plot([-0.5, 0.5], [qR_unbraced, qR_unbraced], 'r--', 'LineWidth', 2.5, 'DisplayName', '$q_R$ unbraced');

% Asintoto verticale
xline(xB_crit_R, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Gestione del salto per il plot
[~, split_idx] = min(abs(xB_vec_R - xB_crit_R));
qR_vec(split_idx) = NaN; 

% Ala con puntone (Braced)
plot(xB_vec_R, qR_vec, '-', 'Color', [0 0.447 0.741], 'LineWidth', 3, 'DisplayName', '$q_R$ braced');

% Linea dello zero
plot([-0.5, 0.5], [0, 0], 'k', 'LineWidth', 1, 'HandleVisibility', 'off'); 

% Titoli e label in nero
xlabel('Strut chordwise position $x_B$ [m]', 'FontSize', 12, 'Interpreter', 'latex', 'Color', 'k');
ylabel('$q_R$ [Pa]', 'FontSize', 12, 'Interpreter', 'latex', 'Color', 'k');
title('\textbf{4: Control Reversal vs Strut Position ($N=1$)}', 'FontSize', 14, 'Interpreter', 'latex', 'Color', 'k');

xlim([-0.5, 0.5]); 
ylim([-1e6, 1e6]); 
legend('Location', 'northeast', 'FontSize', 11, 'Interpreter', 'latex');

% =========================================================================
% Punto 6: Studio di Convergenza Ritz-Galerkin
% =========================================================================
%% 1. Parametri di input (Dati dal problema)
EJ = 1e7;              % [Nm^2] Rigidezza flessionale
GJ = 1e6;              % [Nm^2] Rigidezza torsionale
EA = 1e8;              % [N] Rigidezza assiale puntone
C_Lalpha = 2*pi;       % [-] Pendenza retta portanza
b = 14;                % [m] Apertura alare
c = 1;                 % [m] Corda
e = c/4;               % [m] Distanza AC-EA (Positivo se AC avanti a EA)
gamma = 15 * pi/180;   % [rad] Angolo del puntone
xB = -0;                % [m] Posizione della cerniera B (Da adattare se diverso)

% Calcolo rigidezza equivalente del puntone lungo l'asse z
L_strut = (b/2) / cos(gamma);            % Lunghezza fisica del puntone
k_eq_z = (EA / L_strut) * sin(gamma)^2;  % Rigidezza trasversale equivalente

%% 2. Impostazioni dello studio di convergenza
N_max = 8; % Numero massimo di funzioni di forma da testare (per w e per theta)
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
            % K_ww: Integrale( EJ * phi_wi'' * phi_wj'' dy )
            % phi_wi'' = (i+1)*i/b^2 * eta^(i-1)
            pot_ww = (i-1) + (j-1); 
            K_ww(i,j) = EJ * ((i+1)*i/b^2) * ((j+1)*j/b^2) * (b / (pot_ww + 1));
            
            % K_tt: Integrale( GJ * phi_ti' * phi_tj' dy )
            % phi_ti' = i/b * eta^(i-1)
            pot_tt_strut = (i-1) + (j-1);
            K_tt(i,j) = GJ * (i/b) * (j/b) * (b / (pot_tt_strut + 1));
            
            % -- MATRICI AERODINAMICHE --
            % A_wt: Lavoro Portanza su w -> Integrale( c * CLalpha * phi_wj * phi_ti dy )
            pot_wt = (i+1) + j; 
            A_wt(i,j) = c * C_Lalpha * (b / (pot_wt + 1));
            
            % A_tt: Lavoro Momento su theta -> Integrale( c * e * CLalpha * phi_ti * phi_tj dy )
            pot_tt_aero = i + j;
            A_tt(i,j) = c * e * C_Lalpha * (b / (pot_tt_aero + 1));
            
            % -- CONTRIBUTO DEL PUNTONE (Valutato in y = b/2 -> eta = 0.5) --
            phi_wi_b2 = (0.5)^(i+1);
            phi_wj_b2 = (0.5)^(j+1);
            phi_ti_b2 = (0.5)^i;
            phi_tj_b2 = (0.5)^j;
            
            % Aggiunta dei termini discreti del puntone alle matrici
            K_ww(i,j) = K_ww(i,j) + k_eq_z * phi_wi_b2 * phi_wj_b2;
            K_tt(i,j) = K_tt(i,j) + k_eq_z * (xB^2) * phi_ti_b2 * phi_tj_b2;
            
            % Termine di accoppiamento flesso-torsionale (K_wt e K_tw)
            K_wt(i,j) = -k_eq_z * xB * phi_wi_b2 * phi_tj_b2;
        end
    end
    
    % 4. Assemblaggio del Sistema Globale (Matrici 2N x 2N)
    % Ordine dei gradi di libertà: [qw1..qwN, qtheta1..qthetaN]'
    K_global = [K_ww,       K_wt; 
                K_wt',      K_tt];  % La trasposta assicura la simmetria
            
    A_global = [zeros(N,N), A_wt; 
                zeros(N,N), A_tt];
            
    % 5. Risoluzione del Problema agli Autovalori (Divergenza)
    % (K_global - qD * A_global) * {q} = 0
    [V, D] = eig(K_global, A_global);
    autovalori = diag(D);
    
    % Filtra solo gli autovalori fisici (reali e strettamente positivi)
    valori_fisici = autovalori(imag(autovalori) == 0 & real(autovalori) > 0);
    
    if isempty(autovalori)
        qD_history(N) = NaN; % Nessuna instabilità trovata
    else
        qD_history(N) = min(valori_fisici); % q_D è il valore critico più basso
    end
end

%% 6. Plot dei Risultati
figure('Color', 'w');
plot(N_values, qD_history, '-o', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
set(gca, 'FontSize', 12);
xlabel('Numero di funzioni di forma (N) per DoF', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Pressione di Divergenza q_D [Pa]', 'FontSize', 14, 'FontWeight', 'bold');
title('Studio di Convergenza (Metodo di Ritz-Galerkin)', 'FontSize', 16);

% Stampa in console i valori
disp('--- Risultati dello Studio di Convergenza ---');
for i = 1:N_max
    fprintf('N = %d -> q_D = %.2f Pa\n', i, qD_history(i));
end


%% ============================================================
%                      POINT 5
%    MULTIDISCIPLINARY ANALYSIS: DIVERGENCE vs GEOMETRY
%% ============================================================
% Analyze the trade-off between q_D and strut length L_s.

% --- 1. Parameters for this Study ---
xB_local = -0.1;                 % Fixed chordwise position
N_galerkin = 5;                  % Ritz-Galerkin order

% --- 2. Geometric Constraints (Fixed Point A) ---
gamma_ref_rad = deg2rad(15); 
yB_ref = param.b / 2;
zA_fixed = yB_ref * tan(gamma_ref_rad); 

% --- 3. Iteration Vectors ---
gamma_vec_deg = 15:1:40;         
q_div_history = zeros(size(gamma_vec_deg));
strut_length = zeros(size(gamma_vec_deg));

% --- 4. Main Loop over Gamma ---
for step = 1:length(gamma_vec_deg)
    g_rad = deg2rad(gamma_vec_deg(step));
    
    % Geometry
    yB_curr = zA_fixed / tan(g_rad);
    strut_length(step) = zA_fixed / sin(g_rad);
    
    if yB_curr > (param.b/2 + 1e-6) || yB_curr <= 0
        q_div_history(step) = NaN;
        strut_length(step) = NaN;
        continue;
    end
    
    etaB_curr = yB_curr / param.b;
    keq_z_curr = (param.EA * sin(g_rad)^3) / zA_fixed;
    
    % Galerkin matrices
    K_ww = zeros(N_galerkin, N_galerkin); K_tt = zeros(N_galerkin, N_galerkin); 
    K_wt = zeros(N_galerkin, N_galerkin);
    A_wt = zeros(N_galerkin, N_galerkin); A_tt = zeros(N_galerkin, N_galerkin);
    
    for i = 1:N_galerkin
        for j = 1:N_galerkin
            pow_ww = (i-1) + (j-1);
            K_ww(i,j) = param.EJ * ((i+1)*i/param.b^2) * ((j+1)*j/param.b^2) * (param.b / (pow_ww + 1));
            pow_tt = (i-1) + (j-1);
            K_tt(i,j) = param.GJ * (i/param.b) * (j/param.b) * (param.b / (pow_tt + 1));
            pow_awt = (i+1) + j;
            A_wt(i,j) = param.c * param.CL_a * (param.b / (pow_awt + 1));
            pow_att = i + j;
            A_tt(i,j) = param.c * param.e * param.CL_a * (param.b / (pow_att + 1));
            
            phi_wi = etaB_curr^(i+1); phi_wj = etaB_curr^(j+1);
            phi_ti = etaB_curr^i;     phi_tj = etaB_curr^j;
            
            K_ww(i,j) = K_ww(i,j) + keq_z_curr * phi_wi * phi_wj;
            % FIXED: Use xB_local instead of param.xB to avoid Unrecognized field error
            K_tt(i,j) = K_tt(i,j) + keq_z_curr * (xB_local^2) * phi_ti * phi_tj;
            K_wt(i,j) = -keq_z_curr * xB_local * phi_wi * phi_tj;
        end
    end
    
    K_global = [K_ww, K_wt; K_wt', K_tt];
    A_global = [zeros(N_galerkin, N_galerkin), A_wt; zeros(N_galerkin, N_galerkin), A_tt];
    
    [~, D_eig] = eig(K_global, A_global);
    eigs_val = diag(D_eig);
    physical_roots = eigs_val(imag(eigs_val) == 0 & real(eigs_val) > 0 & ~isinf(eigs_val));
    
    if isempty(physical_roots)
        q_div_history(step) = NaN;
    else
        q_div_history(step) = min(physical_roots);
    end
end

% --- 5. Plotting (Light Mode) ---
figure('Color','w','Name','Point 5: Multidisciplinary Analysis');
grid on; hold on;

yyaxis left
plot(gamma_vec_deg, q_div_history, '-ob', 'LineWidth', 2, 'MarkerFaceColor', 'b');
ylabel('Divergence Pressure q_D [Pa]', 'FontWeight', 'bold');
set(gca, 'YColor', 'b');

yyaxis right
plot(gamma_vec_deg, strut_length, '--or', 'LineWidth', 2, 'MarkerFaceColor', 'r');
ylabel('Strut Length L_s [m]', 'FontWeight', 'bold');
set(gca, 'YColor', 'r');

xlabel('Strut Angle \gamma [deg]', 'FontWeight', 'bold');

% FIXED: Simplified Title syntax to avoid Interpreter errors
title_str = sprintf('5: Divergence and Geometry vs \\gamma (x_B = %.1f, N = %d)', xB_local, N_galerkin);
title(title_str, 'FontSize', 12);

limit_angle = rad2deg(atan(zA_fixed / (param.b/2)));
xline(limit_angle, ':k', 'Span Limit', 'LineWidth', 1.5);
legend({'q_D', 'L_s'}, 'Location', 'north');