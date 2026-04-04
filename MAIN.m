%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%            STRUCTURAL DYNAMICS AND AEROELASTICITY       %%%%%%%%%
%%%%%%%%%                       A.A. 2025-2026                    %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%                        ASSIGNMENT 1                     %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
clear
clc
close all


% CODE
%% ============================================================
%                    PARAMETERS DEFINITION
%% ============================================================
% All physical and geometrical parameters are grouped here
% to keep the code clean and easily modifiable.

wing.EJ = 1e7;          % Bending stiffness of the wing [Nm^2]
wing.GJ = 1e6;          % Torsional stiffness of the wing [Nm^2]
strut.EA = 1e8;          % Axial stiffness of the strut [N]

wing.CL_a = 2*pi;       % Lift curve slope [1/rad]
wing.CL_beta = pi;      % Lift variation due to control deflection
wing.CM_beta = -0.1;    % Moment coefficient due to control deflection

wing.b = 14;            % Wing span [m]
wing.c = 1;             % Chord length [m]
wing.e = wing.c/4;     % Distance between Aerodynamic Center and Elastic Axis [m]

strut.gamma_deg = 15;      % Reference strut angle [deg]
strut.gamma = deg2rad(15); % Reference strut angle [rad]
wing.ba = 2.5;           % Aileron span [m]
wing.beta0= deg2rad(10);     % Input for Point c


%% ============================================================
%                      POINT 1.a
%             AEROELASTIC DIVERGENCE ANALYSIS
%% ============================================================
% Compute the divergence dynamic pressure q_D as a function
% of the chordwise position of the strut attachment (x_B).

% ------------------------------------------------------------
% Equivalent stiffness of the strut in vertical direction
% ------------------------------------------------------------
% The strut provides an additional elastic contribution.
% Only the vertical component contributes to bending-torsion coupling.

strut.k_eqz = (strut.EA * cos(strut.gamma) / (wing.b/2)) * (sin(strut.gamma)^2);

% ------------------------------------------------------------
% Sweep of strut position along the chord
% ------------------------------------------------------------
% xB = -0.5 → Leading Edge
% xB =  0   → Elastic Axis
% xB =  0.5 → Trailing Edge

A.xB = linspace(-wing.c/2, wing.c/2, 1800);

% ------------------------------------------------------------
% Structural stiffness terms (independent of q)
% ------------------------------------------------------------
A.K11 = 4*wing.EJ/wing.b^3 + strut.k_eqz/16;   % bending stiffness term
A.K22_base = wing.GJ/wing.b;             % torsional stiffness (base)
A.K12_coeff = -strut.k_eqz/8;                    % coupling coefficient

% ------------------------------------------------------------
% Aerodynamic influence coefficients
% ------------------------------------------------------------
% These terms come from the aerodynamic matrix
% and represent how lift and moment depend on deformation.

A.A22 = wing.e * wing.b * wing.c * wing.CL_a / 3;
A.A12 = wing.b * wing.c * wing.CL_a / 4;

% ------------------------------------------------------------
% Divergence condition
% ------------------------------------------------------------
% Divergence occurs when det(K - q*A) = 0
% Here we compute q explicitly as a ratio:
%   q_D = (structural terms) / (aerodynamic terms)

A.Num = A.K11*(A.K22_base + strut.k_eqz*A.xB.^2/4) - (A.K12_coeff*A.xB).^2;
A.D   = A.K11*(A.A22) - (A.K12_coeff*A.xB).*(A.A12);

A.q_div = A.Num ./ A.D;

% ------------------------------------------------------------
% Plot results
% ------------------------------------------------------------
figure('Color','w','Name','Divergence Pressure');
plot(A.xB, A.q_div, 'b','LineWidth',2);
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

A.x_critic = (A.K11 * A.A22) / (A.K12_coeff * A.A12);

%% ============================================================
%                      POINT 1.b
%        CONTROL EFFECTIVENESS (FIXED POINT A)
%% ============================================================
% Study how control effectiveness E_C varies with dynamic pressure
% and with the strut angle gamma.

% ------------------------------------------------------------
% Fixed point A condition
% ------------------------------------------------------------
% The geometry is constrained so that point A is fixed.
% This introduces a relation between zA and gamma.

B.gamma_ref = strut.gamma;
B.yB_ref = wing.b/2;
wing.zA = B.yB_ref * tan(B.gamma_ref);

% ------------------------------------------------------------
% Choose a specific chordwise position for the strut
% ------------------------------------------------------------
wing.xB = -0.1; 

% ------------------------------------------------------------
% Define parameter ranges
% ------------------------------------------------------------
B.gamma_deg = [10, 15, 20];       % different strut inclinations
B.gamma_deg = linspace(10,20,5); 
B.qd_vec = linspace(0, 3e5, 1000); % dynamic pressure range

figure('Color','w'); hold on; grid on;
colors = {'r','b','g', 'c', 'm'};
% ------------------------------------------------------------
% Print Header for Control Reversal Pressures
% ------------------------------------------------------------
fprintf('\n============================================================\n');
fprintf('   Dynamic Pressures for Control Reversal (E_C = 0)\n');
fprintf('============================================================\n');

for j = 1:length(B.gamma_deg)

    B.g = deg2rad(B.gamma_deg(j));

    % --------------------------------------------------------
    % Geometry update due to gamma variation
    % --------------------------------------------------------
    % Since A is fixed, changing gamma changes yB

    B.yB = wing.zA / tan(B.g);
    B.etaB = B.yB / wing.b;   % nondimensional position

    % --------------------------------------------------------
    % Equivalent stiffness (nonlinear with gamma)
    % --------------------------------------------------------
    B.keq_z = (strut.EA * sin(B.g)^3) / wing.zA;

    % --------------------------------------------------------
    % Structural stiffness matrix K
    % --------------------------------------------------------
    % Includes coupling between bending and torsion

    B.K11 = 4*wing.EJ/wing.b^3 + B.keq_z * B.etaB^4;
    B.K22 = wing.GJ/wing.b + B.keq_z * wing.xB^2 * B.etaB^2;
    B.K12 = -B.keq_z * wing.xB * B.etaB^3;

    B.K = [B.K11, B.K12;
         B.K12, B.K22];

    % --------------------------------------------------------
    % Aerodynamic matrix
    % --------------------------------------------------------
    B.A = [0, wing.b*wing.c*wing.CL_a/4;
         0, wing.e*wing.b*wing.c*wing.CL_a/3];

    % --------------------------------------------------------
    % Aileron forcing vector
    % --------------------------------------------------------
    % Represents aerodynamic loads due to control deflection

    B.F1b = (wing.c * wing.CL_beta / (3*wing.b^2)) * ...
          (wing.b^3 - (wing.b-wing.ba)^3);

    B.F2b = (wing.c^2 * wing.CM_beta / (2*wing.b) + ...
           wing.e*wing.c*wing.CL_beta / (2*wing.b)) * ...
          (wing.b^2 - (wing.b-wing.ba)^2);

    B.F_beta = [B.F1b; B.F2b];

    % --------------------------------------------------------
    % Solve system for each dynamic pressure
    % --------------------------------------------------------
    % Equation:
    %   (K - q*A) * q_vec = q * F_beta

    B.Ec = zeros(size(B.qd_vec));

    for i = 1:length(B.qd_vec)

        B.q = B.qd_vec(i);

        % Check stability using eigenvalues of the aeroelastic system
        B.eig_vals = eig(B.K - B.q*B.A);

        B.sol = (B.K - B.q*B.A) \ (B.q * B.F_beta);
            B.theta = B.sol(2); % torsional response

            
            % Control Effectiveness expression
            B.Ec(i) = 1 + (wing.b * wing.CL_a * B.theta) / ...
            (2 * wing.CL_beta * wing.ba);         
    end
    % --------------------------------------------------------
    % Find the control reversal point (E_C = 0)
    % --------------------------------------------------------
    % Find the index where E_C crosses from positive to negative
    B.idx_zero = find(B.Ec(1:end-1) > 0 & B.Ec(2:end) <= 0, 1);
    
    if ~isempty(B.idx_zero)
        % Linear interpolation to find the exact q value between the two points
        B.q1 = B.qd_vec(B.idx_zero);
        B.q2 = B.qd_vec(B.idx_zero + 1);
        B.Ec1 = B.Ec(B.idx_zero);
        B.Ec2 = B.Ec(B.idx_zero + 1);
        
        B.q_R_exact = B.q1 - B.Ec1 * (B.q2 - B.q1) / (B.Ec2 - B.Ec1);
        fprintf('Strut Angle \\gamma = %4.1f [deg]  -->  q_R = %8.2f [Pa]\n', B.gamma_deg(j), B.q_R_exact);
    else
        fprintf('Strut Angle \\gamma = %4.1f [deg]  -->  q_R = Out of bounds\n', B.gamma_deg(j));
    end
    % Plot curve
    plot(B.qd_vec, B.Ec, 'Color', colors{j}, 'LineWidth',2, ...
        'DisplayName',['\gamma = ', num2str(B.gamma_deg(j)), '°']);
end

% ------------------------------------------------------------
% Final plot formatting
% ------------------------------------------------------------
yline(0,'k','HandleVisibility','off');
yline(1,'k--','Rigid Wing','HandleVisibility','off');

xlabel('Dynamic Pressure q [Pa]');
ylabel('Control Effectiveness E_C');
title(['Control Effectiveness (x_B = ', num2str(wing.xB), ')']);

legend('Location','best');


%% ============================================================
%                      POINT 1.c
%             DEFORMED SHAPE COMPARISON
%% ============================================================
% 1. Reference: No Strut Divergence
wing.qD_no_strut = (3 * wing.GJ) / (wing.e * wing.b^2 * wing.c * wing.CL_a); %formula ricavata a mano 
C.q_test = 0.4 * wing.qD_no_strut; 

C.F_beta=B.F_beta;
C.A = B.A;

% 2. System WITHOUT Strut
C.K_noStrut = [4*wing.EJ/wing.b^3, 0; 0, wing.GJ/wing.b];
C.sol_noStrut = (C.K_noStrut - C.q_test*C.A) \ (C.q_test * C.F_beta); % A rimane invariata

% 3. System WITH Strut (using gamma = 15 deg)
C.g15 = deg2rad(15); % per cambiare gamma cambia qua! 
C.yB15 = wing.zA / tan(C.g15); 
C.etaB15 = C.yB15 / wing.b;
C.keq15 = (strut.EA * sin(C.g15)^3) / wing.zA;
C.K_strut = [4*wing.EJ/wing.b^3 + C.keq15*C.etaB15^4, -C.keq15*wing.xB*C.etaB15^3;
           -C.keq15*wing.xB*C.etaB15^3, wing.GJ/wing.b + C.keq15*wing.xB^2*C.etaB15^2];
C.sol_strut = (C.K_strut - C.q_test*C.A) \ (C.q_test * C.F_beta);

% 4. Plotting Deformed Shapes
C.y_span = linspace(0, wing.b, 100);
C.w_noStrut = (C.y_span/wing.b).^2 * C.sol_noStrut(1) * wing.beta0;
C.t_noStrut = (C.y_span/wing.b) * C.sol_noStrut(2) * wing.beta0;
C.w_strut = (C.y_span/wing.b).^2 * C.sol_strut(1) * wing.beta0;
C.t_strut = (C.y_span/wing.b) * C.sol_strut(2) * wing.beta0;

figure('Color','w','Name','Deformed Shape Comparison');
subplot(2,1,1);
plot(C.y_span, C.w_noStrut, 'r--', 'LineWidth', 2); hold on;
plot(C.y_span, C.w_strut, 'b', 'LineWidth', 2);
ylabel('Bending w(y) [m]'); grid on; legend('No Strut','With Strut');
title(['Wing Deformation at q = ', num2str(round(C.q_test)), ' Pa']);

subplot(2,1,2);
plot(C.y_span, rad2deg(C.t_noStrut), 'r--', 'LineWidth', 2); hold on;
plot(C.y_span, rad2deg(C.t_strut), 'b', 'LineWidth', 2);
ylabel('Twist \theta(y) [deg]'); xlabel('y [m]'); grid on;


%% ============================================================
%                      POINT 1.d
%       CONTROL REVERSAL vs STRUT POSITION (N=1)
%% ============================================================ 
% Find the dynamic pressure q_R at which control effectiveness 
% becomes zero, as a function of the strut position xB.

% 1. Setup Base Matrices (N=1)
D.Kw_clean = [4*wing.EJ/wing.b^3, 0; 0, wing.GJ/wing.b];
D.Ka = [0, wing.b*wing.c*wing.CL_a/4; 
      0, wing.e*wing.b*wing.c*wing.CL_a/3];

% Aileron forcing vector (Normalized per unit beta)
D.F1b_unit = (wing.c * wing.CL_beta / (3*wing.b^2)) * ...
           (wing.b^3 - (wing.b-wing.ba)^3);
D.F2b_unit = (wing.c^2 * wing.CM_beta / (2*wing.b) + ...
           wing.e*wing.c*wing.CL_beta / (2*wing.b)) * ...
           (wing.b^2 - (wing.b-wing.ba)^2);
D.Fa_unit = [D.F1b_unit; D.F2b_unit];

% 2. Reversal Aerodynamic Matrix (Ka_rev)
D.v_Lalpha = [0; wing.c * wing.CL_a * wing.b / 2];
D.L_beta_unit = wing.c * wing.CL_beta * wing.ba;
D.Ka_rev = D.Ka - (D.Fa_unit * D.v_Lalpha') / D.L_beta_unit;

% 3. Geometry (Fixed Point A logic)
D.etaB = wing.b/2 / wing.b; % yB/total span
D.keq_z = (strut.EA * sin(strut.gamma)^3) / wing.zA;

% 4. Analytic Asymptote (N=1 logic)
D.Ks_11 = D.Kw_clean(1,1) + D.keq_z * D.etaB^4;
D.K_slope = -D.keq_z * D.etaB^3; 
D.xB_crit_R = (D.Ks_11 * D.Ka_rev(2,2)) / (D.K_slope * D.Ka_rev(1,2));

% 5. Baseline: Unbraced Wing Reversal
[~, D.D_unbraced] = eig(D.Kw_clean, D.Ka_rev);
D.q_eig_un = diag(D.D_unbraced);
D.valid_q_un = real(D.q_eig_un(abs(imag(D.q_eig_un)) < 1e-5 & D.q_eig_un > 0));
D.qR_unbraced = min(D.valid_q_un);

% 6. Sweep over Strut Position xB
D.xB_vec_R = linspace(-wing.c/2, wing.c/2, 1000);
D.qR_vec = zeros(size(D.xB_vec_R));
for i = 1:length(D.xB_vec_R)
    D.xB_curr = D.xB_vec_R(i);
    D.K_strut = [D.keq_z * D.etaB^4,            -D.keq_z * D.xB_curr * D.etaB^3;
              -D.keq_z * D.xB_curr * D.etaB^3,   D.keq_z * D.xB_curr^2 * D.etaB^2];
    D.Ks = D.Kw_clean + D.K_strut; 
    [~, D.D_braced] = eig(D.Ks, D.Ka_rev);
    D.q_eig_br = diag(D.D_braced);
    D.valid_q_br = real(D.q_eig_br(abs(imag(D.q_eig_br)) < 1e-5));
    if isempty(D.valid_q_br)
        D.qR_vec(i) = NaN; 
    else
        [~, D.idx_br] = min(abs(D.valid_q_br));
        D.qR_vec(i) = D.valid_q_br(D.idx_br);
    end
end

figure('Name', 'Point 4 - Control Reversal (N=1)'); 
hold on; 
grid on;

% Baseline (Unbraced) - Red dashed line
plot([-0.5, 0.5], [D.qR_unbraced, D.qR_unbraced], 'r--', 'LineWidth', 1.5, 'DisplayName', 'q_R unbraced');

% Vertical asymptote
xline(D.xB_crit_R, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');

% Jump handling for plot (prevents diagonal line crossing the graph)
[~, D.split_idx] = min(abs(D.xB_vec_R - D.xB_crit_R));
D.temp_qR_vec = D.qR_vec; % Using a temporary variable to not pollute the struct
D.temp_qR_vec(D.split_idx) = NaN; 

% Braced wing - Standard blue
plot(D.xB_vec_R, D.temp_qR_vec, 'b-', 'LineWidth', 2, 'DisplayName', 'q_R braced');

% Zero line
yline(0, 'k', 'HandleVisibility', 'off'); 

% Titles and labels (Standard text)
xlabel('Strut chordwise position x_B [m]');
ylabel('q_R [Pa]');
title('4: Control Reversal vs Strut Position (N=1)');

% Axis limits
xlim([-0.5, 0.5]); 
% Note: removed extreme ylim (-1e6, 1e6) to let MATLAB choose a readable scale.
% If needed, uncomment the line below:
% ylim([-5e5, 5e5]); 
legend('Location', 'best');
%% ============================================================
%                      POINT 1.e
%           RITZ-GALERKIN CONVERGENCE STUDY
%% ============================================================ 
% Convergence study settings
E.N_max = 8; % Maximum number of shape functions to test (for w and theta)
E.qD_history = zeros(1, E.N_max); % Vector to store divergence pressures
E.N_values = 1:E.N_max;

% Iterative loop over shape functions
for N = E.N_values
    
    % Initialization of submatrices (N x N)
    E.K_ww = zeros(N, N); % Bending stiffness
    E.K_tt = zeros(N, N); % Torsional stiffness
    E.K_wt = zeros(N, N); % Strut coupling
    
    E.A_wt = zeros(N, N); % Lift work on bending
    E.A_tt = zeros(N, N); % Moment work on torsion
    
    % Fill submatrices
    for i = 1:N
        for j = 1:N
            % Note on integrals: we use non-dimensional coordinate eta = y/b
            % Integral from 0 to b of (y/b)^A dy equals b / (A + 1)
            
            % -- STRUCTURAL MATRICES (Clean Wing) --
            % K_ww: Integral( EJ * phi_wi'' * phi_wj'' dy )
            % phi_wi'' = (i+1)*i/b^2 * eta^(i-1)
            E.pot_ww = (i-1) + (j-1); 
            E.K_ww(i,j) = wing.EJ * ((i+1)*i/wing.b^2) * ((j+1)*j/wing.b^2) * (wing.b / (E.pot_ww + 1));
            
            % K_tt: Integral( GJ * phi_ti' * phi_tj' dy )
            % phi_ti' = i/b * eta^(i-1)
            E.pot_tt_strut = (i-1) + (j-1);
            E.K_tt(i,j) = wing.GJ * (i/wing.b) * (j/wing.b) * (wing.b / (E.pot_tt_strut + 1));
            
            % -- AERODYNAMIC MATRICES --
            % A_wt: Lift work on w -> Integral( c * CLalpha * phi_wj * phi_ti dy )
            E.pot_wt = (i+1) + j; 
            E.A_wt(i,j) = wing.c * wing.CL_a * (wing.b / (E.pot_wt + 1));
            
            % A_tt: Moment work on theta -> Integral( c * e * CLalpha * phi_ti * phi_tj dy )
            E.pot_tt_aero = i + j;
            E.A_tt(i,j) = wing.c * wing.e * wing.CL_a * (wing.b / (E.pot_tt_aero + 1));
            
            % -- STRUT CONTRIBUTION (Evaluated at y = b/2 -> eta = 0.5) --
            E.phi_wi_b2 = (0.5)^(i+1);
            E.phi_wj_b2 = (0.5)^(j+1);
            E.phi_ti_b2 = (0.5)^i;
            E.phi_tj_b2 = (0.5)^j;
            
            % Addition of discrete strut terms to matrices
            E.K_ww(i,j) = E.K_ww(i,j) + B.keq_z * E.phi_wi_b2 * E.phi_wj_b2;
            E.K_tt(i,j) = E.K_tt(i,j) + B.keq_z * (wing.xB^2) * E.phi_ti_b2 * E.phi_tj_b2;
      
            % Flexural-torsional coupling term (K_wt and K_tw)
            E.K_wt(i,j) = -B.keq_z * wing.xB * E.phi_wi_b2 * E.phi_tj_b2;
        end
    end
    
    % 4. Global System Assembly (2N x 2N Matrices)
    % Order of degrees of freedom: [qw1..qwN, qtheta1..qthetaN]'
    E.K_global = [E.K_ww,       E.K_wt; 
                E.K_wt',      E.K_tt];  % Transpose ensures symmetry
            
    E.A_global = [zeros(N,N), E.A_wt; 
                zeros(N,N), E.A_tt];
            
    % 5. Resolution of Eigenvalue Problem (Divergence)
    % (K_global - qD * A_global) * {q} = 0
    [E.V, E.D] = eig(E.K_global, E.A_global);
    E.eig = diag(E.D);
    
    % Filter only physical eigenvalues (real and strictly positive)
    E.qd = E.eig(imag(E.eig) == 0 & real(E.eig) > 0);
    
    if isempty(E.eig)
        E.qD_history(N) = NaN; % No instability found
    else
        E.qD_history(N) = min(E.qd); % q_D is the lowest critical value
    end
end

% 6. Plot Results
figure('Color', 'w');
plot(E.N_values, E.qD_history, '-o', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
set(gca, 'FontSize', 12);
xlabel('Number of shape functions (N) per DoF', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Divergence Pressure q_D [Pa]', 'FontSize', 14, 'FontWeight', 'bold');
title('Convergence Study (Ritz-Galerkin Method)', 'FontSize', 16);

% Print values to console
disp('--- Convergence Study Results ---');
for i = 1:E.N_max
    fprintf('N = %d -> q_D = %.2f Pa\n', i, E.qD_history(i));
end


%% ============================================================
%                      POINT 1.f
%    MULTIDISCIPLINARY ANALYSIS: DIVERGENCE vs GEOMETRY
%% ============================================================
% Analyze the trade-off between q_D and strut length L_s.

% --- Parameters for this Study ---
F.N = 5;                  % Ritz-Galerkin order
F.toll = 1e-6;
% --- Iteration Vectors ---
F.gamma_vec_deg = strut.gamma_deg:1:40;         
F.q_div_history = zeros(size(F.gamma_vec_deg));
F.strut_length = zeros(size(F.gamma_vec_deg));

% --- Main Loop over Gamma ---
for step = 1:length(F.gamma_vec_deg)
    F.g_rad = deg2rad(F.gamma_vec_deg(step));
    
    % Geometry
    F.yB_curr = wing.zA / tan(F.g_rad);
    F.strut_length(step) = wing.zA / sin(F.g_rad);
    
    if F.yB_curr > (wing.b/2 + F.toll) || F.yB_curr <= 0
        F.q_div_history(step) = NaN;
        F.strut_length(step) = NaN;
        continue;
    end
    
    F.etaB_curr = F.yB_curr / wing.b;
    F.keq_z_curr = (strut.EA * sin(F.g_rad)^3) / wing.zA;
    
    % Galerkin matrices
    F.K_ww = zeros(F.N, F.N); 
    F.K_tt = zeros(F.N, F.N); 
    F.K_wt = zeros(F.N, F.N);
    F.A_wt = zeros(F.N, F.N); 
    F.A_tt = zeros(F.N, F.N);
    
    for i = 1:F.N
        for j = 1:F.N
            F.pow_ww = (i-1) + (j-1);
            F.K_ww(i,j) = wing.EJ * ((i+1)*i/wing.b^2) * ((j+1)*j/wing.b^2) * (wing.b / (F.pow_ww + 1));
            F.pow_tt = (i-1) + (j-1);
            F.K_tt(i,j) = wing.GJ * (i/wing.b) * (j/wing.b) * (wing.b / (F.pow_tt + 1));
            F.pow_awt = (i+1) + j;
            F.A_wt(i,j) = wing.c * wing.CL_a * (wing.b / (F.pow_awt + 1));
            F.pow_att = i + j;
            F.A_tt(i,j) = wing.c * wing.e * wing.CL_a * (wing.b / (F.pow_att + 1));
            
            F.phi_wi = F.etaB_curr^(i+1); 
            F.phi_wj = F.etaB_curr^(j+1);
            F.phi_ti = F.etaB_curr^i;     
            F.phi_tj = F.etaB_curr^j;
            
            F.K_ww(i,j) = F.K_ww(i,j) + F.keq_z_curr * F.phi_wi * F.phi_wj;
            F.K_tt(i,j) = F.K_tt(i,j) + F.keq_z_curr * (wing.xB^2) * F.phi_ti * F.phi_tj;
            F.K_wt(i,j) = -F.keq_z_curr * wing.xB * F.phi_wi * F.phi_tj;
        end
    end
    
    F.K_global = [F.K_ww, F.K_wt; F.K_wt', F.K_tt];
    F.A_global = [zeros(F.N, F.N), F.A_wt; zeros(F.N, F.N), F.A_tt];
    
    [~, F.D_eig] = eig(F.K_global, F.A_global);
    F.eigs_val = diag(F.D_eig);
    %F.physical_roots = F.eigs_val(imag(F.eigs_val) == 0 & real(F.eigs_val) > 0 & ~isinf(F.eigs_val));
    F.physical_roots = F.eigs_val(abs(imag(F.eigs_val)) < F.toll & real(F.eigs_val) > 0 & ~isinf(F.eigs_val));
    if isempty(F.physical_roots)
        F.q_div_history(step) = NaN;
    else
        F.q_div_history(step) = min(F.physical_roots);
    end
end

% --- Plotting (Light Mode) ---
figure('Color','w','Name','Point 5: Multidisciplinary Analysis');
grid on; hold on;

yyaxis left
plot(F.gamma_vec_deg, F.q_div_history, '-ob', 'LineWidth', 2, 'MarkerFaceColor', 'b');
ylabel('Divergence Pressure q_D [Pa]', 'FontWeight', 'bold');
set(gca, 'YColor', 'b');

yyaxis right
plot(F.gamma_vec_deg, F.strut_length, '--or', 'LineWidth', 2, 'MarkerFaceColor', 'r');
ylabel('Strut Length L_s [m]', 'FontWeight', 'bold');
set(gca, 'YColor', 'r');

xlabel('Strut Angle \gamma [deg]', 'FontWeight', 'bold');

F.title_str = sprintf('5: Divergence and Geometry vs \\gamma (x_B = %.1f, N = %d)', wing.xB, F.N);
title(F.title_str, 'FontSize', 12);

F.limit_angle = rad2deg(atan(wing.zA / (wing.b/2)));
xline(F.limit_angle, ':k', 'Span Limit', 'LineWidth', 1.5);
legend({'q_D', 'L_s'}, 'Location', 'north');

%% ============================================================
% --- Point 1.f: Design Optimization Assessment (+40% q_D) ---
% ============================================================

% Material Properties Assumption (Aerospace Aluminum 2024 or 7075)
F.E_alu = 70e9;      % Young's Modulus [Pa]
F.rho_alu = 2780;    % Density [kg/m^3]

% Extract Baseline conditions robustly (Find the first non-NaN value)
F.valid_indices = find(~isnan(F.q_div_history));
F.idx_base = F.valid_indices(1); % Prende il primo indice valido (es. 15° o 16°)

F.qD_baseline = real(F.q_div_history(F.idx_base));
F.Ls_baseline = F.strut_length(F.idx_base);
F.gamma_baseline = F.gamma_vec_deg(F.idx_base);

% Define target performance
F.target_multiplier = 1.40; 
F.qD_target = F.qD_baseline * F.target_multiplier;

% Find the first integer angle that satisfies the condition
F.idx_optimal = find(real(F.q_div_history) >= F.qD_target, 1);

%  --- DEBUG ---
% fprintf('--- DEBUG INFO ---\n');
% fprintf('qD Baseline (15 gradi): %.2f Pa\n', F.qD_baseline);
% fprintf('Target da raggiungere (+40%%): %.2f Pa\n', F.qD_target);
% fprintf('Massimo qD raggiunto nel vettore: %.2f Pa\n', max(real(F.q_div_history)));
% fprintf('------------------\n');

if ~isempty(F.idx_optimal)
    
    % Retrieve optimal geometry
    F.gamma_new = F.gamma_vec_deg(F.idx_optimal);
    F.qD_new = real(F.q_div_history(F.idx_optimal));
    F.Ls_new = F.strut_length(F.idx_optimal);
    
    % Structural sizing (Cross-sectional Area A = EA / E)
    F.A_strut = strut.EA / F.E_alu;
    
    % Mass calculations (m = rho * A * L)
    F.mass_baseline = F.rho_alu * F.A_strut * F.Ls_baseline;
    F.mass_new      = F.rho_alu * F.A_strut * F.Ls_new;
    
    % Variations
    F.delta_mass_kg = F.mass_new - F.mass_baseline;
    F.delta_mass_perc = (F.delta_mass_kg / F.mass_baseline) * 100;
    F.delta_Ls = F.Ls_new - F.Ls_baseline;
    
    % Print comprehensive results to Command Window
    fprintf('\n============================================================\n');
    fprintf('   Point 1.f: Structural Redesign for Minimum 40%% q_D Increase\n');
    fprintf('============================================================\n');
    fprintf('Baseline Configuration (\\gamma = %d [deg]):\n', F.gamma_baseline);
    fprintf('  - Divergence q_D : %10.2f [Pa]\n', F.qD_baseline);
    fprintf('  - Strut Length   : %10.3f [m]\n', F.Ls_baseline);
    fprintf('  - Strut Mass     : %10.2f [kg]\n\n', F.mass_baseline);
    
    fprintf('Optimal Target Configuration (Target q_D >= %.2f [Pa]):\n', F.qD_target);
    fprintf('  - Optimal Angle  : %d [deg]\n', F.gamma_new);
    fprintf('  - Achieved q_D   : %10.2f [Pa]\n', F.qD_new);
    fprintf('  - New Length     : %10.3f [m] (Variation: %+.3f [m])\n', F.Ls_new, F.delta_Ls);
    fprintf('  - New Mass       : %10.2f [kg]\n', F.mass_new);
    fprintf('  - Mass Variation : %+10.2f [kg] (%+.2f %%)\n', F.delta_mass_kg, F.delta_mass_perc);
    fprintf('============================================================\n');
    
    % Add visual marker to the existing plot
    yyaxis left;
    plot(F.gamma_new, F.qD_new, 'p', 'MarkerSize', 14, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
    text(F.gamma_new + 0.5, F.qD_new, ' Optimal Design', 'FontSize', 10, 'FontWeight', 'bold');
    
else
    fprintf('\n[!] No integer angle found that satisfies the +40%% requirement within the checked range.\n');
    fprintf('    (Make sure the numerical tolerance is enabled in the physical_roots check!)\n');
end