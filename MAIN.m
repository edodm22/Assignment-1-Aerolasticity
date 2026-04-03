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

%% punto 1.d
function data = solve_and_plot_pt1d(data)


    v_Lalpha = zeros(2, 1);
    v_Lalpha(2) = data.wing.c * data.aero.CL_alpha * data.wing.b / 2; 
    L_beta_unit = data.wing.c * data.aero.CL_beta * data.aero.ba;
    % 2. Matrice Aerodinamica per l'inversione (Ka_rev) - Dimensione 2x2
    Ka_rev = data.sys_mat.Ka - (data.sys_mat.Fa_unit * v_Lalpha') / L_beta_unit;
    % 3. Calcolo dell'asintoto analitico (xB_crit_R)
    gamma_curr = data.strut.gamma_rad; 
    K_strut_0 = compute_strut_matrix_N1(data, data.strut.yB_base, data.strut.xB_base, gamma_curr);
    Ks_11 = data.sys_mat.Kw(1,1) + K_strut_0(1,1);
    K_strut_1 = compute_strut_matrix_N1(data, data.strut.yB_base, 1, gamma_curr);
    K_slope = K_strut_1(1,2); % Pendenza del termine di accoppiamento rispetto a xB
    % Posizione esatta dell'asintoto verticale per l'inversione (valida solo per N=1)
    xB_crit_R = (Ks_11 * Ka_rev(2,2)) / (K_slope * Ka_rev(1,2));
    data.results.pt1d.xB_crit_R = xB_crit_R;
    % 4. Baseline Ala senza puntone (Unbraced Wing) - Sistema 2x2
    [~, D_unbraced] = eig(data.sys_mat.Kw, Ka_rev);
    q_eig_un = diag(D_unbraced);
    valid_q_un = real(q_eig_un(abs(imag(q_eig_un)) < 1e-5));
    [~, idx] = min(abs(valid_q_un(valid_q_un > 0))); 
    data.results.pt1d.qR_unbraced = valid_q_un(idx);
    % 5. Sweep ad alta risoluzione sulla posizione xB del puntone
    xB_vec = linspace(-data.wing.c/2, data.wing.c/2, 1000);
    qR_vec = zeros(size(xB_vec));
    for i = 1:length(xB_vec)
        xB_curr = xB_vec(i);
        K_strut = compute_strut_matrix_N1(data, data.strut.yB_base, xB_curr, gamma_curr);
        Ks = data.sys_mat.Kw + K_strut; 
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
    data.results.pt1d.xB_vec = xB_vec;
    data.results.pt1d.qR_vec = qR_vec;
 
    % --- Routine di Plot ---
    figure('Name', 'Point 1.d - Control Reversal (N=1)', 'Position', [200, 200, 800, 550], 'Color', [0.1 0.1 0.1]);
    ax = gca; hold on; grid on;
    set(ax, 'Color', [0.1 0.1 0.1], 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'GridColor', [0.4 0.4 0.4]);
    % Baseline (Ala senza puntone)
    plot([-0.5, 0.5], [data.results.pt1d.qR_unbraced, data.results.pt1d.qR_unbraced], ...
         'r--', 'LineWidth', 2.5, 'DisplayName', 'q_R senza puntone');
    % Linea Asintoto Verticale
    xline(xB_crit_R, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    % Gestione salto iperbolico (singolarità)
    [~, split_idx] = min(abs(xB_vec - xB_crit_R));
    qR_vec(split_idx) = NaN; 
    % Plot Ala con puntone
    plot(xB_vec, qR_vec, '-', 'Color', [0 0.447 0.741], 'LineWidth', 3, 'DisplayName', 'q_R con puntone');
    plot([-0.5, 0.5], [0, 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off'); 
    xlabel('Posizione cordale puntone $x_B$ [m]', 'FontSize', 12, 'Interpreter', 'latex', 'Color', 'w');
    ylabel('$q_R$ [Pa]', 'FontSize', 12, 'Interpreter', 'latex', 'Color', 'w');
    title('\textbf{1.d: Inversione del Comando vs Posizione Puntone ($N=1$)}', 'FontSize', 14, 'Interpreter', 'latex', 'Color', 'w');
    xlim([-0.5, 0.5]);
    ylim([-1e6, 1e6]); 
    legend('Location', 'northeast', 'FontSize', 11, 'Interpreter', 'latex');
end