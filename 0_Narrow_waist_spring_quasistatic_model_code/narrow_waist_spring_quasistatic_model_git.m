%% ========================================================================
%  Displacement-driven quasi-static model of a narrow-waist spring
%
%  Key model equations:
%  ------------------------------------------------------------------------
%  Narrow-waist geometry:r(u) = r_w + (r_e - r_w) · (2u - 1)²
%
%  Equivalent bending stiffnes:
%  EI_i = (E π d_w⁴ / 64) · 2 sin(β_i) / [2 + ν cos²(β_i)]
%
%  Local axial stiffness：
%  EA_i = G · d_w⁴ · p_i / [8 · (2 r_i)³]
%  ------------------------------------------------------------------------
%  Script outputs:
%  1. A comparison table of the tip bending angle θtip predicted by the
%     model and the Abaqus simulation reference.
%  2. Curve of relative error in tip bending angle vs. number of segments N
%  3. 3D bending shape visualization
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  Section 1: Spring parameter definition
% =========================================================================
 
% ------------------------ 1.1 Geometric parameters ----------------------------------
geom.r_e    = 50e-3;         % Maximum end radius r_e of the spring (m)
geom.r_w    = 30e-3;         % Minimum waist radius r_w of the spring (m)
geom.delta  = 0e-3;          % Offset from the cable guide hole to the
                             % neutral axis of the spring wire (m).
geom.n_turn = 8;             % Total number of helical turns
geom.d_w    = 2e-3;          % Wire diameter (m)

% Segmented pitch (per turn, starting from one end, unit: m)
geom.p_per_turn = [30, 30, 30, 30, 30, 30, 30, 30] * 1e-3;
geom.L0    = sum(geom.p_per_turn);    % Natural axial length of the spring
geom.p_avg = mean(geom.p_per_turn);   % Average pitch

% ------------------------ 1.2 Material parameters ----------------------------------
mat.E  = 201.1e9;                  % Young's modulus E (Pa)
mat.nu = 0.30;                     % Poisson's ratio ν
mat.G  = mat.E/(2*(1 + mat.nu));   % Shear modulus G (Pa)

% ------------------------ 1.3 Circumferential cable angles ------------------------------
alpha = [0; 2*pi/3; 4*pi/3];       %  Three cables uniformly spaced by 120 deg

% ------------------------ 1.4 Model/simulation cases -----------------------------
% theta_E is the Abaqus simulation reference value.
% Negative sign in dL_input = motor tightens (cable shortens)
% Positive sign = cable relaxes (cable lengthens)
cases(1).label    = 'Case 1: [-20,  0, 0] mm';
cases(1).dL_input = [-20;  0; 0] * 1e-3;
cases(1).theta_E  = 0.3391*57.3;    

cases(2).label    = 'Case 2: [-40,  0, 0] mm';
cases(2).dL_input = [-40;  0; 0] * 1e-3;
cases(2).theta_E  = 0.677*57.3;    


cases(3).label    = 'Case 3: [-60,  0, 0] mm';
cases(3).dL_input = [-60;  0; 0] * 1e-3;
cases(3).theta_E  = 1.013*57.3;    


cases(4).label    = 'Case 4: [-80,  0, 0] mm';
cases(4).dL_input = [-80;  0; 0] * 1e-3;
cases(4).theta_E  = 1.347*57.3;   


cases(5).label    = 'Case 5: [-100,  0, 0] mm';
cases(5).dL_input = [-100;  0; 0] * 1e-3;
cases(5).theta_E  = 1.676*57.3;   

cases(6).label    = 'Case 6: [-130, -10, -10] mm';
cases(6).dL_input = [-130; -10; -10] * 1e-3;
cases(6).theta_E  = 1.993*57.3;   
% ------------------------ 1.5 Solver options --------------------------------
solver_opts = optimoptions('fsolve', ...
    'Display',                'off', ...
    'Algorithm',              'levenberg-marquardt', ...
    'FunctionTolerance',      1e-13, ...
    'StepTolerance',          1e-13, ...
    'MaxIterations',          3000, ...
    'MaxFunctionEvaluations', 100000);

% ------------------------ 1.6 Number of Spring Segments --------------------
N_list = [12, 18, 24, 30, 60, 90, 120, 180, 240];

% ------------------------ 1.7 Print current configuration --------------------------------
fprintf('\n=================================================================\n');
fprintf(' Narrow-waist spring parameters\n');
fprintf('=================================================================\n');
fprintf('Geometric parameters：\n');
fprintf('  r_e = %5.2f mm,  r_w = %5.2f mm,  δ = %+5.2f mm\n', ...
    geom.r_e*1000, geom.r_w*1000, geom.delta*1000);
fprintf('  n_turn = %d,  d_w = %.2f mm\n', geom.n_turn, geom.d_w*1000);
fprintf('  pitch (mm) = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]\n', ...
    geom.p_per_turn*1000);
fprintf('  L_0 = %.2f mm  (= sum of per-turn pitches)\n', geom.L0*1000);
fprintf('  p_avg = %.2f mm\n', geom.p_avg*1000);
fprintf('Material parameters：\n');
fprintf('  E = %.0f GPa,  ν = %.2f,  G = %.0f GPa\n\n', ...
    mat.E/1e9, mat.nu, mat.G/1e9);

%% ========================================================================
%  Section 2: Main solution loop
% =========================================================================
nC = numel(cases);
nN = numel(N_list);

results = repmat(struct(...
    'theta_S',      NaN, ...   % Model-predicted tip bending angle (deg)
    'rel_err_th',   NaN, ...   % Relative error in θ vs. Abaqus (%) 相对 Abaqus 误差（%）
    'kappa',        [], ...    % Segment curvature κ_i（1/m）
    'phi',          NaN, ...   % Common bending-plane angle φ（rad）
    'eps',          [], ...    % Segment axial strain
    'lambda',       [], ...    % Lagrange multipliers of the three cable length constraints
    'centerline',   [], ...    % Discrete centerline points
    'axial_change', NaN, ...   % Overall axial length change of the sprinG（m）
    'res_norm',     NaN, ...   % Norm of the KKT residual
    'exitflag',     NaN), nC, nN);

fprintf('=== Numerical Solving Begins ===\n');
for iC = 1:nC
    for iN = 1:nN
        N = N_list(iN);
        out = solveSpringConfig(cases(iC).dL_input, alpha, geom, mat, N, solver_opts);
        F1 = abs(out.lambda(1));

        results(iC,iN).theta_S      = out.theta_tip_deg;
        results(iC,iN).F1           = F1;
        results(iC,iN).rel_err_th   = abs(out.theta_tip_deg - cases(iC).theta_E) ...
                                      / cases(iC).theta_E * 100;

        results(iC,iN).kappa        = out.kappa;
        results(iC,iN).phi          = out.phi;
        results(iC,iN).eps          = out.eps;
        results(iC,iN).lambda       = out.lambda;
        results(iC,iN).centerline   = out.centerline;
        results(iC,iN).axial_change = out.axial_change;
        results(iC,iN).res_norm     = out.res_norm;
        results(iC,iN).exitflag     = out.exitflag;

        fprintf('  %-26s | N=%3d | θ=%7.2f° | err_θ=%5.2f%% |res=%.1e | flag=%d |\n' , ...
            cases(iC).label, N, results(iC,iN).theta_S, results(iC,iN).rel_err_th, ...
            out.res_norm, out.exitflag);
    end
    fprintf('  ---\n');
end

%% ========================================================================
%  Section 3:  Print tip-angle comparison table (Model vs. Abaqus - Max N)
% =========================================================================
N_final = N_list(end);
fprintf('\n=== Comparison with Abaqus (N = %d) ===\n', N_final);
fprintf('  Case                       | θ_calc   θ_Abq    err_θ    |\n');
fprintf('  ---------------------------|----------------------------|\n');
for iC = 1:nC
    sgn_th = sign(results(iC,end).theta_S - cases(iC).theta_E);

    fprintf('  %-26s | %6.2f°  %6.2f°  %+5.1f%%   |\n', ...
        cases(iC).label, ...
        results(iC,end).theta_S, cases(iC).theta_E, sgn_th * results(iC,end).rel_err_th);
end


%% ========================================================================
%  Section 4: Plot results
%  Fig.1: Variation of tip bending angle with number of segments N
% =========================================================================
figure('Name', 'Mesh convergence', 'Color', 'w', 'Position', [80, 80, 500, 500]);
markers = {'o-', 's-', '^-', 'd-', 'v-','^-',':'};
colors  = lines(nC);

% 图 1： Relative error in θtip versus N
hold on; grid on; box on;
for iC = 1:nC
    err = arrayfun(@(k) results(iC,k).rel_err_th, 1:nN);
    plot(N_list, err, markers{iC}, 'LineWidth', 1.7, 'MarkerSize', 8, ...
         'Color', colors(iC,:), 'DisplayName', cases(iC).label);
end
xlabel('Number of segments N', 'FontSize', 12);
ylabel('Relative error in tip bending angle θtip (%)', 'FontSize', 12);
title('Mesh convergence: RE_θtip vs. N', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 9); set(gca, 'FontSize', 11);

%% ========================================================================
%   Fig.2: 3D Bending Shape Visualization
% =========================================================================
[~, idx_plot] = min(abs(N_list - 120));
N_plot = N_list(idx_plot);

figure('Name', 'Spring bending 3D', 'Color', 'w', 'Position', [120, 120, 1700, 480]);
for iC = 1:nC
    subplot(1, nC, iC);
    cl_mm = results(iC, idx_plot).centerline * 1000;
    plot3(cl_mm(1,:), cl_mm(2,:), cl_mm(3,:), '-', ...
          'LineWidth', 2.0, 'Color', colors(iC,:));
    hold on; grid on; box on;
    plot3(cl_mm(1,1),   cl_mm(2,1),   cl_mm(3,1),   'ko', ...
          'MarkerFaceColor','k','MarkerSize',9);
    plot3(cl_mm(1,end), cl_mm(2,end), cl_mm(3,end), 'rp', ...
          'MarkerFaceColor','r','MarkerSize',14);
    text(cl_mm(1,1),   cl_mm(2,1),   cl_mm(3,1)-8,   '  Base', ...
         'FontSize', 9, 'FontWeight', 'bold');
    text(cl_mm(1,end), cl_mm(2,end), cl_mm(3,end)+8, '  Tip', ...
         'FontSize', 9, 'FontWeight', 'bold', 'Color', 'r');
    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
end




%% ========================================================================
%                            Model functions
% ==========================================================================
function out = solveSpringConfig(dL_input, alpha, geom, mat, N, solver_opts)
%SOLVESPRINGCONFIG Solve the spring deformation under prescribed
%cable displacements.

    %% --- Step 1: Normalized midpoint coordinate u_i = (i-0.5)/N ---
    % The geometry and stiffness values at the segment midpoint are used to
    % represent the whole segment
    u_mid = ((1:N) - 0.5) / N;

    %% --- Step 2: Local radius r_i at each segment midpoint r_i ---
    % r_ji is the radial distance from the cable guide hole to the neutral
    % axis at the current segment
    r_i  = geom.r_w + (geom.r_e - geom.r_w) .* (2*u_mid - 1).^2;
    r_ji = r_i + geom.delta;

    %% --- Step 3: Pitch assignment ---
    % The turn index for each segment midpoint is
    % idx_turn = floor(u_mid*n_turn) + 1.
    % The local pitch p_i of this segment is assigned from geom.p_per_turn.
    p_i = zeros(1, N);
    for k = 1:N
        idx_turn = min(floor(u_mid(k) * geom.n_turn) + 1, geom.n_turn);
        p_i(k) = geom.p_per_turn(idx_turn);
    end
    n_per_turn = N / geom.n_turn;
    dL_seg_raw = p_i / n_per_turn;
    dL_seg     = dL_seg_raw * (geom.L0 / sum(dL_seg_raw));

    %% --- Step 4: Local helix angle and equivalent bending stiffness ---
     % Local helix angle: beta_i = atan(p_i/(2*pi*r_i)).
    beta_i    = atan(p_i ./ (2*pi*r_i));
   
    EI_factor = (mat.E * pi * geom.d_w^4) / 64;
    EI_i      = EI_factor .* 2 .* sin(beta_i) ./ ...
                ( 2 + mat.nu .* cos(beta_i).^2 );

    %% --- Step 5: Local equivalent axial stiffness of each segment ---
    EA_i = mat.G * geom.d_w^4 * p_i ./ (8 * (2*r_i).^3);

    %% --- Step 6: Geometry-stiffness integration constants K_a, K_b ---
    K_b = sum( r_ji.^2 ./ EI_i .* dL_seg );    % Bending-related
    K_a = sum( dL_seg ./ EA_i );               % Axial-related

    %% --- Step 7: Construct the four-dimensional KKT residual function ---
    function res = kkt4(x)
        phi = x(1); lam = x(2:4);
        S_c = lam(1)*cos(phi-alpha(1)) + lam(2)*cos(phi-alpha(2)) + lam(3)*cos(phi-alpha(3));
        S_s = lam(1)*sin(phi-alpha(1)) + lam(2)*sin(phi-alpha(2)) + lam(3)*sin(phi-alpha(3));
        S_l = sum(lam);
        
        eq_a = S_s;
         
        eq_b = zeros(3, 1);
        for j = 1:3
            eq_b(j) = dL_input(j) - ( -S_l*K_a - S_c*K_b*cos(phi - alpha(j)) );
        end
        res = [eq_a; eq_b];
    end

    %% --- Step 8: Estimate the initial guess ---
    % A simplified single-arc constant-curvature model, together with the
    % dominant axial contraction effect, is used to estimate the initial solution.
   
    A = sum(r_ji .* dL_seg);
    M = [cos(alpha(:)), sin(alpha(:))];
    PQ = M \ (dL_input(:) / A);
    phi_0 = atan2(PQ(2), PQ(1));
    lam_0 = (dL_input(:) / (-3*K_a))';
    x0 = [phi_0; lam_0(:)];

    %% --- Step 9: Solve the nonlinear KKT system ---
    [x_sol, ~, exitflag] = fsolve(@kkt4, x0, solver_opts);
    res_norm = norm(kkt4(x_sol));

    %% --- Step 10: Recover segment variables  ---
    phi    = x_sol(1);
    lambda = x_sol(2:4);
    S_c = lambda(1)*cos(phi-alpha(1)) + lambda(2)*cos(phi-alpha(2)) + ...
          lambda(3)*cos(phi-alpha(3));
    S_l = sum(lambda);
     
    kappa = -r_ji(:) ./ EI_i(:) * S_c;
     
    eps_v = -S_l ./ EA_i(:);

    %% --- Step 11:Compute the end pose through homogeneous transforms---
    T = eye(4);
    centerline = zeros(3, N+1);
    centerline(:, 1) = T(1:3, 4);
    for i = 1:N
        theta_i = abs(kappa(i)) * dL_seg(i) ;                   
        dL_pos  = dL_seg(i) * (1 + eps_v(i));               
        T = T * segment_transform(theta_i, phi, dL_pos);
        centerline(:, i+1) = T(1:3, 4);
    end
    % Tip bending angle: angle between the base z-axis and tip z-axis.
    z_tip  = T(1:3, 3);
    cos_th = max(-1, min(1, dot(z_tip, [0; 0; 1])));
    theta_tip_rad = acos(cos_th);

    % Overall axial length change of the spring
    axial_change = sum(eps_v(:) .* dL_seg(:));

    %% --- Step 12: Output variables  ---
    out.theta_tip_deg = theta_tip_rad * 180/pi;         % 末端弯曲角（°）
    out.kappa         = kappa(:);                       % N×1
    out.phi           = phi;                            % 标量
    out.eps           = eps_v(:);                       % N×1
    out.lambda        = lambda(:);                      % 3×1
    out.centerline    = centerline;                     % 3×(N+1)
    out.axial_change  = axial_change;                   % 标量（m）
    out.res_norm      = res_norm;                       % KKT 残差范数
    out.exitflag      = exitflag;                       % fsolve 标志
    out.EI_i          = EI_i(:);
    out.EA_i          = EA_i(:);
    out.r_i           = r_i(:);
    out.beta_i        = beta_i(:);
    out.p_i           = p_i(:);
    out.dL_seg        = dL_seg(:);
end


function T = segment_transform(theta, phi, dL_pos)
% SEGMENT_TRANSFORM Homogeneous transformation matrix of one segment.
% T = R_z(phi) * T_y(theta) * R_z(-phi)
% theta  : total bending angle of the segment (rad)
% phi    : bending-plane angle (rad)
% dL_pos : actual axial contribution length of the segment

    if abs(theta) < 1e-12
        T = eye(4);
        T(1:3, 4) = [0; 0; dL_pos];
        return;
    end
    c  = cos(theta); s  = sin(theta);
    cp = cos(phi);   sp = sin(phi);
    radius = dL_pos / theta;
    p = radius * [(1 - c)*cp; (1 - c)*sp; s];
    R = [ cp^2*(c-1)+1,  cp*sp*(c-1),    cp*s;
          cp*sp*(c-1),   sp^2*(c-1)+1,   sp*s;
         -cp*s,         -sp*s,           c    ];
    T = [R, p; 0, 0, 0, 1];
end
