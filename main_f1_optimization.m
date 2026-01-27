%% F1 TRAJECTORY OPTIMIZATION
% % 1. Setup & Data (Common)
% 
% % 2. V1: Geometric (Shortest Path)
% 
% % 3. V2: Dynamic Mechanical (No Aero - Static Grip Limit)
% 
% % 4. V3: Dynamic Aero (With Downforce - Variable Grip Limit)
% 
% % 5. Comparison & Plots

clear; clc; close all;

% 1. SETUP & DATA LOADING

fprintf('=== 1. SETUP & TRACK DATA ===\n');

%% 1.1 Load Data

try
    track_data = readmatrix('Circuits_Data/Nuerburgring_track.csv');
catch
    error('Track file not found. Please check the file path.');
end

x_center = track_data(:,1);
y_center = track_data(:,2);
w_right  = track_data(:,3);
w_left   = track_data(:,4);
N = length(x_center);
%--- CRITICAL FIX: DATA CLEANING ---
% Remove duplicate points or points too close together (ds < 0.1m)
% This prevents "Division by Zero" (NaN) in curvature calculations
dist_sq = [1; (diff(x_center).^2 + diff(y_center).^2)];
keep_idx = dist_sq > 0.01; % Keep points at least 10cm apart
x_center = x_center(keep_idx);
y_center = y_center(keep_idx);
w_right  = w_right(keep_idx);
w_left   = w_left(keep_idx);
% Re-check loop closure (if track is circular, ensure start != end in data, we constrain it later)
if sqrt((x_center(1)-x_center(end))^2 + (y_center(1)-y_center(end))^2) < 1.0
    fprintf('   -> Removed duplicate end point for clean looping.\n');
    x_center(end) = []; y_center(end) = [];
    w_right(end) = [];  w_left(end) = [];
end
%% 1.2 Vehicle Parameters

L = 2.5;                % Wheelbase [m]
delta_max = deg2rad(20);
delta_min = -deg2rad(20);

% Physics Parameters for Speed Profile
param.a_lat_max = 35;   % ~3.5 G (Lateral Grip)
param.a_acc_max = 12;   % Engine Power Limit
param.a_brk_max = 45;   % ~4.5 G (Braking Grip)
param.v_max     = 95;   % ~340 km/h (Top Speed)

% Aero Parameters
param.mu   = 1.5;       % Tire Friction
param.rho  = 1.225;     % Air Density
param.Cl_A = 2.5; % Réduire de 3.8 à 2.5 pour éviter l'effet "rail"param.m    = 750;       % Mass [kg]
param.g    = 9.81;
param.m = 798;

% Aero Coefficients
coeff_aero   = (param.m  * 0.5 * param.rho * param.Cl_A) / param.m;
limit_static = param.mu * param.g;
%% 1.3 Pre-compute Track Headings & Boundaries

% Apply cubic smoothing spline to remove numerical noise
% p is the smoothing parameter: 0.999 is close to the original data but smooths jitters
p_smooth = 0.999; 
x_ref_smooth = csaps(1:N, x_center, p_smooth, 1:N)';
y_ref_smooth = csaps(1:N, y_center, p_smooth, 1:N)';

theta_track = zeros(N,1);
for k = 2:N-1
    theta_track(k) = atan2(y_ref_smooth(k+1)-y_ref_smooth(k-1), x_ref_smooth(k+1)-x_ref_smooth(k-1));
end
theta_track(1) = atan2(y_ref_smooth(2)-y_ref_smooth(N), x_ref_smooth(2)-x_ref_smooth(N));
theta_track(N) = atan2(y_ref_smooth(1)-y_ref_smooth(N-1), x_ref_smooth(1)-x_ref_smooth(N-1));
theta_track = unwrap(theta_track);


% Left and right boundaries (absolute coordinates)
x_L = x_center - w_left .* sin(theta_track);
y_L = y_center + w_left .* cos(theta_track);
x_R = x_center + w_right .* sin(theta_track);
y_R = y_center - w_right .* cos(theta_track);

% Reference data package for solver
ref.x = x_center; 
ref.y = y_center; 
ref.th = theta_track;
ref.xL = x_L; 
ref.yL = y_L;
ref.xR = x_R;
ref.yR = y_R;
% Solver Settings
ops = sdpsettings('solver','gurobi','verbose',0,'usex0',1);

fprintf('   -> Setup Complete. N = %d points.\n', N);

%% 2. PHASE V1: GEOMETRIC OPTIMIZATION (Standard)

fprintf('\n=== 2. RUNNING V1 (Geometric) ===\n');
% YALMIP Variables (Decision Variables)
x=sdpvar(N,1); 
y=sdpvar(N,1); 
theta=sdpvar(N,1); 
delta=sdpvar(N-1,1);
z=sdpvar(N,1);
dtheta=sdpvar(N,1);
slack_cycle = sdpvar(4,1);
% --- INITIALIZATION ---
assign(x, x_center); assign(y, y_center);
assign(theta, theta_track); assign(z, 0.5*ones(N,1));
assign(slack_cycle, zeros(4,1));
%% --- Constraints V1 ---

cons = [];
%--- ROBUST CYCLIC CONSTRAINTS ---
% Allow tiny slack instead of hard equality to prevent crashing
cons = [cons, x(1) - x(N) == slack_cycle(1)];
cons = [cons, y(1) - y(N) == slack_cycle(2)];
cons = [cons, z(1) - z(N) == slack_cycle(3)];
cons = [cons, delta(1) - delta(end) == slack_cycle(4)];; 
% RELATIVE Heading Constraint (Handles 2pi wrap correctly)
cons = [cons, (theta(1) - ref.th(1)) == (theta(N) - ref.th(N))];
% Boundaries & Dynamics
cons = [cons, x == x_L + z.*(x_R-x_L), y == y_L + z.*(y_R-y_L), 0 <= z <= 1];
cons = [cons, dtheta == theta - ref.th];
% Removed hard +/- 90 constraint to avoid feasibility issues
cons = [cons, delta_min <= delta <= delta_max];

% Geometric Kinematics
for k = 1:N-1
    d = sqrt((ref.x(k+1)-ref.x(k))^2 + (ref.y(k+1)-ref.y(k))^2);
    cons = [cons, ((x(k+1)-x(k))*(-sin(ref.th(k))) + (y(k+1)-y(k))*cos(ref.th(k))) == d * dtheta(k)];
    cons = [cons, theta(k+1) == theta(k) + (d/L)*delta(k)];
end
% --- Objective V1 (Minimize Curvature & Length) ---
dx_seg = diff(x); 
dy_seg = diff(y); 
kappa=sdpvar(N-1,1);

cons_obj = [];
for k = 1:N-1
    d = sqrt((ref.x(k+1)-ref.x(k))^2 + (ref.y(k+1)-ref.y(k))^2);
    cons_obj = [cons_obj, kappa(k) == (theta(k+1)-theta(k))/d];
end
% Weights (V1 Standard)
w_len = 10; w_smooth = 1000; w_curv = 1e5;

obj_v1 = w_len * sum(dx_seg.^2 + dy_seg.^2) + ... % Approx length
      w_smooth * sum(diff(delta).^2) + ...        % Smooth inputs
      w_curv * (kappa'*kappa);                 % Minimize curvature

sol_v1 = optimize([cons, cons_obj], obj_v1, ops);


% --- CHECK LINEARIZATION VALIDITY (Insert at end of Phase 2) ---
if sol_v1.problem == 0
    % Calculate heading error relative to reference
    th_error = value(theta) - ref.th;
    
    % Wrap error to [-pi, pi] to handle the loop reset correctly
    th_error = atan2(sin(th_error), cos(th_error));
    
    max_err_deg = max(abs(rad2deg(th_error)));
    fprintf('   -> Max Heading Deviation: %.2f deg\n', max_err_deg);
    
    if max_err_deg > 35
        warning('LINEARIZATION ISSUE: Deviation > 30 deg. The model assumption sin(x)=x is failing.');
        fprintf('      Fix: Increase w_len weight or smooth the reference centerline.\n');
    else
        fprintf('   -> Linearization Assumption: VALID (Model is accurate)\n');
    end

else
    warning('V1 Optimization failed! Check constraints.');
end
% Store V1 Results (Keep this line)
res_v1.x = value(x); res_v1.y = value(y); res_v1.th = value(theta);
% NaN Check before Phase 3
if any(isnan(res_v1.x))
    error('Optimization returned NaNs. Check track data for duplicates.');
end

%--- Calculate Speed for V1 (to get a Lap Time) ---
% We calculate the speed the car could assume on this path using V2
% physics.
% This serves as our baseline time.
[res_v1.v, res_v1.t] = solve_speed_profile(res_v1.x, res_v1.y, res_v1.th, limit_static, 0, param, ops, 'diamond');
fprintf('   -> V1 Lap Time (Geom path, Static Grip): %.3f s\n', res_v1.t);
%% 3. PHASE V2: DYNAMIC MECHANICAL (No Aero)

fprintf('\n=== 3. RUNNING V2 (Dynamic - NO Aero) ===\n');
% Initialize with V1 path
path_curr = res_v1;
v_profile = res_v1.v; 
max_iter = 3;
for iter=1:max_iter
    fprintf('   -> V2 Iteration %d/%d... ', iter, max_iter);
    
    % A. Path Optimization (Weighted by speed)
    % Note: Even without aero, faster sections benefit from straighter lines
    w_vec = 1e5 * (v_profile(1:N-1)/param.v_max).^1.5 + 1000;
    path_curr = solve_path(path_curr, w_vec, ref, L, delta_min, delta_max, ops);
    
    % B. Speed Optimization (NO AERO TERM)
    % We pass '0' as the aero coefficient
    [v_profile, t_curr] = solve_speed_profile(path_curr.x, path_curr.y, path_curr.th, limit_static, 0, param, ops,'diamond');
    
    fprintf('Time: %.3f s\n', t_curr);
end
res_v2 = path_curr; res_v2.v = v_profile; res_v2.t = t_curr;
98.546
%% === 4. PHASE V3: DYNAMIC AERO (With Downforce) ===

fprintf('\n=== 4. RUNNING V3 (Dynamic - WITH Aero) ===\n');

% Initialize with V1 path (to compare fairly how V3 deviates)
path_curr = res_v1;
v_profile = res_v1.v; 
max_iter = 3;

for iter=1:max_iter
    fprintf('   -> V3 Iteration %d/%d... ', iter, max_iter);
    % A. Path Optimization
    w_vec = 1e5 * (v_profile(1:N-1)/param.v_max).^2 + 1000; % Stronger weight for Aero
    path_curr = solve_path(path_curr, w_vec, ref, L, delta_min, delta_max, ops);
    
    % B. Speed Optimization (WITH AERO TERM)
    % We pass 'coeff_aero' here
    [v_profile, t_curr] = solve_speed_profile(path_curr.x, path_curr.y, path_curr.th, limit_static, coeff_aero, param, ops,'ellipse');
    
    fprintf('Time: %.3f s\n', t_curr);
end
res_v3 = path_curr; res_v3.v = v_profile; res_v3.t = t_curr;
%%
%% === 5. COMPARISON & PLOTS ===
fprintf('\n=== 5. FINAL ANALYSIS ===\n');

% Calculate Distance Axis
dx=diff(x_center); dy=diff(y_center); ds_ref=[sqrt(dx.^2+dy.^2); 0]; d_cum = cumsum(ds_ref);

% --- Figure 1: Trajectories ---
figure('Name','Trajectory Comparison','Color','w');
plot(ref.xL, ref.yL, 'w', 'LineWidth',0.5); hold on;
plot(ref.xR, ref.yR, 'w', 'LineWidth',0.5);
% Plot V1, V2, V3
h1 = plot(res_v1.x, res_v1.y, 'g--', 'LineWidth', 1.5, 'DisplayName', 'V1: Geometric');
h2 = plot(res_v2.x, res_v2.y, 'b-',  'LineWidth', 2,   'DisplayName', 'V2: Mechanical (No Aero)');
h3 = plot(res_v3.x, res_v3.y, 'r-',  'LineWidth', 2,   'DisplayName', 'V3: Aero (Downforce)');
legend([h1 h2 h3]); axis equal; grid on; title('Trajectory Deviation');

% --- Figure 2: Speed Profiles ---
figure('Name','Speed Comparison','Color','w');
plot(d_cum, res_v1.v*3.6, 'g--', 'LineWidth',1.5, 'DisplayName',['V1 (Static) - ' num2str(res_v1.t,'%.2f') 's']); hold on;
plot(d_cum, res_v2.v*3.6, 'b-',  'LineWidth',2,   'DisplayName',['V2 (Mech) - ' num2str(res_v2.t,'%.2f') 's']);
plot(d_cum, res_v3.v*3.6, 'r-',  'LineWidth',2,   'DisplayName',['V3 (Aero) - ' num2str(res_v3.t,'%.2f') 's']);
xlabel('Dist [m]'); ylabel('Speed [km/h]'); legend; grid on; title('Speed Profiles');

% --- Figure 3: GG Diagram (V2 vs V3) ---
figure('Name','GG Diagram','Color','w');

% --- CALC ACCELS FOR V2 (NO AERO) ---
dx = diff(res_v2.x); dy = diff(res_v2.y); 
ds = sqrt(dx.^2+dy.^2); ds = [ds; ds(end)]; % ds devient taille N

dth = diff(res_v2.th); 
dth(dth>pi) = dth(dth>pi)-2*pi; 
dth(dth<-pi) = dth(dth<-pi)+2*pi;
dth = [dth; 0];

k2 = abs(dth) ./ (ds + 1e-3); % Maintenant N ./ N fonctionne !

lat2 = res_v2.v.^2 .* k2; 
% Pour lon, on utilise ds(1:end-1) car diff réduit la taille de 1
lon2 = [diff(res_v2.v.^2)./(2*ds(1:end-1) + 1e-3); 0]; 

% --- CALC ACCELS FOR V3 (AERO) ---
dx = diff(res_v3.x); dy = diff(res_v3.y); 
ds = sqrt(dx.^2+dy.^2); ds = [ds; ds(end)];

dth = diff(res_v3.th); 
dth(dth>pi) = dth(dth>pi)-2*pi; 
dth(dth<-pi) = dth(dth<-pi)+2*pi;
dth = [dth; 0]; 

k3 = abs(dth) ./ (ds + 1e-3);

lat3 = res_v3.v.^2 .* k3; 
lon3 = [diff(res_v3.v.^2)./(2*ds(1:end-1) + 1e-3); 0];

% --- PLOT ---
scatter(lat2, lon2, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'V2 (No Aero)'); hold on;
scatter(lat3, lon3, 15, 'r', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'V3 (Aero)');

% Draw Limits
th_c = linspace(0,2*pi,100);
plot(limit_static*cos(th_c), limit_static*sin(th_c), 'w--', 'LineWidth',2, 'DisplayName','Static Limit');

legend('Location','best'); axis equal; grid on; 
title('GG Diagram: Grip Usage Comparison');
xlabel('Lateral Acc [m/s^2]'); ylabel('Longitudinal Acc [m/s^2]');

% --- Save Data ---
save('Comparison_Data.mat', 'res_v1', 'res_v2', 'res_v3');
fprintf('\nDone. All results saved.\n');
%% 
% Helper funtion


% 2. Helper for Speed Optimization (Unified V2 & V3)
function [v_prof, t_lap] = solve_speed_profile(x, y, th, lim_stat, coeff_aero, p, ops, mode)
    % Inputs:
    %   ...
    %   mode : 'diamond' (V2 - QP - Kamm Diamond) ou 'ellipse' (V3 - SOCP - Friction Circle)
    
    N = length(x);
    
    % 1. Calculs géométriques (Courbure, distances)
    dx = diff(x); dy = diff(y); 
    ds = sqrt(dx.^2+dy.^2); ds = [ds; ds(end)]; % Pad end
    
    dth = diff(th); 
    dth(dth>pi) = dth(dth>pi)-2*pi; 
    dth(dth<-pi) = dth(dth<-pi)+2*pi;
    
    % Courbure (k = dtheta/ds)
    k_opt = abs([dth; 0]) ./ (ds + 1e-3);
    
    % 2. Variables d'optimisation (E = v^2)
    E = sdpvar(N,1);
    
    % Contraintes de base (Positivité, Départ=Arrivée, Vmax)
    cons = [0 <= E <= p.v_max^2, E(1) == E(end)];
    
    % 3. Boucle sur les points
    for k = 1:N-1
        % Définition des accélérations requises
        lat_acc = E(k) * k_opt(k);             % Latérale (Centrifuge)
        lon_acc = (E(k+1)-E(k)) / (2*ds(k));   % Longitudinale (Moteur/Frein)
        
        % Limite d'adhérence disponible (Statique + Aéro)
        grip_available = lim_stat + coeff_aero * E(k);
        
        % --- LE CHOIX DU MODELE PHYSIQUE ---
        if strcmp(mode, 'diamond')
            % MODE V2: Kamm Diamond (Somme des valeurs absolues)
            % Compatible avec les solveurs QP (rapides)
            % |lat| + |lon| <= Grip
            cons = [cons, norm([lat_acc; lon_acc], 1) <= grip_available];
            
        elseif strcmp(mode, 'ellipse')
            % MODE V3: Friction Circle/Ellipse (Norme Euclidienne)
            % Nécessite un solveur SOCP (plus précis pour la physique réelle)
            % sqrt(lat^2 + lon^2) <= Grip
            cons = [cons, norm([lat_acc; lon_acc], 2) <= grip_available];
        end
        %G constraint
        max_lateral_g = 5.5 * 9.81; 
        cons = [cons, lat_acc <= max_lateral_g];
        % Contrainte Moteur (Uniquement positive, le frein est géré par le grip)
        cons = [cons, lon_acc <= p.a_acc_max];
    end
    
    % 4. Résolution (Maximiser l'énergie cinétique -> Minimiser temps)
    optimize(cons, -sum(E), ops);
    
    % 5. Reconstruction des résultats
    v_prof = sqrt(value(E));
    
    % Calcul du temps (t = d / v_moyen)
    v_avg = 0.5 * (v_prof(1:end-1) + v_prof(2:end)) + 1e-3;
    t_lap = sum(ds(1:end-1) ./ v_avg);
end

function new_path = solve_path(prev_path, w_vec , ref, L, d_min, d_max, ops)
    % Inputs: w_vec (vector), w_len (scalar), w_smooth (scalar)
    
    N = length(ref.x);
    x=sdpvar(N,1); y=sdpvar(N,1); theta=sdpvar(N,1); delta=sdpvar(N-1,1);
    z=sdpvar(N,1); dtheta=sdpvar(N,1); kappa=sdpvar(N-1,1); slack=sdpvar(4,1);
    
    % Init
    assign(x, prev_path.x); assign(y, prev_path.y); assign(theta, prev_path.th); assign(z, 0.5*ones(N,1));
    
    cons = [x(1)-x(N)==slack(1), y(1)-y(N)==slack(2), z(1)-z(N)==slack(3), delta(1)-delta(end)==slack(4)];
    cons = [cons, (theta(1)-ref.th(1)) == (theta(N)-ref.th(N))];
    cons = [cons, x == ref.xL + z.*(ref.xR-ref.xL), y == ref.yL + z.*(ref.yR-ref.yL), 0<=z<=1];
    cons = [cons, dtheta == theta - ref.th];
    cons = [cons, d_min <= delta <= d_max];
    
    % Rate of change delta constraint
    d_delta_max = deg2rad(5);
    for k=1:N-2, cons=[cons, -d_delta_max <= delta(k+1)-delta(k) <= d_delta_max]; end

    for k = 1:N-1
        d = sqrt((ref.x(k+1)-ref.x(k))^2 + (ref.y(k+1)-ref.y(k))^2);
        cons = [cons, ((x(k+1)-x(k))*(-sin(ref.th(k))) + (y(k+1)-y(k))*cos(ref.th(k))) == d * dtheta(k)];
        cons = [cons, theta(k+1) == theta(k) + (d/L)*delta(k)];
        cons = [cons, kappa(k) == (theta(k+1)-theta(k))/d];
    end
    % --- OBJECTIVE FUNCTION WITH DYNAMIC WEIGHTS ---
    obj = 0.20 * sum(diff(x).^2 + diff(y).^2) + ...   % Poids Longueur (Grid)
          200 * sum(diff(delta).^2) + ...         % Poids Volant (Grid)
          sum(w_vec .* (kappa.^2)) + ...               % Poids Courbure Adaptatif (Grid)
          1e9*(slack'*slack);
          
    optimize(cons, obj, ops);
    new_path.x = value(x); new_path.y = value(y); new_path.th = value(theta);
end
%%
%% === 4.2 WEIGHT GRID SEARCH (EXPLICIT) ===
fprintf('\n=== 4.2 WEIGHT GRID SEARCH (Explicit Loop) ===\n');

% 1. Settings
ops_grid = sdpsettings('solver','gurobi','verbose',0,'usex0',1,...
    'gurobi.TimeLimit',8, 'gurobi.Method',2); 

% 2. Weights to test
w_len_base = 0.1; 
w_smooth_base = 100;
w_scale_base = 2.5e5;

w_len_grid    = [w_len_base*0.5, w_len_base*2];
w_smooth_grid = [w_smooth_base*0.5, w_smooth_base*2];
w_scale_grid  = [w_scale_base*0.5, w_scale_base*2];

results = struct('w_len',{},'w_smooth',{},'w_scale',{},'t_final',{},'path_v2',{},'v_profile',{});
trial = 0;

% 3. Warm Start Setup (Using V1 Result)
path_warm.x = res_v1.x; 
path_warm.y = res_v1.y; 
path_warm.th = res_v1.th;
v_ref = res_v1.v; % Speed reference for adaptive weights

fprintf('Starting 8 explicit trials...\n');

for i = 1:2
    for j = 1:2
        for k = 1:2
            trial = trial + 1;
            
            % Current Weights
            w_len = w_len_grid(i);
            w_smooth = w_smooth_grid(j);
            w_scale = w_scale_grid(k);
            
            fprintf('Trial %d/8: Len=%.2f Smth=%.0f Scale=%.0e ... ', trial, w_len, w_smooth, w_scale);
            drawnow; % Force display before calculation
            
            % ---------------------------------------------------------
            % STEP A: EXPLICIT PATH OPTIMIZATION (No Helper Function)
            % ---------------------------------------------------------
            
            % Define Variables
            x=sdpvar(N,1); y=sdpvar(N,1); theta=sdpvar(N,1); delta=sdpvar(N-1,1);
            z=sdpvar(N,1); dtheta=sdpvar(N,1); kappa=sdpvar(N-1,1); slack=sdpvar(4,1);
            
            % Warm Start Assignments
            assign(x, path_warm.x); assign(y, path_warm.y); 
            assign(theta, path_warm.th); assign(z, 0.5*ones(N,1));
            
            % Constraints
            cons = [x(1)-x(N)==slack(1), y(1)-y(N)==slack(2), z(1)-z(N)==slack(3), delta(1)-delta(end)==slack(4)];
            cons = [cons, (theta(1)-ref.th(1)) == (theta(N)-ref.th(N))];
            cons = [cons, x == ref.xL + z.*(ref.xR-ref.xL), y == ref.yL + z.*(ref.yR-ref.yL), 0<=z<=1];
            cons = [cons, dtheta == theta - ref.th];
            cons = [cons, delta_min <= delta <= delta_max];
            
            % Smoothness
            d_delta_max = deg2rad(5);
            for m=1:N-2, cons=[cons, -d_delta_max <= delta(m+1)-delta(m) <= d_delta_max]; end
            
            % Kinematics
            for m = 1:N-1
                d = sqrt((ref.x(m+1)-ref.x(m))^2 + (ref.y(m+1)-ref.y(m))^2);
                cons = [cons, ((x(m+1)-x(m))*(-sin(ref.th(m))) + (y(m+1)-y(m))*cos(ref.th(m))) == d * dtheta(m)];
                cons = [cons, theta(m+1) == theta(m) + (d/L)*delta(m)];
                cons = [cons, kappa(m) == (theta(m+1)-theta(m))/d];
            end
            
            % Adaptive Weights Calculation
            speed_ratio = (v_ref(1:N-1) / param.v_max);
            w_vec = w_scale * (speed_ratio.^2) + 1e4;
            
            % Objective
            obj = w_len * sum(diff(x).^2 + diff(y).^2) + ...
                  w_smooth * sum(diff(delta).^2) + ...
                  sum(w_vec .* (kappa.^2)) + ...
                  1e9 * (slack'*slack);
            
            % Solve Path
            sol = optimize(cons, obj, ops_grid);
            
            if sol.problem ~= 0
                fprintf('PATH FAILED\n');
                results(trial).t_final = Inf;
                continue; % Skip to next trial
            end
            
            % Save Path
            path_tmp.x = value(x); path_tmp.y = value(y); path_tmp.th = value(theta);
            
            % ---------------------------------------------------------
            % STEP B: EXPLICIT SPEED OPTIMIZATION (V3 AERO PHYSICS)
            % ---------------------------------------------------------
            
            % Geometry calc on new path
            dx = diff(path_tmp.x); dy = diff(path_tmp.y); 
            ds = sqrt(dx.^2+dy.^2); ds = [ds; ds(end)];
            dth = diff(path_tmp.th);
            dth(dth>pi) = dth(dth>pi)-2*pi; dth(dth<-pi) = dth(dth<-pi)+2*pi;
            k_opt = abs([dth; 0]) ./ (ds+1e-3);
            
            % Speed Variables
            E = sdpvar(N,1);
            cons_v = [0 <= E <= param.v_max^2, E(1)==E(end)];
            
            for m = 1:N-1
                lat_acc = E(m) * k_opt(m);
                lon_acc = (E(m+1)-E(m))/(2*ds(m));
                
                % AERO CONSTRAINT (V3 Physics: Ellipse + Downforce)
                grip_avail = limit_static + coeff_aero * E(m);
                cons_v = [cons_v, norm([lat_acc; lon_acc], 2) <= grip_avail];
                cons_v = [cons_v, lon_acc <= param.a_acc_max];
            end
            
            % Solve Speed
            sol_v = optimize(cons_v, -sum(E), ops_grid);
            
            if sol_v.problem ~= 0
                fprintf('SPEED FAILED\n');
                results(trial).t_final = Inf;
            else
                v_tmp = sqrt(value(E));
                t_tmp = sum(ds(1:end-1)./(0.5*(v_tmp(1:end-1)+v_tmp(2:end))+1e-3));
                fprintf('-> %.3f s\n', t_tmp);
                
                % Store
                results(trial).w_len = w_len;
                results(trial).w_smooth = w_smooth;
                results(trial).w_scale = w_scale;
                results(trial).t_final = t_tmp;
                results(trial).path_v2 = path_tmp;
                results(trial).v_profile = v_tmp;
            end
        end
    end
end

% --- RESULT SELECTION ---
[t_best, idx_best] = min([results.t_final]);
best = results(idx_best);

res_v4 = best.path_v2; 
res_v4.v = best.v_profile; 
res_v4.t = best.t_final;

fprintf('\n🏆 BEST WEIGHTS: w_len=%.2f w_smooth=%.0f w_scale=%.0e\n', ...
    best.w_len, best.w_smooth, best.w_scale);
fprintf('⏱️  V4 Time: %.3f s (Gain vs V3: %.3f s)\n', t_best, res_v3.t - t_best);

%% === PLOT GRID SEARCH RESULTS ===
figure('Name','BEST WEIGHT PATH','Position',[100 100 1200 700],'Color','w');

% 1. Map Comparison
subplot(1,2,1);
plot(ref.xL, ref.yL, 'k-', 'LineWidth',1.2,'HandleVisibility','off'); hold on;
plot(ref.xR, ref.yR, 'k-', 'LineWidth',1.2,'HandleVisibility','off');
plot(ref.x, ref.y, 'k--', 'LineWidth',0.5, 'DisplayName','Centerline');

% Plot all trials (thin lines)
for t=1:8
    if results(t).t_final < Inf
        plot(results(t).path_v2.x, results(t).path_v2.y, '-', ...
             'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility','off');
    end
end

% Plot Winner
plot(res_v4.x, res_v4.y, 'r-', 'LineWidth',2.0, 'DisplayName',sprintf('V4 BEST (%.2fs)', t_best));
plot(res_v1.x, res_v1.y, 'g--', 'LineWidth',1.5, 'DisplayName','V1 Geometric');

title('Grid Search: All Paths vs Winner', 'FontSize',14);
legend('Location','best'); axis equal; grid on;
xlabel('X [m]'); ylabel('Y [m]');

% 2. Telemetry of Winner
subplot(1,2,2);
dx = diff(res_v4.x); dy = diff(res_v4.y);
d_cum = [0; cumsum(sqrt(dx.^2 + dy.^2))];

yyaxis left;
plot(d_cum, res_v4.v*3.6, 'r-', 'LineWidth',2);
ylabel('Speed [km/h]','FontSize',12,'Color','r');

yyaxis right;
ds_seg = sqrt(dx.^2 + dy.^2); 
dth = diff(res_v4.th);
dth(dth>pi)=dth(dth>pi)-2*pi; dth(dth<-pi)=dth(dth<-pi)+2*pi;
kappa_profile = zeros(N,1); 
kappa_profile(1:end-1) = abs(dth)./(ds_seg+1e-3);
kappa_profile(end) = kappa_profile(end-1);

plot(d_cum, kappa_profile, 'b-', 'LineWidth',1);
ylabel('Curvature [1/m]','FontSize',12,'Color','b');

title(['V4 Telemetry (Time: ' num2str(t_best, '%.2f') 's)'], 'FontSize',14);
xlabel('Distance [m]'); grid on;
legend('Speed', 'Curvature');

sgtitle('GRID SEARCH RESULTS (V4)', 'FontSize', 16, 'FontWeight', 'bold');