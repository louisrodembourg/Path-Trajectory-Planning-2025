%% F1 TRAJECTORY OPTIMIZATION
% Project: Optimal Racing Line Generation
% Description: Multi-stage optimization (Geometric -> Mechanical -> Aero)
%              for a generic F1 car on a variation of the Nürburgring.

clear; clc; close all;

% ---------------------------------------------------------
% 1. SETUP & DATA LOADING
% ---------------------------------------------------------
fprintf('=== 1. SETUP & TRACK DATA ===\n');

% Load Track
try
    track_data = readmatrix('Circuits_Data/Nuerburgring_track.csv');
catch
    error('Track file not found in Circuits_Data/.');
end

x_center = track_data(:,1);
y_center = track_data(:,2);
w_right  = track_data(:,3);
w_left   = track_data(:,4);
N_raw = length(x_center);

% Data Cleaning (Remove points too close together to prevent instability)
dist_sq = [1; (diff(x_center).^2 + diff(y_center).^2)];
keep_idx = dist_sq > 0.01; 
x_center = x_center(keep_idx);
y_center = y_center(keep_idx);
w_right  = w_right(keep_idx);
w_left   = w_left(keep_idx);
N = length(x_center);

% 1.2 Vehicle Parameters
L = 2.5;                % Wheelbase [m]
delta_max = deg2rad(20);
delta_min = -deg2rad(20);

% Physics Parameters 
param.m    = 798;       % Mass [kg] include driver
param.g    = 9.81;
param.v_max     = 95;   % Top speed [m/s] (~340 km/h)
param.a_acc_max = 12;   % Max longitudinal acc [m/s^2]

% Aero Parameters
param.mu   = 1.5;       % Tire Friction Coefficient using Soft Compound
param.rho  = 1.225;     % Air Densitry
param.Cl_A = 2.5;       % Downforce Area Product (Adjusted)

% Derived Physics
coeff_aero   = (0.5 * param.rho * param.Cl_A) / param.m;
limit_static = param.mu * param.g;

% 1.3 Pre-compute Track Headings & Boundaries (Spline Smoothing)
p_smooth = 0.999; 
x_ref_smooth = csaps(1:N, x_center, p_smooth, 1:N)';
y_ref_smooth = csaps(1:N, y_center, p_smooth, 1:N)';

theta_track = zeros(N,1);
for k = 2:N-1
    theta_track(k) = atan2(y_ref_smooth(k+1)-y_ref_smooth(k-1), x_ref_smooth(k+1)-x_ref_smooth(k-1));
end
% Handle endpoints (Cyclic)
theta_track(1) = atan2(y_ref_smooth(2)-y_ref_smooth(N), x_ref_smooth(2)-x_ref_smooth(N));
theta_track(N) = atan2(y_ref_smooth(1)-y_ref_smooth(N-1), x_ref_smooth(1)-x_ref_smooth(N-1));
theta_track = unwrap(theta_track);

% Boundaries
x_L = x_center - w_left .* sin(theta_track);
y_L = y_center + w_left .* cos(theta_track);
x_R = x_center + w_right .* sin(theta_track);
y_R = y_center - w_right .* cos(theta_track);

% Structs for Functions
ref.x = x_center; ref.y = y_center; ref.th = theta_track;
ref.xL = x_L; ref.yL = y_L; ref.xR = x_R; ref.yR = y_R;
ops = sdpsettings('solver','gurobi','verbose',0,'usex0',1);

fprintf('   -> Track Loaded: %d points.\n', N);

% ---------------------------------------------------------
% 2. PHASE V1: GEOMETRIC OPTIMIZATION
% ---------------------------------------------------------
fprintf('\n=== 2. RUNNING V1 (Geometric) ===\n');

% Weights for V1
w_v1.len = 10; 
w_v1.smooth = 1000; 
w_v1.curv = 1e5;

% Run Helper Function (V1 is essentially solve_path with uniform Speed)
% We create a dummy speed profile for V1 weighting (just 1.0)
dummy_speed_ratio = zeros(N,1); 
path_init.x = x_center; path_init.y = y_center; path_init.th = theta_track;

res_v1 = solve_path(path_init, dummy_speed_ratio, w_v1, ref, L, delta_min, delta_max, ops);

% Calculate V1 Speed & Time (Baseline)
[res_v1.v, res_v1.t] = solve_speed_profile(res_v1, limit_static, 0, param, ops, 'diamond');
fprintf('   -> V1 Lap Time: %.3f s\n', res_v1.t);


% ---------------------------------------------------------
% 3. PHASE V2: DYNAMIC MECHANICAL (No Aero)
% ---------------------------------------------------------
fprintf('\n=== 3. RUNNING V2 (Dynamic - No Aero) ===\n');

res_v2 = res_v1; 
w_v2.len = 0.20; w_v2.smooth = 200; w_v2.curv_base = 1000; w_v2.curv_scale = 1e5;

for iter = 1:3
    fprintf('   -> Iter %d... ', iter);
    
    % Path Optimization
    speed_ratio = (res_v2.v / param.v_max).^1.5; % Less aggressive than Aero
    res_v2 = solve_path(res_v2, speed_ratio, w_v2, ref, L, delta_min, delta_max, ops);
    
    % Speed Optimization (Diamond Constraint for speed)
    [res_v2.v, res_v2.t] = solve_speed_profile(res_v2, limit_static, 0, param, ops, 'diamond');
    fprintf('Time: %.3f s\n', res_v2.t);
end


% ---------------------------------------------------------
% 4. PHASE V3: DYNAMIC AERO (With Downforce)
% ---------------------------------------------------------
fprintf('\n=== 4. RUNNING V3 (Dynamic - WITH Aero) ===\n');

res_v3 = res_v1; % Start fresh from V1 to compare divergence
w_v3 = w_v2;     % Use same weights base

for iter = 1:3
    fprintf('   -> Iter %d... ', iter);
    
    % Path Optimization (Curvature penalty increases with speed)
    speed_ratio = (res_v3.v / param.v_max).^2; 
    res_v3 = solve_path(res_v3, speed_ratio, w_v3, ref, L, delta_min, delta_max, ops);
    
    % Speed Optimization (Ellipse Constraint + Aero)
    [res_v3.v, res_v3.t] = solve_speed_profile(res_v3, limit_static, coeff_aero, param, ops, 'ellipse');
    fprintf('Time: %.3f s\n', res_v3.t);
end


% ---------------------------------------------------------
% 5. PHASE V4: GRID SEARCH (Hyperparameter Tuning)
% ---------------------------------------------------------
fprintf('\n=== 5. GRID SEARCH (Weights Tuning) ===\n');

% Search Space
w_len_grid    = [0.05, 0.2];
w_smooth_grid = [100, 400];
best_t = Inf;
res_v4 = res_v3;

fprintf('   -> Running simplified grid search...\n');
count = 0;
for wl = w_len_grid
    for ws = w_smooth_grid
        count = count + 1;
        
        % Setup Weights
        w_tmp = w_v3;
        w_tmp.len = wl;
        w_tmp.smooth = ws;
        
        % Run Single Pass Optimization
        p_tmp = solve_path(res_v3, (res_v3.v/param.v_max).^2, w_tmp, ref, L, delta_min, delta_max, ops);
        [v_tmp, t_tmp] = solve_speed_profile(p_tmp, limit_static, coeff_aero, param, ops, 'ellipse');
        
        fprintf('      [%d] W_len:%.2f W_sm:%.0f -> T: %.3f s\n', count, wl, ws, t_tmp);
        
        if t_tmp < best_t
            best_t = t_tmp;
            res_v4 = p_tmp;
            res_v4.v = v_tmp;
            res_v4.t = t_tmp;
            best_w = w_tmp;
        end
    end
end
fprintf('   -> Best Time: %.3f s (Gain: %.3f s)\n', best_t, res_v3.t - best_t);


% ---------------------------------------------------------
% 6. ANALYSIS & PLOTS
% ---------------------------------------------------------
save('Comparison_Data.mat', 'res_v1', 'res_v2', 'res_v3', 'res_v4');

% Calculate Distance Axis
dx=diff(x_center); dy=diff(y_center); ds_ref=[sqrt(dx.^2+dy.^2); 0]; d_cum = cumsum(ds_ref);

% --- Figure 1: Trajectories ---
figure('Name','Trajectory Comparison','Color','w');
plot(ref.xL, ref.yL, 'k', 'LineWidth',0.5); hold on;
plot(ref.xR, ref.yR, 'k', 'LineWidth',0.5);
h1 = plot(res_v1.x, res_v1.y, 'g--', 'LineWidth', 1.5, 'DisplayName', 'V1: Geometric');
h2 = plot(res_v2.x, res_v2.y, 'b-',  'LineWidth', 2,   'DisplayName', 'V2: Mechanical (No Aero)');
h4 = plot(res_v4.x, res_v4.y, 'r-',  'LineWidth', 2,   'DisplayName', 'V4: Best Aero');
legend([h1 h2 h4], 'Location', 'best'); axis equal; grid on; title('Trajectory Deviation');
xlabel('X [m]'); ylabel('Y [m]');

% --- Figure 2: Speed Profiles ---
figure('Name','Speed Comparison','Color','w');
plot(d_cum, res_v1.v*3.6, 'g--', 'LineWidth',1.5, 'DisplayName',['V1 (Static) - ' num2str(res_v1.t,'%.2f') 's']); hold on;
plot(d_cum, res_v2.v*3.6, 'b-',  'LineWidth',2,   'DisplayName',['V2 (Mech) - ' num2str(res_v2.t,'%.2f') 's']);
plot(d_cum, res_v4.v*3.6, 'r-',  'LineWidth',2,   'DisplayName',['V4 (Aero) - ' num2str(res_v4.t,'%.2f') 's']);
xlabel('Dist [m]'); ylabel('Speed [km/h]'); legend; grid on; title('Speed Profiles');

% --- Figure 3: GG Diagram (V2 vs V4) ---
figure('Name','GG Diagram','Color','w');

% Helper for GG (anonymous)
calc_acc = @(res) deal(...
    res.v.^2 .* (abs([diff(res.th);0]) ./ ([sqrt(diff(res.x).^2+diff(res.y).^2);1e-3] + 1e-3)), ... % Lat
    [diff(res.v.^2)./(2*[sqrt(diff(res.x).^2+diff(res.y).^2);1e-3] + 1e-3); 0] ... % Lon
);

[lat2, lon2] = calc_acc(res_v2);
[lat4, lon4] = calc_acc(res_v4);

scatter(lat2, lon2, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'V2 (No Aero)'); hold on;
scatter(lat4, lon4, 15, 'r', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'V4 (Aero)');

% Draw Static Limit Circle
th_c = linspace(0,2*pi,100);
plot(limit_static*cos(th_c), limit_static*sin(th_c), 'k--', 'LineWidth',2, 'DisplayName','Static Friction Limit');
legend('Location','best'); axis equal; grid on; 
title('GG Diagram: Grip Usage Comparison');
xlabel('Lateral Acc [m/s^2]'); ylabel('Longitudinal Acc [m/s^2]');

fprintf('\nDone. All results saved.\n');


% ---------------------------------------------------------
% 7. HELPER FUNCTIONS
% ---------------------------------------------------------

function new_path = solve_path(prev_path, speed_ratio, w, ref, L, d_min, d_max, ops)
    % Unpack Weights
    if isfield(w, 'curv')
        % Static weighting (V1)
        w_vec = w.curv;
    else
        % Dynamic weighting (V2/V3)
        w_vec = w.curv_scale * speed_ratio + w.curv_base;
    end
    
    N = length(ref.x);
    
    % YALMIP Variables
    x=sdpvar(N,1); y=sdpvar(N,1); theta=sdpvar(N,1); delta=sdpvar(N-1,1);
    z=sdpvar(N,1); dtheta=sdpvar(N,1); kappa=sdpvar(N-1,1); slack=sdpvar(4,1);
    
    % Warm Start
    assign(x, prev_path.x); assign(y, prev_path.y); 
    assign(theta, prev_path.th); assign(z, 0.5*ones(N,1));
    assign(slack, zeros(4,1));
    
    cons = [];
    % Cyclic Constraints
    cons = [cons, x(1)-x(N)==slack(1), y(1)-y(N)==slack(2)];
    cons = [cons, z(1)-z(N)==slack(3), delta(1)-delta(end)==slack(4)];
    cons = [cons, (theta(1)-ref.th(1)) == (theta(N)-ref.th(N))];
    
    % Boundaries
    cons = [cons, x == ref.xL + z.*(ref.xR-ref.xL), y == ref.yL + z.*(ref.yR-ref.yL), 0<=z<=1];
    
    % Kinematics
    cons = [cons, dtheta == theta - ref.th];
    cons = [cons, d_min <= delta <= d_max];
    
    % Delta Rate (Smoothness)
    d_delta_max = deg2rad(5);
    for k=1:N-2
        cons=[cons, -d_delta_max <= delta(k+1)-delta(k) <= d_delta_max]; 
    end

    for k = 1:N-1
        d = sqrt((ref.x(k+1)-ref.x(k))^2 + (ref.y(k+1)-ref.y(k))^2);
        % Linearized Bicycle Model
        cons = [cons, ((x(k+1)-x(k))*(-sin(ref.th(k))) + (y(k+1)-y(k))*cos(ref.th(k))) == d * dtheta(k)];
        cons = [cons, theta(k+1) == theta(k) + (d/L)*delta(k)];
        cons = [cons, kappa(k) == (theta(k+1)-theta(k))/d];
    end
    
    % Objective
    obj = w.len * sum(diff(x).^2 + diff(y).^2) + ...
          w.smooth * sum(diff(delta).^2) + ...
          sum(w_vec .* (kappa.^2)) + ...
          1e9*(slack'*slack);
          
    optimize(cons, obj, ops);
    
    new_path.x = value(x); 
    new_path.y = value(y); 
    new_path.th = value(theta);
end

function [v_prof, t_lap] = solve_speed_profile(path, lim_stat, coeff_aero, p, ops, mode)
    N = length(path.x);
    
    % Geometry
    dx = diff(path.x); dy = diff(path.y); 
    ds = sqrt(dx.^2+dy.^2); ds = [ds; ds(end)];
    
    dth = diff(path.th); 
    dth(dth>pi) = dth(dth>pi)-2*pi; 
    dth(dth<-pi) = dth(dth<-pi)+2*pi;
    
    k_opt = abs([dth; 0]) ./ (ds + 1e-3);
    
    % Variables (E = v^2)
    E = sdpvar(N,1);
    
    cons = [0 <= E <= p.v_max^2, E(1) == E(end)];
    
    for k = 1:N-1
        lat_acc = E(k) * k_opt(k);
        lon_acc = (E(k+1)-E(k)) / (2*ds(k));
        
        grip = lim_stat + coeff_aero * E(k);
        
        if strcmp(mode, 'diamond')
            cons = [cons, norm([lat_acc; lon_acc], 1) <= grip];
        else
             % SOCP
            cons = [cons, norm([lat_acc; lon_acc], 2) <= grip];
        end
        
        cons = [cons, lat_acc <= 5.5 * 9.81]; % Max Lateral G
        cons = [cons, lon_acc <= p.a_acc_max];
    end
    
    optimize(cons, -sum(E), ops);
    
    v_prof = sqrt(value(E));
    v_avg = 0.5 * (v_prof(1:end-1) + v_prof(2:end)) + 1e-3;
    t_lap = sum(ds(1:end-1) ./ v_avg);
end
