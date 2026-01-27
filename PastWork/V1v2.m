%% F1 TRAJECTORY OPTIMIZATION - V1 (GEOMETRIC SHORTEST PATH)
% Ce script calcule la trajectoire la plus courte possible (Géométrique)
% Il utilise le moteur robuste de la V3 pour garantir la comparaison.

clear; clc; close all;

%% 1. SETUP & DATA
fprintf('=== 1. SETUP V1 (Shortest Path) ===\n');

try
    track_data = readmatrix('Circuits_Data/Nuerburgring_track.csv');
catch
    error('Fichier track introuvable. Vérifie le chemin.');
end

x_center = track_data(:,1);
y_center = track_data(:,2);
w_right  = track_data(:,3);
w_left   = track_data(:,4);

% Nettoyage des données (comme dans V3 pour la cohérence)
dist_sq = [1; (diff(x_center).^2 + diff(y_center).^2)];
keep_idx = dist_sq > 0.01;
x_center = x_center(keep_idx);
y_center = y_center(keep_idx);
w_right  = w_right(keep_idx);
w_left   = w_left(keep_idx);

N = length(x_center);

% Paramètres Véhicule
L = 2.5;
delta_max = deg2rad(35); % Augmenté pour permettre les épingles
delta_min = -deg2rad(35);

%% 2. PRE-COMPUTE TRACK
% Lissage (Spline) pour aider le solveur
p_smooth = 0.999; 
x_ref = csaps(1:N, x_center, p_smooth, 1:N)';
y_ref = csaps(1:N, y_center, p_smooth, 1:N)';

theta_track = zeros(N,1);
for k = 2:N-1
    theta_track(k) = atan2(y_ref(k+1)-y_ref(k-1), x_ref(k+1)-x_ref(k-1));
end
theta_track(1) = atan2(y_ref(2)-y_ref(N), x_ref(2)-x_ref(N));
theta_track(N) = atan2(y_ref(1)-y_ref(N-1), x_ref(1)-x_ref(N-1));
theta_track = unwrap(theta_track);

% Limites gauche/droite
x_L = x_center - w_left .* sin(theta_track);
y_L = y_center + w_left .* cos(theta_track);
x_R = x_center + w_right .* sin(theta_track);
y_R = y_center - w_right .* cos(theta_track);

ref.th = theta_track;

%% 3. OPTIMIZATION
fprintf('=== 2. SOLVING GEOMETRIC PATH ===\n');

% Variables
x = sdpvar(N,1); 
y = sdpvar(N,1); 
theta = sdpvar(N,1); 
delta = sdpvar(N-1,1);
z = sdpvar(N,1);        % Position latérale (0=Gauche, 1=Droite)
dtheta = sdpvar(N,1);   % Erreur de cap
slack_cycle = sdpvar(4,1); % Pour boucler le circuit proprement

% Initialisation pour aider le solveur
assign(x, x_center); assign(y, y_center); assign(theta, theta_track); assign(z, 0.5*ones(N,1));

cons = [];

% 3.1 Contraintes de Piste (Convexes)
cons = [cons, x == x_L + z.*(x_R-x_L), y == y_L + z.*(y_R-y_L), 0 <= z <= 1];

% 3.2 Contraintes Cycliques (Bouclage)
cons = [cons, x(1)-x(N)==slack_cycle(1), y(1)-y(N)==slack_cycle(2), ...
              z(1)-z(N)==slack_cycle(3), delta(1)-delta(end)==slack_cycle(4)];
cons = [cons, theta(1) - ref.th(1) == theta(N) - ref.th(N)]; % Cap relatif

% 3.3 Cinématique linéarisée (Robuste)
cons = [cons, dtheta == theta - ref.th];
cons = [cons, delta_min <= delta <= delta_max];

for k = 1:N-1
    d = sqrt((x_ref(k+1)-x_ref(k))^2 + (y_ref(k+1)-y_ref(k))^2);
    
    % Mouvement
    cons = [cons, ((x(k+1)-x(k))*(-sin(ref.th(k))) + (y(k+1)-y(k))*cos(ref.th(k))) == d * dtheta(k)];
    
    % Rotation (Empattement L)
    cons = [cons, theta(k+1) == theta(k) + (d/L)*delta(k)];
end

% 3.4 OBJECTIF SPÉCIFIQUE V1 (Shortest Path)
% C'est ici que ça change par rapport à V3
dx_path = diff(x);
dy_path = diff(y);

% POIDS V1 : On veut minimiser la longueur, on s'en fiche de la courbure (un peu)
w_length = 1000;   % Priorité absolue à la distance la plus courte
w_steer  = 1;      % Juste pour éviter les vibrations
w_cycle  = 1e9;    % Fermeture stricte

obj = w_length * sum(dx_path.^2 + dy_path.^2) + ... % Minimise la distance (Pythagore approx)
      w_steer  * sum(diff(delta).^2) + ...          % Douceur minimale
      w_cycle  * (slack_cycle'*slack_cycle);

% Résolution
ops = sdpsettings('solver','gurobi','verbose',1, 'usex0',1);
sol = optimize(cons, obj, ops);

if sol.problem ~= 0
    warning('Problème solver : %s', sol.info);
end

%% 4. EXPORT & PLOT
x_opt = value(x);
y_opt = value(y);
th_opt = value(theta);

% Sauvegarde pour le comparateur
fprintf('Sauvegarde des données V1...\n');
data_v1.x = x_opt;
data_v1.y = y_opt;
data_v1.lbl = 'V1 (Geometric - Shortest Path)';
save('DATA_V1_GEO.mat', '-struct', 'data_v1');
fprintf('Fichier DATA_V1_GEO.mat créé !\n');

% Plot rapide
figure;
plot(x_L, y_L, 'k'); hold on; plot(x_R, y_R, 'k');
plot(x_center, y_center, 'g:');
plot(x_opt, y_opt, 'b-', 'LineWidth', 2);
title('V1 Result: Shortest Path');
axis equal; grid on;