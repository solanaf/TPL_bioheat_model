%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
% 1D hyperbolic bioheat w/ Gaussian tumor source, single k=0.3, finer mesh
% and reduced dt for stability; snapshots at 2.5, 5, 7.5, 10 s in one plot.
%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

clear; clc;

%% 1) CONSTANTS AND PARAMETERS
rho    = 1000;    % [kg/m^3]  Tissue density (uniform)
c      = 4000;    % [J/(kg·°C)]  Tissue specific heat

k      = 0.3;     % [W/(m·°C)]  Thermal conductivity (single value)
k_star = 0.1;     % [W/(m·°C·s)]  Hyperbolic‐term coefficient

h      = 4.5;     % [W/(m^2·°C)]  Convective coefficient at boundaries
wb     = 0.0098;  % [1/s]           Blood perfusion
rho_b  = 1056;    % [kg/m^3]        Blood density
cb     = 4000;    % [J/(kg·°C)]     Blood specific heat

Qm0    = 50.65;   % [W/m^3]         Metabolic heat (uniform)
Tb     = 37;      % [°C]            Arterial blood temperature
Tl     = 37;      % [°C]            Ambient/tissue reference
T0     = 37;      % [°C]            Initial tissue temperature

% ------------ Gaussian‐source parameters ------------
% Q_r(i) = rho * S * P * exp( -a0*( x(i) - x_ast )^2 )
S      = 1;       % “per‐kg” scaling factor
P      = 30;      % “power” factor (tune as needed)
a0     = 1e5;     % [1/m^2]        Gaussian width control
x_ast  = 0.025;   % [m]            Tumor center (midpoint of 0.05 m)

Lx     = 0.05;    % [m]   Domain length (50 mm)
dx     = 0.0025;  % [m]   Spatial step → now 21 nodes from 0 to 0.05
dt     = 0.01;    % [s]   Time step (smaller for stability)
time_steps = round(10/dt) + 1; % ensures we reach exactly 10 s: (n−1)*dt = 10

% Hyperbolic relaxation times
tau_q = 600;   % [s]
tau_T = 300;   % [s]
tau_v = 100;   % [s]

%% 2) DISCRETIZE 1D DOMAIN
x  = 0:dx:Lx;          % x = [0, 0.0025, 0.005, …, 0.05]
nx = numel(x);

% Find index of tumor center (closest grid point to x_ast)
[~, iTumor] = min(abs(x - x_ast));

% Snapshot times (in seconds) and their corresponding loop‐indices
snapshot_times = [2.5, 5.0, 7.5, 10.0];
snap_idx = round(snapshot_times/dt) + 1;
% e.g., 2.5/0.01 = 250 → +1 = 251 ⇒ t = (251−1)*0.01 = 2.50 s

%% 3) INITIALIZE STORAGE FOR “SNAPSHOTS”
T_snapshots = zeros(nx, numel(snapshot_times));  
% Each column j holds T(x) at t = snapshot_times(j)

% Initialize temperature fields at t=0
T     = ones(nx,1) * T0;  
T_new = T;   
dTdt   = zeros(nx,1);
d2Tdt2 = zeros(nx,1);

%% 4) TIME‐MARCHING LOOP
for t = 1:time_steps
    % (a) Uniform volumetric terms at each timestep:
    Qm = Qm0 * (1 + (Tl - T0)/10);         % [W/m^3] metabolic
    Qb = wb * rho_b * cb * (Tb - Tl);      % [W/m^3] perfusion (constant)

    % (b) Update interior nodes i=2:(nx−1)
    for i = 2:(nx-1)
        % Gaussian tumor heat source at node i:
        Qr = rho * S * P * exp( -a0 * (x(i) - x_ast)^2 );  

        % 2nd‐spatial derivative (finite difference)
        d2Tdx2 = (T(i+1) - 2*T(i) + T(i-1)) / dx^2;

        % Hyperbolic bioheat terms:
        dTdt(i)   = (k*d2Tdx2 + Qb + Qm + Qr) / (rho * c);
        d2Tdt2(i) = (k_star*d2Tdx2) / (rho * c);

        % Update temperature at node i (eq. 7):
        T_new(i) = T(i) + dt * ( ...
           dTdt(i) + tau_q*dTdt(i) ...
         - tau_T*d2Tdt2(i) ...
         + (k + k_star*tau_v)*dTdt(i) );
    end

    % (c) Convective (Robin) BC at left boundary (i = 1):
    T_new(1) = (h*dx*Tl + k * T_new(2)) / (h*dx + k);

    % (d) Convective (Robin) BC at right boundary (i = nx):
    T_new(nx) = (h*dx*Tl + k * T_new(nx-1)) / (h*dx + k);

    % (e) Advance to next timestep
    T = T_new;

    % (f) If current “t” matches any snapshot index, store T(x)
    idx_this = find(snap_idx == t);
    if ~isempty(idx_this)
        T_snapshots(:, idx_this) = T;
    end
end

%% 5) PLOT ALL FOUR PROFILES IN A SINGLE FIGURE
figure('Position',[200,200,800,500]);
hold on;

colors = lines(numel(snapshot_times));  
for j = 1:numel(snapshot_times)
    plot(x*1000, T_snapshots(:,j), 'Color', colors(j,:), 'LineWidth',1.5, ...
         'DisplayName', sprintf('t = %.1f s', snapshot_times(j)));
end

xlabel('Distance (mm)');
ylabel('Temperature (°C)');
title('Temperature Profile T(x) at 2.5, 5, 7.5, 10 s (k = 0.3)');
legend('Location','best');
grid on;
hold off;
