%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
% 1D hyperbolic bioheat w/ Gaussian tumor source, 3 skin layers, finer mesh
% and reduced dt for stability; snapshots at 2.5, 5, 7.5, 10 s in one plot.
%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

clear; clc;

%% 1) DISCRETIZE TIME & SPATIAL DOMAINS
% Skin Layer thickness
L_epi = 0.0015;   % Epidermis (1.5 mm)
L_derm = 0.0035;  % Dermis (3.5 mm)
L_subq = 0.01;    % Subcutaneous (10 mm)
Lx = L_epi + L_derm + L_subq; % Domain Length

x_ast  = 0.008;   % [m]   Tumor center (midpoint of 0.05 m)

dx     = 0.0005;  % [m]   Spatial step → now 21 nodes from 0 to 0.05
dt     = 0.015;    % [s]   Time step (smaller for stability)
max_time = 10;   % plot up to 10s
time_steps = round(max_time/dt) + 1; % ensures we reach exactly 10 s: (n−1)*dt = 10
time = 0:dt:max_time+dt;

x  = 0:dx:Lx;          % x = [0, 0.0025, 0.005, …, 0.05]
nx = numel(x);

% Define Layers
layer = zeros(1, nx);
layer(x <= L_epi) = 1;                                  % Epidermis
layer(x > L_epi & x <= L_epi + L_derm) = 2;             % Dermis
layer(x > L_epi + L_derm) = 3;                          % Subcutaneous

% Find index of tumor center (closest grid point to x_ast)
[~, iTumor] = min(abs(x - x_ast));

%% 2) CONSTANTS AND PARAMETERS
% ------------ Physical Parameters ------------
% Skin Layer params
rho    = [1150 1116 900];   % [kg/m^3]  Tissue density (layers)
c      = [3590 3300 2500];  % [J/(kg·°C)]  Tissue specific heat
k      = [0.2 0.45 0.3];    % [W/(m·°C)]  Thermal conductivity
k_star = [0.1 0.1 0.1];     % [W/(m·°C·s)]  Hyperbolic‐term coefficient

% Blood params
h      = 4.5;     % [W/(m^2·°C)]    Convective coefficient at boundaries
wb     = 0.0098;  % [1/s]           Blood perfusion
rho_b  = 1056;    % [kg/m^3]        Blood density
cb     = 4000;    % [J/(kg·°C)]     Blood specific heat

Qm0    = 50.65;   % [W/m^3]         Metabolic heat (uniform)
Tb     = 37;      % [°C]            Arterial blood temperature
Tl     = [37 37 37];      % [°C]    Ambient/tissue reference
T0     = 37;      % [°C]            Initial tissue temperature

% ------------ Gaussian‐source parameters ------------
% Q_r(i) = rho * S * P * exp( -a0*( x(i) - x_ast )^2 )
S      = 15;       % “per‐kg” scaling factor
P      = 30;      % “power” factor (tune as needed)
a0     = 1e6;   % [1/m]  Gaussian width control

% ------------ Water vaporization and diffusion ------------
Da = 2.5e-5;    % m^2/s (air)
Mw = 0.018;     % kg/mol
Ra = 8.314;     % J/mol·K
Tw = 306;       % K (~33°C)
Pw = 5600;      % Pa (sat. vapor pressure at 33°C)
RH = 0.5;       % Relative humidity (fraction)
delta = 1e-4;   % m
c_air = 1005;   % J/kg·K
Df = [2e-9 2e-9 2e-9];  % m^2/s
cw = 4180;      % J/kg°C
rho_s = 1100;   % kg/m^3
rho_c = 1000;   % kg/m^3
nabla_r2 = (Lx/3)^2;

% ---------- Hyperbolic relaxation times ---------------
tau_q = 10;   % [s]
tau_T = 2;   % [s]
tau_v = 1;   % [s]

%% 3) INITIALIZE STORAGE FOR “SNAPSHOTS”
% Snapshot times (in seconds) and their corresponding loop‐indices
snapshot_times = linspace(2,max_time,5); % [2.5, 5.0, 7.5, 10.0];
snap_idx = round(snapshot_times/dt) + 1;
% e.g., 2.5/0.01 = 250 → +1 = 251 ⇒ t = (251−1)*0.01 = 2.50 s

for P = 10:5:50
    T_snapshots = zeros(nx, numel(snapshot_times));  
    % Each column j holds T(x) at t = snapshot_times(j)
    
    % Tumor temp array
    T_tumor = zeros(1, time_steps);
    
    % --------- Initialize temperature fields at t=0 -------------
    T     = ones(nx,1) * T0;  
    T_new = T;   
    dTdt   = zeros(nx,1);
    d2Tdt2 = zeros(nx,1);


    %% 4) TIME‐MARCHING LOOP
    for n = 1:time_steps
    
        % (a) Update interior nodes i=2:(nx−1)
        for i = 2:(nx-1)
            L = layer(i);
    
            % (b) Calculate Heat Source values for given spatial step
            % Metabolic Heat Source:
            Qm = Qm0 * (1 + (Tl(L) - T0)/10);         % [W/m^3] metabolic
    
            % Water diffusion:
            Qd = (Df(L) * cw * (rho_s - rho_c) / nabla_r2) * (T(i) - Tl(L));
            
            % Blood Perfusion (dermis, subq):
            if L > 1
                Qb = wb * rho_b * cb * (Tb - Tl(L));      % [W/m^3] perfusion (constant)
            else
                Qb = 0;
            end
    
            % Water vaporization (epidermis only)
            if L == 1
                Delta_m = (Da*Mw/(Ra*Tw)) * (Pw/Tw) * RH/(delta*c_air);
                Delta_Hvap = 2400e3; % J/kg
                Qv = Delta_m * Delta_Hvap / (delta * c_air);
            else 
                Qv = 0;
            end
    
            % Gaussian tumor heat source at node i:
            Qr = rho(L) * S * P * exp( -a0 * (x(i) - x_ast)^2 );  
    
            % Add up all heat sources to consolidate
            Q_tot = Qm + Qd + Qb + Qv + Qr;
    
            % 2nd‐spatial derivative (finite difference)
            d2Tdx2 = (T(i+1) - 2*T(i) + T(i-1)) / dx^2;
    
            % Hyperbolic bioheat terms:
            dTdt(i)   = (k(L)*d2Tdx2 + Q_tot) / (rho(L) * c(L));
            d2Tdt2(i) = (k_star(L)*d2Tdx2) / (rho(L) * c(L));
    
            % Update temperature at node i (eq. 7):
            T_new(i) = T(i) + dt * ( ...
               dTdt(i) + tau_q*dTdt(i) ...
             - tau_T*d2Tdt2(i) ...
             + (k(L) + k_star(L)*tau_v)*dTdt(i) );
        end
    
        % (c) Convective (Robin) BC at left boundary (i = 1):
        T_new(1) = (h*dx*Tl(1) + k(L) * T_new(2)) / (h*dx + k(L));
    
        % (d) Convective (Robin) BC at right boundary (i = nx):
        T_new(nx) = (h*dx*Tl(end) + k(L) * T_new(nx-1)) / (h*dx + k(L));
    
        % (e) Advance to next timestep
        T = T_new;
    
        % (f) If current “t” matches any snapshot index, store T(x)
        idx_this = find(snap_idx == n);
        if ~isempty(idx_this)
            T_snapshots(:, idx_this) = T;
        end
    
        % (g) Store T at tumor location
        T_tumor(n) = T_new(iTumor);
    end

    %% 5) PLOT ALL FOUR PROFILES IN A SINGLE FIGURE
    figure('Position',[200,200,800,500]);
    hold on;
    
    colors = lines(numel(snapshot_times));  
    for j = 1:numel(snapshot_times)
        plot(x*1000, T_snapshots(:,j), 'Color', colors(j,:), 'LineWidth',1.5, ...
             'DisplayName', sprintf('t = %.1f s', snapshot_times(j)));
    end
    
    ylim([35 45])
    xline(0, '-',{'Epiderm'},'HandleVisibility','off')
    xline((L_epi)*10^3,'-',{'Derm'},'HandleVisibility','off')
    xline((L_epi+L_derm)*10^3,'-',{'Subcutaneous'},'HandleVisibility','off')
    xline(x_ast*10^3,'-r',{'Tumor'},'LineWidth',10,'HandleVisibility','off')
    xlabel('Distance (mm)');
    ylabel('Temperature (°C)');
    title(['Temperature Profile T(x); P =',num2str(P)]);
    legend('Location','best');
    grid on;
    hold off;

end 
%% 6. PLOT: TUMOR TEMPERATURE VS TIME
figure;
plot(time(1:length(T_tumor)), T_tumor, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Temperature at Tumor (°C)');
title('Temperature at Tumor Location vs. Time');
grid on;