%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
% 1D hyperbolic bioheat w/ Gaussian tumor source, 3 skin layers, finer mesh
% Computes T(x) at t=10 s for three different tau_T / tau_v / tau_q values.
% Then plots each triplet in a 1×3 subplot.
%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

clear; clc;

%% 1) DISCRETIZE TIME & SPATIAL DOMAINS
% Skin Layer thicknesses
L_epi  = 0.0015;   % Epidermis (1.5 mm)
L_derm = 0.0035;   % Dermis   (3.5 mm)
L_subq = 0.01;     % Subcutaneous (10 mm)
Lx     = L_epi + L_derm + L_subq;  % Total slab thickness (m)

x_ast  = 0.008;    % [m] Tumor center (8 mm from surface)

dx         = 0.0005;         % [m] spatial step →  (Lx/dx + 1) nodes
dt         = 0.015;          % [s] time step (for stability)
max_time   = 10;             % [s] we want solution at t=10 s
time_steps = round(max_time/dt) + 1;  % ensures (time_steps–1)*dt ≈ 10
time       = (0:(time_steps-1)) * dt;

x     = 0:dx:Lx;    % grid from 0 to Lx in steps of dx
nx    = numel(x);

% Build a “layer” index for each x(i):
layer = zeros(1,nx);
layer(x <= L_epi) = 1;                               % Epidermis
layer(x > L_epi & x <= L_epi + L_derm) = 2;          % Dermis
layer(x > L_epi + L_derm) = 3;                       % Subcutaneous

% Find the grid index closest to x_ast for tumor‐location tracking:
[~, iTumor] = min(abs(x - x_ast));

%% 2) CONSTANTS AND PARAMETERS (unchanged except τ’s)
% Physical parameters per layer:
rho    = [1150, 1116, 900];    % [kg/m^3] densities for [epi, derm, subq]
c      = [3590, 3300, 2500];   % [J/(kg·°C)] specific heats
k      = [0.2, 0.45, 0.3];     % [W/(m·°C)] thermal conductivities
k_star = [0.1, 0.1, 0.1];       % [W/(m·°C·s)] hyperbolic‐term coefficient

% Blood and convection:
h     = 4.5;       % [W/(m^2·°C)] convective‐BC coefficient
wb    = 0.0098;    % [1/s] blood perfusion rate
rho_b = 1056;      % [kg/m^3] blood density
cb    = 4000;      % [J/(kg·°C)] blood specific heat

Qm0 = 50.65;       % [W/m^3] metabolic heat
Tb  = 37;          % [°C]   blood inlet temperature
Tl  = [37, 37, 37];% [°C]   ambient/tissue reference for each layer
T0  = 37;          % [°C]   initial tissue temperature

% Gaussian‐source parameters (same for all runs):
S     = 15;        % [unitless] “per‐kg” scaling
P     = 35;        % [unitless] “power” factor
a0    = 1e6;       % [1/m^2] Gaussian width control
                   % Q_r(i) = rho(layer(i)) * S * P * exp(-a0*(x(i)-x_ast)^2)

%% 3) PREALLOCATE “SNAPSHOT” MATRICES
% Each final‐time profile has 3 columns, one for each tau‐value in that set:
T_tauT_snap = zeros(nx, 3);  
T_tauV_snap = zeros(nx, 3);
T_tauQ_snap = zeros(nx, 3);

% We will hold the “default” values for the taus not being varied:
tauT_default = 20;    % [s]
tauV_default = 1;    % [s]
tauQ_default = 10;   % [s]

%% 4) VARY τ_T ∈ {2,3,4}  (fix τ_v=1, τ_q=10)
tauT_list = [1, 10, 20]; % values over 25 begin to look unstable
for j = 1:3
    tau_T = tauT_list(j);
    tau_v = tauV_default;
    tau_q = tauQ_default;
    
    % Initialize temperature fields at t=0:
    T     = ones(nx,1) * T0;
    T_new = T;
    
    % Time‐marching loop to t=10 s:
    for n = 1:time_steps
        for i = 2:(nx-1)
            L = layer(i);
            
            % 1) Metabolic heat (layer-dependent Tl):
            Qm = Qm0 * (1 + (T(i) - Tl(L))/10);
            % 2) Water diffusion (not varying here, same as original):
            Da     = 2.5e-5;  Mw = 0.018;  Ra  = 8.314;
            Tw     = 306;     Pw = 5600;  RH  = 0.5;
            delta  = 1e-4;    c_air = 1005;
            Df     = [2e-9, 2e-9, 2e-9];  cw = 4180;
            rho_s  = 1100;    rho_c = 1000; nabla_r2 = (Lx/3)^2;
            Qd     = (Df(L)*cw*(rho_s-rho_c)/nabla_r2)*(T(i)-Tl(L));
            
            % 3) Perfusion (dermis+subq only):
            if L > 1
                Qb = wb * rho_b * cb * (Tb - Tl(L));
            else
                Qb = 0;
            end
            
            % 4) Vaporization (epidermis only):
            if L == 1
                Delta_m    = (Da*Mw/(Ra*Tw)) * (Pw/Tw) * RH/(delta*c_air);
                Delta_Hvap = 2400e3; % J/kg
                Qv         = Delta_m * Delta_Hvap / (delta*c_air);
            else
                Qv = 0;
            end
            
            % 5) Gaussian tumor source:
            Qr = rho(L) * S * P * exp(-a0*(x(i)-x_ast)^2);
            
            % 6) Sum all sources:
            Q_tot = Qm + Qd + Qb + Qv + Qr;
            
            % 7) Second‐derivative in x:
            d2Tdx2 = (T(i+1) - 2*T(i) + T(i-1)) / dx^2;
            
            % 8) Hyperbolic bioheat:
            dTdt(i)   = (k(L)*d2Tdx2 + Q_tot) / (rho(L)*c(L));
            d2Tdt2(i) = (k_star(L)*d2Tdx2) / (rho(L)*c(L));
            
            T_new(i) = T(i) + dt * ( ...
               dTdt(i) + tau_q*dTdt(i) ...
             - tau_T*d2Tdt2(i) ...
             + (k(L) + k_star(L)*tau_v) * dTdt(i) );
        end
        
        % Convective BC at i=1 (left):
        L = layer(1);
        T_new(1)  = (h*dx*Tl(L) + k(L)*T_new(2)) / (h*dx + k(L));
        % Convective BC at i=nx (right):
        L = layer(nx);
        T_new(nx) = (h*dx*Tl(L) + k(L)*T_new(nx-1)) / (h*dx + k(L));
        
        % Advance:
        T = T_new;
    end
    
    % After 10 s, store the final profile in column j:
    T_tauT_snap(:,j) = T;
end

%% 5) VARY τ_v ∈ {1,2,3}  (fix τ_T=2, τ_q=10)
tauV_list = [1, 10, 20]; % minimizing works best
for j = 1:3
    tau_T = tauT_default;
    tau_v = tauV_list(j);
    tau_q = tauQ_default;
    
    % Initialize:
    T     = ones(nx,1) * T0;
    T_new = T;
    
    for n = 1:time_steps
        for i = 2:(nx-1)
            L = layer(i);
            
            % (Same source calculations as above)
            Qm = Qm0 * (1 + (T(i) - Tl(L))/10);
            Da     = 2.5e-5;  Mw = 0.018;  Ra  = 8.314;
            Tw     = 306;     Pw = 5600;  RH  = 0.5;
            delta  = 1e-4;    c_air = 1005;
            Df     = [2e-9,2e-9,2e-9];  cw = 4180;
            rho_s  = 1100;   rho_c = 1000; nabla_r2 = (Lx/3)^2;
            Qd     = (Df(L)*cw*(rho_s - rho_c)/nabla_r2)*(T(i) - Tl(L));
            
            if L > 1
                Qb = wb * rho_b * cb * (Tb - Tl(L));
            else
                Qb = 0;
            end
            
            if L == 1
                Delta_m    = (Da*Mw/(Ra*Tw)) * (Pw/Tw) * RH/(delta*c_air);
                Delta_Hvap = 2400e3; 
                Qv         = Delta_m * Delta_Hvap / (delta*c_air);
            else
                Qv = 0;
            end
            
            Qr = rho(L) * S * P * exp(-a0*(x(i) - x_ast)^2);
            Q_tot = Qm + Qd + Qb + Qv + Qr;
            
            d2Tdx2 = (T(i+1) - 2*T(i) + T(i-1)) / dx^2;
            dTdt(i)   = (k(L)*d2Tdx2 + Q_tot) / (rho(L)*c(L));
            d2Tdt2(i) = (k_star(L)*d2Tdx2) / (rho(L)*c(L));
            
            T_new(i) = T(i) + dt * ( ...
               dTdt(i) + tau_q*dTdt(i) ...
             - tau_T*d2Tdt2(i) ...
             + (k(L) + k_star(L)*tau_v) * dTdt(i) );
        end
        
        % BCs at i=1 and i=nx:
        L = layer(1);
        T_new(1)  = (h*dx*Tl(L) + k(L)*T_new(2))   / (h*dx + k(L));
        L = layer(nx);
        T_new(nx) = (h*dx*Tl(L) + k(L)*T_new(nx-1)) / (h*dx + k(L));
        
        T = T_new;
    end
    
    T_tauV_snap(:,j) = T;
end

%% 6) VARY τ_q ∈ {5,10,15}  (fix τ_T=2, τ_v=1)
tauQ_list = [7, 10, 15]; % below 8, unstable -- above 15, localization decreases
for j = 1:3
    tau_T = tauT_default;
    tau_v = tauV_default;
    tau_q = tauQ_list(j);
    
    T     = ones(nx,1) * T0;
    T_new = T;
    
    for n = 1:time_steps
        for i = 2:(nx-1)
            L = layer(i);
            
            Qm = Qm0 * (1 + (T(i) - Tl(L))/10);
            Da     = 2.5e-5;  Mw = 0.018;  Ra  = 8.314;
            Tw     = 306;     Pw = 5600;  RH  = 0.5;
            delta  = 1e-4;    c_air = 1005;
            Df     = [2e-9,2e-9,2e-9];  cw = 4180;
            rho_s  = 1100;    rho_c = 1000; nabla_r2 = (Lx/3)^2;
            Qd     = (Df(L)*cw*(rho_s - rho_c)/nabla_r2)*(T(i)-Tl(L));
            
            if L > 1
                Qb = wb * rho_b * cb * (Tb - Tl(L));
            else
                Qb = 0;
            end
            
            if L == 1
                Delta_m    = (Da*Mw/(Ra*Tw)) * (Pw/Tw) * RH/(delta*c_air);
                Delta_Hvap = 2400e3;
                Qv         = Delta_m * Delta_Hvap / (delta*c_air);
            else
                Qv = 0;
            end
            
            Qr = rho(L) * S * P * exp(-a0*(x(i)-x_ast)^2);
            Q_tot = Qm + Qd + Qb + Qv + Qr;
            
            d2Tdx2 = (T(i+1) - 2*T(i) + T(i-1)) / dx^2;
            dTdt(i)   = (k(L)*d2Tdx2 + Q_tot) / (rho(L)*c(L));
            d2Tdt2(i) = (k_star(L)*d2Tdx2) / (rho(L)*c(L));
            
            T_new(i) = T(i) + dt * ( ...
               dTdt(i) + tau_q*dTdt(i) ...
             - tau_T*d2Tdt2(i) ...
             + (k(L) + k_star(L)*tau_v) * dTdt(i) );
        end
        
        % BCs:
        L = layer(1);
        T_new(1)  = (h*dx*Tl(L) + k(L)*T_new(2))   / (h*dx + k(L));
        L = layer(nx);
        T_new(nx) = (h*dx*Tl(L) + k(L)*T_new(nx-1)) / (h*dx + k(L));
        
        T = T_new;
    end
    
    T_tauQ_snap(:,j) = T;
end

%% 7) PLOT ALL THREE SETS IN A SINGLE 1×3 SUBPLOT FIGURE
figure('Position',[200,200,1200,400]);

% 7.1) Subplot #1: varying tau_T
subplot(1,3,1);
plot(x*1000, T_tauT_snap(:,1), 'b', 'LineWidth',1.5);  hold on;
plot(x*1000, T_tauT_snap(:,2), 'r','LineWidth',1.5);
plot(x*1000, T_tauT_snap(:,3), 'g','LineWidth',1.5);
ylim([36 50]);
xlabel('Distance (mm)');
ylabel('Temperature (°C)');
title('Temperature profile T (t = 10 s), varying \tau_T');
legend('\tau_T = 1','\tau_T = 10','\tau_T = 20','Location','best');
grid on;
hold off;

% 7.2) Subplot #2: varying tau_v
subplot(1,3,2);
plot(x*1000, T_tauV_snap(:,1), 'b', 'LineWidth',1.5);  hold on;
plot(x*1000, T_tauV_snap(:,2), 'r','LineWidth',1.5);
plot(x*1000, T_tauV_snap(:,3), 'g','LineWidth',1.5);
ylim([36 50]);
xlabel('Distance (mm)');
title('Temperature profile T (t = 10 s), varying \tau_v');
legend('\tau_v = 1','\tau_v = 10','\tau_v = 20','Location','best');
grid on;
hold off;

% 7.3) Subplot #3: varying tau_q
subplot(1,3,3);
plot(x*1000, T_tauQ_snap(:,1), 'b', 'LineWidth',1.5);  hold on;
plot(x*1000, T_tauQ_snap(:,2), 'r','LineWidth',1.5);
plot(x*1000, T_tauQ_snap(:,3), 'g','LineWidth',1.5);
ylim([36 50]);
xlabel('Distance (mm)');
title('Temperature profile T (t = 10 s), varying \tau_q');
legend('\tau_q = 7','\tau_q = 10','\tau_q = 15','Location','best');
grid on;
hold off;
