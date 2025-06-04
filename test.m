%% THREE-LAYER 1D TPL BIOHEAT MODEL WITH HEAT SOURCES AND INTERFACE BCs

clear; clc;

%% 1. PHYSICAL PARAMETERS AND LAYER PROPERTIES

% Skin layer thicknesses (m)
L_epi = 0.0015;   % Epidermis (1.5 mm)
L_derm = 0.0035;  % Dermis (3.5 mm)
L_subq = 0.01;    % Subcutaneous (10 mm)
L_total = L_epi + L_derm + L_subq;

% Spatial discretization
Nx = 300; dx = L_total/(Nx-1);
x = linspace(0, L_total, Nx);

% Assign layer indices
layer = zeros(1, Nx);
layer(x <= L_epi) = 1;                                  % Epidermis
layer(x > L_epi & x <= L_epi + L_derm) = 2;             % Dermis
layer(x > L_epi + L_derm) = 3;                          % Subcutaneous

% Layers' Interface locations
idx1 = find(layer(1:end-1)==1 & layer(2:end)==2, 1, 'last')+1;
idx2 = find(layer(1:end-1)==2 & layer(2:end)==3, 1, 'last')+1;

% TPL and material parameters (per layer)
rho     = [1150 1116 900];         % kg/m^3
c       = [3590 3300 2500];         % J/kg°C
k       = [0.2 0.45 0.3];           % W/m°C
k_star  = [0.1 0.1 0.1];         % W/m°C/s (dynamic, small)
tau_q   = [28 28 28];            % s
tau_T   = [2 2 20];            % s
tau_v   = [30 30 30];            % s

% Metabolic heat (baseline, W/m^3)
Qm0 = [50.65 50.65 50.65];

% Blood perfusion (1/s)
wb = [0 0.008 0.008];
rho_b = 1056; cb = 4000; Tb = 37; % blood

% Water vaporization and diffusion (all constants)
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
nabla_r2 = (L_total/3)^2;

% External heat source (tumor device)
P = 30; % Watts
S = 12.5; % 1/kg
a0 = -127; % (1/m) tune a0 for spot size
% Qr0 will get calculated on a per-layer basis
x_star = 0.008; % m (tumor location)

% Simulation settings
T0 = [34 36 37];    % °C, initial
% dTdt0 = 0; % IC
% d2Tdt20 = 0; % IC
Nt = 2000;  % time steps
dt = 0.01;  % s
time = (0:Nt-1)*dt;

%% 2. INITIALIZATION

T = ones(Nx,1);          % initial temp
T(1:idx1-1) = T(idx1-1)*T0(1);
T(idx1:idx2-1) = T(idx1:idx2-1)*T0(2);
T(idx2:end) = T(idx2:end)*T0(3);
T_old = T;               % for second time-level
T_new = T;               % updated at each step

% For 2nd order time derivative
T_history = zeros(Nx, 3); % [n+1, n, n-1]
T_history(:,1) = T;       % n+1
T_history(:,2) = T;       % n
T_history(:,3) = T;       % n-1

tumor_idx = find(abs(x-x_star) == min(abs(x-x_star)),1); % nearest grid to tumor
T_tumor = zeros(1, Nt);

%% 3. MAIN TIME LOOP

for t = 1:2
    % --- Heat sources (at current T) ---
    Qm = zeros(Nx,1); Qd = zeros(Nx,1); Qb = zeros(Nx,1); Qv = zeros(Nx,1); Qr = zeros(Nx,1);

    for i = 1:Nx
        L = layer(i);

        %% CALCULATE HEAT SOURCES %%

        % Metabolic
        Qm(i) = Qm0(L) * (1 + (T(i) - T0(L))/10);

        % Water diffusion
        Qd(i) = (Df(L) * cw * (rho_s - rho_c) / nabla_r2) * (T(i) - T0(L));

        % Blood perfusion (dermis, subq)
        if L > 1
            Qb(i) = wb(L) * rho_b * cb * (Tb - T(i));
        end

        % Water vaporization (epidermis only)
        if L == 1
            Delta_m = (Da*Mw/(Ra*Tw)) * (Pw/Tw) * RH/(delta*c_air);
            Delta_Hvap = 2400e3; % J/kg
            Qv(i) = Delta_m * Delta_Hvap / (delta * c_air);
        end

        % External device (all layers, spatial Gaussian)
        Qr0 = rho(L)*S*P;
        Qr(i) = Qr0 * exp(-a0 * (x(i)-x_star)^2);

        Q_total = Qm + Qd + Qb + Qv + Qr;

        %% FINITE DIFFERENCE METHODS - STEP UPDATE %%
            % --- Update T_new using explicit Euler (illustrative) with TPL memory (2nd order) ---
            % For pedagogical clarity, below is a simplified TPL scheme using 2-level time stepping:
        L = layer(i);
        rho_i = rho(L); c_i = c(L); k_i = k(L);
        tq = tau_q(L); tT = tau_T(L); tv = tau_v(L);

        % Discrete spatial 2nd derivative
        if i == 1
            d2Tdx2 = (T(i+1) - 2*T(i) + 22) / dx^2; %% Assuming room temp outside of skin - BC
        elseif i == Nx
            d2Tdx2 = (37 - 2*T(i) + 22) / dx^2; %% Assuming body temp at end of tissue - BC
        else
            d2Tdx2 = (T(i+1) - 2*T(i) + T(i-1)) / dx^2;
        end

        % Temporal 1st and 2nd derivatives
        dTdt = (k(L) * (d2Tdx2) + Q_total(i)) / (rho(L) * c(L));
        d2Tdt2 = (k_star(L) * (d2Tdx2)) / (rho(L) * c(L));

        % % Simple TPL finite difference (needs refinement for strict TPL, see note)
        % % Central difference in time for 2nd order
        % if t > 2
        %     Tnm1 = T_history(i,3); % n-1
        % else
        %     Tnm1 = T0(L);
        % end

        % Main update (simplified; true TPL requires implicit or staggered multi-step)
        T_new(i) = T(i) + dt * (dTdt + tau_q(L) .* dTdt - tau_T(L) .* d2Tdt2 + (k(L) + k_star(L) .* tau_v(L)) .* dTdt);

    end

    % --- INTERFACE ENFORCEMENT ---

    % Epidermis-dermis
    kL = k(layer(idx1-1)); kR = k(layer(idx1));
    T_new(idx1) = (kL * T_new(idx1-1) + kR * T_new(idx1+1)) / (kL + kR);
    T_new(idx1-1) = T_new(idx1);
    % Dermis-subcutaneous
    kL = k(layer(idx2-1)); kR = k(layer(idx2));
    T_new(idx2) = (kL * T_new(idx2-1) + kR * T_new(idx2+1)) / (kL + kR);
    T_new(idx2-1) = T_new(idx2);

    % --- Store T at tumor location ---
    T_tumor(t) = T_new(tumor_idx);

    % --- Shift time levels ---
    T_history(:,3) = T_history(:,2); % n-1 <= n
    T_history(:,2) = T;              % n <= T^n
    T_history(:,1) = T_new;          % n+1 <= T^{n+1}

    % --- Advance solution ---
    T = T_new;

    % (Optional: plot temperature profiles in real time)
    % if mod(t,200)==0
    %     plot(x*1000,T,'LineWidth',2); drawnow;
    % end
end

%% 4. PLOT: TUMOR TEMPERATURE VS TIME

figure;
plot(time, T_tumor, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Temperature at Tumor (°C)');
title('Temperature at Tumor Location vs. Time');
grid on;