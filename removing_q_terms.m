%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
% 1D hyperbolic bioheat with Gaussian tumor source, 3 skin layers, finer mesh
% Produces three T(x) profiles at t=10 s:
%   1) FULL model (Qm + Qd + Qb + Qv + Qr)
%   2) NO Qv      (Qm + Qd + Qb     + Qr)
%   3) NO Qv,Qd   (Qm       + Qb     + Qr)
% Then plots these three curves together, with a zoomed‐in inset.
%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

clear; clc;

%% 1) DISCRETIZE TIME & SPACE
L_epi  = 0.0015;   % [m] Epidermis thickness (1.5 mm)
L_derm = 0.0035;   % [m] Dermis thickness   (3.5 mm)
L_subq = 0.01;     % [m] Subcutaneous      (10 mm)
Lx     = L_epi + L_derm + L_subq;   % [m] total thickness (15 mm)

x_ast  = 0.008;    % [m] Tumor center (8 mm from surface)

dx         = 0.0005;         % [m] spatial step
dt         = 0.015;          % [s] time step
max_time   = 10;             % [s] final time
time_steps = round(max_time/dt) + 1;
time       = (0:(time_steps-1)) * dt;

x  = 0:dx:Lx;     
nx = numel(x);

% Assign layer index: 1=epi, 2=derm, 3=subq
layer = zeros(1,nx);
layer(x <= L_epi)                          = 1;
layer(x > L_epi & x <= (L_epi + L_derm))   = 2;
layer(x > (L_epi + L_derm))                = 3;

% Find index for tumor location
[~, iTumor] = min(abs(x - x_ast));

%% 2) CONSTANTS AND PARAMETERS
% ––– Physical parameters per layer
rho    = [1150, 1116, 900];    % [kg/m^3]
c      = [3590, 3300, 2500];   % [J/(kg·°C)]
k      = [0.2, 0.45, 0.3];     % [W/(m·°C)]
k_star = [0.1, 0.1, 0.1];       % [W/(m·°C·s)]

% ––– Blood‐ and convection
h     = 4.5;       % [W/(m^2·°C)]
wb    = 0.0098;    % [1/s]
rho_b = 1056;      % [kg/m^3]
cb    = 4000;      % [J/(kg·°C)]

Qm0 = 50.65;       % [W/m^3]
Tb  = 37;          % [°C]
Tl  = [37, 37, 37];% [°C] ambient‐reference per layer
T0  = 37;          % [°C]

% ––– Gaussian heat‐source parameters (same for all runs)
S     = 15;       
P     = 35;       
a0    = 1e6;      
% Qr(i) = rho(layer(i))*S*P * exp(-a0*(x(i)-x_ast)^2)

% ––– Water diffusion / vaporization constants
Da      = 2.5e-5;    % [m^2/s]
Mw      = 0.018;     
Ra      = 8.314;     
Tw      = 306;       % [K] (~33 °C)
Pw      = 5600;      % [Pa]
RH      = 0.5;       
delta   = 1e-4;      
c_air   = 1005;      
Df      = [2e-9, 2e-9, 2e-9];  % [m^2/s] per layer
cw      = 4180;      
rho_s   = 1100;      
rho_c   = 1000;      
nabla_r2 = (Lx/3)^2;          % [m^2]

% ––– Hyperbolic relaxation times (overridden per scenario)
tauQ_default = 10;  % [s]
tauT_default = 20;  % [s]
tauV_default = 1;   % [s]

%% 3) PREALLOCATE final‐time profiles for 3 scenarios
T_full   = zeros(nx,1);   % Full model at t = 10 s
T_noQv   = zeros(nx,1);   % No Qv at t = 10 s
T_noQvQd = zeros(nx,1);   % No Qv & Qd at t = 10 s

%% 4) RUN SCENARIOS
for scenario = 1:3
    % Reset initial temperature
    T     = ones(nx,1) * T0;
    T_new = T; 
    
    % Decide which terms to zero out:
    switch scenario
        case 1
            include_Qv = true;
            include_Qd = true;
        case 2
            include_Qv = false;
            include_Qd = true;
        case 3
            include_Qv = false;
            include_Qd = false;
    end
    
    % Set tau’s to defaults (only Qv/Qd vary here)
    tau_q = tauQ_default;
    tau_T = tauT_default;
    tau_v = tauV_default;
    
    % Time‐march to t = 10 s
    for n = 1:time_steps
        for i = 2:(nx-1)
            L = layer(i);
            
            % 1) Metabolic heat
            Qm = Qm0 * (1 + (Tl(L) - T0)/10);
            
            % 2) Water diffusion (only if include_Qd == true)
            if include_Qd
                Qd = (Df(L) * cw * (rho_s - rho_c) / nabla_r2) * (T(i) - Tl(L));
            else
                Qd = 0;
            end
            
            % 3) Blood perfusion (dermis & subq)
            if L > 1
                Qb = wb * rho_b * cb * (Tb - Tl(L));
            else
                Qb = 0;
            end
            
            % 4) Vaporization (only if include_Qv == true AND L == 1)
            if include_Qv && (L == 1)
                Delta_m    = (Da*Mw/(Ra*Tw)) * (Pw/Tw) * RH/(delta*c_air);
                Delta_Hvap = 2400e3;   % J/kg
                Qv = Delta_m * Delta_Hvap / (delta * c_air);
            else
                Qv = 0;
            end
            
            % 5) Gaussian tumor source
            Qr = rho(L) * S * P * exp( -a0 * (x(i) - x_ast)^2 );
            
            % 6) Total volumetric source
            Q_tot = Qm + Qd + Qb + Qv + Qr;
            
            % 7) Second spatial derivative
            d2Tdx2 = (T(i+1) - 2*T(i) + T(i-1)) / dx^2;
            
            % 8) Hyperbolic bioheat update
            dTdt(i)   = (k(L) * d2Tdx2 + Q_tot) / (rho(L) * c(L));
            d2Tdt2(i) = (k_star(L) * d2Tdx2) / (rho(L) * c(L));
            
            T_new(i) = T(i) + dt * ( ...
               dTdt(i) + tau_q * dTdt(i) ...
             - tau_T * d2Tdt2(i) ...
             + (k(L) + k_star(L)*tau_v) * dTdt(i) );
        end
        
        % Convective BC at i=1 (left) & i=nx (right)
        L1 = layer(1);
        T_new(1)  = (h * dx * Tl(L1) + k(L1) * T_new(2))   / (h*dx + k(L1));
        L2 = layer(nx);
        T_new(nx) = (h * dx * Tl(L2) + k(L2) * T_new(nx-1)) / (h*dx + k(L2));
        
        % Advance
        T = T_new;
    end
    
    % Store final profile at t=10 s
    switch scenario
        case 1
            T_full = T;
        case 2
            T_noQv = T;
        case 3
            T_noQvQd = T;
    end
end

%% 5) PLOT THE THREE FINAL PROFILES AT t=10 s, WITH ZOOM‐INSET
hFig = figure('Position',[200,200,800,500]);
hMainAx = axes('Parent',hFig);
hold(hMainAx,'on');

plot(hMainAx, x*1000, T_full,   'r',  'LineWidth',1.5, 'DisplayName','Full model');
plot(hMainAx, x*1000, T_noQv,   'b--', 'LineWidth',1.5, 'DisplayName','No Q_v');
plot(hMainAx, x*1000, T_noQvQd, 'g',  'LineWidth',1.5, 'DisplayName','No Q_v & Q_d');

% Mark layer boundaries and tumor on main plot
xline(hMainAx, 0,                '-', {'Epidermis'},        'HandleVisibility','off','FontSize',12);
xline(hMainAx, L_epi*1e3,        '-', {'Dermis'},           'HandleVisibility','off','FontSize',12);
xline(hMainAx, (L_epi+L_derm)*1e3,'-', {'Subcutaneous'},    'HandleVisibility','off','FontSize',12);
xline(hMainAx, x_ast*1e3,        '-r', {'Tumor'},          'LineWidth',3, 'HandleVisibility','off','FontSize',14);

xlabel(hMainAx,'Distance (mm)','FontSize',14);
ylabel(hMainAx,'Temperature (°C)','FontSize',14);
title(hMainAx,'T(x) at t = 10 s: full vs. no Q_v vs. no Q_v & Q_d','FontSize',16);
legend(hMainAx,'Location','best','FontSize',12);
grid(hMainAx,'on');
hold(hMainAx,'off');

% --------------- ZOOM-IN INSET ----------------
% Zoom window: x ∈ [7.8, 7.9] mm, T ∈ [44.82, 44.87] °C
x_zoom_min = 7.8;   x_zoom_max = 7.9;
y_zoom_min = 44.82; y_zoom_max = 44.87;

% 1) Draw dashed rectangle on main axes to denote zoom region
hold(hMainAx,'on');
rectangle(hMainAx, ...
    'Position',[x_zoom_min, y_zoom_min, x_zoom_max - x_zoom_min, y_zoom_max - y_zoom_min], ...
    'EdgeColor','k', 'LineStyle','-', 'LineWidth',1.2);
hold(hMainAx,'off');

% 2) Create inset axes in upper-right portion of the figure
insetPos = [0.75, 0.6, 0.1, 0.1];  % [left bottom width height] in normalized units
hInsetAx = axes('Parent',hFig, 'Position',insetPos);
box(hInsetAx,'on');
hold(hInsetAx,'on');

% 3) Re-plot the same three curves within the inset
plot(hInsetAx, x*1000, T_full,   'r-',  'LineWidth',1);
plot(hInsetAx, x*1000, T_noQv,   'b--', 'LineWidth',1);
plot(hInsetAx, x*1000, T_noQvQd, 'g',  'LineWidth',1);

% 4) Set zoom‐in limits on inset
xlim(hInsetAx, [x_zoom_min, x_zoom_max]);
ylim(hInsetAx, [y_zoom_min, y_zoom_max]);
xlabel('Zoomed in region');

% 5) Optionally remove tick labels for a clean inset look
set(hInsetAx, 'XTick', [], 'YTick', []);
hold(hInsetAx,'off');

%% 6) PLOT: TUMOR TEMPERATURE VS. TIME (FULL MODEL)
% figure;
% plot(time, T_full(iTumor)*ones(size(time)), 'r', 'LineWidth',2);
% xlim([0, max_time]);
% xlabel('Time (s)',            'FontSize',14);
% ylabel('Tumor Temperature (°C)','FontSize',14);
% title('Tumor Temperature vs. Time (full model final)','FontSize',16);
% grid on;
