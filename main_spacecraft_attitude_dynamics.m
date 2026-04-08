%% Spacecraft Attitude Dynamics Simulation
%
% Authors:
% Pasquale Marzaioli
% Paolo Iaccarino
% Federica Pirozzi
% Silvia Preziosi

%%

close all
clear
clc

%% --- Orbit characterization

%--- Earth's Costants
mu_earth = 3.986*10^14;                     % Standard gravitational parameter [m^3 s^-2]
R_earth  = 6378.1363;                       % [km]
H = 700;                                    % [km]
r_earth = (R_earth + H)*10^3;               % [m]

%--- Orbital Parameters
T = 98.6 * 60;                              % [s]
N = 2;                                      % [days]
LAT = 18.00;                                % [hours]
a = (((T/(2*pi))^2*mu_earth)^(1/3));        % [m]
e = 0.001;                                  % [-]
rp = a * (1-e);                             % [km]
ra = a * (1+e);                             % [km]
i = deg2rad(98.18);                         % [rad]
om = deg2rad(0);                            % [rad]
th_0 = deg2rad(0);                          % [rad]
theta_G0 = deg2rad(0);                      % [rad]
RAAN0 = 0;                                  % [°]
precession_rate = 0.9856;                   % [°/s]

% Typically RAAN0 is defined in reference to a fixed direction 
% such the vernal equinox and must account for the inclination of the orbit 
% (generally high for most Earth observation as Sentinel 1)

for j = 1:N
    RAAN = RAAN0 + precession_rate*j;
end
RAAN = deg2rad(RAAN);                       % [rad]

%% --- Geometry Data

% Upload S/C Model
[EXTRACTED_DATA, MODEL_SC] = SCModel;

% Extract Structure's data 
N_B = EXTRACTED_DATA.N_B;           % Surface normals
r_F = EXTRACTED_DATA.r_F;           % Force vector radii
surface = EXTRACTED_DATA.surface;   % Surface areas
rho_S = EXTRACTED_DATA.rho_S;       % Reflection coefficients
rho_D = EXTRACTED_DATA.rho_D;       % Dissipation coefficients
cm_shift = EXTRACTED_DATA.cm_shift; % Center of mass shift

% Assign the extracted data to the required variables
% Main body
N_B1 = N_B(1, :)'; r_f1 = r_F(1, :)'; A_body1 = surface(1); ps_body1 = rho_S(1); pd_body1 = rho_D(1);
N_B2 = N_B(2, :)'; r_f2 = r_F(2, :)'; A_body2 = surface(2); ps_body2 = rho_S(2); pd_body2 = rho_D(2);
N_B3 = N_B(3, :)'; r_f3 = r_F(3, :)'; A_body3 = surface(3); ps_body3 = rho_S(3); pd_body3 = rho_D(3);
N_B4 = N_B(4, :)'; r_f4 = r_F(4, :)'; A_body4 = surface(4); ps_body4 = rho_S(4); pd_body4 = rho_D(4);
N_B5 = N_B(5, :)'; r_f5 = r_F(5, :)'; A_body5 = surface(5); ps_body5 = rho_S(5); pd_body5 = rho_D(5);
N_B6 = N_B(6, :)'; r_f6 = r_F(6, :)'; A_body6 = surface(6); ps_body6 = rho_S(6); pd_body6 = rho_D(6);

% Solar panels
N_P1 = N_B(7, :)'; r_f7 = r_F(7, :)'; A_panel1 = surface(7); ps_panel1 = rho_S(7); pd_panel1 = rho_D(7);
N_P2 = N_B(8, :)'; r_f8 = r_F(8, :)'; A_panel2 = surface(8); ps_panel2 = rho_S(8); pd_panel2 = rho_D(8);
N_P3 = N_B(9, :)'; r_f9 = r_F(9, :)'; A_panel3 = surface(9); ps_panel3 = rho_S(9); pd_panel3 = rho_D(9);
N_P4 = N_B(10, :)'; r_f10 = r_F(10, :)'; A_panel4 = surface(10); ps_panel4 = rho_S(10); pd_panel4 = rho_D(10);

% Total moments of inertia
J = EXTRACTED_DATA.I_total;  % Inertia matrix

%% --- Perturbations

% --- SOLAR RADIATION PRESSURE
T_S = 365*24*60*60;                      % [s]
n_S = (2*pi)/T_S;                        % [rad/s]
eps = deg2rad(23.45);                    % [rad]
r_Sun = 149597870;                       % [km]
Fe = 1358;                               % [W/m^2]
c = 299792458;                           % [m/s]
P = Fe/c;                                % [N/m^2]
w_E = (2*pi)/(24*60*60);                 % [rad/s]

% --- MAGNETIC FIELD
DGRF2020 = IGRF2020MagneticField;
alpha_M = deg2rad(11.5);                 % [rad]
g1_0 = DGRF2020(1).g_nm;                 % [T]
g1_1 = DGRF2020(2).g_nm;                 % [T]
h1_1 = DGRF2020(2).h_nm;                 % [T]
j_B = [0.01, 0.05, 0.01]';               % [A m^2]

% --- GRAVITY GRADIENT
G  = 6.67430e-11;                        % [m^3/kg^1*s^2]
M_E = 5.972e24;                          % [kg]  
n  = sqrt((G*M_E)/(r_earth^3));          % [rad]/[s]

% --- AIR DRAG
CD = 2.6;   
v0 = [0 0 0];
AtmosphereModel = AtmosphereData;

%% --- Definition of sensor values and Attitude determination Parameters

A_epss = eye(3);

% Gyroscope
a_gyro = 43.69;                               % 0.01 [arcsec/LSB]
sample_f = 1000 ;                             % [Hz]
Ts = 1/sample_f;                              % [s]
ARW = (0.015);                                % [deg/h^1/2]
sigma_n = ARW*(pi/180)*(1/60)*(1/sqrt(Ts));   % [rad/s]
bias_instab = 0.08;
RRW = bias_instab*(sqrt(Ts));                 % [deg/h^3/2]
sigma_b = RRW*(pi/180)*(1/3600)*(1/sqrt(Ts)); % [rad/s^2]
f_cut = 150;                                  % [Hz]
zetaf = 0.707;
Wfn = 1/f_cut;

% --- From datasheet: mini fine sun sensor BRADFORD (SUN SENSOR)
% We are using 6 Sun sensors
FOV = 128 * pi/180;
a_sun = 1.8;
variance_Ss = (a_sun)^2;
R1 = @(orienatation_Ss) [1 0 0
    0 cos(orienatation_Ss) -sin(orienatation_Ss)
    0 sin(orienatation_Ss) cos(orienatation_Ss)];
R2 = @(orienatation_Ss) [cos(orienatation_Ss) 0 -sin(orienatation_Ss)
    0 1 0
    sin(orienatation_Ss) 0 cos(orienatation_Ss)];
R3 = @(orienatation_Ss) [cos(orienatation_Ss) sin(orienatation_Ss) 0
    -sin(orienatation_Ss) cos(orienatation_Ss) 0
    0 0 1];
A_Ss_orienation(:,:,1) = eye(3);
A_Ss_orienation(:,:,2) = R2(pi/2);
A_Ss_orienation(:,:,3) = R1(-pi/2);
A_Ss_orienation(:,:,4) = R2(pi);
A_Ss_orienation(:,:,5) = R2(-pi/2);
A_Ss_orienation(:,:,6) = R1(pi/2);

% --- From datasheet: MAG3 (MAGNETOMETER)
a_mag = 20.97; % in [arcminutes] from [0.75% of FS]
sigma_mag = 10^-6;

% --- Normalisation of weight factor
alpha_sun_hat = 1/a_sun;
alpha_mag_hat =  1/a_mag;
alpha_gyro_hat = 1/a_gyro;

% come gestirli in questo caso?
alpha1 = alpha_sun_hat/(alpha_sun_hat+alpha_gyro_hat+alpha_mag_hat); 
alpha2 = alpha_mag_hat/ (alpha_sun_hat+alpha_gyro_hat+alpha_mag_hat);
alpha3 = alpha_gyro_hat/(alpha_sun_hat+alpha_gyro_hat+alpha_mag_hat);
alpha = [alpha1,alpha2,alpha3];

%% LQR with Normal State for Detumbling and Slew Maneuver

% ---- Satellite Parameters ----
omega_0 = 2 * pi / T;                      % Orbital angular velocity [rad/s]
Jx = J(1,1);                               % Principal moment of inertia along x-axis
Jy = J(2,2);                               % Principal moment of inertia along y-axis
Jz = J(3,3);                               % Principal moment of inertia along z-axis

% Coefficients derived from rotational dynamics
f41 = 8 * (Jz - Jy) * omega_0^2 - 2 * omega_0;
f46 = (- Jx - Jz + Jy) * omega_0;
f52 = 6 * (Jz - Jx) * omega_0^2; 
f63 = 2 * (Jx - Jy) * omega_0^2 - 2 * omega_0;
f64 = -f46; 

% ---- State-Space Matrices ----
% A-matrix: State transition matrix
A1 = [ 0 0 0 1/2  0   0
       0 0 0  0  1/2  0
       0 0 0  0   0  1/2];                % Quaternion dynamics

A2 = [f41  0  0   0  0 f46
        0  f52 0   0  0  0
        0   0 f63 f64 0  0];              % Angular velocity dynamics

A = [A1; A2];

% B-matrix: Input matrix (magnetic torquer control)
B_MAG = [ 0 0 0
          0 0 0
          0 0 0
          1 0 0
          0 1 0 
          0 0 1];

% ---- Controllability Check ----
n_A = size(A);                              % Size of A-matrix
Co = ctrb(A, B_MAG);                      % Controllability matri

% ---- LQR Design for Different Phases ----
% Phase 1: Detumbling
Q_detumbling = diag([50, 50, 50, 50, 50, 50]);  % State weighting matrix
R_detumbling = diag([100, 100, 100]);                    % Control weighting matrix

% Phase 2: Slew Maneuver
Q_slew = diag([0.1, 0.1, 0.1, 0.1, 0.1, 0.1]);                 % State weighting matrix
R_slew = diag([100, 100, 100]);                       % Control weighting matrix

% ---- Compute LQR Gains ----
% Detumbling Phase
[K_detumbling, ~, ~] = lqr(A, B_MAG, Q_detumbling, R_detumbling);
kp_detumbling = K_detumbling(:, 1:3) ;                % Proportional gain for detumbling (angular velocities)

% Slew Maneuver Phase
[K_slew, ~, ~] = lqr(A, B_MAG, Q_slew, R_slew);
kp_slew = K_slew(:, 1:3)     ;                        % Proportional gain for slew (quaternions)
kd_slew = K_slew(:, 4:6)    ;                         % Derivative gain for slew (angular velocities)

% ---- Assign Results ----
kp = struct('detumbling', kp_detumbling, 'slew', kp_slew); % Proportional gains
kd = struct('slew', kd_slew);                              % Derivative gains

Kw = max(max(kp_detumbling));
k2 = max(max(kp_slew));
k1 = max(max(kd_slew));
Detumbling_time = 5000;

%% --- Actuators 

% Magnetorquer
MagnetoTorquerDipoleError =  0.01;
MagnetoTorquerSampleRate = 0.01;
%MagnetoTorquerMaxDipole = 0.3045;
D_max = 50; % [A*m^2]
%T_coil_max = 50 * 0.0141;
T_coil_max = 50 * 45*10^-6;

% Inertia Wheel
% Ir = 3.4e-4;
Ir = 0.1;
Jr = Ir.*eye(3);                
Ar = [8; 1; 1];                    % Configuration matrix
Ap = pinv(Ar);
Ar = Ar./norm(Ar);
h0 = 0;                           % Nominal condition of inertia wheel
w_wheel_max = 250; 
h_r_max = w_wheel_max * Ir;
t_saturation = h_r_max/T_coil_max;
w_r0 = [0,0,0];

%% --- Simulation

theta_0 = 0;                       % [rad] Initial condition for dynamics integration
wx0 = 2.5;                         % [rad/s] Initial angular velocity x-comp
wy0 = 2.6;                         % [rad/s] Initial angular velocity y-comp
wz0 = -2.8;                        % [rad/s] Initial angular velocity z-comp
w_0 = [wx0; wy0; wz0];             % Initial angular velocity [rad s^-1]
th0 = 0;                           % [rad]
InvJ = inv(J);
Ix = J(1,1);
Iy = J(2,2);
Iz = J(3,3);
A0 = eye(3);
A_BN_0 = eye(3);
w_LN = [0 0 0]';
q4_0 = 0.5*sqrt(1+A_BN_0(1,1)+A_BN_0(2,2)+A_BN_0(3,3)); % we consider the positive case 
q1_0 =(1/(4*q4_0))*(A_BN_0(2,3)-A_BN_0(3,2));
q2_0 =(1/(4*q4_0))*(A_BN_0(3,1)-A_BN_0(1,3));
q3_0 =(1/(4*q4_0))*(A_BN_0(1,2)-A_BN_0(2,1));
q0 = [q1_0 q2_0 q3_0 q4_0];

% Time
n_orbit = 5;
t0 = 0;
tf = n_orbit * T;

%% Results
out = sim("spacecraft_attitude_model.slx","StartTime","t0","StopTime","tf");

Slew_time = 25000;       % Tempo di fine slew (s)
Tracking_time = 40000;   % Tempo finale di tracking (s)
time      = out.time;
% Trova l'indice corrispondente al tempo di detumbling
index_detumbling = find(time >= Detumbling_time, 1);
% Crea il nuovo vettore time_new
time_new = time(index_detumbling:end);
T_GG      = squeeze(out.T_GG);
T_SRP     = out.T_sun_tot;
T_B       = squeeze(out.T_M);
T_DRAG    = squeeze(out.T_Drag);
Point_err = out.pointingerror;
w_bn      = out.w; 
w_meas    = squeeze(out.w_measured); 
q_bn      = out.q;
q_err     = out.q_error;  
q_mes     = out.q_measured; 
T_act    = squeeze(out.T_act);

norm_T_GG     = vecnorm(T_GG, 2, 2);
norm_T_SRP    = vecnorm(T_SRP, 2, 2);
norm_T_B      = vecnorm(T_B, 2, 2);
norm_T_DRAG   = vecnorm(T_DRAG, 2, 2);


%% First plot with logarithmic axis

figure;
semilogy(time, norm_T_GG, 'LineWidth', 1); hold on;
semilogy(time, norm_T_SRP, 'LineWidth', 1.5);
semilogy(time, norm_T_B,  'LineWidth', 1.5);
semilogy(time, norm_T_DRAG, 'LineWidth', 1.5);
grid on;
legend('T_{GG}', 'T_{SRP}', 'T_{B}','T_{DRAG}' );
xlabel('Time [s]');
ylabel('Torque Magnitude [Nm]');
title('Comparison of Torques (Log Scale)');

avg_T_GG    = mean(norm_T_GG);
avg_T_SRP   = mean(norm_T_SRP);
avg_T_B     = mean(norm_T_B);
avg_T_DRAG  = mean(norm_T_DRAG);

figure;
bar([avg_T_GG, avg_T_SRP, avg_T_B, avg_T_DRAG], 'FaceColor', 'flat');
xticklabels({'T_{GG}', 'T_{SRP}', 'T_{B}','T_{DRAG}'});
ylabel('Average Torque [Nm]');
title('Average Torque Comparison');
grid on

% Bar plot with logarithmic scale
figure;
bar([avg_T_GG, avg_T_SRP, avg_T_B, avg_T_DRAG], 'FaceColor', 'flat');
set(gca, 'YScale', 'log'); % Set Y-axis to logarithmic scale
xticklabels({'T_{GG}', 'T_{SRP}', 'T_{B}','T_{DRAG}'});
ylabel('Average Torque [Nm] (Log Scale)');
title('Average Torque Comparison (Log Scale)');
grid on;

%% Plot 1: Pointing Error
figure;
hold on;
plot(time, Point_err, 'r-', 'LineWidth', 1.5);

% Dashed lines for transition times
xline(Detumbling_time, 'k--', 'LineWidth', 1.2);
xline(Slew_time, 'k--', 'LineWidth', 1.2);

% Phase annotations
text(Detumbling_time / 2, max(Point_err) * 0.9, 'Detumbling', 'HorizontalAlignment', 'center', 'FontSize', 10);
text((Detumbling_time + Slew_time) / 2, max(Point_err) * 0.9, 'Slew', 'HorizontalAlignment', 'center', 'FontSize', 10);
text((12000 + Tracking_time) / 2, max(Point_err) * 0.9, 'Tracking', 'HorizontalAlignment', 'center', 'FontSize', 10);

% Axis configuration
xlabel('Time [s]', 'FontSize', 12);
ylabel('Pointing Error [deg]', 'FontSize', 12);
ylim([0 max(Point_err) * 1.1]);
xlim([0 tf]);
title('Pointing Error', 'FontSize', 14);
grid on;
hold off;

%% Plot 2: Angular Velocity
figure;
hold on;
plot(time, w_meas(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x');
plot(time, w_meas(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y');
plot(time, w_meas(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z');

% Dashed lines for transition times
xline(Detumbling_time, 'k--', 'LineWidth', 1.2);
xline(Slew_time, 'k--', 'LineWidth', 1.2);

% Phase annotations
text(Detumbling_time / 2, max(max(w_meas)) * 0.9, 'Detumbling', 'HorizontalAlignment', 'center', 'FontSize', 10);
text((Detumbling_time + Slew_time) / 2, max(max(w_meas)) * 0.9, 'Slew', 'HorizontalAlignment', 'center', 'FontSize', 10);
text((12000 + Tracking_time) / 2, max(max(w_meas)) * 0.9, 'Tracking', 'HorizontalAlignment', 'center', 'FontSize', 10);

% Axis configuration
xlabel('Time [s]', 'FontSize', 12);
ylabel('Angular Velocity [rad/s]', 'FontSize', 12);
ylim([min(min(w_meas)) * 1.1, max(max(w_meas)) * 1.1]);
xlim([0 tf]);
title('Angular velocity measured', 'FontSize', 14);
grid on;
hold off;

%% Plot 3: Attitude Error
figure;
hold on;
q_err = q_err';
plot(time_new, q_err(1, :), 'LineWidth', 1.5, 'DisplayName', 'q1');
plot(time_new, q_err(2, :), 'LineWidth', 1.5, 'DisplayName', 'q2');
plot(time_new, q_err(3, :), 'LineWidth', 1.5, 'DisplayName', 'q3');
plot(time_new, q_err(4, :), 'LineWidth', 1.5, 'DisplayName', 'q0');

% Dashed lines for transition times
xline(Detumbling_time, 'k--', 'LineWidth', 1.2);
xline(Slew_time, 'k--', 'LineWidth', 1.2);

% Phase annotations
text(Detumbling_time / 2, max(max(w_meas)) * 0.9, 'Detumbling', 'HorizontalAlignment', 'center', 'FontSize', 10);
text((Detumbling_time + Slew_time) / 2, max(max(w_meas)) * 0.9, 'Slew', 'HorizontalAlignment', 'center', 'FontSize', 10);
text((12000 + Tracking_time) / 2, max(max(w_meas)) * 0.9, 'Tracking', 'HorizontalAlignment', 'center', 'FontSize', 10);

% Axis configuration
xlabel('Time [s]', 'FontSize', 12);
ylabel('Quaternion error', 'FontSize', 12);
ylim([min(min(q_err)) * 1.1, max(max(q_err)) * 1.1]);
legend('q1','q2','q3','q4','Location','Best')
xlim([0 tf]);
title('q error', 'FontSize', 14);
grid on;
hold off;
 


%% Random initial angular velocities

w_0 = [0.5; 0.5; -0.5];
out05 = sim("spacecraft_attitude_model.slx","StartTime","t0","StopTime","tf");
w05 = squeeze(out05.w);

w_0 = [1, -1, 1]';
out1 = sim("spacecraft_attitude_model.slx","StartTime","t0","StopTime","tf");
w1 = squeeze(out1.w);

w_0 = [1.5, 1.5, -1.5]';
out15 = sim("spacecraft_attitude_model.slx","StartTime","t0","StopTime","tf");
w15 = squeeze(out15.w);

w_0 = [-2, 2, 2]';
out2 = sim("spacecraft_attitude_model.slx","StartTime","t0","StopTime","tf");
w2 = squeeze(out2.w);

w_0 = [2.5, -2.5, -2.5]';
out25 = sim("spacecraft_attitude_model.slx","StartTime","t0","StopTime","tf");
w25 = squeeze(out25.w);

w_0 = [-3, 3, 3]';
out3 = sim("spacecraft_attitude_model.slx","StartTime","t0","StopTime","tf");
w3 = squeeze(out3.w);


figure;
hold on;
plot(time, w05(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x', 'Color', 'k');
plot(time, w05(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y', 'Color', 'k');
plot(time, w05(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z', 'Color', 'k');
plot(time, w1(1,:), 'LineWidth', 1.5, 'DisplayName', 'w_x','Color', 'g');
plot(time, w1(2,:), 'LineWidth', 1.5, 'DisplayName', 'w_y','Color', 'g');
plot(time, w1(3,:), 'LineWidth', 1.5, 'DisplayName', 'w_z','Color', 'g');
plot(time, w15(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x','Color', 'r');
plot(time, w15(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y','Color', 'r');
plot(time, w15(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z','Color', 'r');
plot(time, w2(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x','Color', 'b');
plot(time, w2(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y','Color', 'b');
plot(time, w2(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z','Color', 'b');
plot(time, w25(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x','Color', 'y');
plot(time, w25(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_y','Color', 'y');
plot(time, w25(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z','Color', 'y');
plot(time, w3(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x','Color','c');
plot(time, w3(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y','Color','c');
plot(time, w3(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z','Color','c');
xlabel('Time');
ylabel('w0');
legend;
hold off;

% Create a new figure
figure;
% First subplot
subplot(2, 3, 1); % 2 rows, 3 columns, position 1
hold on;
plot(time, w05(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x', 'Color', 'k');
plot(time, w05(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y', 'Color', 'k');
plot(time, w05(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z', 'Color', 'k');
title('w0 = [0.5; 0.5; -0.5]');
xlabel('Time');
ylabel('w0');
legend;
hold off;
% Second subplot
subplot(2, 3, 2); % 2 rows, 3 columns, position 2
hold on;
plot(time, w1(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x', 'Color', 'g');
plot(time, w1(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y', 'Color', 'g');
plot(time, w1(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z', 'Color', 'g');
title('w0 = [1, -1, 1]');
xlabel('Time');
ylabel('w0');
legend;
hold off;
% Third subplot
subplot(2, 3, 3); % 2 rows, 3 columns, position 3
hold on;
plot(time, w15(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x', 'Color', 'r');
plot(time, w15(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y', 'Color', 'r');
plot(time, w15(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z', 'Color', 'r');
title('w0 = [1.5, 1.5, -1.5]');
xlabel('Time');
ylabel('w0');
legend;
hold off;
% Fourth subplot
subplot(2, 3, 4); % 2 rows, 3 columns, position 4
hold on;
plot(time, w2(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x', 'Color', 'b');
plot(time, w2(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y', 'Color', 'b');
plot(time, w2(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z', 'Color', 'b');
title('w_0 = [-2, 2, 2]');
xlabel('Time');
ylabel('w0');
legend;
hold off;
% Fifth subplot
subplot(2, 3, 5); % 2 rows, 3 columns, position 5
hold on;
plot(time, w25(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x', 'Color', 'y');
plot(time, w25(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y', 'Color', 'y');
plot(time, w25(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z', 'Color', 'y');
title('w_0 = [2.5, -2.5, -2.5]');
xlabel('Time');
ylabel('w0');
legend;
hold off;
% Sixth subplot
subplot(2, 3, 6); % 2 rows, 3 columns, position 6
hold on;
plot(time, w3(1, :), 'LineWidth', 1.5, 'DisplayName', 'w_x', 'Color', 'c');
plot(time, w3(2, :), 'LineWidth', 1.5, 'DisplayName', 'w_y', 'Color', 'c');
plot(time, w3(3, :), 'LineWidth', 1.5, 'DisplayName', 'w_z', 'Color', 'c');
title('w_0 = [-3, 3, 3]');
xlabel('Time');
ylabel('w0');
legend;
hold off;

%% --- Functions

function AtmosphereModel = AtmosphereData
    
    altitude_values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, ...
                       110, 120, 130, 140, 150, 160, 170, 180, 190, 200, ...
                       210, 220, 230, 240, 250, 260, 270, 280, 290, 300, ...
                       310, 320, 330, 340, 350, 360, 370, 380, 390, 400, ...
                       410, 420, 430, 440, 450, 460, 470, 480, 490, 500, ...
                       510, 520, 530, 540, 550, 560, 570, 580, 590, 600, ...
                       610, 620, 630, 640, 650, 660, 670, 680, 690, 700, ...
                       710, 720, 730, 740, 750, 760, 770, 780, 790, 800, ...
                       810, 820, 830, 840, 850, 860, 870, 880, 890, 900, ...
                       910, 920, 930, 940, 950, 960, 970, 980, 990, 1000] ; % Altitude [m]
    
    density_values = [4.02e-04, 8.34e-05, 1.57e-05, 3.18e-06, 8.37e-07, ...
                      2.33e-07, 5.86e-08, 1.40e-08, 2.99e-09, 5.17e-10, ...
                      8.42e-11, 1.84e-11, 7.36e-12, 3.78e-12, 2.19e-12, ...
                      1.37e-12, 9.00e-13, 6.15e-13, 4.32e-13, 3.10e-13, ...
                      2.27e-13, 1.68e-13, 1.26e-13, 9.58e-14, 7.35e-14, ...
                      5.68e-14, 4.43e-14, 3.48e-14, 2.75e-14, 2.18e-14, ...
                      1.74e-14, 1.40e-14, 1.13e-14, 9.10e-15, 7.39e-15, ...
                      6.02e-15, 4.92e-15, 4.03e-15, 3.31e-15, 2.72e-15, ...
                      2.25e-15, 1.86e-15, 1.54e-15, 1.28e-15, 1.07e-15, ...
                      8.89e-16, 7.43e-16, 6.22e-16, 5.22e-16, 4.39e-16, ...
                      3.71e-16, 3.13e-16, 2.66e-16, 2.26e-16, 1.93e-16, ...
                      1.65e-16, 1.41e-16, 1.22e-16, 1.05e-16, 9.14e-17, ...
                      7.96e-17, 6.97e-17, 6.12e-17, 5.41e-17, 4.79e-17, ...
                      4.27e-17, 3.82e-17, 3.44e-17, 3.10e-17, 2.82e-17, ...
                      2.57e-17, 2.35e-17, 2.16e-17, 1.99e-17, 1.84e-17, ...
                      1.71e-17, 1.59e-17, 1.49e-17, 1.39e-17, 1.31e-17, ...
                      1.23e-17, 1.16e-17, 1.10e-17, 1.04e-17, 9.86e-18, ...
                      9.37e-18, 8.91e-18, 8.48e-18, 8.09e-18, 7.72e-18, ...
                      7.37e-18, 7.04e-18, 6.74e-18, 6.45e-18, 6.18e-18, ...
                      5.92e-18, 5.68e-18, 5.44e-18, 5.22e-18, 5.02e-18] * 10^3; % Density [kg/m^3]
    
    % Structure Creation 
    AtmosphereModel = struct();
    
    for i = 1:length(altitude_values)
        AtmosphereModel(i).Altitude_m = altitude_values(i);
        AtmosphereModel(i).Density_kg_m3 = density_values(i); 
        AtmosphereModel(i).Description = sprintf('Altitude Data %d m', altitude_values(i));
    end
end

function DGRF2020 = IGRF2020MagneticField
    
    % --- Data IGRF2020
    n_values = [1, 1, 2, 2, 2, 3, 3, 3]; 
    m_values = [0, 1, 0, 1, 2, 0, 1, 3]; 
    g_nm = [-29404.8, -1450.9, -2499.6, 2982.0, 1677.0, 1363.2, -2381.2, 525.7]; % [nT]
    h_nm = [0, 4652.5, 0, -2991.6, -734.6, 0, -82.1, -543.4]; % [nT]
    
    g_nm = g_nm * 10^-9; 
    h_nm = h_nm * 10^-9; 
    
    DGRF2020 = struct();
    
    for i = 1:length(n_values)
        DGRF2020(i).n = n_values(i);              
        DGRF2020(i).m = m_values(i);              
        DGRF2020(i).g_nm = g_nm(i);               
        DGRF2020(i).h_nm = h_nm(i);               
        DGRF2020(i).description = sprintf('Harmonic data for n=%d, m=%d', ...
                                           n_values(i), m_values(i)); 
    end
end 

function [EXTRACTED_DATA, MODEL_SC] = SCModel
    
    % -------------------------------
    % Initial parameters
    % -------------------------------
    clc;
    num_body  = 6;                     % Number of main body surfaces
    num_panel = 4;                     % Number of solar panels
    n_sur = num_body + num_panel;      % Total number of surfaces
    
    mass_body = 500;                   % Mass of the main body [kg]
    mass_panel = 3;                    % Mass of a solar panel [kg]
    
    % Dimensions of the main body (cuboid)
    length_body = 1;   % Length [m]
    width_body  = 1;   % Width [m]
    height_body = 1;   % Height [m]
    surface_value_body = width_body * height_body; % Area of one surface of the main body
    
    % Dimensions of the panels
    length_panel = 2;  % Length [m]
    width_panel  = 2;  % Width [m]
    surface_value_panel = length_panel * width_panel; % Area of a solar panel
    
    % -------------------------------
    % Surface area vector
    % -------------------------------
    surface = [ones(num_body, 1) * surface_value_body; ... % Areas of the main body
               ones(num_panel, 1) * surface_value_panel];  % Areas of the panels
    
    % -------------------------------
    % Moment of inertia calculations
    % -------------------------------
    % Moments of inertia for the main body
    I_body = diag([
        (mass_body / 12) * (width_body^2 + height_body^2), % Ixx
        (mass_body / 12) * (length_body^2 + height_body^2), % Iyy
        (mass_body / 12) * (length_body^2 + width_body^2)   % Izz
    ]);
    
    % Moments of inertia for a solar panel
    I_panel = diag([
        (mass_panel / 12) * (width_panel^2 + 0),  % Ixx (width)
        (mass_panel / 12) * (length_panel^2 + 0), % Iyy (length)
        (mass_panel / 12) * 0                     % Izz negligible
    ]);
    
    % Total moment of inertia (body + panels)
    I = I_body + num_panel * I_panel;
    InvI = inv(I);  % Inverse inertia matrix
    
    % -------------------------------
    % Preparation of vectors and matrices
    % -------------------------------
    N_B = [  1  0  0  1
             0  1  0  1
            -1  0  0  1
             0 -1  0  1
             0  0  1  1
             0  0 -1  1
             1  0  0  0
            -1  0  0  0
             1  0  0  0
            -1  0  0  0];  % Surface normals
    
    cm_shift = [0.0, -0.1, 0.0]; % Center of mass shift vector [m]
    % Initialization of vectors
    rho_S = zeros(n_sur, 1);          % Reflection coefficients
    rho_D = zeros(n_sur, 1);          % Dissipation coefficients
    r_F = zeros(n_sur, 3);            % Force vector radii
    type = cell(n_sur, 1);            % Surface type ('body' or 'panel')
    
    % -------------------------------
    % Calculation of force radii and types
    % -------------------------------
    for i = 1:n_sur
        if N_B(i, 4) == 1 % Main body
            type{i}  = 'body';
            rho_S(i) = 0.5;
            rho_D(i) = 0.3;
    
            if N_B(i, 1) == 1
                r_F(i, :) = [length_body / 2, 0, 0] + cm_shift;
            elseif N_B(i, 1) == -1
                r_F(i, :) = [-length_body / 2, 0, 0] + cm_shift;
            elseif N_B(i, 2) == 1
                r_F(i, :) = [0, width_body / 2, 0] + cm_shift;
            elseif N_B(i, 2) == -1
                r_F(i, :) = [0, -width_body / 2, 0] + cm_shift;
            elseif N_B(i, 3) == 1
                r_F(i, :) = [0, 0, height_body / 2] + cm_shift;
            elseif N_B(i, 3) == -1
                r_F(i, :) = [0, 0, -height_body / 2] + cm_shift;
            end
        else % Solar panels
            type{i}  = 'panel';
            rho_S(i) = 0.1;
            rho_D(i) = 0.05;
    
            % Calculation of relative positions of panels
            % panel_offset = width_body / 2 + width_panel / 2;
            panel_offset = 0;
            if i == 7 % First panel
                r_F(i, :) = [0, length_panel / 2, panel_offset] + cm_shift;
            elseif i == 8 % Second panel
                r_F(i, :) = [0, length_panel / 2, panel_offset] + cm_shift;
            elseif i == 9 % Third panel
                r_F(i, :) = [0, -length_panel / 2, panel_offset] + cm_shift;
            elseif i == 10 % Fourth panel
                r_F(i, :) = [0, -length_panel / 2, panel_offset] + cm_shift;
            end
        end
    end
    
    % -------------------------------
    % Updating the MODEL_SC structure
    % -------------------------------
    MODEL_SC = struct();
    for i = 1:n_sur
        MODEL_SC(i).N_B     = N_B(i, 1:3);
        MODEL_SC(i).rho_S   = rho_S(i);
        MODEL_SC(i).rho_D   = rho_D(i);
        MODEL_SC(i).surface = surface(i);
        MODEL_SC(i).r_F     = r_F(i, :);
        MODEL_SC(i).type    = type{i};
    end
    
    % -------------------------------
    % Consolidating data
    % -------------------------------
    EXTRACTED_DATA = struct();
    EXTRACTED_DATA.N_B       = N_B(:, 1:3);
    EXTRACTED_DATA.type      = type;
    EXTRACTED_DATA.rho_S     = rho_S;
    EXTRACTED_DATA.rho_D     = rho_D;
    EXTRACTED_DATA.surface   = surface;
    EXTRACTED_DATA.r_F       = r_F;
    EXTRACTED_DATA.I_total   = I;
    EXTRACTED_DATA.cm_shift  = cm_shift;
    
    % -------------------------------
    % Display
    % -------------------------------
    disp('Total Moments of Inertia [kg·m²]:');
    disp(I);
    
    disp('Center of Mass Shift [m]:');
    disp(cm_shift);
    
    disp('Force Vectors with Shift [m]:');
    disp(r_F);
    
    disp('Updated MODEL_SC structure:');
    disp(MODEL_SC);
end