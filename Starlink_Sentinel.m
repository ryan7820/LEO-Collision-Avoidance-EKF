% =========================================================================
% THE STARLINK SENTINEL: Autonomous COLA & EKF Estimation
% =========================================================================
% Simulates a LEO satellite and a debris object with J2 perturbations.
% Uses an Extended Kalman Filter (EKF) to estimate debris trajectory from
% noisy radar data, and executes a continuous low-thrust maneuver to avoid.
% =========================================================================
clear; clc; close all;
% 1. SYSTEM PARAMETERS & CONSTANTS
mu = 3.986004418e14;    
Re = 6378137.0;         
J2 = 1.08262668e-3;     
% Simulation Settings
dt = 5;                 
t_end = 3000;           
time = 0:dt:t_end;
N = length(time);
% 2. INITIAL CONDITIONS
% Satellite (Starlink-ish: ~550km altitude)
r0_sat = [Re + 550000; 0; 0];
v0_sat = [0; sqrt(mu/norm(r0_sat))*cosd(53); sqrt(mu/norm(r0_sat))*sind(53)];
state_sat = zeros(6, N);
state_sat(:,1) = [r0_sat; v0_sat];
% Debris (Truth Model - Intersecting trajectory)
r0_deb = [Re + 540000; -20000; 10000]; 
v0_deb = [10; sqrt(mu/norm(r0_deb))*cosd(54); sqrt(mu/norm(r0_deb))*sind(54)];
state_deb_truth = zeros(6, N);
state_deb_truth(:,1) = [r0_deb; v0_deb];
%% 3. EKF INITIALIZATION (Debris Estimation)
% We only get noisy position data from ground radar
meas_noise_std = 50; 
R = eye(3) * meas_noise_std^2; 
Q = eye(6) * 1e-4;             
H = [eye(3), zeros(3,3)];      
state_deb_est = zeros(6, N);
state_deb_est(:,1) = [r0_deb + randn(3,1)*meas_noise_std; v0_deb]; 
P = eye(6) * 1e3; 
%% 4. CONTROL SETTINGS (Low-Thrust Hall Effect)
maneuver_active = false;
thrust_accel = 0.005; 
safe_distance = 2500; 
Kp = 1e-5;
Kd = 5e-4;
%% 5. MAIN SIMULATION LOOP
disp('Simulating Constellation Dynamics & EKF...');
for k = 1:N-1
    
    % --- A. TRUTH DYNAMICS (RK4 with J2) ---
    state_deb_truth(:,k+1) = rk4_step(@(t,x) orbit_dynamics_J2(x, mu, Re, J2), ...
                                      time(k), state_deb_truth(:,k), dt);
                                  
    % --- B. SENSOR MEASUREMENT (Noisy Radar) ---
    % Ground station measures debris position with noise
    z_k = H * state_deb_truth(:,k+1) + randn(3,1) * meas_noise_std;
    
    % --- C. EXTENDED KALMAN FILTER (Predict & Update) ---
    % 1. Predict state
    x_pred = rk4_step(@(t,x) orbit_dynamics_J2(x, mu, Re, J2), ...
                      time(k), state_deb_est(:,k), dt);
    
    % 2. Predict Covariance (Using numerical Jacobian Phi)
    % Computes Numerical Jacobian for the EKF State Transition Matrix (Phi)
    r_est = state_deb_est(1:3,k);
    r_norm = norm(r_est);
    G = (3*mu / r_norm^5) * (r_est * r_est') - (mu / r_norm^3) * eye(3);
    A = [zeros(3,3), eye(3); G, zeros(3,3)];
    Phi = eye(6) + A * dt + 0.5 * A^2 * dt^2;
    P_pred = Phi * P * Phi' + Q;
    
    % 3. Update
    y_tilde = z_k - H * x_pred;               
    S = H * P_pred * H' + R;                  
    K = P_pred * H' / S;                      
    
    state_deb_est(:,k+1) = x_pred + K * y_tilde;
    P = (eye(6) - K * H) * P_pred;
    
    % --- D. COLLISION AVOIDANCE CONTROL (LQR / PD concept) ---
    % Calculate estimated distance
    dist_est = norm(state_sat(1:3, k) - state_deb_est(1:3, k+1));
    v_rel = norm(state_sat(4:6, k) - state_deb_est(4:6, k+1));
    
    % Control Logic: If debris is too close, fire Hall thrusters
    u_thrust = [0;0;0]; 
    if dist_est < safe_distance || maneuver_active
        maneuver_active = true; 
        
        % Calculate along-track vector for efficient orbit raising
        v_vec = state_sat(4:6, k);
        along_track_dir = v_vec / norm(v_vec);
        
        error_dist = (safe_distance * 1.5) - dist_est;
        pd_mag = Kp * error_dist + Kd * v_rel;
        u_thrust = along_track_dir * min(max(pd_mag, 0), thrust_accel);
        
        % Turn off thrust if distance opens up safely
        if dist_est > safe_distance * 2
            maneuver_active = false;
        end
    end
    
    % Propagate Satellite Truth with Control Input
    state_sat(:,k+1) = rk4_step(@(t,x) orbit_dynamics_J2(x, mu, Re, J2, u_thrust), ...
                                time(k), state_sat(:,k), dt);
end
disp('Simulation Complete. Generating telemetry visualizers...');
% 6. VISUALIZATION & TELEMETRY
figure('Name', 'Starlink Sentinel Telemetry', 'Color', 'w', 'Position', [100, 100, 1200, 500]);
% Plot 1: 3D Trajectory
subplot(1,2,1); hold on; grid on; view(3);
plot3(state_sat(1,:)/1000, state_sat(2,:)/1000, state_sat(3,:)/1000, 'b', 'LineWidth', 2);
plot3(state_deb_truth(1,:)/1000, state_deb_truth(2,:)/1000, state_deb_truth(3,:)/1000, 'r--', 'LineWidth', 1.5);
plot3(state_deb_est(1,:)/1000, state_deb_est(2,:)/1000, state_deb_est(3,:)/1000, 'g:', 'LineWidth', 1.5);
% Draw a simple Earth sphere
[Ex, Ey, Ez] = sphere(30);
surf(Ex*Re/1000, Ey*Re/1000, Ez*Re/1000, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
title('Orbital Trajectories (km)');
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
legend('Starlink Sat', 'Debris (Truth)', 'Debris (EKF Est)', 'Earth');
axis equal;
% Plot 2: Miss Distance & EKF Error
subplot(1,2,2); hold on; grid on;
true_dist = vecnorm(state_sat(1:3,:) - state_deb_truth(1:3,:));
est_dist = vecnorm(state_sat(1:3,:) - state_deb_est(1:3,:));
ekf_error = vecnorm(state_deb_truth(1:3,:) - state_deb_est(1:3,:));
plot(time/60, true_dist/1000, 'b', 'LineWidth', 2);
plot(time/60, est_dist/1000, 'g--', 'LineWidth', 1.5);
plot(time/60, ekf_error/1000, 'r', 'LineWidth', 1.5);
yline(safe_distance/1000, 'k:', 'Safe Threshold', 'LineWidth', 1.5);
title('Telemetry & Filter Performance');
xlabel('Time (Minutes)'); ylabel('Distance / Error (km)');
legend('True Sat-Debris Dist', 'Estimated Sat-Debris Dist', 'EKF Tracking Error');
% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function dx = orbit_dynamics_J2(x, mu, Re, J2, u)
    % Calculates state derivatives with J2 perturbation and optional thrust
    if nargin < 5
        u = [0;0;0]; 
    end
    
    r_vec = x(1:3);
    v_vec = x(4:6);
    r = norm(r_vec);
    
    % Spherical J2 factor
    factor = (3/2) * J2 * (Re/r)^2;
    z2_r2 = (r_vec(3)/r)^2;
    
    % Accelerations
    ax = -(mu*r_vec(1)/r^3) * (1 - factor * (5*z2_r2 - 1)) + u(1);
    ay = -(mu*r_vec(2)/r^3) * (1 - factor * (5*z2_r2 - 1)) + u(2);
    az = -(mu*r_vec(3)/r^3) * (1 - factor * (5*z2_r2 - 3)) + u(3);
    
    dx = [v_vec; ax; ay; az];
end
function x_next = rk4_step(dyn_func, t, x, dt)
    % Runge-Kutta 4th Order Integrator
    k1 = dyn_func(t, x);
    k2 = dyn_func(t + dt/2, x + k1*dt/2);
    k3 = dyn_func(t + dt/2, x + k2*dt/2);
    k4 = dyn_func(t + dt, x + k3*dt);
    x_next = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end