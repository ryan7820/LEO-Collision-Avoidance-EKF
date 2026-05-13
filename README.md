# LEO-Collision-Avoidance-EKF
A 6-DOF orbital dynamics simulator in MATLAB featuring J2 perturbations, an Extended Kalman Filter (EKF) for state estimation, and autonomous low-thrust collision avoidance (COLA) maneuvers.


Overview
Starlink Sentinel is a high-fidelity MATLAB simulation designed to model autonomous Collision Avoidance (COLA) for Low Earth Orbit (LEO) mega-constellations. The project bridges the gap between theoretical astrodynamics and practical flight software logic. It features a custom Runge-Kutta 4th Order (RK4) integrator that accounts for Earth's oblateness (J2 perturbations), an Extended Kalman Filter (EKF) to predict debris trajectories from noisy ground radar data, and a closed-loop Proportional-Derivative (PD) control system to execute continuous-thrust evasion maneuvers.

This project was designed by Ryan Epperly in May 2026 as an undergraduate Aerospace Engineering Student at Virginia Tech
