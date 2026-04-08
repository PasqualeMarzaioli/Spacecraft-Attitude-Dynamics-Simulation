# Spacecraft Attitude Dynamics Project

## Overview
This project models the attitude dynamics of a large spacecraft using MATLAB and Simulink.
The folder contains the simulation model, main execution script, support functions, and the project report.

## Authors
- Pasquale Marzaioli
- Paolo Iaccarino
- Federica Pirozzi
- Silvia Preziosi

## Files in this folder
- `spacecraft_attitude_model.slx`: main Simulink model.
- `spacecraft_attitude_model.slxc`: Simulink cache file.
- `main_spacecraft_attitude_dynamics.m`: MATLAB script that configures and runs the simulation and plots results.
- `spacecraft_attitude_report.pdf`: project report with mission definition, dynamic model, sensors, actuators, and results.
- `SCModel.m`: spacecraft geometry and inertia extraction.
- `IGRF2020MagneticField.m`: simplified geomagnetic field coefficients.
- `AtmosphereData.m`: altitude-density atmospheric model.
- `Model_Images/`: folder containing images of the Simulink model blocks.

## How to run the project
1. Open MATLAB.
2. Set the working folder to `.../Spacecraft Attitude Dynamics/22`.
3. Make sure the folder is on the MATLAB path.
4. Run `main_spacecraft_attitude_dynamics.m`.

The script runs the Simulink model `spacecraft_attitude_model.slx` and produces plots of:
- disturbance torques
- pointing error
- measured angular velocity
- quaternion attitude errors
- comparison between random initial angular velocity conditions

## Support functions

### `SCModel.m`
Creates the spacecraft structure used by the simulation:
- spacecraft geometry and surface areas
- surface normals and panel orientation
- diffuse and specular optical coefficients
- mass and inertia properties
- center of mass shift and force radii

### `IGRF2020MagneticField.m`
Provides a simplified IGRF 2020 magnetic field data set:
- Gauss spherical harmonic coefficients `g_nm` and `h_nm`
- associated degree and order values
- output structure used for magnetic torque calculation

### `AtmosphereData.m`
Returns a table of atmospheric densities vs altitude:
- altitude grid from 10 m to 1000 m
- density values for drag modeling
- returns a structure array with density and description fields

## Model Images
The `Model_Images/` folder contains screenshots of the Simulink model subsystems:
- `Dynamics_Kinematics_Enviroment.png`: dynamics, kinematics, and environment blocks.
- `Disturbances_Errors.png`: disturbance torques and attitude error calculations.
- `Sensors_Actuators.png`: sensor and actuator modeling blocks.
- `Manouvers.png`: control and maneuver execution blocks.

## Simulink model structure
The main Simulink model is organized into subsystems that reflect the project architecture.
Relevant blocks include:

- `Envoirment`
  - computes environmental inputs for the spacecraft.
- `Disturbances Torque`
  - calculates external torques: gravity gradient, solar radiation pressure, drag, and magnetic torque.
- `Euler_Dynamics`
  - spacecraft rotational dynamics in the body frame.
- `Kinematics`
  - quaternion integration and attitude update.
- `Attitude_determination_Sun`
  - attitude estimation using sun sensors and magnetometer data.
- `Gyroscope`
  - gyroscope measurement model including noise and bias.
- `Sun Sensor `
  - sun sensor measurement model.
- `Torque to w_dot_control`
  - converts actuator torque commands into angular acceleration inputs.
- `Detumbling`
  - detumbling control phase using magnetorquers.
- `Slew & Pointing`
  - pointing control phase for attitude acquisition and tracking.
- `Attitude errors`
  - evaluates attitude error between current and desired orientation.
- `Pointing error`
  - computes pointing error output.
- `Filter`
  - sensor and state estimation filtering.
- `A_LN`
  - auxiliary navigation matrix computations.
- `Inetial_to_Body`
  - inertial-to-body frame transformation.
- `Norm`
  - vector and error normalization modules.
- `Flag`
  - logic and phase management for mode transitions.

## Key elements

### Spacecraft and mission
- large spacecraft with a symmetric cubic body and deployable solar panels
- total mass around 510 kg
- mission scenario derived from Earth observation applications in Sun-synchronous LEO

### Sensors
- gyroscope: high frequency sampling with ARW noise and bias modeling
- sun sensors: six sun sensors used for attitude determination
- magnetometer: measurements of the Earth magnetic field for attitude estimation

### Actuators
- one reaction wheel for torque generation by momentum exchange
- three magnetorquers for torque generation via the geomagnetic field

### Control strategy
- LQR-based control design for the detumbling phase
- LQR-based control design for slew and pointing
- gains and dynamics tuned for the spacecraft inertia properties

## Notes
- `spacecraft_attitude_model.slxc` is a Simulink cache file and is not required for model editing.
- The project assumes the support functions `SCModel.m`, `IGRF2020MagneticField.m`, and `AtmosphereData.m` are available in the MATLAB path.
