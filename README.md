# Optimal Racing Line Optimization (F1)

This project computes the optimal racing line for a generic F1 car on a variation of the Nürburgring track using convex optimization (Bicycle Model).

## Overview

The optimization process is divided into three logical phases to ensure convergence and physical accuracy:

1.  **Phase V1 (Geometric):** Shortest path optimization minimizing curvature.
2.  **Phase V2 (Mechanical):** Optimization including vehicle dynamics (friction circle) but without aerodynamic downforce.
3.  **Phase V3 (Aerodynamic):** Full optimization including speed-dependent downforce (higher grip at high speeds).
4.  **Phase V4 (Grid Search):** Refinement of weight parameters for marginal gains.

## Analysis

The project compares the three phases to demonstrate:
- The difference between the geometric shortest path and the time-optimal path.
- The impact of aerodynamic downforce on cornering speeds and braking points.
- The "G-G Diagram" usage (Friction Circle vs Diamond constraints).

## Requirements

- MATLAB (R2020b or later recommended)
- **YALMIP** Toolbox
- **Gurobi** (or any QP/SOCP solver compatible with YALMIP)

## Usage

1. Ensure the `Circuits_Data` folder contains the track CSV file (`Nuerburgring_track.csv`).
2. Run `main_f1_optimization.m`.
3. Results are saved to `Comparison_Data.mat`.

## Authors

Project developed for Path Planning & Trajectory Analysis.
