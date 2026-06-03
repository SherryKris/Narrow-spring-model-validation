# Abaqus Simulation Files for Narrow-waist Spring Validation

This folder contains the Abaqus finite-element simulation files used to validate the quasi-static deformation model of the cable-driven variable-radius narrow-waist spring backbone.

## Purpose

The simulations are used to compare the tip bending angle predicted by the theoretical model with finite-element results under the same prescribed cable displacements. The validation focuses on deformation prediction.


## Model description

- The narrow-waist spring is built with the same geometric and material parameters used in the theoretical model.
- Cable-guide locations are represented by reference points located at the centers of selected spring-wire cross sections.
- The corresponding reference points are coupled to the local wire cross sections to approximate guide rings fixed to the spring.
- The driving cables are routed through the guide locations using sliding-ring constraints.
- The bottom end of the spring is fixed.
- Prescribed displacement boundary conditions are applied to the driving cables so that the cable-length changes are consistent with the theoretical inputs.
- Cable-guide friction is not included, in order to approximate the ideal cable-length constraint assumed in the theoretical model.

## Included spring configurations

The model uses a linear-elastic material with Young's modulus `E = 201.1 GPa` and Poisson's ratio `nu = 0.3`.

## How to run the simulations

1. Open the corresponding `.cae` file in Abaqus/CAE.
2. Check the model parameters, boundary conditions, and prescribed cable displacements for the selected case.
3. Submit the job in Abaqus/CAE.
4. After the job is completed, open the `.odb` file.
5. Extract the deformed configuration and calculate the tip bending angle. The tip bending angle is defined as the angle between the base `z`-axis and the deformed tip `z`-axis.
6. Compare the extracted tip bending angle with the theoretical model prediction reported in the Supplementary Materials.

## Output used for comparison

The main validation metric is: tip bending angle, theta_tip

## Notes

- These simulations are intended to validate the deformation trend of the narrow-waist spring under prescribed cable-length inputs.
- The models do not include cable-guide friction, cable creep, manufacturing imperfections, or detailed local contact effects.
- For compact springs under large deformation, end-induced geometric locking or local contact may occur. In such cases, higher-fidelity contact modeling or experimental calibration may be required.
- Mesh density and solver settings can be refined if higher numerical accuracy is needed.

## Citation

If you use these files, please cite the associated paper:

```text
Shasha Wang, Jialin Zang, Shengyu He, et al. Adaptive soft robot for complex multiple scenes: navigating pipelines, valves, and pressure vessels. Authorea. 19 September 2024.
DOI: https://doi.org/10.22541/au.172674437.72728445/v1
```
