# Narrow-waist Spring Quasi-static Model

This repository contains a MATLAB implementation of a displacement-driven quasi-static model for a cable-driven variable-radius narrow-waist spring. The script compares the predicted tip bending angle with Abaqus finite-element simulation reference values.

## File

- `narrow_waist_spring_quasistatic_model_git.m`  
  Main MATLAB script, including:
  - geometric and material parameter definition;
  - cable-length input cases;
  - nonlinear KKT-system solution;
  - tip bending angle comparison with Abaqus reference data;
  - Curve of relative error in tip bending angle vs. number of segments N;
  - 3D centerline visualization.
  
## Usage Instructions

1. **Clone the repository** and open `narrow_waist_spring_quasistatic_model_git.m` in MATLAB.
2. **Configure Parameters (Optional):**
   * **Geometry:** Modify `geom.r_e`, `geom.r_w`, `geom.delta`, `geom.n_turn`, and `geom.p_per_turn` to match your specific spring design.
   * **Material Properties:** Adjust Young's Modulus (`mat.E`) and Poisson's ratio (`mat.nu`) .
   * **Load Cases:** Define your input cable displacement configurations in the `cases` struct inside Negative values imply the motor is retracting the cable (shortening).
3. **Run the Script:** Execute the script directly from the MATLAB editor or type `narrow_waist_spring_quasistatic_model_git` in the Command Window.

## Citation

If this code is used in a publication, please cite the associated paper:

```text
Shasha Wang, Jialin Zang, Shengyu He, et al. Adaptive soft robot for complex multiple scenes: navigating pipelines, valves, and pressure vessels. Authorea. 19 September 2024.
DOI: https://doi.org/10.22541/au.172674437.72728445/v1
```