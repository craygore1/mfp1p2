# mfp1p2
MATLAB implementation of the mfp1p2 multifractal model.
# README

## Overview
This repository contains two pieces of MATLAB code that generate multifractal patterns using a probabilistic subdivision approach based on two probabilities, `p1` and `p2`. These patterns are created by iteratively subdividing a matrix and applying different probabilities to each region. The first script also includes a theoretical multifractal analysis of the generated pattern. The second code is a function that allows for flexible inputs to generate similar patterns but does not include the analysis.

### Files:
1. **Scripts:** `mfp_script.m` (first code block)
2. **Function:** `mfp1p2gen.m` (second code block)

---

## 1. Scripts: `mfp1p2.m`, `mfp1p2circle.m`, `mfp1p2radial.m`

### Description
These scripts generate a binary matrix representing a multifractal pattern using the mfp1p2 probabilistic model. The matrix is subdivided recursively into smaller blocks, and each block is assigned values based on the probabilities `p1` and `p2`. After generating the pattern, the script performs a **multifractal analysis** by calculating the generalized fractal dimensions \( D(q) \), the multifractal spectrum \( \tau(q) \), the singularity strength \( \alpha \), and the multifractal spectrum \( f(\alpha) \).

The variations `mfp1p2circle.m` and `mfp1p2radial.m` generate the multifractal pattern on a circular support. `mfp1p2circle.m` only subdivides by angle while `mfp1p2radial.m` subdivides by both radius and angle.

---

## 2. Function: `mfp1p2gen.m`

### Description
This is a more generalized version of the multifractal pattern generation code. It is encapsulated in a function that allows the user to specify the matrix size (`n`), the probabilities (`prob1` and `prob2`), and the number of iterations (`iterations`). The function generates a binary matrix based on the same probabilistic subdivision method on a rectangular support as the script but does **not** perform any multifractal analysis or visualization.

### Key Features:
- **Flexible Input Parameters:**
  - `n`: Size of the matrix (nxn).
  - `prob1`: First probability used in the model.
  - `prob2`: Second probability used in the model.
  - `iterations`: Number of subdivision iterations to apply.
- **Output:** Returns a binary matrix `A` representing the generated multifractal pattern.
- **No Analysis/Visualization:** The function focuses solely on matrix generation, leaving visualization and analysis to the user.

### Usage:
To use the function, call it with the desired input parameters:
```matlab
A = mfp1p2gen(n, prob1, prob2, iterations);
```
For example:
```matlab
A = mfp1p2gen(2048, 1, 0.75, 10);
imshow(A)
```
This will generate a 2048x2048 matrix with 10 iterations using probabilities 1 and 0.75, and display the result.

## Requirements
- **MATLAB**: Both the scripts and function require MATLAB to run. No additional toolboxes are needed.

---

## License
These scripts are provided as-is for educational purposes. Feel free to modify and extend them for your own use.

---

## Contact
If you have any questions or encounter any issues, feel free to reach out via the repository's contact information or open an issue.
