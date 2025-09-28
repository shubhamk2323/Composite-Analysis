This repository contains a Java program for macromechanical analysis of composite laminates, 
including computation of laminate stiffness matrices (A, B, D), transformation of stresses and strains, 
and evaluation of ply-level stress and strain under applied loads.

## Features

- Compute **compliance (S)** and **stiffness (Q) matrices** for orthotropic materials.
- Transform material properties based on **fiber orientation**.
- Generate **laminate stiffness matrices**:  
  - **A matrix**: In-plane stiffness  
  - **B matrix**: Coupling stiffness  
  - **D matrix**: Bending stiffness  
- Compute **ABD matrix** and its inverse for full laminate analysis.
- Calculate **apparent laminate properties**: Ex, Ey, Gxy, bending stiffnesses.
- Perform **ply-level stress-strain calculations**.
- Transform stresses and strains between **material axes (1-2)** and **global axes (X-Y)**.
- Evaluate **allowable strengths** in global coordinates.

---

## How to Run

1. Clone the repository:
2. Run the Macromechanics.java file and Micromechanics will  be calculated from this file only
