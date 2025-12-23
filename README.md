AGS (Analog Ground State) Physics Framework: Cosmological Inverse Solvers
Project Overview
The AGS Physics Framework is a suite of numerical tools designed to bridge the gap between theoretical high-energy physics and observational cosmology.
This framework specializes in the AGS Potential—a scalar field model designed to explain cosmic inflation while remaining consistent with
Cosmic Microwave Background (CMB) data from the Planck Satellite.

The framework uses a "closed-loop" logic:
Dynamic Simulation: Solving the equations of motion for the early universe.

Inverse Optimization: Back-calculating the fundamental constants of the universe from observed data.

Core Components1. AGS_TESTIMONY_FRAMEWORK (AGSCMB.py)
This is the Dynamic Engine. It simulates the physical evolution of the scalar field (the Inflaton) over time using a memory-safe iterative ODE solver.

Differential Equation Engine: Utilizes scipy.integrate.solve_ivp (Runge-Kutta 45) to track field trajectory ($\phi$) and velocity ($\dot{\phi}$).

Automated Termination: Implements a high-precision event handler to detect the exact moment inflation ends ($\epsilon \approx 1$).

Observable Extraction: Calculates the Spectral Index ($n_s$) and Tensor-to-Scalar ratio ($r$) specifically at the $N=60$ e-folds horizon exit point.

:: 2. AGS Inverter (AGSInverter.py) ::
This is the Optimization Engine. It treats the fundamental constants of the AGS model ($C_{AGS}$ and $\alpha_{NL}$) as unknowns in an inverse problem.

L-BFGS-B Optimization: Employs a multi-dimensional minimization algorithm to find the specific constants that yield a universe matching our actual
cosmological observations.

Weighted Cost Function: Minimizes the squared error between model predictions and the Planck 2018 targets ($n_s = 0.965$).

Verification Module: Includes a post-optimization check that re-calculates slow-roll parameters to ensure physical consistency.

Technical Specifications
Language: Python 3.8
Numerical Libraries: NumPy for vector operations, SciPy for ODE integration and non-linear optimization.
Physics Framework: Standard Slow-Roll Inflationary Theory in Reduced Planck Units ($M_{Pl} = 1$).

-- How to Run --

MUST HAVE THESE LIBRARIES INSTALLED! :: numpy, scipy.
Bash command line:
pip install numpy
pip install scipy

Derive Constants: Run AGSInverter.py to find the optimized $C_{AGS}$ and $\alpha_{NL}$ constants for a specific field exit value ($\phi_k$).

Simulate Universe: Plug those derived constants into AGSCMB.py to run a full dynamic simulation and verify the e-fold accumulation ($N \ge 60$).
Research ImplicationsThis framework demonstrates that the AGS Potential can successfully produce a "Flat" universe with the correct spectral tilt,
resolving the fine-tuning problem by deriving constants directly from observational constraints.
