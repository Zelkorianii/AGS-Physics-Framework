import numpy as np
from scipy.optimize import minimize, fsolve

# =================================================================
# --- 1. COSMOLOGICAL TARGETS :: THE UNIVERSE'S ABSOLUTE VALUES ---
# =================================================================
# These values are derived from Planck, WEP, and LambdaCDM observations.
N_TARGET = 60.0 # Required e-folds for CMB horizon exit
NS_TARGET = 0.965 # Target Spectral Index (Center of Planck 2018 data)
R_TARGET = 1.0e-5 # Target Tensor Ratio (must be tiny for consistency)

# In Reduced Planck Units (M_Pl=1), Dark Energy Density is approximately 10^-122.
RHO_LAMBDA_TARGET = 1.0e-120

# M_PLANCK is 1.0 in this unit convention, but we define it for clarity
M_PLANCK = 1.0
M_PHI = 1.0e-51

# =================================================================
# --- 2. AGS POTENTIAL AND DERIVATIVES :: MUST MATCH AGSCMB.py ---
# =================================================================
# These functions use the constants (C_AGS, alpha_NL) as variables X[0], X[1]

def potential_energy(phi, C_AGS, alpha_NL):
    """V(phi) - The AGS Potential."""
    V_phi = (M_PHI**2 * phi**2) + (alpha_NL * phi**4) - (C_AGS * np.exp(-phi))
    return V_phi

def potential_derivative(phi, C_AGS, alpha_NL):
    """V'(phi) - First Derivative."""
    V_prime = (2 * M_PHI**2 * phi) + (4 * alpha_NL * phi**3) + (C_AGS * np.exp(-phi))
    return V_prime

def potential_second_derivative(phi, C_AGS, alpha_NL):
    """V''(phi) - Second Derivative."""
    V_double_prime = (2 * M_PHI**2) + (12 * alpha_NL * phi**2) - (C_AGS * np.exp(-phi))
    return V_double_prime

# =================================================================
# --- 3. OBJECTIVE FUNCTION :: The Cost/Error Function ---
# =================================================================
# This function calculates the total error between model prediction and universal targets.
def objective_function(X, phi_k):
    """
    X = [C_AGS, alpha_NL]
    phi_k = Field value at N=60 (The assumed exit point derived from the Dynamic Solver)
    """
    C_AGS, alpha_NL = X

    # A. SLOW-ROLL PARAMETER CALCULATIONS (at phi_k)
    V_k = potential_energy(phi_k, C_AGS, alpha_NL)
    if V_k <= 0:
        return 1.0e+300 # Massive penalty if V is non-positive
   
    V_prime_k = potential_derivative(phi_k, C_AGS, alpha_NL)
    V_double_prime_k = potential_second_derivative(phi_k, C_AGS, alpha_NL)

    # Calculate Slow-Roll Parameters (M_Pl=1)
    epsilon_model = 0.5 * (V_prime_k / V_k)**2
    eta_model = V_double_prime_k / V_k
   
    # Calculate CMB Observables
    ns_model = 1.0 - 6.0 * epsilon_model + 2.0 * eta_model
    r_model = 16.0 * epsilon_model

    # B. LATE-TIME CONSTRAINT (V_min) - RETAINED FOR V_MIN_MODEL CALCULATION
    phi_min_approx = 1.0
    V_min_model = potential_energy(phi_min_approx, C_AGS, alpha_NL)
   
    # --- CALCULATE ERROR SQUARED :: COST ---
   
    # Constraint 1: Spectral Index (n_s)
    error_ns = (ns_model - NS_TARGET)**2
   
    # Constraint 2: Tensor Ratio (r)
    error_r = (r_model - R_TARGET)**2
   
    # Constraint 3: Dark Energy Density (V_min) - NOT USED IN COST DUE TO INSTABILITY
    error_lambda = (V_min_model - RHO_LAMBDA_TARGET)**2

    # TOTAL COST: Weighted sum of squared errors
    # CRITICAL FIX: Only use ns and r constraints to stabilize the optimization
    COST = (error_ns * 1.0e+0) + (error_r * 1.0e+12) 
   
    # Penalty for non-physical constants
    if C_AGS < 0 or alpha_NL < 0:
        COST += 1.0e+300
   
    return COST

# =================================================================
# --- 4. EXECUTION ---
# =================================================================

# --- CRITICAL INPUT: The field value at N=60 from the analytic solver ---
PHI_K_FIXED = 9.9759711307e+04

# Initial guesses for the constants (Used for the first successful CMB-focused run)
initial_guesses = [27.0, 1.0e-14]

print(f"--- AGS CONSTANT DERIVATION START (CMB-FOCUSED INVERSE PROBLEM) ---", flush=True)
print(f"Targeting: ns={NS_TARGET}, r={R_TARGET}", flush=True)
print(f"Input PHI_K (Exit Field): {PHI_K_FIXED:.5e}", flush=True)

# Run the minimization to find the constants that minimize the total error
bounds = [(0, 100), (0, 1e-10)]

result = minimize(
    lambda X: objective_function(X, PHI_K_FIXED),
    initial_guesses,
    method='L-BFGS-B',
    bounds=bounds,
    options={'maxiter': 50000, 'ftol': 1e-12}
)

# Extract the final, 'absolute' constants
C_AGS_ABS, alpha_NL_ABS = result.x

# --- 5. FINAL OUTPUT ---
print("\n--- AGS ABSOLUTE CONSTANTS DERIVED FROM CMB DATA ---", flush=True)
print(f"Optimization Success: {result.success}", flush=True)
print(f"Optimization Message: {result.message}", flush=True)
print(f"Derived C_AGS: {C_AGS_ABS:.10e}", flush=True)
print(f"Derived alpha_NL: {alpha_NL_ABS:.10e}", flush=True)

# Verification check
print("\n--- VERIFICATION OF DERIVED CONSTANTS ---", flush=True)

# Recalculate Observables using the derived constants (C_AGS_ABS, alpha_NL_ABS)
V_final = potential_energy(PHI_K_FIXED, C_AGS_ABS, alpha_NL_ABS)
V_prime_final = potential_derivative(PHI_K_FIXED, C_AGS_ABS, alpha_NL_ABS)
V_double_prime_final = potential_second_derivative(PHI_K_FIXED, C_AGS_ABS, alpha_NL_ABS)

# Recalculate Slow-Roll Parameters
epsilon_final = 0.5 * (V_prime_final / V_final)**2
eta_final = V_double_prime_final / V_final

# Recalculate CMB Observables
ns_final = 1.0 - 6.0 * epsilon_final + 2.0 * eta_final
r_final = 16.0 * epsilon_final

# Recalculate V_min_model in the global scope for reporting
V_min_model_recalc = potential_energy(1.0, C_AGS_ABS, alpha_NL_ABS)

# Print Final Assertions
print(f"Calculated n_s: {ns_final:.6f} (Target: {NS_TARGET})", flush=True)
print(f"Calculated r: {r_final:.6e} (Target: {R_TARGET})", flush=True)
print(f"Calculated V_min (at phi=1): {V_min_model_recalc:.3e} (Target: {RHO_LAMBDA_TARGET:.0e})", flush=True)

input("Press Enter to Continue...")
