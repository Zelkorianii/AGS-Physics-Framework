# AGS_TESTIMONY_FRAMEWORK_INTEGRATED.py :: Version 1.1.0
import numpy as np
from scipy.integrate import solve_ivp

# --- Global logging variable ---
LAST_LOGGED_TIME = 0.0

# =================================================================
# --- 1. DERIVED INTRINSIC CONSTANTS :: CALIBRATED FOR CMB ---
# =================================================================
# To get n_s = 0.965 and r > 0, we use alpha-attractor scaling
C_AGS_SOLVED = 1.2e-10    
alpha_NL_SOLVED = 1.0     # This value drives 'r' up to ~0.03
PHI_REF_SOLVED = 7.043  # Your high starting field

# REAL-WORLD SCALARS
M_PHI = 1.0e-6
M_PLANCK = 1.0

# =================================================================
# --- 2. THE AGS POTENTIAL ENGINE :: STABLE PLATEAU ---
# =================================================================

def potential_energy(phi, C_AGS, alpha_NL):
    """V(phi) - Merged AGS Plateau Potential."""
    # We use tanh to prevent the exp(-75000) overflow while keeping the AGS tilt
    s6a = np.sqrt(6 * alpha_NL)
    return C_AGS * np.tanh(phi / s6a)**2

def potential_derivative(phi, C_AGS, alpha_NL):
    """V'(phi) - The Force."""
    s6a = np.sqrt(6 * alpha_NL)
    return 2 * C_AGS * np.tanh(phi / s6a) * (1 / (s6a * np.cosh(phi / s6a)**2))

def potential_second_derivative(phi, C_AGS, alpha_NL):
    """V''(phi) - The Curvature."""
    s6a = np.sqrt(6 * alpha_NL)
    sech = 1.0 / np.cosh(phi / s6a)
    return (2 * C_AGS / (6 * alpha_NL)) * (sech**4 - 2 * np.tanh(phi / s6a)**2 * sech**2)

def slow_roll_epsilon(phi, C_AGS, alpha_NL):
    V = potential_energy(phi, C_AGS, alpha_NL)
    Vp = potential_derivative(phi, C_AGS, alpha_NL)
    return 0.5 * (Vp / (V + 1e-18))**2

def slow_roll_eta(phi, C_AGS, alpha_NL):
    V = potential_energy(phi, C_AGS, alpha_NL)
    Vpp = potential_second_derivative(phi, C_AGS, alpha_NL)
    return Vpp / (V + 1e-18)

# =================================================================
# --- 3. FULL DYNAMIC SIMULATION :: YOUR CHUNKING ENGINE ---
# =================================================================

def full_dynamic_simulation(C_AGS, alpha_NL, PHI_REF):
    global LAST_LOGGED_TIME
    LAST_LOGGED_TIME = 0.0 

    # --- INITIAL VELOCITY :: Set to 0 to allow Hubble Friction to take over ---
    PHI_DOT_START = 0.0 
    print(f"DEBUG: Initial Phi_dot set to {PHI_DOT_START:.10e}")

    def system_of_odes(t, state):
        phi, phi_dot = state
        V = potential_energy(phi, C_AGS, alpha_NL)
        Vp = potential_derivative(phi, C_AGS, alpha_NL)
        H = np.sqrt((0.5 * phi_dot**2 + V) / 3.0)
        phi_double_dot = -3 * H * phi_dot - Vp
        return [phi_dot, phi_double_dot]

    def end_of_inflation_event(t, state):
        phi, phi_dot = state
        V = potential_energy(phi, C_AGS, alpha_NL)
        H_sq = (0.5 * phi_dot**2 + V) / 3.0
        epsilon_H = (phi_dot**2 / (2 * H_sq + 1e-18))
        return epsilon_H - 1.0

    end_of_inflation_event.terminal = True
    end_of_inflation_event.direction = 1

    # --- CHUNKING LOGIC ---
    N_TARGET = 60.0
    t_start = 0.0
    t_chunk_duration = 1000.0 
    phi_current, phi_dot_current = PHI_REF, PHI_DOT_START 
    TOTAL_E_FOLDS_ACCUMULATED = 0.0
    INFLATION_TERMINATED = False
    phi_k = None 

    print(f"DEBUG: Starting iterative solution with PHI={PHI_REF}")

    while t_start < 1e7 and not INFLATION_TERMINATED:
        chunk_span = [0, t_chunk_duration]
        sol = solve_ivp(system_of_odes, chunk_span, [phi_current, phi_dot_current],
                        events=end_of_inflation_event, rtol=1e-8, atol=1e-10)
        
        # Calculate E-folds: Delta N = H * Delta t
        # (Using a more precise dN = -V/V' dPhi for the accumulator)
        V_avg = potential_energy(sol.y[0].mean(), C_AGS, alpha_NL)
        Vp_avg = potential_derivative(sol.y[0].mean(), C_AGS, alpha_NL)
        delta_phi = sol.y[0][-1] - sol.y[0][0]
        dN = -(V_avg / (Vp_avg + 1e-18)) * delta_phi
        TOTAL_E_FOLDS_ACCUMULATED += dN

        # Capture Horizon Exit (N=60)
        if phi_k is None and TOTAL_E_FOLDS_ACCUMULATED >= N_TARGET:
            phi_k = sol.y[0][-1]
            eps = slow_roll_epsilon(phi_k, C_AGS, alpha_NL)
            eta = slow_roll_eta(phi_k, C_AGS, alpha_NL)
            print(f"\n!!! CMB HORIZON EXIT AT N=60 !!!")
            print(f"Field Phi_k: {phi_k:.4f} | n_s: {1-6*eps+2*eta:.5f} | r: {16*eps:.5e}\n")

        if len(sol.t_events[0]) > 0:
            INFLATION_TERMINATED = True
            break
        
        t_start += t_chunk_duration
        phi_current, phi_dot_current = sol.y[0][-1], sol.y[1][-1]

    print("--- SIMULATION COMPLETE ---")
    return "Success"

# --- EXECUTION ---
full_dynamic_simulation(C_AGS_SOLVED, alpha_NL_SOLVED, PHI_REF_SOLVED)
input("Press Enter to Continue...")
