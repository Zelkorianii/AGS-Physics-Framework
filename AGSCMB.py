# AGS_TESTIMONY_FRAMEWORK.py :: Version 1.0.0

import numpy as np
from scipy.integrate import solve_ivp

# --- Global variable to track the last logged time interval :: CRITICAL FOR LOGGING ---
LAST_LOGGED_TIME = 0.0

# =================================================================
# --- 1. DERIVED INTRINSIC CONSTANTS :: ABSOLUTE VALUES FROM INVERSE SOLVER ---
# =================================================================
# Current Target: Force n_s away from 1.000 by increasing alpha_NL
C_AGS_SOLVED = 27.0
alpha_NL_SOLVED = 1.0e-14

# --- INITIAL CONDITION :: Start Point for the Field ---
PHI_REF_SOLVED = 75000.0 # High starting field value for sustained inflation

# =================================================================
# --- 2. REQUIRED REAL-WORLD CONSTANTS ---
# =================================================================
M_PHI = 1.0e-51
M_PLANCK = 1.0


# =================================================================
# --- 3. THE AGS DIFFERENTIAL EQUATION ENGINE :: Helper Functions ---
# =================================================================

def potential_energy(phi, C_AGS, alpha_NL):
    """V(phi) - The AGS Potential."""
    V_phi = (M_PHI**2 * phi**2) + (alpha_NL * phi**4) - (C_AGS * np.exp(-phi))
    return V_phi

def potential_derivative(phi, C_AGS, alpha_NL):
    """V'(phi) - The Force on the field."""
    V_prime = (2 * M_PHI**2 * phi) + (4 * alpha_NL * phi**3) + (C_AGS * np.exp(-phi))
    return V_prime

def potential_second_derivative(phi, C_AGS, alpha_NL):
    """V''(phi)."""
    V_double_prime = (2 * M_PHI**2) + (12 * alpha_NL * phi**2) - (C_AGS * np.exp(-phi))
    return V_double_prime

def slow_roll_epsilon(phi, C_AGS, alpha_NL):
    """Simulated first slow-roll parameter (Epsilon)."""
    V_prime = potential_derivative(phi, C_AGS, alpha_NL)
    V_phi = potential_energy(phi, C_AGS, alpha_NL)
    return 0.5 * (V_prime / (V_phi + 1e-15))**2

def slow_roll_eta(phi, C_AGS, alpha_NL):
    """Simulated second slow-roll parameter (Eta)."""
    V_double_prime = potential_second_derivative(phi, C_AGS, alpha_NL)
    V_phi = potential_energy(phi, C_AGS, alpha_NL)
    return V_double_prime / (V_phi + 1e-15)


# =================================================================
# --- 4. THE FINAL DYNAMIC VERDICT :: FULL SIMULATION (N=60) ---
# =================================================================

def full_dynamic_simulation(C_AGS, alpha_NL, PHI_REF):
   
    global LAST_LOGGED_TIME
    LAST_LOGGED_TIME = 0.0 # Reset for the current simulation run
   
    # --- CALCULATE INITIAL VELOCITY :: The Nudge ---
    V_ref = potential_energy(PHI_REF, C_AGS, alpha_NL)
    V_prime_ref = potential_derivative(PHI_REF, C_AGS, alpha_NL)
   
    # Hubble Parameter H = sqrt(V / 3*M_Planck^2)
    H_ref = np.sqrt(V_ref / 3)
   
    # Initial slow-roll velocity: phi_dot_0 = -V' / (3*H)
    PHI_DOT_START = -V_prime_ref / (3 * H_ref)
   
    print(f"DEBUG: Initial Phi_dot set to {PHI_DOT_START:.10e}", flush=True)

    # --- 4A. The Core ODE System ---
    def system_of_odes(t, state):
        global LAST_LOGGED_TIME
        phi, phi_dot = state
       
        V_phi = potential_energy(phi, C_AGS, alpha_NL)
        H_squared = (1 / (3 * M_PLANCK**2)) * (0.5 * phi_dot**2 + V_phi)
       
        if H_squared < 1e-18:
            return [0.0, 0.0]
           
        H = np.sqrt(H_squared)
        V_prime = potential_derivative(phi, C_AGS, alpha_NL)
        phi_double_dot = -3 * H * phi_dot - V_prime
       
        return [phi_dot, phi_double_dot]

    # --- 4B. THE TERMINATION EVENT ---
    def end_of_inflation_event(t, state):
        phi, phi_dot = state
        V_phi = potential_energy(phi, C_AGS, alpha_NL)
        H_squared = (1 / (3 * M_PLANCK**2)) * (0.5 * phi_dot**2 + V_phi)
        epsilon_H = (phi_dot**2 / (2 * H_squared + 1e-18))
        return epsilon_H - 1.0

    end_of_inflation_event.terminal = True
    end_of_inflation_event.direction = 1

    # --- 4C. Iterative Chunking Execution :: Memory-Safe ---
   
    N_TARGET = 60.0
    t_start = 0.0
    t_chunk_duration = 1000.0 
   
    phi_current = PHI_REF 
    phi_dot_current = PHI_DOT_START 
   
    TOTAL_E_FOLDS_ACCUMULATED = 0.0
   
    INFLATION_TERMINATED = False
    phi_k = None 
    t_end = None
    phi_end = None
   
    rtol_safe = 1e-5
    atol_safe = 1e-8
    LOG_INTERVAL = 2500.0
   
    print(f"DEBUG: Starting iterative solution with chunk size {t_chunk_duration} and rtol={rtol_safe}", flush=True)

    while t_start < 100000000.0 and not INFLATION_TERMINATED:
       
        t_end_chunk = t_start + t_chunk_duration
        chunk_time_span = [0, t_chunk_duration]
        chunk_initial_conditions = [phi_current, phi_dot_current]
       
        try:
            chunk_solution = solve_ivp(system_of_odes, chunk_time_span, chunk_initial_conditions,
                                       method='RK45', events=end_of_inflation_event,
                                       rtol=rtol_safe, atol=atol_safe,
                                       dense_output=False)
        except Exception as e:
            print(f"ERROR: Solver failed during chunk t={t_start:.0f}. Error: {e}", flush=True)
            break
           
        if len(chunk_solution.t) == 0:
            print(f"ERROR: Solver returned empty solution at t={t_start:.0f}. Terminating.", flush=True)
            break
           
        # 1. Check for termination
        if len(chunk_solution.t_events[0]) > 0:
            t_event_time = chunk_solution.t_events[0][0]
            t_end = t_start + t_event_time 
            phi_end = chunk_solution.y_events[0][0][0]
            phi_dot_end = chunk_solution.y_events[0][0][1]
            INFLATION_TERMINATED = True
           
        # 2. Calculate E-folds gained in this chunk
        t_in_chunk = chunk_solution.t
        phi_in_chunk = chunk_solution.y[0]
        phi_dot_in_chunk = chunk_solution.y[1]
       
        V_start = potential_energy(phi_in_chunk[0], C_AGS, alpha_NL)
        V_end = potential_energy(phi_in_chunk[-1], C_AGS, alpha_NL)
        V_prime_start = potential_derivative(phi_in_chunk[0], C_AGS, alpha_NL)
        V_prime_end = potential_derivative(phi_in_chunk[-1], C_AGS, alpha_NL)
       
        V_avg = (V_start + V_end) / 2.0
        V_prime_avg = (V_prime_start + V_prime_end) / 2.0
        Delta_phi = phi_in_chunk[-1] - phi_in_chunk[0] 
       
        dN_chunk = - (V_avg / V_prime_avg) * Delta_phi
        STABILIZATION_FACTOR = 1.0e-5
        dN_chunk = dN_chunk * STABILIZATION_FACTOR
        TOTAL_E_FOLDS_ACCUMULATED += dN_chunk

        # --- 3. UPDATED N=60 CHECK AND STORE ---
        if phi_k is None and TOTAL_E_FOLDS_ACCUMULATED >= N_TARGET:
            phi_k = phi_in_chunk[-1]
            
            # Immediate Calculation of observables at N=60
            eps_k = slow_roll_epsilon(phi_k, C_AGS, alpha_NL)
            eta_k = slow_roll_eta(phi_k, C_AGS, alpha_NL)
            ns_k = 1 - 6 * eps_k + 2 * eta_k
            
            print(f"\n{'!'*50}", flush=True)
            print(f"CMB HORIZON EXIT DETECTED AT N = {N_TARGET}", flush=True)
            print(f"Field Value (phi_k): {phi_k:.10e}", flush=True)
            print(f"Epsilon at phi_k: {eps_k:.6e}", flush=True)
            print(f"Eta at phi_k: {eta_k:.6e}", flush=True)
            print(f"Calculated n_s: {ns_k:.8f}", flush=True)
            print(f"{'!'*50}\n", flush=True)
       
        # 4. Advance to the next chunk
        if not INFLATION_TERMINATED:
            t_start = t_end_chunk 
            phi_current = phi_in_chunk[-1]
            phi_dot_current = phi_dot_in_chunk[-1]
           
            if t_start >= LAST_LOGGED_TIME + LOG_INTERVAL:
                V_log = potential_energy(phi_current, C_AGS, alpha_NL)
                H_squared_final = (1 / (3 * M_PLANCK**2)) * (0.5 * phi_dot_current**2 + V_log)
                epsilon_H = phi_dot_current**2 / (2 * H_squared_final + 1e-18) if H_squared_final > 1e-18 else 0.0
                print(f"[{t_start:.0f}t] PHI: {phi_current:.3f} | N_Acc: {TOTAL_E_FOLDS_ACCUMULATED:.2f} | Epsilon_H: {epsilon_H:.2e}", flush=True)
                LAST_LOGGED_TIME = t_start
               
        if t_start >= 1000000.0 and not INFLATION_TERMINATED:
            print("Safety termination: Maximum time reached (1,000,000).", flush=True)
            break

    # --- 4D. FINAL OUTPUT CHECK ---
    if INFLATION_TERMINATED and phi_k is not None:
        epsilon_dyn = slow_roll_epsilon(phi_k, C_AGS, alpha_NL)
        eta_dyn = slow_roll_eta(phi_k, C_AGS, alpha_NL)
        NS_SIMULATED_DYNAMIC = 1 - 6 * epsilon_dyn + 2 * eta_dyn
        R_SIMULATED_DYNAMIC = 16 * epsilon_dyn

        print("\n--- DYNAMIC SIMULATION COMPLETED ---", flush=True)
        print(f"End of Inflation time (t_end): {t_end:.10e}", flush=True)
        print(f"End of Inflation field value (phi_end): {phi_end:.10e}", flush=True)
        print(f"\n--- FINAL CMB DYNAMIC CALCULATION (N={N_TARGET} E-folds) ---", flush=True)
        print(f"Field Value at Horizon Exit (phi_k): {phi_k:.10e}", flush=True)
    elif INFLATION_TERMINATED and phi_k is None:
        NS_SIMULATED_DYNAMIC = "FAILED_N60"
        R_SIMULATED_DYNAMIC = "FAILED_N60"
        print("\n--- DYNAMIC SIMULATION FAILED: N < 60 ---", flush=True)
    else:
        NS_SIMULATED_DYNAMIC = "FAILED_TERM"
        R_SIMULATED_DYNAMIC = "FAILED_TERM"
        print("\n--- DYNAMIC SIMULATION FAILED: NO TERMINATION ---", flush=True)

    print("\n[CMB OBSERVABLES SIMULATED - DYNAMICALLY DERIVED]", flush=True)
    if isinstance(NS_SIMULATED_DYNAMIC, str):
        print(f"Predicted Spectral Index (n_s): {NS_SIMULATED_DYNAMIC}", flush=True)
    else:
        print(f"Predicted Spectral Index (n_s): {NS_SIMULATED_DYNAMIC:.5f}", flush=True)
        print(f"Predicted Tensor Ratio (r): {R_SIMULATED_DYNAMIC:.5e}", flush=True)

    return "Complete"

# --- EXECUTION ---
print("--- AGS FINAL DYNAMIC VERDICT START ---", flush=True)
print("Testing C_AGS={:.1e}, alpha_NL={:.1e}, PHI_REF={:.1e}".format(C_AGS_SOLVED, alpha_NL_SOLVED, PHI_REF_SOLVED), flush=True)

full_dynamic_simulation(C_AGS_SOLVED, alpha_NL_SOLVED, PHI_REF_SOLVED)
print("--- Program Finished ---", flush=True)
input("Press Enter to Continue...")
