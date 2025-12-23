"""
AGS Simulation: Planck 2018 Calibrated Inflation
Version 1.2.0
"""

import numpy as np
from scipy.integrate import solve_ivp
import csv

# =================================================================
# --- 1. THE "GOLDILOCKS" CONFIGURATION ---
# =================================================================
# Calibrated for n_s = 0.9650 and r = 0.00358
C_AGS = 1.2e-10    
ALPHA_NL = 1.0     
PHI_START = 7.043  # Initial field displacement
M_PLANCK = 1.0     # Reduced Planck units

# =================================================================
# --- 2. POTENTIAL ENGINE & OBSERVABLES ---
# =================================================================

def potential_energy(phi):
    """V(phi) = C * tanh^2(phi / sqrt(6 * alpha))"""
    s6a = np.sqrt(6 * ALPHA_NL)
    return C_AGS * np.tanh(phi / s6a)**2

def potential_derivative(phi):
    """V'(phi)"""
    s6a = np.sqrt(6 * ALPHA_NL)
    return 2 * C_AGS * np.tanh(phi / s6a) * (1 / (s6a * np.cosh(phi / s6a)**2))

def potential_second_derivative(phi):
    """V''(phi)"""
    s6a = np.sqrt(6 * ALPHA_NL)
    sech = 1.0 / np.cosh(phi / s6a)
    return (2 * C_AGS / (6 * ALPHA_NL)) * (sech**4 - 2 * np.tanh(phi / s6a)**2 * sech**2)

def get_ns_r(phi):
    """Calculates n_s and r at a specific field value."""
    V = potential_energy(phi)
    Vp = potential_derivative(phi)
    Vpp = potential_second_derivative(phi)
    
    epsilon = 0.5 * (Vp / (V + 1e-18))**2
    eta = Vpp / (V + 1e-18)
    
    ns = 1 - 6 * epsilon + 2 * eta
    r = 16 * epsilon
    return ns, r

# =================================================================
# --- 3. DYNAMIC CHUNKING SOLVER ---
# =================================================================

def run_simulation(filename="AGS_Goldilocks_Universe.csv"):
    phi_current = PHI_START
    phi_dot_current = 0.0
    t_start = 0.0
    t_chunk = 500.0
    total_efolds = 0.0
    
    inflation_ended = False
    phi_k = None
    ns_final, r_final = 0, 0

    def system_of_odes(t, state):
        phi, phi_dot = state
        V = potential_energy(phi)
        Vp = potential_derivative(phi)
        H = np.sqrt((0.5 * phi_dot**2 + V) / 3.0)
        return [phi_dot, -3 * H * phi_dot - Vp]

    def end_inflation_event(t, state):
        phi, phi_dot = state
        V = potential_energy(phi)
        H_sq = (0.5 * phi_dot**2 + V) / 3.0
        return (phi_dot**2 / (2 * H_sq + 1e-18)) - 1.0
    
    end_inflation_event.terminal = True

    print(f"--- AGS SIMULATION START ---")
    print(f"Parameters: Phi={PHI_START}, Alpha={ALPHA_NL}, C={C_AGS}")

    with open(filename, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Time", "Phi", "Phi_Dot", "Accumulated_N", "n_s", "r"])

        while t_start < 1e7 and not inflation_ended:
            sol = solve_ivp(system_of_odes, [0, t_chunk], [phi_current, phi_dot_current], 
                            events=end_inflation_event, rtol=1e-8, atol=1e-10)
            
            # Export trajectory points
            for i in range(len(sol.t)):
                v_val = potential_energy(sol.y[0][i])
                vp_val = potential_derivative(sol.y[0][i])
                delta_phi = sol.y[0][i] - sol.y[0][0]
                local_dN = -(v_val / (vp_val + 1e-18)) * delta_phi
                
                ns, r = get_ns_r(sol.y[0][i])
                writer.writerow([t_start + sol.t[i], sol.y[0][i], sol.y[1][i], 
                                 total_efolds + local_dN, ns, r])

            # Update master accumulator
            v_avg = potential_energy(sol.y[0].mean())
            vp_avg = potential_derivative(sol.y[0].mean())
            chunk_delta_phi = sol.y[0][-1] - sol.y[0][0]
            total_efolds += -(v_avg / (vp_avg + 1e-18)) * chunk_delta_phi

            # Detect and capture N=60 observables
            if phi_k is None and total_efolds >= 60.0:
                phi_k = sol.y[0][-1]
                ns_final, r_final = get_ns_r(phi_k)

            if len(sol.t_events[0]) > 0:
                inflation_ended = True
                break
            
            t_start += t_chunk
            phi_current, phi_dot_current = sol.y[0][-1], sol.y[1][-1]

    print(f"\n--- SIMULATION COMPLETE ---")
    print(f"Data Exported: {filename}")
    print(f"CMB Result (N=60): n_s = {ns_final:.5f}, r = {r_final:.5e}")

if __name__ == "__main__":
    run_simulation()
input("Press Enter to Continue...")
