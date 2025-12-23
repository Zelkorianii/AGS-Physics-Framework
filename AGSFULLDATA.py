"""
AGS Full Trajectory Simulation: Plateau to Bottom of the Well
Version 1.3.0 - High-Resolution Export
"""

import numpy as np
from scipy.integrate import solve_ivp
import csv

# =================================================================
# --- 1. CONFIGURATION (PLANCK CALIBRATED) ---
# =================================================================
C_AGS = 1.2e-10    
ALPHA_NL = 1.0     
PHI_START = 7.043  
M_PLANCK = 1.0

# =================================================================
# --- 2. POTENTIAL ENGINE ---
# =================================================================

def potential_energy(phi):
    s6a = np.sqrt(6 * ALPHA_NL)
    return C_AGS * np.tanh(phi / s6a)**2

def potential_derivative(phi):
    s6a = np.sqrt(6 * ALPHA_NL)
    # Hyperbolic secant squared for derivative stability
    return 2 * C_AGS * np.tanh(phi / s6a) * (1 / (s6a * np.cosh(phi / s6a)**2))

def get_observables(phi, phi_dot):
    """Calculates physical state of the universe."""
    V = potential_energy(phi)
    Vp = potential_derivative(phi)
    H_sq = (0.5 * phi_dot**2 + V) / 3.0
    
    # Hubble Epsilon (Inflation ends at epsilon_H = 1)
    epsilon_H = (phi_dot**2 / (2 * H_sq + 1e-18))
    return V, epsilon_H

# =================================================================
# --- 3. DYNAMIC HIGH-RES SOLVER ---
# =================================================================

def run_full_trajectory(filename="AGS_Full_Trajectory.csv"):
    phi_val = PHI_START
    phi_dot_val = 0.0
    t_start = 0.0
    
    # Dynamic settings
    # Inflation lasts millions of Planck times; Reheating lasts thousands.
    t_max = 50_000_000  
    chunk_size = 50_000 
    
    print(f"--- STARTING FULL TRAJECTORY SIMULATION ---")
    print(f"Tracking from Plateau (Phi={PHI_START}) to Bottom...")

    def system_of_odes(t, state):
        phi, phi_dot = state
        V = potential_energy(phi)
        Vp = potential_derivative(phi)
        H = np.sqrt((0.5 * phi_dot**2 + V) / 3.0)
        return [phi_dot, -3 * H * phi_dot - Vp]

    with open(filename, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Time", "Phi", "Phi_Dot", "V_phi", "Epsilon_H"])

        # Loop until field settles or time runs out
        while t_start < t_max:
            # Check current state to determine recording density
            V, eps = get_observables(phi_val, phi_dot_val)
            
            # During the "Drop" (eps > 0.01), we increase resolution 10x
            points = 200 if eps < 0.01 else 2000
            t_eval = np.linspace(0, chunk_size, points)
            
            sol = solve_ivp(system_of_odes, [0, chunk_size], [phi_val, phi_dot_val], 
                            t_eval=t_eval, rtol=1e-8, atol=1e-10)
            
            # Write points to CSV
            for i in range(len(sol.t)):
                p, pd = sol.y[0][i], sol.y[1][i]
                v, e_h = get_observables(p, pd)
                writer.writerow([t_start + sol.t[i], p, pd, v, e_h])
            
            phi_val, phi_dot_val = sol.y[0][-1], sol.y[1][-1]
            t_start += chunk_size
            
            # Termination: The field has reached the bottom and stopped oscillating
            if t_start > 1e6 and abs(phi_val) < 1e-5 and abs(phi_dot_val) < 1e-9:
                print(f"Target Reached: Field settled at t = {t_start:.1f}")
                break
                
            # Progress update
            if int(t_start / chunk_size) % 50 == 0:
                print(f"T: {t_start:,.0f} | Phi: {phi_val:.4f} | Epsilon_H: {eps:.2e}")

    print(f"\n--- SIMULATION COMPLETE ---")
    print(f"Total points recorded: (Check {filename})")

if __name__ == "__main__":
    run_full_trajectory()
input("Press Enter to Continue...")
