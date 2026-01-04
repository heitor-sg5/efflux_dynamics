import numpy as np
from scipy.integrate import solve_ivp

def growth_rate(A, p):
    return p['G_max'] * (1 / (1 + (A / p['IC50'])**p['h']))

def ode_system(t, y, p):
    A, P, M, Q = y

    E = p['E_h'] + p['beta'] * Q
    G = growth_rate(A, p)

    dA = (p["k_in"] * p["A_ext"]) - (p["k_out"] * E * A)
    dP = (p["r_P"] * G * P) - (p["c_P"] * P) - (p["gamma"] * P**2)
    dM = (p["k_m"] * P) - (p["delta_m"] * M)
    dQ = (p["k_q"] * M) - (p["delta_q"] * Q)

    return [dA, dP, dM, dQ]

def run_odes(p, TMAX, sc):
    results = {}

    y0 = [0.0, 1.0, 0.0, 0.0]
    t_eval = np.linspace(0, TMAX, 1000)

    sol = solve_ivp(lambda t, y: ode_system(t, y, p), [0, TMAX], y0, t_eval=t_eval)
    results[sc["name"]] = {"solution": sol, "params": sc}

    return results