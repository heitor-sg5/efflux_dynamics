import numpy as np
from scipy.integrate import solve_ivp

y0 = [0.0, 1.0, 0.0, 0.0]

def plasmid_free_equilibrium(p):
    a = (p["k_in"] * p["A_ext"])
    b = (p["k_out"] * p["E_h"])
    A_star = a / b if b != 0 else 0 
    return A_star

def invasion_eigenvalue(p):
    A_star = plasmid_free_equilibrium(p)
    G_star = growth_rate(A_star, p)
    return (p['r_P'] * G_star) - p['c_P']

def critical_Eh(p):
    a = (p['k_in'] * p['A_ext']) / (p['k_out'] * p['IC50'])
    b = (p['c_P'] / ((p['r_P'] * p['G_max']) - p['c_P']))**(1 / p['h'])
    return a * b

def critical_Aext(p):
    a = (p['k_out'] * p['IC50'] * p['E_h']) / p['k_in']
    b = ((p['r_P'] * p['G_max']) - p['c_P']) / (p['c_P'])**(1 / p['h'])
    return a * b

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

def run_time_series(base_params, TMAX):
    results = {}
    A_E_values=[(0.0, 0.0), (15.0, 0.0), (15.0, 2.0)]
    t_eval = np.linspace(0, TMAX, 500)
    for a, e in A_E_values:
        params = base_params.copy()
        params['E_h'] = e
        params['A_ext'] = a
        sol = solve_ivp(lambda t, y: ode_system(t, y, params), [0, TMAX], y0, t_eval=t_eval)
        results[(a, e)] = sol
    return results

def run_invasion(base_params):
    Eh_vals = np.linspace(0.01, 5.0, 100)
    Aext_vals = np.linspace(0.1, 30.0, 80)
    lambda_grid = np.zeros((len(Eh_vals), len(Aext_vals)))
    for i, e in enumerate(Eh_vals):
        for j, a in enumerate(Aext_vals):
            params = base_params.copy()
            params["E_h"] = e
            params["A_ext"] = a
            lambda_grid[i, j] = invasion_eigenvalue(params)
    return lambda_grid

def run_efflux_sweep(base_params):
    results = {}
    Eh_vals = np.linspace(0.01, 5.0, 100)
    Aext_vals = [5.0, 15.0, 25.0]
    for a in Aext_vals:
        G_star_list = []
        for e in Eh_vals:
            params = base_params.copy()
            params["E_h"] = e
            params["A_ext"] = a
            A_star = plasmid_free_equilibrium(params)
            G_star = growth_rate(A_star, params)
            G_star_list.append(G_star)
        results[a] = {"E_h": np.array(Eh_vals), "G_star": np.array(G_star_list)}
    return results

def analyse_ode_results(base_params):
    A_E_values = [(0.0, 0.0), (15.0, 0.0), (15.0, 2.0)]
    for a, e in A_E_values:
        params = base_params.copy()
        params['E_h'] = e
        params['A_ext'] = a
        A_star = plasmid_free_equilibrium(params)
        lambda_P = invasion_eigenvalue(params)
        Eh_crit = critical_Eh(params)
        Aext_crit = critical_Aext(params)
        print(f"\nScenario: Eh={e}, Aext={a}")
        print(f"Plasmid-free equilibrium (A*, P*, M*, Q*): ({A_star:.1f}, 0, 0, 0)")
        print(f"Plasmid invasion eigenvalue (λP): {lambda_P:.4f}")
        print(f"Critical intrinsic efflux (Eh) = {Eh_crit:.4f}")
        print(f"Critical extracellular drug (Aext) = {Aext_crit:.4f}")