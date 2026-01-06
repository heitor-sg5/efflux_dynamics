import numpy as np
from scipy.integrate import solve_ivp

def a1_(p, A, P):
    if P <= 0:
        return 0.0
    G = p["G_max"] * (1 / (1 + (A / p["IC50"]) ** p["h"]))
    return p["r_P"] * G * P

def a2_(p, P):
    if P <= 0:
        return 0.0
    return (p["c_P"] * P) + (p["gamma"] * P ** 2)

def gillespie_step(p, state):
    A, P, M, Q = state

    a1 = a1_(p, A, P)
    a2 = a2_(p, P)
    a0 = a1 + a2

    if a0 <= 0:
        return None, None

    r1, r2 = np.random.rand(), np.random.rand()
    tau = -np.log(r1) / a0

    if r2 < a1 / a0:
        event = "replication"
    else:
        event = "loss"

    return tau, event

def run_pdmp(p, TMAX, ode_model):
    t = 0.0
    A, P, M, Q = 0.0, 1, 0.0, 0.0

    times = [t]
    traj = [[A], [P], [M], [Q]]

    while t < TMAX and P > 0:
        tau, event = gillespie_step(p, (A, P, M, Q))

        if tau is None:
            break

        t_next = min(t + tau, TMAX)
        sol = solve_ivp(lambda tt, yy: ode_model.ode_system(tt, yy, p), [t, t_next], [A, P, M, Q], t_eval=[t_next])

        A = sol.y[0, -1]
        M = sol.y[2, -1]
        Q = sol.y[3, -1]

        if event == "replication":
            P += 1
        elif event == "loss":
            P -= 1

        t = t_next

        times.append(t)
        traj[0].append(A)
        traj[1].append(P)
        traj[2].append(M)
        traj[3].append(Q)

    return {
        "t": np.array(times),
        "A": np.array(traj[0]),
        "P": np.array(traj[1]),
        "M": np.array(traj[2]),
        "Q": np.array(traj[3]),
    }

def run_pdmp_scenarios(base_params, TMAX, steps, ode_model):
    results = {}
    E_A_values = [(0.0, 0.0), (15.0, 0.0), (15.0, 2.0)]
    t_grid = np.linspace(0, TMAX, 500)

    for a, e in E_A_values:
        A_runs = np.zeros((steps, len(t_grid)))
        P_runs = np.zeros((steps, len(t_grid)))
        M_runs = np.zeros((steps, len(t_grid)))
        Q_runs = np.zeros((steps, len(t_grid)))

        params = base_params.copy()
        params['E_h'] = e
        params['A_ext'] = a

        for run_idx in range(steps):
            traj = run_pdmp(params, TMAX, ode_model)
            t = traj['t']
            A_runs[run_idx, :] = np.interp(t_grid, t, traj['A'])
            P_runs[run_idx, :] = np.interp(t_grid, t, traj['P'])
            M_runs[run_idx, :] = np.interp(t_grid, t, traj['M'])
            Q_runs[run_idx, :] = np.interp(t_grid, t, traj['Q'])

        results[(a, e)] = {
            't': t_grid,
            'A': A_runs,
            'P': P_runs,
            'M': M_runs,
            'Q': Q_runs
        }

    return results

def extinction_probability_pdmp(scenario_data):
    P_final = scenario_data['P'][:, -1]
    extinct = (P_final == 0)
    return extinct.mean()

def rescue_probability_pdmp(scenario_data):
    return 1.0 - extinction_probability_pdmp(scenario_data)

def analyse_ssa_results(trajectories):
    E_A_values = [(0.0, 0.0), (15.0, 0.0), (15.0, 2.0)]
    for a, e in E_A_values:
        scenario = trajectories[(a, e)]
        P_ext = extinction_probability_pdmp(scenario)
        P_rescue = rescue_probability_pdmp(scenario)
        print(f"\nScenario: Eh={e}, Aext={a}")
        print(f"Probability of plasmid extinction (Pext): {P_ext:.3f}")
        print(f"Probability of resistance rescue (Prescue): {P_rescue:.3f}")