import numpy as np
from scipy.integrate import solve_ivp

from model.deterministic.odes import ode_system
from model.stochastic.gillespie import gillespie_step

def run_pdmp(p, TMAX):
    t = 0.0

    A, P, M, Q = 0.0, 1, 0.0, 0.0

    times = [t]
    traj = [[A], [P], [M], [Q]]

    while t < TMAX and P > 0:
        tau, event = gillespie_step(p, (A, P, M, Q))

        if tau is None:
            break

        t_next = min(t + tau, TMAX)

        sol = solve_ivp(lambda tt, yy: ode_system(tt, yy, p), [t, t_next], [A, P, M, Q], t_eval=[t_next])

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

def run_pdmp_scenarios(p, TMAX, sc, steps):
    results = {}
    trajectories = []
    for _ in range(steps):
        traj = run_pdmp(p, TMAX)
        trajectories.append(traj)
    results[sc["name"]] = { "trajectories": trajectories,"params": sc}
    return results