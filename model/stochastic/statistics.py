import numpy as np

def is_extinct(traj):
    return 1 if traj["P"][-1] == 0 else 0

def time_to_resistance(traj, Q_res):
    Q = traj["Q"]
    t = traj["t"]
    idx = np.where(Q >= Q_res)[0]
    if len(idx) == 0:
        return np.inf
    return t[idx[0]]

def extinction_probability(trajectories):
    N = len(trajectories)
    I_ext = sum(is_extinct(traj) for traj in trajectories)
    return I_ext / N

def rescue_probability(trajectories):
    return 1.0 - extinction_probability(trajectories)

def time_to_resistance_distribution(trajectories, Q_res):
    times = []
    for traj in trajectories:
        T = time_to_resistance(traj, Q_res)
        if np.isfinite(T):
            times.append(T)
    return np.array(times)

def mean_time_to_resistance(trajectories, Q_res):
    times = time_to_resistance_distribution(trajectories, Q_res)
    if len(times) == 0:
        return np.inf
    return np.mean(times)