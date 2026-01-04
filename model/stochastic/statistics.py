def is_extinct(traj):
    return 1 if traj["P"][-1] == 0 else 0

def extinction_probability(trajectories):
    N = len(trajectories)
    I_ext = sum(is_extinct(traj) for traj in trajectories)
    return I_ext / N

def rescue_probability(trajectories):
    return 1.0 - extinction_probability(trajectories)