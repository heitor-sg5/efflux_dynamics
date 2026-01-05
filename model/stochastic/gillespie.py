import numpy as np
from model.stochastic.propensities import a1_, a2_

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