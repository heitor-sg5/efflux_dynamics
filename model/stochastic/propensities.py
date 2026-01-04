def a1_(p, A, P):
    if P <= 0:
        return 0.0
    G = p["G_max"] * (1 / (1 + (A / p["IC50"]) ** p["h"]))
    return p["r_P"] * G * P

def a2_(p, P):
    if P <= 0:
        return 0.0
    return (p["c_P"] * P) + (p["gamma"] * P ** 2)