from model.deterministic.odes import growth_rate

def plasmid_free_equilibrium(p):
    a = (p["k_in"] * p["A_ext"])
    b = (p["k_out"] * p["E_h"])
    A_star = a / b if b != 0 else 0 
    return A_star

def invasion_eigenvalue(p):
    A_star = plasmid_free_equilibrium(p)
    G_star = growth_rate(A_star, p)
    return (p['r_P'] * G_star) - p['c_P']