def critical_Eh(p):
    a = (p['k_in'] * p['A_ext']) / (p['k_out'] * p['IC50'])
    b = (p['c_P'] / ((p['r_P'] * p['G_max']) - p['c_P']))**(1 / p['h'])
    return a * b

def critical_Aext(p):
    a = (p['k_out'] * p['IC50'] * p['E_h']) / p['k_in']
    b = ((p['r_P'] * p['G_max']) - p['c_P']) / (p['c_P'])**(1 / p['h'])
    return a * b