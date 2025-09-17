import numpy as np
import testlablib as lab


def activity_coefficient_truesdell_jones(LiquidStream, id):

    # Truesdell Jones...
    I = LiquidStream.get_solution_ionic_strength_mol_kg()
    sqrtI = np.sqrt(I)
    T = LiquidStream.get_solution_temp_K()
    A = 0.51
    B = 0.328
    # A = 2.74 * 10**(-6) * T**2 - 7.60*10**(-4) * T + 0.4916
    # B = 1.62 * 10**(-4) * T + 0.2799
    z = LiquidStream.get_specie_charge(id)
    a, b = 4.57, 0
    if z == 0:
        log10_gamma = b * I
    else:
        log10_gamma = - A * z ** 2 * sqrtI / (1 + a * B * sqrtI) + b * I
    gamma = 10 ** log10_gamma
    return gamma

def activity_coefficient_ideal(LiquidStream, id):
    return 1


