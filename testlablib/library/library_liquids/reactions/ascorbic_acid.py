import numpy as np
import testlablib as lab


library = lab.Library()


def ascorbic_acid_rxn_1_info(LiquidStream):
    info = {}
    info["Reactants"] = {"H2C6H6O6": 1}
    info["Products"] = {"HC6H6O6-": 1, "H+": 1}
    info["Reactants Type"] = {"H2C6H6O6": "c"}
    info["Products Type"] = {"HC6H6O6-": "c", "H+": "c"}
    return info

def ascorbic_acid_rxn_1_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-4.1) * np.exp(0 * (1 / T - 1 / 298.15))
    return K

def ascorbic_acid_rxn_2_info(LiquidStream):
    info = {}
    info["Reactants"] = {"HC6H6O6-": 1}
    info["Products"] = {"C6H6O6-2": 1, "H+": 1}
    info["Reactants Type"] = {"HC6H6O6-": "c"}
    info["Products Type"] = {"C6H6O6-2": "c", "H+": "c"}
    return info

def ascorbic_acid_rxn_2_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-11.8) * np.exp(0 * (1 / T - 1 / 298.15))
    return K



