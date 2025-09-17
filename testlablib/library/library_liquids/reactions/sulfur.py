import numpy as np
import testlablib as lab


"""""""""
Notes:
- Equilibrium Constants are Calculated from Enthalpy, Entropy and Gibbs Free Energy at 25C.
- Rate Law for Oxidation of Sulfite to Sulfate is probably inaccurate and need to be Reviewed (!!!)

"""""""""


library = lab.Library()


def sulfurous_acid_rxn_10_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-1.86) * np.exp(2140 * (1 / T - 1 / 298.15))
    return K

def sulfurous_acid_rxn_11_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-7.17) * np.exp(439 * (1 / T - 1 / 298.15))
    return K

def sulfurous_acid_rxn_12_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-12.13) * np.exp(-8863 * (1 / T - 1 / 298.15))
    return K

def sulfurous_acid_rxn_13_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-6.815) * np.exp(-7161 * (1 / T - 1 / 298.15))
    return K

#def rate_law_kmol_m3s_SO3_O2_rxn_1(LiquidStream):
#    epsilon = 10 ** (-15)
#    a_O2 = LiquidStream.get_specie_activity_kmol_m3(id="O2")
#    a_SO3 = LiquidStream.get_specie_activity_kmol_m3(id="SO3-2")
#    a_HSO3 = LiquidStream.get_specie_activity_kmol_m3(id="HSO3-")
#    rate = 5.77 * 10 ** (-5) * (a_SO3 + a_HSO3) ** 0.65 * (a_O2 > epsilon)
#    return rate


def sulfuric_acid_rxn_1_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-1.97689) * np.exp(2637 * (1 / T - 1 / 298.15))
    return K


library.add_LiquidStream_rxn_insta(id="SO2 + H2O = HSO3- + H+",
                                   stoch={"SO2": -1, "H2O": -1, "HSO3-": 1, "H+": 1},
                                   unit={"SO2": "c", "H2O": "x", "HSO3-": "c", "H+": "c"},
                                   equilibrium_constant=sulfurous_acid_rxn_10_eq_const)

library.add_LiquidStream_rxn_insta(id="HSO3- = SO3-2 + H+",
                                   stoch={"HSO3-": -1, "SO3-2": 1, "H+": 1},
                                   unit={"HSO3-": "c", "SO3-2": "c", "H+": "c"},
                                   equilibrium_constant=sulfurous_acid_rxn_11_eq_const)

library.add_LiquidStream_rxn_insta(id="HSO3- = SO2 + OH-",
                                   stoch={"HSO3-": -1, "SO2": 1, "OH-": 1},
                                   unit={"HSO3-": "c", "SO2": "c", "OH-": "c"},
                                   equilibrium_constant=sulfurous_acid_rxn_12_eq_const)

library.add_LiquidStream_rxn_insta(id="SO3-2 + H2O = HSO3- + OH-",
                                   stoch={"SO3-2": -1, "H2O": -1, "HSO3-": 1, "OH-": 1},
                                   unit={"SO3-2": "c", "H2O": "x", "HSO3-": "c", "OH-": "c"},
                                   equilibrium_constant=sulfurous_acid_rxn_13_eq_const)

#library.add_LiquidStream_rxn_irreversible(id="2SO3-2 + O2 -> 2SO4-2",
#                                          stoch={"SO3-2": -2, "O2": -1, "SO4-2": 2},
#                                          rate_kmol_m3s=rate_law_kmol_m3s_SO3_O2_rxn_1,
#                                          exothermic_heat_kJ_kmol=None)


library.add_LiquidStream_rxn_insta(id="HSO4- = SO4-2 + H+",
                                   stoch={"HSO4-": -1, "SO4-2": 1, "H+": 1},
                                   unit={"HSO4-": "c", "SO4-2": "c", "H+": 1},
                                   equilibrium_constant=sulfuric_acid_rxn_1_eq_const)


