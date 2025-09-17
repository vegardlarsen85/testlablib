import numpy as np
import testlablib as lab


library = lab.Library()



def citric_acid_rxn_1_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 100.446 - 52.159 * (100 / T) - 15.823 * np.log(T)
    return K

def citric_acid_rxn_2_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 134.23 - 67.107 * (100 / T) - 21.532 * np.log(T)
    return K

def citric_acid_rxn_3_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 188.717 - 87.047 * (100 / T) - 30.583 * np.log(T)
    return K


library.add_LiquidStream_rxn_insta(id="Ci = Ci- + H+",
                                   stoch={"Ci": -1, "Ci-": 1, "H+": 1},
                                   unit={"Ci": "c", "Ci-": "c", "H+": "c"},
                                   equilibrium_constant=citric_acid_rxn_1_eq_const)

library.add_LiquidStream_rxn_insta(id="Ci- = Ci-2 + H+",
                                   stoch={"Ci-": -1, "Ci-2": 1, "H+": 1},
                                   unit={"Ci-": "c", "Ci-2": "c", "H+": "c"},
                                   equilibrium_constant=citric_acid_rxn_2_eq_const)

library.add_LiquidStream_rxn_insta(id="Ci-2 = Ci-3 + H+",
                                   stoch={"Ci-2": -1, "Ci-3": 1, "H+": 1},
                                   unit={"Ci-2": "c", "Ci-3": "c", "H+": "c"},
                                   equilibrium_constant=citric_acid_rxn_3_eq_const)


