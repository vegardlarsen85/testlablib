import numpy as np
import testlablib as lab


library = lab.Library()



def nh3_rxn_1_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-9.25) * np.exp(-6205 * (1 / T - 1 / 313))
    return K

def nh3_rxn_2_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (0.0436) * np.exp(1479 * (1 / T - 1 / 313))
    return K


library.add_LiquidStream_rxn_insta(id="NH4+ = NH3 + H+",
                                   stoch={"NH4+": -1, "NH3": 1, "H+": 1},
                                   unit={"NH4+": "c", "NH3": "c", "H+": "c"},
                                   equilibrium_constant=nh3_rxn_1_eq_const)

library.add_LiquidStream_rxn_insta(id="HCO3- + NH3 = NH2COO- + H2O",
                                   stoch={"HCO3-": -1, "NH3": -1, "NH2COO-": 1, "H2O": 1},
                                   unit={"HCO3-": "c", "NH3": "c", "NH2COO-": "c", "H2O": "x"},
                                   equilibrium_constant=nh3_rxn_2_eq_const)



