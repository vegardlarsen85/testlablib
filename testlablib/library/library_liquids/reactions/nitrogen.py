import numpy as np
import testlablib as lab


library = lab.Library()





def nitric_acid_rxn_1_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-0.44) * np.exp(-284 * (1 / T - 1 / 298.15))
    return K


library.add_LiquidStream_rxn_insta(id="HNO3 = NO3- + H+",
                                   stoch={"HNO3": -1, "NO3-": 1, "H+": 1},
                                   unit={"HNO3": "c", "NO3-": "c", "H+": "c"},
                                   equilibrium_constant=nitric_acid_rxn_1_eq_const)



def nitrous_acid_rxn_1_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-3.39) * np.exp(0 * (1 / T - 1 / 298.15))  # Eq.Coeff TBD
    return K


library.add_LiquidStream_rxn_insta(id="HNO2 = NO2- + H+",
                                   stoch={"HNO2": -1, "NO2-": 1, "H+": 1},
                                   unit={"HNO2": "c", "NO2-": "c", "H+": "c"},
                                   equilibrium_constant=nitrous_acid_rxn_1_eq_const)


