import numpy as np
import testlablib as lab


library = lab.Library()


def water_autoprotolysis_eq_constant(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    # Kw = 10 ** (-14) * np.exp(-6716 * (1/T - 1/298.15))
    Kw = 10 ** ((61.2 - 9.761 * np.log(T)) - 5839 / T)
    return Kw

library.add_LiquidStream_rxn_insta(id="H2O = H+ + OH-",
                                   stoch={"H2O": -1, "H+": 1, "OH-": 1},
                                   unit={"H2O": "x", "H+": "c", "OH-": "c"},
                                   equilibrium_constant=water_autoprotolysis_eq_constant)







