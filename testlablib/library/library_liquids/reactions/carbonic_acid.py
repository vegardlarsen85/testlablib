import numpy as np
import testlablib as lab
from testlablib.library.library_liquids.reactions.water import water_autoprotolysis_eq_constant



library = lab.Library()


def carbonic_acid_rxn_1_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** ((T * (102.2678 - 15.9739 * np.log(T)) - 5251.064) / T)
    return K

def carbonic_acid_rxn_2_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (95.573 - 15.4095 * np.log(T) - 5398.98 / T)
    return K

def carbonic_acid_rxn_5_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    Kw = water_autoprotolysis_eq_constant(LiquidStream)
    K = K1 / Kw
    return K

def carbonic_acid_rxn_6_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    K2 = carbonic_acid_rxn_2_eq_const(LiquidStream)
    Kw = water_autoprotolysis_eq_constant(LiquidStream)
    K = K1 * K2 / Kw **2
    return K

def carbonic_acid_rxn_1_rate_law_kmol_m3s(LiquidStream):

    # Temperature
    T = LiquidStream.temp_K

    # Concentrations
    a_CO2 = LiquidStream.get_specie_activity_kmol_m3(id="CO2")
    x_H2O = LiquidStream.get_specie_activity_coefficient(id="H2O") * LiquidStream.get_specie_molar_fraction(id="H2O")
    a_HCO3 = LiquidStream.get_specie_activity_kmol_m3(id="HCO3-")
    a_H = LiquidStream.get_specie_activity_kmol_m3(id="H+")
    a_OH = LiquidStream.get_specie_activity_kmol_m3(id="OH-")

    # Rate Coefficients
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    k1_forward = 0.0328 * np.exp(-8590 * ( 1 /T - 1/ 298.15))
    k1_backward = k1_forward / K1

    K5 = carbonic_acid_rxn_5_eq_const(LiquidStream)
    k5_forward = 8416 * np.exp(-6667 * (1 / T - 1 / 298.15))
    k5_backward = k5_forward / K5

    # Rate Law
    r_forward = (k1_forward * x_H2O + k5_forward * a_OH) * a_CO2
    r_backward = (k1_backward * a_H + k5_backward) * a_HCO3
    return r_forward, r_backward

def carbonic_acid_rxn_1_exothermic_heat_kJ_kmol(LiquidStream):

    K0 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    T0 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K + 0.1

    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    T1 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K - 0.1

    return 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))



library.add_LiquidStream_rxn_insta(id="CO2 + H2O = HCO3- + H+",
                                   stoch={"H2O": -1, "CO2": -1, "H+": 1, "HCO3-": 1},
                                   unit={"H2O": "x", "CO2": "c", "H+": "c", "HCO3-": "c"},
                                   equilibrium_constant=carbonic_acid_rxn_1_eq_const)

library.add_LiquidStream_rxn_insta(id="HCO3- = CO3-2 + H+",
                                   stoch={"HCO3-": -1, "CO3-2": 1, "H+": 1},
                                   unit={"HCO3-": "c", "CO3-2": "c", "H+": "c"},
                                   equilibrium_constant=carbonic_acid_rxn_2_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + OH- = HCO3-",
                                   stoch={"CO2": -1, "OH-": -1, "HCO3-": 1},
                                   unit={"CO2": "c", "OH-": "c", "HCO3-": "c"},
                                   equilibrium_constant=carbonic_acid_rxn_5_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + 2OH- = CO3-2 + H2O",
                                   stoch={"CO2": -1, "OH-": -2, "CO3-2": 1, "H2O": 1},
                                   unit={"CO2": "c", "OH-": "c", "CO3-2": "c", "H2O": "x"},
                                   equilibrium_constant=carbonic_acid_rxn_6_eq_const)

library.add_LiquidStream_rxn_reversible(id="CO2 + H2O -> HCO3- + H+",
                                        stoch={"CO2": -1, "H2O": -1, "HCO3-": 1, "H+": 1},
                                        dependencies=["CO2", "H2O", "HCO3-", "H+", "OH-"],
                                        rate_kmol_m3s=carbonic_acid_rxn_1_rate_law_kmol_m3s,
                                        exothermic_heat_kJ_kmol=carbonic_acid_rxn_1_exothermic_heat_kJ_kmol)



