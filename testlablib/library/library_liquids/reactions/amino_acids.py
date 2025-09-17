import numpy as np
import testlablib as lab
from testlablib.library.library_liquids.reactions.carbonic_acid import carbonic_acid_rxn_1_eq_const
from testlablib.library.library_liquids.reactions.water import water_autoprotolysis_eq_constant


library = lab.Library()


def lysine_carboxyl_group_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    ln_K = np.log(10 ** (-2.15)) + (4.54 * 1000 / 8.314) * (1 / 298.15 - 1 / T) - (120.9 / 8.314) * (
            np.log(T / 298.15) + 298.15 / T - 1)
    K = np.exp(ln_K)
    return K

def lysine_amino_group_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    ln_K = np.log(10 ** (-9.16)) + (44.4 * 1000 / 8.314) * (1 / 298.15 - 1 / T) - (172.9 / 8.314) * (
            np.log(T / 298.15) + 298.15 / T - 1)
    K = np.exp(ln_K)
    return K

def lysine_carbamate_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    ln_K = -2.49116083 - 1.52403685 * (1 / 298.15 - 1 / T) - 0.00237651031 * (
            np.log(T / 298.15) + 298.15 / T - 1)
    K = np.exp(ln_K)
    return K

def lysine_side_chain_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    ln_K = np.log(10 ** (-10.67)) + 9634.44530 * (1 / 298.15 - 1 / T) - 17.4833263 * (
            np.log(T / 298.15) + 298.15 / T - 1)
    K = np.exp(ln_K)
    return K

def lysine_side_chain_carbamate_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    ln_K = 2.05546021 + 1.08171196 * (1 / 298.15 - 1 / T) - 0.00166878193 * (np.log(T / 298.15) + 298.15 / T - 1)
    K = np.exp(ln_K)
    return K

def lysine_carbon_dioxide_eq_const(LiquidStream):
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    KC2 = lysine_side_chain_carbamate_eq_const(LiquidStream)
    KCO2 = K1 / KC2
    return KCO2

def lysine_carbamate_to_carbamate_eq_const(LiquidStream):
    KC1 = lysine_carbamate_eq_const(LiquidStream)
    KR = lysine_side_chain_eq_const(LiquidStream)
    KC2 = lysine_side_chain_carbamate_eq_const(LiquidStream)
    KC3 = KC1 * KR / KC2
    return KC3

# Reaction Rate Coeff. TBD
def lysine_carbon_dioxide_rate_kmol_m3s(LiquidStream):
    T = LiquidStream.temp_K

    a_CO2 = LiquidStream.get_specie_activity_kmol_m3(id="CO2")
    a_Lys = LiquidStream.get_specie_activity_kmol_m3(id="Lys-")
    a_LysCOO = LiquidStream.get_specie_activity_kmol_m3(id="LysCOO-2")
    a_H = LiquidStream.get_specie_activity_kmol_m3(id="H+")

    K = lysine_carbon_dioxide_eq_const(LiquidStream)
    k_forward = 15000 * np.exp(-3127 * (1 / T - 1 / 298.15))
    k_backward = k_forward / K

    r_forward = k_forward * a_CO2 * a_Lys
    r_backward = k_backward * a_LysCOO * a_H

    return r_forward, r_backward

def lysine_carbon_dioxide_exothermic_heat_kJ_kmol(LiquidStream):
    K0 = lysine_carbon_dioxide_eq_const(LiquidStream)
    T0 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K + 0.1

    K1 = lysine_carbon_dioxide_eq_const(LiquidStream)
    T1 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K - 0.1

    q = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))

    return q


library.add_LiquidStream_rxn_insta(id="Lys+2 = Lys+ + H+",
                                   stoch={"Lys+2": -1, "Lys+": 1, "H+": 1},
                                   unit={"Lys+2": "c", "Lys+": "c", "H+": "c"},
                                   equilibrium_constant=lysine_carboxyl_group_eq_const)

library.add_LiquidStream_rxn_insta(id="Lys+ = Lys + H+",
                                   stoch={"Lys+": -1, "Lys": 1, "H+": 1},
                                   unit={"Lys+": "c", "Lys": "c", "H+": "c"},
                                   equilibrium_constant=lysine_amino_group_eq_const)

library.add_LiquidStream_rxn_insta(id="Lys = Lys- + H+",
                                   stoch={"Lys": -1, "Lys-": 1, "H+": 1},
                                   unit={"Lys": "c", "Lys-": "c", "H+": "c"},
                                   equilibrium_constant=lysine_side_chain_eq_const)

library.add_LiquidStream_rxn_insta(id="LysCOO- + H2O = Lys + HCO3-",
                                   stoch={"LysCOO-": -1, "H2O": -1, "Lys": 1, "HCO3-": 1},
                                   unit={"LysCOO-": "c", "H2O": "x", "Lys": "c", "HCO3-": "c"},
                                   equilibrium_constant=lysine_carbamate_eq_const)

library.add_LiquidStream_rxn_insta(id="LysCOO-2 + H2O = Lys- + HCO3-",
                                   stoch={"LysCOO-2": -1, "H2O": -1, "Lys-": 1, "HCO3-": 1},
                                   unit={"LysCOO-2": "c", "H2O": "x", "Lys-": "c", "HCO3-": "c"},
                                   equilibrium_constant=lysine_side_chain_carbamate_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + Lys- = LysCOO-2 + H+",
                                   stoch={"CO2": -1, "Lys-": -1, "LysCOO-2": 1, "H+": 1},
                                   unit={"CO2": "c", "Lys-": "c", "LysCOO-2": "c", "H+": "c"},
                                   equilibrium_constant=lysine_carbon_dioxide_eq_const)

library.add_LiquidStream_rxn_insta(id="LysCOO- = LysCOO-2 + H+",
                                   stoch={"LysCOO-": -1, "LysCOO-2": 1, "H+": 1},
                                   unit={"LysCOO-": "c", "LysCOO-2": "c", "H+": "c"},
                                   equilibrium_constant=lysine_carbamate_to_carbamate_eq_const)

library.add_LiquidStream_rxn_reversible(id="CO2 + Lys- -> H+ + LysCOO-2",
                                        stoch={"CO2": -1, "Lys-": -1, "H+": 1, "LysCOO-2": 1},
                                        dependencies=["CO2", "Lys-", "H+", "LysCOO-2"],
                                        rate_kmol_m3s=lysine_carbon_dioxide_rate_kmol_m3s,
                                        exothermic_heat_kJ_kmol=lysine_carbon_dioxide_exothermic_heat_kJ_kmol)


def sarcosine_carboxyl_group_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    ln_K = np.log(10 ** (-2.36)) + (4.54 * 1000 / 8.314) * (
            1 / 298.15 - 1 / T) - (120.9 / 8.314) * (np.log(T / 298.15) + 298.15 / T - 1)
    K = np.exp(ln_K)
    return K

def sarcosine_amino_group_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    ln_K = np.log(10 ** (-11.64)) + (44.4 * 1000 / 8.314) * (1 / 298.15 - 1 / T) - (
            172.9 / 8.314) * (np.log(T / 298.15) + 298.15 / T - 1)
    K = np.exp(ln_K)
    return K

def sarcosine_carbamate_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    ln_K = 3.5 + 1000 * 2.406 * (1 / T - 1 / 298)
    K = np.exp(ln_K)
    return K

def sarcosine_carbon_dioxide_eq_const(LiquidStream):
    Kc = sarcosine_carbamate_eq_const(LiquidStream)
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    K = K1 * Kc
    return K

def sarcosine_carbon_dioxide_rate_kmol_m3s(LiquidStream):
    T = LiquidStream.temp_K

    a_CO2 = LiquidStream.get_specie_activity_kmol_m3(id="CO2")
    a_Sar = LiquidStream.get_specie_activity_kmol_m3(id="Sar-")
    a_SarCOO = LiquidStream.get_specie_activity_kmol_m3(id="SarCOO-2")
    a_H = LiquidStream.get_specie_activity_kmol_m3(id="H+")

    K = sarcosine_carbon_dioxide_eq_const(LiquidStream)
    k_forward = 15000 * np.exp(-3127 * (1 / T - 1 / 298.15))
    k_backward = k_forward / K

    r_forward = k_forward * a_CO2 * a_Sar
    r_backward = k_backward * a_SarCOO * a_H

    return r_forward, r_backward

def sarcosine_carbon_dioxide_exothermic_heat_kJ_kmol(LiquidStream):
    K0 = sarcosine_carbon_dioxide_eq_const(LiquidStream)
    T0 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K + 0.1

    K1 = sarcosine_carbon_dioxide_eq_const(LiquidStream)
    T1 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K - 0.1

    q = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))

    return q



library.add_LiquidStream_rxn_insta(id="Sar+ = Sar + H+",
                                   stoch={"Sar+": -1, "Sar": 1, "H+": 1},
                                   unit={"Sar+": "c", "Sar": "c", "H+": "c"},
                                   equilibrium_constant=sarcosine_carboxyl_group_eq_const)

library.add_LiquidStream_rxn_insta(id="Sar = Sar- + H+",
                                   stoch={"Sar": -1, "Sar-": 1, "H+": 1},
                                   unit={"Sar": "c", "Sar-": "c", "H+": "c"},
                                   equilibrium_constant=sarcosine_amino_group_eq_const)

library.add_LiquidStream_rxn_insta(id="Sar- + HCO3- = SarCOO-2 + H2O",
                                   stoch={"Sar-": -1, "HCO3-": -1, "SarCOO-2": 1, "H2O": 1},
                                   unit={"Sar-": "c", "HCO3-": "c", "SarCOO-2": "c", "H2O": "x"},
                                   equilibrium_constant=sarcosine_carbamate_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + Sar- = H+ + SarCOO-2",
                                   stoch={"CO2": -1, "Sar-": -1, "H+": 1, "SarCOO-2": 1},
                                   unit={"CO2": "c", "Sar-": "c", "H+": "c", "SarCOO-2": "c"},
                                   equilibrium_constant=sarcosine_carbon_dioxide_eq_const)

library.add_LiquidStream_rxn_reversible(id="CO2 + Sar- -> H+ + SarCOO-2",
                                        stoch={"CO2": -1, "Sar-": -1, "H+": 1, "SarCOO-2": 1},
                                        dependencies=["CO2", "Sar-", "H+", "SarCOO-2"],
                                        rate_kmol_m3s=sarcosine_carbon_dioxide_rate_kmol_m3s,
                                        exothermic_heat_kJ_kmol=sarcosine_carbon_dioxide_exothermic_heat_kJ_kmol)



def glycine_carboxyl_group_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    ln_K = np.log(10 ** (-2.34)) + (4.54 * 1000 / 8.314) * (1 / 298.15 - 1 / T) - (120.9 / 8.314) * (
            np.log(T / 298.15) + 298.15 / T - 1)
    K = np.exp(ln_K)
    return K

def glycine_amino_group_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    ln_K = np.log(10 ** (-9.78)) + (44.4 * 1000 / 8.314) * (1 / 298.15 - 1 / T) - (172.9 / 8.314) * (
            np.log(T / 298.15) + 298.15 / T - 1)
    K = np.exp(ln_K)
    return K

def glycine_carbamate_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    # ln_K = 1.9633155 + 1000 * 1.97422924 * (1 / T - 1 / 298)
    ln_K = 1.13969548 + 8.86625640e-04 * (1 / 298.15 - 1 / T) + 4.87705989e-04 * (
            np.log(T / 298.15) + 298.15 / T - 1)
    K = np.exp(ln_K)
    return K

def glycine_carbon_dioxide_eq_const(LiquidStream):
    Kc = glycine_carbamate_eq_const(LiquidStream)
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    K = K1 * Kc
    return K

def glycine_carbon_dioxide_rate_kmol_m3s(LiquidStream):
    T = LiquidStream.temp_K

    a_CO2 = LiquidStream.get_specie_activity_kmol_m3(id="CO2")
    a_Gly = LiquidStream.get_specie_activity_kmol_m3(id="Gly-")
    a_GlyCOO = LiquidStream.get_specie_activity_kmol_m3(id="GlyCOO-2")
    a_H = LiquidStream.get_specie_activity_kmol_m3(id="H+")

    K = glycine_carbon_dioxide_eq_const(LiquidStream)
    k_forward = 13900 * np.exp(-5459 * (1 / T - 1 / 298.15))
    k_backward = k_forward / K

    r_forward = k_forward * a_CO2 * a_Gly
    r_backward = k_backward * a_GlyCOO * a_H
    return r_forward, r_backward

def glycine_carbon_dioxide_exothermic_heat_kJ_kmol(LiquidStream):
    K0 = glycine_carbon_dioxide_eq_const(LiquidStream)
    T0 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K + 0.1

    K1 = glycine_carbon_dioxide_eq_const(LiquidStream)
    T1 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K - 0.1

    q = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))

    return q


library.add_LiquidStream_rxn_insta(id="Gly+ = Gly + H+",
                                   stoch={"Gly+": -1, "Gly": 1, "H+": 1},
                                   unit={"Gly+": "c", "Gly": "c", "H+": "c"},
                                   equilibrium_constant=glycine_carboxyl_group_eq_const)

library.add_LiquidStream_rxn_insta(id="Gly = Gly- + H+",
                                   stoch={"Gly": -1, "Gly-": 1, "H+": 1},
                                   unit={"Gly": "c", "Gly-": "c", "H+": "c"},
                                   equilibrium_constant=glycine_amino_group_eq_const)

library.add_LiquidStream_rxn_insta(id="Gly- + HCO3- = GlyCOO-2 + H2O",
                                   stoch={"Gly-": -1, "HCO3-": -1, "GlyCOO-2": 1, "H2O": 1},
                                   unit={"Gly-": "c", "HCO3-": "c", "GlyCOO-2": "c", "H2O": "x"},
                                   equilibrium_constant=glycine_carbamate_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + Gly- = H+ + GlyCOO-2",
                                   stoch={"CO2": -1, "Gly-": -1, "H+": 1, "GlyCOO-2": 1},
                                   unit={"CO2": "c", "Gly-": "c", "H+": "c", "GlyCOO-2": "c"},
                                   equilibrium_constant=glycine_carbon_dioxide_eq_const)

library.add_LiquidStream_rxn_reversible(id="CO2 + Gly- -> H+ + GlyCOO-2",
                                        stoch={"CO2": -1, "Gly-": -1, "H+": 1, "GlyCOO-2": 1},
                                        dependencies=["CO2", "Gly-", "H+", "GlyCOO-2"],
                                        rate_kmol_m3s=glycine_carbon_dioxide_rate_kmol_m3s,
                                        exothermic_heat_kJ_kmol=glycine_carbon_dioxide_exothermic_heat_kJ_kmol)



