import numpy as np
import testlablib as lab
from testlablib.library.library_liquids.reactions.carbonic_acid import carbonic_acid_rxn_1_eq_const
from testlablib.library.library_liquids.reactions.water import water_autoprotolysis_eq_constant


library = lab.Library()


def mea_rxn_20_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-9.27) * np.exp(-6375 * (1 / T - 1 / 313))
    return K

def mea_rxn_21_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (1.44) * np.exp(2406 * (1 / T - 1 / 313))
    return K

def mea_rxn_23_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K20 = mea_rxn_20_eq_const(LiquidStream)
    K21 = mea_rxn_21_eq_const(LiquidStream)
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    K = K21 * K1 / K20
    return K

def mea_rxn_24_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K21 = mea_rxn_21_eq_const(LiquidStream)
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    K = K21 * K1
    return K

def mea_rxn_24_rate_law_kmol_m3s(LiquidStream):
    T = LiquidStream.temp_K
    K = mea_rxn_24_eq_const(LiquidStream)

    a_CO2 = LiquidStream.get_specie_activity_kmol_m3(id="CO2")
    a_MEA = LiquidStream.get_specie_activity_kmol_m3(id="MEA")
    a_MEACOO = LiquidStream.get_specie_activity_kmol_m3(id="MEACOO-")
    a_H = LiquidStream.get_specie_activity_kmol_m3(id="H+")

    k_forward = 5993 * np.exp(-5400 * (1 / T - 1 / 298.15))
    k_backward = k_forward / K

    r_forward = k_forward * a_CO2 * a_MEA
    r_backward = k_backward * a_MEACOO * a_H
    return r_forward, r_backward

def mea_rxn_24_exothermic_heat_kJ_kmol(LiquidStream):
    K0 = mea_rxn_24_eq_const(LiquidStream)
    T0 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K + 0.1

    K1 = mea_rxn_24_eq_const(LiquidStream)
    T1 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K - 0.1

    q = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))

    return q


library.add_LiquidStream_rxn_insta(id="MEAH+ = MEA + H+",
                                   stoch={"MEAH+": -1, "MEA": 1, "H+": 1},
                                   unit={"MEAH+": "c", "MEA": "c", "H+": "c"},
                                   equilibrium_constant=mea_rxn_20_eq_const)

library.add_LiquidStream_rxn_insta(id="HCO3- + MEA = MEACOO- + H2O",
                                   stoch={"HCO3-": -1, "MEA": -1, "MEACOO-": 1, "H2O": 1},
                                   unit={"HCO3-": "c", "MEA": "c", "MEACOO-": "c", "H2O": "x"},
                                   equilibrium_constant=mea_rxn_21_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + 2MEA = MEACOO- + MEAH+",
                                   stoch={"CO2": -1, "MEA": -2, "MEACOO-": 1, "MEAH+": 1},
                                   unit={"CO2": "c", "MEA": "c", "MEACOO-": "c", "MEAH+": "c"},
                                   equilibrium_constant=mea_rxn_23_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + MEA = MEACOO- + H+",
                                   stoch={"CO2": -1, "MEA": -1, "MEACOO-": 1, "H+": 1},
                                   unit={"CO2": "c", "MEA": "c", "MEACOO-": "c", "H+": "c"},
                                   equilibrium_constant=mea_rxn_24_eq_const)

library.add_LiquidStream_rxn_reversible(id="CO2 + MEA -> MEACOO- + H+",
                                        stoch={"CO2": -1, "MEA": -1, "MEACOO-": 1, "H+": 1},
                                        dependencies=["CO2", "MEA", "MEACOO-", "H+"],
                                        rate_kmol_m3s=mea_rxn_24_rate_law_kmol_m3s,
                                        exothermic_heat_kJ_kmol=mea_rxn_24_exothermic_heat_kJ_kmol)


def mdea_rxn_40_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    # K = 10**(-14.01 + 0.0184 * T)                                  # Older
    # K = 10 ** (-7.75) * np.exp(-1500 * (1 / T - 1 / 298.15))       # Old
    K = 10 ** (-8.35) * np.exp(-4900 * (1 / T - 1 / 298.15))
    return K

library.add_LiquidStream_rxn_insta(id="MDEAH+ = MDEA + H+",
                                   stoch={"MDEAH+": -1, "H+": 1, "MDEA": 1},
                                   unit={"MDEAH+": "c", "H+": "c", "MDEA": "c"},
                                   equilibrium_constant=mdea_rxn_40_eq_const)



def amp_rxn_70_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    # K = 10**(-9.39) * np.exp(-6976 * (1/T - 1/313))
    K = 10 ** (-9.37) * np.exp(-4691 * (1 / T - 1 / 313))
    return K

def amp_rxn_71_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-1.3) * np.exp(5292 * (1 / T - 1 / 313))
    return K

def amp_rxn_73_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K70 = amp_rxn_70_eq_const(LiquidStream)
    K71 = amp_rxn_71_eq_const(LiquidStream)
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    K = K71 * K1 / K70
    return K

def amp_rxn_74_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K71 = amp_rxn_71_eq_const(LiquidStream)
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    K = K71 * K1
    return K

def amp_rxn_74_rate_law_kmol_m3s(LiquidStream):
    T = LiquidStream.temp_K

    a_CO2 = LiquidStream.get_specie_activity_kmol_m3(id="CO2")
    a_AMP = LiquidStream.get_specie_activity_kmol_m3(id="AMP")
    a_AMPCOO = LiquidStream.get_specie_activity_kmol_m3(id="AMPCOO-")
    a_H = LiquidStream.get_specie_activity_kmol_m3(id="H+")

    K = amp_rxn_74_eq_const(LiquidStream)
    k_forward = 1500 * np.exp(-5176 * (1 / T - 1 / 298.15))
    k_backward = k_forward / K

    r_forward = k_forward * a_CO2 * a_AMP
    r_backward = k_backward * a_AMPCOO * a_H
    return r_forward, r_backward

def amp_rxn_74_exothermic_heat_kJ_kmol(LiquidStream):
    K0 = amp_rxn_74_eq_const(LiquidStream)
    T0 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K + 0.1

    K1 = amp_rxn_74_eq_const(LiquidStream)
    T1 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K - 0.1

    q = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))

    return q


library.add_LiquidStream_rxn_insta(id="AMPH+ = AMP + H+",
                                   stoch={"AMPH+": -1, "H+": 1, "AMP": 1},
                                   unit={"AMPH+": "c", "H+": "c", "AMP": "c"},
                                   equilibrium_constant=amp_rxn_70_eq_const)

library.add_LiquidStream_rxn_insta(id="HCO3- + AMP = AMPCOO- + H2O",
                                   stoch={"HCO3-": -1, "AMP": -1, "AMPCOO-": 1, "H2O": 1},
                                   unit={"HCO3-": "c", "AMP": "c", "AMPCOO-": "c", "H2O": "x"},
                                   equilibrium_constant=amp_rxn_71_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + 2AMP = AMPCOO- + AMPH+",
                                   stoch={"CO2": -1, "AMP": -2, "AMPCOO-": 1, "AMPH+": 1},
                                   unit={"CO2": "c", "AMP": "c", "AMPCOO-": "c", "AMPH+": "c"},
                                   equilibrium_constant=amp_rxn_73_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + AMP = AMPCOO- + H+",
                                   stoch={"CO2": -1, "AMP": -1, "AMPCOO-": 1, "H+": 1},
                                   unit={"CO2": "c", "AMP": "c", "AMPCOO-": "c", "H+": "c"},
                                   equilibrium_constant=amp_rxn_74_eq_const)

library.add_LiquidStream_rxn_reversible(id="CO2 + AMP -> AMPCOO- + H+",
                                        stoch={"CO2": -1, "AMP": -1, "AMPCOO-": 1, "H+": 1},
                                        dependencies=["CO2", "AMP", "AMPCOO-", "H+"],
                                        rate_kmol_m3s=amp_rxn_74_rate_law_kmol_m3s,
                                        exothermic_heat_kJ_kmol=amp_rxn_74_exothermic_heat_kJ_kmol)



def pz_rxn_50_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-9.30) * np.exp(-4570 * (1 / T - 1 / 313))
    return K

def pz_rxn_51_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (1.2) * np.exp(3849 * (1 / T - 1 / 313))
    return K

def pz_rxn_52_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-8.9) * np.exp(-2165 * (1 / T - 1 / 313))
    return K

def pz_rxn_53_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    # K = 10**(0.8) * np.exp(2887 * (1/T - 1/313))
    K = 10 ** (0.6) * np.exp(6014 * (1 / T - 1 / 313))
    return K

def pz_rxn_54_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    # K = 10**(-8.9) * np.exp(-2887 * (1/T - 1/313))
    K = 10 ** (-8.44) * np.exp(-7096 * (1 / T - 1 / 313))
    return K

def pz_rxn_55_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-5.07) * np.exp(-3850 * (1 / T - 1 / 313))
    return K

def pz_rxn_56_eq_const(LiquidStream):
    K51 = pz_rxn_51_eq_const(LiquidStream)
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    K56 = K51 * K1
    return K56

def pz_rxn_57_eq_const(LiquidStream):
    K53 = pz_rxn_53_eq_const(LiquidStream)
    K1 = carbonic_acid_rxn_1_eq_const(LiquidStream)
    K57 = K53 * K1
    return K57

def pz_rxn_56_rate_law_kmol_m3s(LiquidStream):

    T = LiquidStream.temp_K
    a_CO2 = LiquidStream.get_specie_activity_kmol_m3(id="CO2")
    a_PZ = LiquidStream.get_specie_activity_kmol_m3(id="PZ")
    a_PZCOO = LiquidStream.get_specie_activity_kmol_m3(id="PZCOO-")
    a_H = LiquidStream.get_specie_activity_kmol_m3(id="H+")

    k_forward = 58000 * np.exp(-4209 * (1 / T - 1 / 298.15))
    r_forward = k_forward * a_CO2 * a_PZ

    k_backward = k_forward / pz_rxn_56_eq_const(LiquidStream)
    r_backward = k_backward * a_PZCOO * a_H

    return r_forward, r_backward

def pz_rxn_56_exothermic_heat_kJ_kmol(LiquidStream):
    K0 = pz_rxn_56_eq_const(LiquidStream)
    T0 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K + 0.1

    K1 = pz_rxn_56_eq_const(LiquidStream)
    T1 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K - 0.1

    q = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))

    return q

def pz_rxn_57_rate_law_kmol_m3s(LiquidStream):

    T = LiquidStream.temp_K

    a_CO2 = LiquidStream.get_specie_activity_kmol_m3(id="CO2")
    a_PZCOO = LiquidStream.get_specie_activity_kmol_m3(id="PZCOO-")
    a_PZCOO2 = LiquidStream.get_specie_activity_kmol_m3(id="PZ(COO)2-2")
    a_H = LiquidStream.get_specie_activity_kmol_m3(id="H+")

    k_forward = 14500 * np.exp(-4209 * (1 / T - 1 / 298.15))
    r_forward = k_forward * a_CO2 * a_PZCOO

    k_backward = k_forward / pz_rxn_57_eq_const(LiquidStream)
    r_backward = k_backward * a_PZCOO2 * a_H

    return r_forward, r_backward

def pz_rxn_57_exothermic_heat_kJ_kmol(LiquidStream):
    K0 = pz_rxn_57_eq_const(LiquidStream)
    T0 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K + 0.1

    K1 = pz_rxn_57_eq_const(LiquidStream)
    T1 = LiquidStream.temp_K
    LiquidStream.temp_K = LiquidStream.temp_K - 0.1

    q = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))

    return q



library.add_LiquidStream_rxn_insta(id="PZH+ = PZ + H+",
                                   stoch={"PZH+": -1, "PZ": 1, "H+": 1},
                                   unit={"PZH+": "c", "PZ": "c", "H+": "c"},
                                   equilibrium_constant=pz_rxn_50_eq_const)

library.add_LiquidStream_rxn_insta(id="HCO3- + PZ = PZCOO- + H2O",
                                   stoch={"HCO3-": -1, "PZ": -1, "PZCOO-": 1, "H2O": 1},
                                   unit={"HCO3-": "c", "PZ": "c", "PZCOO-": "c", "H2O": "x"},
                                   equilibrium_constant=pz_rxn_51_eq_const)

library.add_LiquidStream_rxn_insta(id="HPZCOO = PZCOO- + H+",
                                   stoch={"HPZCOO": -1, "PZCOO-": 1, "H+": 1},
                                   unit={"HPZCOO": "c", "PZCOO-": "c", "H+": "c"},
                                   equilibrium_constant=pz_rxn_52_eq_const)

library.add_LiquidStream_rxn_insta(id="HCO3- + PZCOO- = PZ(COO)2-2 + H2O",
                                   stoch={"HCO3-": -1, "PZCOO-": -1, "PZ(COO)2-2": 1, "H2O": 1},
                                   unit={"HCO3-": "c", "PZCOO-": "c", "PZ(COO)2-2": "c", "H2O": "x"},
                                   equilibrium_constant=pz_rxn_53_eq_const)

library.add_LiquidStream_rxn_insta(id="HOOCPZCOO- = PZ(COO)2-2 + H+",
                                   stoch={"HOOCPZCOO-": -1, "PZ(COO)2-2": 1, "H+": 1},
                                   unit={"HOOCPZCOO-": "c", "PZ(COO)2-2": "c", "H+": "c"},
                                   equilibrium_constant=pz_rxn_54_eq_const)

library.add_LiquidStream_rxn_insta(id="PZ(2H)+ = PZH+ + H+",
                                   stoch={"PZ(2H)+": -1, "PZH+": 1, "H+": 1},
                                   unit={"PZ(2H)+": "c", "PZH+": "c", "H+": "c"},
                                   equilibrium_constant=pz_rxn_55_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + PZ = PZCOO- + H+",
                                   stoch={"CO2": -1, "PZ": -1, "PZCOO-": 1, "H+": 1},
                                   unit={"CO2": "c", "PZ": "c", "PZCOO-": "c", "H+": "c"},
                                   equilibrium_constant=pz_rxn_56_eq_const)

library.add_LiquidStream_rxn_insta(id="CO2 + PZCOO- = PZ(COO)2-2 + H+",
                                   stoch={"CO2": -1, "PZCOO-": -1, "PZ(COO)2-2": 1, "H+": 1},
                                   unit={"CO2": "c", "PZCOO-": "c", "PZ(COO)2-2": "c", "H+": "c"},
                                   equilibrium_constant=pz_rxn_57_eq_const)


library.add_LiquidStream_rxn_reversible(id="CO2 + PZ -> PZCOO- + H+",
                                        stoch={"CO2": -1, "PZ": -1, "PZCOO-": 1, "H+": 1},
                                        dependencies=["CO2", "PZ", "PZCOO-", "H+"],
                                        rate_kmol_m3s=pz_rxn_56_rate_law_kmol_m3s,
                                        exothermic_heat_kJ_kmol=pz_rxn_56_exothermic_heat_kJ_kmol)

library.add_LiquidStream_rxn_reversible(id="CO2 + PZCOO- -> PZ(COO)2-2 + H+",
                                        stoch={"CO2": -1, "PZCOO-": -1, "PZ(COO)2-2": 1, "H+": 1},
                                        dependencies=["CO2", "PZCOO-", "PZ(COO)2-2", "H+"],
                                        rate_kmol_m3s=pz_rxn_57_rate_law_kmol_m3s,
                                        exothermic_heat_kJ_kmol=pz_rxn_57_exothermic_heat_kJ_kmol)





