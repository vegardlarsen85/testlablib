import numpy as np
import testlablib as lab


library = lab.Library()



def nox_rxn_1_eq_const(GasStream):
    T = GasStream.get_gas_temp_K()
    K = 4.7 * 10 ** (-31) * np.exp(-21710 * (1 / T - 1 / 298.15))
    return K

def nox_rxn_2_eq_const(GasStream):
    T = GasStream.get_gas_temp_K()
    K = 2.23 * 10 ** 12 * np.exp(13729 * (1 / T - 1 / 298.15))
    return K

def nox_rxn_3_eq_const(GasStream):
    T = GasStream.get_gas_temp_K()
    K = 6.740719533987337 * np.exp(6879.96 * (1 / T - 1 / 298.15))
    return K

def nox_rxn_4_eq_const(GasStream):
    T = GasStream.get_gas_temp_K()
    K = 0.5349121477524398 * np.exp(4776.059 * (1 / T - 1 / 298.15))
    return K

def nox_rxn_5_eq_const(GasStream):
    T = GasStream.get_gas_temp_K()
    K = 6.054816450642387 * 10 ** (-6) * np.exp(1218.128 * (1 / T - 1 / 298.15))
    return K

def nox_rxn_6_eq_const(GasStream):
    T = GasStream.get_gas_temp_K()
    K = 0.5925728685473731 * np.exp(-93.548231 * (1 / T - 1 / 298.15))
    return K


library.add_GasStream_rxn_insta(id="O2 + N2 = 2NO",
                                stoch={"O2": -1, "N2": -1, "NO": 2},
                                equilibrium_constant=nox_rxn_1_eq_const)

library.add_GasStream_rxn_insta(id="2NO + O2 = 2NO2",
                                stoch={"NO": -2, "O2": -1, "NO2": 2},
                                equilibrium_constant=nox_rxn_2_eq_const)

library.add_GasStream_rxn_insta(id="2NO2 = N2O4",
                                stoch={"NO2": -2, "N2O4": 1},
                                equilibrium_constant=nox_rxn_3_eq_const)

library.add_GasStream_rxn_insta(id="NO2 + NO = N2O3",
                                stoch={"NO2": -1, "NO": -1, "N2O3": 1},
                                equilibrium_constant=nox_rxn_4_eq_const)

library.add_GasStream_rxn_insta(id="2NO = N2O2",
                                stoch={"NO": -2, "N2O2": 1},
                                equilibrium_constant=nox_rxn_5_eq_const)

library.add_GasStream_rxn_insta(id="N2O3 + H2O = 2HNO2",
                                stoch={"N2O3": -1, "H2O": -1, "HNO2": 2},
                                equilibrium_constant=nox_rxn_6_eq_const)





