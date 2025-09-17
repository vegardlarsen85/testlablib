import numpy as np
import testlablib as lab


library = lab.Library()


def CO2_henrys_constant(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    # H_CO2 = 0.0343 * np.exp(2440 * (1 / T - 1 / 298.15))
    # H_CO2 = 1.153 * np.exp((-T*(1713*(1-0.0015453*T)**(1/3) + 3680) + 1198506)/T**2)
    H_CO2 = 0.4572 * np.exp((941290.2 - 3934.4 * T) / T ** 2)
    return H_CO2

def SO2_henrys_constant(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    H_SO2 = 1.23 * np.exp(3241 * (1 / T - 1 / 298.15))
    return H_SO2

def MEA_henrys_constant(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    H_MEA = 2.31 * 10 ** (5) * np.exp(6243 * (1 / T - 1 / 298.15))
    return H_MEA

def MDEA_henrys_constant(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    H_MDEA = 1.89 * 10 ** (6) * np.exp(6191 * (1 / T - 1 / 298.15))
    return H_MDEA

def PZ_henrys_constant(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    H_PZ = 6.91 * 10 ** (5) * np.exp(8901 * (1 / T - 1 / 298.15))
    return H_PZ

def O2_henrys_constant(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    H_O2 = 0.0013 * np.exp(1700 * (1 / T - 1 / 298.15))
    return H_O2

def N2_henrys_constant(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    H_N2 = 0.0006 * np.exp(1300 * (1 / T - 1 / 298.15))
    return H_N2

def H2O_vapor_pressure_bara(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    T = np.minimum(T, 273.15 + 150)
    pc = 220.64
    Tc = 647.096
    tau = 1 - T / Tc
    a1 = -7.85951783
    a2 = 1.84408259
    a3 = -11.7866497
    a4 = 22.6807411
    a5 = -15.9618719
    a6 = 1.80122502
    p = pc * np.exp((Tc / T) * (a1 * tau + a2 * tau ** 1.5 + a3 * tau ** 3 + a4 * tau ** 3.5 + a5 * tau ** 4 + a6 * tau ** 7.5))
    return p


library.add_LiquidStream_vapor_pressure_bara_henry(id="SO2(g) = SO2(aq)",
                                                   gas_id="SO2",
                                                   liq_id="SO2",
                                                   liq_unit="c",
                                                   henrys_coefficient=SO2_henrys_constant)

library.add_LiquidStream_vapor_pressure_bara_henry(id="CO2(g) = CO2(aq)",
                                                   gas_id="CO2",
                                                   liq_id="CO2",
                                                   liq_unit="c",
                                                   henrys_coefficient=CO2_henrys_constant)

library.add_LiquidStream_vapor_pressure_bara_henry(id="MEA(g) = MEA(aq)",
                                                   gas_id="MEA",
                                                   liq_id="MEA",
                                                   liq_unit="m",
                                                   henrys_coefficient=MEA_henrys_constant)

library.add_LiquidStream_vapor_pressure_bara_henry(id="MDEA(g) = MDEA(aq)",
                                                   gas_id="MDEA",
                                                   liq_id="MDEA",
                                                   liq_unit="m",
                                                   henrys_coefficient="MDEA_henrys_constant")

library.add_LiquidStream_vapor_pressure_bara_henry(id="PZ(g) = PZ(aq)",
                                                   gas_id="PZ",
                                                   liq_id="PZ",
                                                   liq_unit="m",
                                                   henrys_coefficient=PZ_henrys_constant)

library.add_LiquidStream_vapor_pressure_bara_henry(id="O2(g) = O2(aq)",
                                                   gas_id="O2",
                                                   liq_id="O2",
                                                   liq_unit="c",
                                                   henrys_coefficient=O2_henrys_constant)

library.add_LiquidStream_vapor_pressure_bara_henry(id="N2(g) = N2(aq)",
                                                   gas_id="N2",
                                                   liq_id="N2",
                                                   liq_unit="c",
                                                   henrys_coefficient=N2_henrys_constant)

library.add_LiquidStream_vapor_pressure_bara_raoult(id="H2O(g) = H2O(l)",
                                                    gas_id="H2O",
                                                    liq_id="H2O",
                                                    pure_vapor_pressure_bara=H2O_vapor_pressure_bara)






