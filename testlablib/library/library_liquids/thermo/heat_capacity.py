import numpy as np
import testlablib as lab



def heat_capacity_H2O_kJ_kgK(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    t = T - 273.15
    t = np.maximum(t, 0.0)
    a = 4.2174356
    b = -0.0056181625
    c = 0.0012992528
    d = -0.00011535353
    e = 4.14964 * 10 ** (-6)
    cp = a + b * t + c * t ** 1.5 + d * t ** 2 + e * t ** 2.5
    return cp

def heat_capacity_MEA_kJ_kgK(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    # t = T - 273.15
    # cp = 2.75 + 3.47 * 10**(-3) * (t - 25)
    cp = (115228 + 99.98 * T + 0.231 * T ** 2) / (1000 * 61)
    return cp

def heat_capacity_MDEA_kJ_kgK(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    t = T - 273.15
    cp = 2.28 + 5.4 * 10 ** (-3) * (t - 30)
    return cp

def heat_capacity_PZ_kJ_kgK(LiquidStream):
    return 2.76

def heat_capacity_AMP_kJ_kgK(LiquidStream):
    return 2.58

def heat_capacity_CO2_kJ_kgK(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    A, B, C, D, E = 29.370, 34.540, 1428, 26.4, 588
    cp_CO2_gas = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
    cp_CO2_gas = cp_CO2_gas / 44
    return cp_CO2_gas

def heat_capacity_amino_kJ_kgK(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    t = T - 273.15
    cp_salt = 1.32 + 8.53 * 10 ** (-3) * t
    # ------------------------------------------------------------------------------------------
    cp = {}
    cp["CO2"] = heat_capacity_CO2_kJ_kgK(LiquidStream)
    cp["H2O"] = heat_capacity_H2O_kJ_kgK(LiquidStream)
    cp["H+"] = cp["H2O"]
    cp["OH-"] = cp["H2O"]
    cp["Na+"] = cp_salt
    cp["K+"] = cp_salt
    cp["K'+"] = cp_salt
    cp["Ca+2"] = cp_salt
    cp["Ci"] = cp_salt
    cp["Ci-"] = cp_salt
    cp["Ci-2"] = cp_salt
    cp["Ci-3"] = cp_salt
    cp["Cl-"] = cp_salt
    cp["H2O2"] = cp_salt
    cp["HSO3-"] = cp_salt
    cp["HSO4-"] = cp_salt

    cp["MEA"] = (115228 + 99.98 * T + 0.231 * T ** 2) / (1000 * 61)
    cp["MDEA"] = 2.28 + 5.4 * 10 ** (-3) * (t - 30)
    cp["PZ"] = 2.76
    cp["AMP"] = 2.58
    # ------------------------------------------------------------------------------------------
    cp["MEAH+"] = (1 / 62) * (61 * cp["MEA"] + 1 * cp["H+"])
    cp["MEACOO-"] = (1 / 104) * (61 * cp["MEA"] + 44 * cp["CO2"] - 1 * cp["H+"])
    # ------------------------------------------------------------------------------------------
    cp["MDEAH+"] = (1 / 120) * (119 * cp["MDEA"] + 1 * cp["H+"])
    # ------------------------------------------------------------------------------------------
    cp["AMPH+"] = (1 / 90) * (89 * cp["AMP"] + 1 * cp["H+"])
    cp["AMPCOO-"] = (1 / 132) * (89 * cp["AMP"] + 44 * cp["CO2"] - 1 * cp["H+"])
    # ------------------------------------------------------------------------------------------
    cp["PZH+"] = (1 / 87) * (86 * cp["PZ"] + 1 * cp["H+"])
    cp["PZ(2H)+"] = (1 / 88) * (86 * cp["PZ"] + 2 * cp["H+"])
    cp["PZCOO-"] = (1 / 129) * (86 * cp["PZ"] + 44 * cp["CO2"] - 1 * cp["H+"])
    cp["PZ(COO)2-2"] = (1 / 172) * (86 * cp["PZ"] + 88 * cp["CO2"] - 2 * cp["H+"])
    cp["HOOCPZCOO-"] = (1 / 173) * (86 * cp["PZ"] + 88 * cp["CO2"] - 1 * cp["H+"])
    # ------------------------------------------------------------------------------------------
    cp["HCO3-"] = (1 / 61) * (44 * cp["CO2"] + 18 * cp["H2O"] - 1 * cp["H+"])
    cp["CO3-2"] = (1 / 60) * (44 * cp["CO2"] + 18 * cp["H2O"] - 2 * cp["H+"])
    # ------------------------------------------------------------------------------------------
    cp["Sar-"] = cp_salt
    cp["Sar"] = (1 / 89) * (88 * cp["Sar-"] + 1 * cp["H+"])
    cp["Sar+"] = (1 / 90) * (88 * cp["Sar-"] + 2 * cp["H+"])
    cp["SarCOO-2"] = (1 / 131) * (88 * cp["Sar-"] + 44 * cp["CO2"] - 1 * cp["H+"])
    # ------------------------------------------------------------------------------------------
    cp["Gly-"] = cp_salt
    cp["Gly"] = (1 / 75) * (74 * cp["Gly-"] + 1 * cp["H+"])
    cp["Gly+"] = (1 / 76) * (74 * cp["Gly-"] + 2 * cp["H+"])
    cp["GlyCOO-2"] = (1 / 117) * (74 * cp["Gly-"] + 44 * cp["CO2"] - 1 * cp["H+"])
    # ------------------------------------------------------------------------------------------
    cp["Lys-"] = cp_salt
    cp["Lys"] = (1 / 146) * (145 * cp["Lys-"] + 1 * cp["H+"])
    cp["Lys+"] = (1 / 147) * (145 * cp["Lys-"] + 2 * cp["H+"])
    cp["Lys+2"] = (1 / 148) * (145 * cp["Lys-"] + 3 * cp["H+"])
    cp["LysCOO-2"] = (1 / 188) * (145 * cp["Lys-"] + 44 * cp["CO2"] - 1 * cp["H+"])
    cp["LysCOO-"] = (1 / 189) * (145 * cp["Lys-"] + 44 * cp["CO2"])
    # ------------------------------------------------------------------------------------------

    cp_tot = 0
    for id in LiquidStream.specie.keys():
        if id in cp.keys():
            cp_tot = cp_tot + cp[id] * LiquidStream.get_specie_mass_fraction(id=id)
        else:
            cp_tot = cp_tot + cp_salt * LiquidStream.get_specie_mass_fraction(id=id)

    return cp_tot



