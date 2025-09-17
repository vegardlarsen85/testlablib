import numpy as np
import testlablib as lab


def density_H2O_kg_m3(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    Tc = 647.096  # Critical Temp [K]
    rhoc = 322  # Critical Density [kg/m3]
    tau = 1 - T / Tc
    b1 = 1.99274064
    b2 = 1.09965342
    b3 = -0.510839303
    b4 = -1.75493479
    b5 = -45.5170352
    b6 = -6.74694450 * 10 ** 5
    rho = rhoc * (1 + b1 * tau ** (1 / 3) + b2 * tau ** (2 / 3) + b3 * tau ** (5 / 3) + b4 * tau ** (16 / 3) + b5 * tau ** (43 / 3) + b6 * tau ** (110 / 3))
    return rho

def density_amino_kg_m3(LiquidStream):

    T = LiquidStream.get_solution_temp_K()
    alpha = LiquidStream.CO2Load(LiquidStream)
    shape = np.ones(shape=T.shape)

    # Unloaded Mass Fractions
    w = {}
    w["MEA"] = LiquidStream.info["MEA Mass Fraction"] if "MEA Mass Fraction" in LiquidStream.info.keys() else 0 * shape
    w["MDEA"] = LiquidStream.info["MDEA Mass Fraction"] if "MDEA Mass Fraction" in LiquidStream.info.keys() else 0 * shape
    w["AMP"] = LiquidStream.info["AMP Mass Fraction"] if "AMP Mass Fraction" in LiquidStream.info.keys() else 0 * shape
    w["PZ"] = LiquidStream.info["PZ Mass Fraction"] if "PZ Mass Fraction" in LiquidStream.info.keys() else 0 * shape
    w["KSer"] = LiquidStream.info["KSer Mass Fraction"] if "KSer Mass Fraction" in LiquidStream.info.keys() else 0 * shape
    w["KSar"] = LiquidStream.info["KSar Mass Fraction"] if "KSar Mass Fraction" in LiquidStream.info.keys() else 0 * shape
    w["KLys"] = LiquidStream.info["KLys Mass Fraction"] if "KLys Mass Fraction" in LiquidStream.info.keys() else 0 * shape
    w["NGly"] = LiquidStream.info["NGly Mass Fraction"] if "NGly Mass Fraction" in LiquidStream.info.keys() else 0 * shape
    w["KGly"] = LiquidStream.info["KGly Mass Fraction"] if "KGly Mass Fraction" in LiquidStream.info.keys() else 0 * shape
    w["K2CO3"] = LiquidStream.info["K2CO3 Mass Fraction"] if "K2CO3 Mass Fraction" in LiquidStream.info.keys() else 0 * shape


    # Total Mass Fractions of Amines and Aminos
    w_Amino = 0
    for id in w.keys():
        w_Amino = w_Amino + w[id]

    # Mass Ratio of the Amines and Aminos w.r.t. Overall Amine/Amino Mass Fraction
    r = {}
    for id in w.keys():
        r[id] = w[id] / w_Amino

    # Mass Fraction of Water
    w["H2O"] = 1 - w_Amino

    # Molar Masses
    M = {}
    M["MEA"] = 61
    M["PZ"] = 86
    M["AMP"] = 89
    M["MDEA"] = 119
    M["KSar"] = 39 + 88
    M["KLys"] = 39 + 145
    M["NGly"] = 23 + 74
    M["KGly"] = 39 + 74
    M["KSer"] = 39 + 104
    M["K2CO3"] = 2 * 39 + 12 + 3 * 16

    # Average Molar Mass of Amine/Aminos
    M_Amino = 0
    for id in M.keys():
        M_Amino = M_Amino + r[id] * M[id]

    # Density in pure form
    rho = {}
    rho["MEA"] = 1000 + ((1000 - 920) / (310 - 410)) * (T - 310)
    rho["MDEA"] = 1035 + ((1035 - 920) / (300 - 450)) * (T - 300)
    rho["AMP"] = 920 + ((920 - 868) / (310 - 370)) * (T - 310)
    rho["PZ"] = 1000 + ((1000 - 950) / (313 - 393)) * (T - 313)
    rho["KSar"] = 1534 - 1.075 * (T - 298.15)
    rho["KLys"] = 1356 - 1.075 * (T - 298.15)
    rho["NGly"] = 1473 - 1.075 * (T - 298.15)
    rho["KGly"] = 1473 - 1.075 * (T - 298.15)       # Assumed to be same as NGly
    rho["KSer"] = 1551 - 1.075 * (T - 298.15)
    rho["K2CO3"] = 1473 - 1.075 * (T - 298.15)      # Assumed to be same as NGly
    rho["H2O"] = density_H2O_kg_m3(LiquidStream)

    # Density of Ideal, Unloaded Solution
    rho_ideal = 0
    for id in rho.keys():
        rho_ideal = rho_ideal + w[id] / rho[id]
    rho_ideal = rho_ideal ** (-1)

    # Density of Non-Ideal, Unloaded Solution
    rho_unloaded = rho_ideal + 198 * w["H2O"] ** 1.955 * w_Amino ** 1.34

    # Density of Loaded Solution
    rho_loaded = rho_unloaded + 37.76518 * 1000 * alpha * w_Amino / M_Amino
    return rho_loaded


