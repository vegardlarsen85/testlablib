import numpy as np
import testlablib as lab



def __viscosity_H2O_Pas__(T):
    t = T - 273.15
    t = np.maximum(t, 1.0)
    a = 557.82468
    b = 19.408782
    c = 0.1360459
    d = -3.1160832 * 10 ** (-4)
    mu = 1 / (a + b * t + c * t ** 2 + d * t ** 3)
    return mu

def viscosity_H2O_Pas(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    mu = __viscosity_H2O_Pas__(T=T)
    return mu

def viscosity_amino_Pas(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    alpha = LiquidStream.CO2Load(LiquidStream)
    shape = np.ones(shape=T.shape)

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

    # Total Weight Fraction of Amines & Amino Acids
    w_Amino = 0
    for id in w.keys():
        w_Amino = w_Amino + w[id]

    # Amine and Amino Acid Weight Fraction (w.r.t Total Amines) in Solution
    r = {}
    for id in w.keys():
        r[id] = w[id] / w_Amino

    # Viscosity Parameters
    k = {}
    k["MEA"] = [7.93907867, 1.87968826, 2.27906571]
    k["MDEA"] = [6.32736485, 1.55625674, 6.09127699]
    k["AMP"] = [5.44103233, 1.0014569, 2.92675613]
    k["PZ"] = [6.40305716, 1.41253454, 3.60395713]
    k["KSer"] = [2.65166181, 1.02116664, 1.60266472]
    k["KSar"] = [4.11521217, 1.28617671, 2.13611497]
    k["KLys"] = [8.02967125, 1.80224515, 2.75306825]
    k["NGly"] = [7.93907867, 1.87968826, 2.27906571]                      # Guess. Assuming same as MEA
    k["KGly"] = [7.93907867, 1.87968826, 2.27906571]                      # Guess. Assuming same as MEA

    # Viscosity of Individual Amines
    mu = {}
    for id in w.keys():
        mu[id] = 0.9 * np.exp(k[id][0] * w_Amino ** k[id][1] + 1000 * k[id][2] * (1 / T - 1 / 298) + 0.80851878 * alpha ** 0.55112124)

    # Viscosity of Mixture
    mu_solution = 0
    for id in w.keys():
        mu_solution = mu_solution + mu[id] * r[id]
    return mu_solution / 1000


