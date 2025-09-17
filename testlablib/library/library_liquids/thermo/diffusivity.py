import numpy as np
import testlablib as lab


def diffusivity_m2_s(LiquidStream, id):
    D_in_H2O_at_298K = {}
    D_in_H2O_at_298K["H2O"] = 2.27 * 10 ** (-9)
    D_in_H2O_at_298K["H+"] = 9.311 * 10 ** (-9)
    D_in_H2O_at_298K["OH-"] = 4.56 * 10 ** (-9)
    D_in_H2O_at_298K["SO2"] = 1.83 * 10 ** (-9)
    D_in_H2O_at_298K["HSO3-"] = 1.33 * 10 ** (-9)
    D_in_H2O_at_298K["SO3-2"] = 0.959 * 10 ** (-9)
    D_in_H2O_at_298K["HSO4-"] = 1.33 * 10 ** (-9)
    D_in_H2O_at_298K["SO4-2"] = 1.07 * 10 ** (-9)
    D_in_H2O_at_298K["CO2"] = 1.92 * 10 ** (-9)
    D_in_H2O_at_298K["HCO3-"] = 1.19 * 10 ** (-9)
    D_in_H2O_at_298K["CO3-2"] = 0.92 * 10 ** (-9)
    D_in_H2O_at_298K["O2"] = 2.01 * 10 ** (-9)
    D_in_H2O_at_298K["Ca+2"] = 0.792 * 10 ** (-9)
    D_in_H2O_at_298K["Mg+2"] = 0.706 * 10 ** (-9)
    D_in_H2O_at_298K["Na+"] = 1.334 * 10 ** (-9)
    D_in_H2O_at_298K["NH4+"] = 1.957 * 10 ** (-9)
    D_in_H2O_at_298K["Cl-"] = 2.03 * 10 ** (-9)
    D_in_H2O_at_298K["MEA"] = 1.078025911 * 10 ** (-9)
    D_in_H2O_at_298K["MEAH+"] = 0.75 * 10 ** (-9)
    D_in_H2O_at_298K["MEACOO-"] = 0.75 * 10 ** (-9)

    # Temperature Compensation
    T = LiquidStream.get_solution_temp_K()
    tempComp = 10 ** (-8.1764) * 10 ** (712.5 / T - 2.591 * 10 ** 5 / T ** 2) / (1.92 * 10 ** (-9))
    D_in_H2O = D_in_H2O_at_298K[id] * tempComp if id in D_in_H2O_at_298K.keys() else 1.0 * 10 ** (-9) * tempComp

    # Viscosity of Pure Water
    t = T - 273.15
    t = np.maximum(t, 1.0)
    a = 557.82468
    b = 19.408782
    c = 0.1360459
    d = -3.1160832 * 10 ** (-4)
    mu_H2O = 1 / (a + b * t + c * t ** 2 + d * t ** 3)

    # Viscosity Compensation
    mu_solution = LiquidStream.get_solution_viscosity_Pas()
    D_in_solution = D_in_H2O * (mu_H2O / mu_solution) ** 0.8

    return D_in_solution


