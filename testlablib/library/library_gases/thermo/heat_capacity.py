import numpy as np
import testlablib as lab


def exhaust_gas_heat_capacity_kJ_kmolK(GasStream, id):

    T = GasStream.get_gas_temp_K()
    if id == "O2":
        A, B, C, D, E = 29.103, 10.040, 2526.5, 9.356, 1153.8
        Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
    elif id == "N2":
        A, B, C, D, E = 29.105, 8.6149, 1701.6, 0.10347, 909.79
        Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
    elif id == "H2O":
        A, B, C, D, E = 33.363, 26.790, 2610.5, 8.896, 1169
        Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
    elif id == "CO2":
        A, B, C, D, E = 29.370, 34.540, 1428, 26.4, 588
        Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
    else:
        Cp = 30 * np.ones(shape=T.shape)
    return Cp



