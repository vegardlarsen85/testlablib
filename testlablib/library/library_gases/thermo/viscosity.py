import numpy as np
import testlablib as lab




def exhaust_gas_viscosity_Pas(GasStream):
    T = GasStream.get_gas_temp_K()
    t = T - 273.15
    mu = (17.36589187 + 4.22323528 * (t / 100) - 0.11854257 * (t / 100) ** 2) * 10 ** (-6)
    return mu


