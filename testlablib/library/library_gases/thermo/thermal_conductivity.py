import numpy as np
import testlablib as lab




def exhaust_gas_thermal_conductivity_kW_mK(GasStream):
    T = GasStream.get_gas_temp_K()
    t = T - 273.15
    k = (23.84429941 + 7.19127866 * (t / 100) - 0.14985262 * (t / 100) ** 2) * 10 ** (-6)
    return k


