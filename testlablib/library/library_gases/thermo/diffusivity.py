import numpy as np
import testlablib as lab



def exhaust_gas_diffusivity_m2_s(GasStream, id):
    T = GasStream.get_gas_temp_K()
    t = T - 273.15
    D_H2O = (20.18437495 + 19.74335108 * (t / 100) + 0.90921174 * (t / 100) ** 2) * 10 ** (-6)
    M = GasStream.get_specie_molar_mass_kg_kmol(id)
    D = D_H2O * np.sqrt(18 / M)
    return D





