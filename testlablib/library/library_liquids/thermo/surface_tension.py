import numpy as np
import testlablib as lab


def surface_tension_H2O_N_m(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    t = T - 273.15
    t = np.maximum(t, 0.0)
    a = 0.075652711
    b = -0.00013936956
    c = -3.0842103 * 10 ** (-7)
    d = 2.7588435 * 10 ** (-10)
    sigma = a + b * t + c * t ** 2 + d * t ** 3
    return sigma



