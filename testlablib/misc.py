import numpy as np


def fuel2gas(FuelFlow):
    # Testlab Motor
    # Exhaustgas Flow as a function of Fuel Flow
    return (FuelFlow < 122) * (26.5 * FuelFlow + 750) + (FuelFlow >= 122) * (27.5 * FuelFlow + 100)


