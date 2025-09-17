import numpy as np
import matplotlib.pyplot as plt
import testlablib.main as lab


class Column_StructuredPacking(lab.GasLiquidContactor_PFR):

    def __init__(self, height_m, cross_sectional_area_m2, void_fraction_m3_m3, packing_area_m2_m3, corrugation_angle_degree):

        position_m = np.array([0, height_m])
        cross_sectional_area_m2 = np.array([cross_sectional_area_m2, cross_sectional_area_m2])
        void_fraction_m3_m3 = np.array([void_fraction_m3_m3, void_fraction_m3_m3])

        super().__init__(position_m=position_m,
                         cross_sectional_area_m2=cross_sectional_area_m2,
                         void_fraction_m3_m3=void_fraction_m3_m3)

        self.packing_area_m2_m3 = packing_area_m2_m3
        self.corrugation_angle_degree = corrugation_angle_degree

    def react(self, GasStreamIn, LiquidStreamIn, step_size_m, max_epochs=np.inf, lr=0.75, flow="countercurrent"):

        if flow == "countercurrent":
            GasStreamOut, LiquidStreamOut = self.__react_countercurrent__(GasStreamIn=GasStreamIn,
                                                                          LiquidStreamIn=LiquidStreamIn,
                                                                          step_size_m=step_size_m,
                                                                          max_epochs=max_epochs,
                                                                          lr=lr)
        else:
            GasStreamOut, LiquidStreamOut = self.__react_cocurrent__(GasStreamIn=GasStreamIn,
                                                                     LiquidStreamIn=LiquidStreamIn,
                                                                     step_size_m=step_size_m,
                                                                     max_epochs=max_epochs,
                                                                     lr=lr)
        return GasStreamOut, LiquidStreamOut

    def get_packing_area_m2_m3(self):
        return self.packing_area_m2_m3

    def get_corrugation_angle_degree(self):
        return self.corrugation_angle_degree


class Column_RandomPacking_CounterCurrent():

    def __init__(self):
        pass


class Column_RandomPacking_CoCurrent():

    def __init__(self):
        pass


class Column_OpenSpray_CounterCurrent():

    def __init__(self):
        pass


class Column_OpenSpray_CoCurrent():

    def __init__(self):
        pass


class RPB_StructuredPacking_CounterCurrent():

    def __init__(self):
        pass


class RPB_StructuredPacking_CoCurrent():

    def __init__(self):
        pass


class WettedWallColumn_CounterCurrent():

    def __init__(self):
        pass


class WettedWallColumn_CoCurrent():

    def __init__(self):
        pass


class WettedWallColumn_GasMixed():

    def __init__(self):
        pass


class StirredCellReactor():

    def __init__(self):
        pass


class Venturi():

    def __init__(self):
        pass


# ---------------------------------------------------------------------------------


class Reboiler_CSTR():

    def __init__(self):
        pass


class Reboiler_LiquidEquilibrium():

    def __init__(self):
        pass


class NozzleFlash():

    def __init__(self):
        pass


# ---------------------------------------------------------------------------------


class PlateHeatExchanger_CounterCurrent():

    def __init__(self):
        pass


# ---------------------------------------------------------------------------------



