import numpy as np
import testlablib as lab
from testlablib.library.library import library
import time



# ---------------------------------------------------------------------------------

class ExhaustGas(lab.GasStream):

    def __init__(self, flow_Nm3_h_dry, pressure_bara, temp_K, CO2_pct_dry, H2O_pct):
        super().__init__()
        self.load_viscosity_Pas(id="Exhaust Gas", library=library)
        self.load_diffusivity_m2_s(id="Exhaust Gas", library=library)
        self.load_thermal_conductivity_kW_mK(id="Exhaust Gas", library=library)
        self.load_heat_capacity_kJ_kmolK(id="Exhaust Gas", library=library)
        self.add_specie(id="CO2", library=library)
        self.add_specie(id="O2", library=library)
        self.add_specie(id="H2O", library=library)
        self.add_specie(id="N2", library=library)

        self.set_gas_temp_K(value=temp_K)
        self.set_gas_pressure_bara(value=pressure_bara)
        shape = np.ones(shape=temp_K.shape)
        for id in self.specie.keys():
            self.set_specie_molar_fraction(id=id, value=0 * shape)
        y_CO2_dry = CO2_pct_dry / 100
        y_O2_dry = 12.3 * shape / 100
        y_N2_dry = 1 - y_CO2_dry - y_O2_dry


        if H2O_pct is None:
            p_H2O = self.H2O_vapor_pressure_bara(temp_K=temp_K)
            y_H2O = p_H2O / pressure_bara
        else:
            y_H2O = H2O_pct / 100

        y_CO2_wet = y_CO2_dry * (1 - y_H2O)
        y_O2_wet = y_O2_dry * (1 - y_H2O)
        y_N2_wet = y_N2_dry * (1 - y_H2O)

        self.set_specie_molar_fraction(id="CO2", value=y_CO2_wet)
        self.set_specie_molar_fraction(id="O2", value=y_O2_wet)
        self.set_specie_molar_fraction(id="N2", value=y_N2_wet)
        self.set_specie_molar_fraction(id="H2O", value=y_H2O)

        flow_kmol_h_dry = flow_Nm3_h_dry / (0.08314 * 273.15)
        flow_kmol_h_H2O = flow_kmol_h_dry * y_H2O
        flow_kmol_h_wet = flow_kmol_h_dry + flow_kmol_h_H2O
        self.set_gas_flow_kmol_h(value=flow_kmol_h_wet)

        self.normalize_molar_fractions()

    def H2O_vapor_pressure_bara(self, temp_K):
        T = np.minimum(temp_K, 273.15 + 150)
        pc = 220.64
        Tc = 647.096
        tau = 1 - T / Tc
        a1 = -7.85951783
        a2 = 1.84408259
        a3 = -11.7866497
        a4 = 22.6807411
        a5 = -15.9618719
        a6 = 1.80122502
        p = pc * np.exp((Tc / T) * (a1 * tau + a2 * tau ** 1.5 + a3 * tau ** 3 + a4 * tau ** 3.5 + a5 * tau ** 4 + a6 * tau ** 7.5))
        return p


class ReboilerVapor(lab.GasStream):

    def __init__(self, flow_kmol_h, pressure_bara, temp_K, CO2_molar_fraction, H2O_molar_fraction):
        super().__init__()
        self.load_viscosity_Pas(id="Exhaust Gas", library=library)
        self.load_diffusivity_m2_s(id="Exhaust Gas", library=library)
        self.load_thermal_conductivity_kW_mK(id="Exhaust Gas", library=library)
        self.load_heat_capacity_kJ_kmolK(id="Exhaust Gas", library=library)
        self.add_specie(id="CO2", library=library)
        self.add_specie(id="H2O", library=library)
        self.set_gas_flow_kmol_h(value=flow_kmol_h)
        self.set_gas_temp_K(value=temp_K)
        self.set_gas_pressure_bara(value=pressure_bara)
        self.set_specie_molar_fraction(id="CO2", value=CO2_molar_fraction)
        self.set_specie_molar_fraction(id="H2O", value=H2O_molar_fraction)


class Solvent(lab.LiquidStream):

    def __init__(self, stream_id, temp_K, flow_kg_h, CO2Load, MEA_mass_fraction=None, MDEA_mass_fraction=None, PZ_mass_fraction=None, AMP_mass_fraction=None, K2CO3_mass_fraction=None, KLys_mass_fraction=None, KSar_mass_fraction=None, KGly_mass_fraction=None, NGly_mass_fraction=None):

        super().__init__(stream_id=stream_id, solvent_id="H2O")
        z = np.zeros(shape=temp_K.shape)

        MEA_mass_fraction = z if MEA_mass_fraction is None else MEA_mass_fraction
        MDEA_mass_fraction = z if MDEA_mass_fraction is None else MDEA_mass_fraction
        PZ_mass_fraction = z if PZ_mass_fraction is None else PZ_mass_fraction
        AMP_mass_fraction = z if AMP_mass_fraction is None else AMP_mass_fraction
        K2CO3_mass_fraction = z if K2CO3_mass_fraction is None else K2CO3_mass_fraction
        KLys_mass_fraction = z if KLys_mass_fraction is None else KLys_mass_fraction
        KSar_mass_fraction = z if KSar_mass_fraction is None else KSar_mass_fraction
        KGly_mass_fraction = z if KGly_mass_fraction is None else KGly_mass_fraction
        NGly_mass_fraction = z if NGly_mass_fraction is None else NGly_mass_fraction

        self.load_viscosity(id="Aqueous Solution w/Amines and Amino Acids", library=library)
        self.load_heat_capacity(id="Aqueous Solution w/Amines and Amino Acids", library=library)
        self.load_density(id="Aqueous Solution w/Amines and Amino Acids", library=library)
        self.load_activity_coefficient(id="Truesdell Jones", library=library)
        self.load_diffusivity(id="All", library=library)
        self.load_surface_tension_N_m(id="Water", library=library)

        self.add_info(key="PZ Mass Fraction", value=PZ_mass_fraction)
        self.add_info(key="AMP Mass Fraction", value=AMP_mass_fraction)
        self.add_info(key="MEA Mass Fraction", value=MEA_mass_fraction)
        self.add_info(key="MDEA Mass Fraction", value=MDEA_mass_fraction)
        self.add_info(key="K2CO3 Mass Fraction", value=K2CO3_mass_fraction)
        self.add_info(key="KLys Mass Fraction", value=KLys_mass_fraction)
        self.add_info(key="KSar Mass Fraction", value=KSar_mass_fraction)
        self.add_info(key="KGly Mass Fraction", value=KGly_mass_fraction)
        self.add_info(key="NGly Mass Fraction", value=NGly_mass_fraction)

        if np.max(self.info["MEA Mass Fraction"]) > 0:
            self.add_specie(id="MEA", library=library)
            self.add_specie(id="MEAH+", library=library)
            self.add_specie(id="MEACOO-", library=library)
            self.add_rxn_insta(id="MEAH+ = MEA + H+", library=library)
            self.add_rxn_insta(id="CO2 + MEA = MEACOO- + H+", library=library)
            #self.add_rxn_reversible(id="CO2 + MEA -> MEACOO- + H+", library=library)
            #self.add_vapor_pressure_bara(id="MEA(g) = MEA(aq)", library=library)

        if np.max(self.info["MDEA Mass Fraction"]) > 0:
            self.add_specie(id="MDEA", library=library)
            self.add_specie(id="MDEAH+", library=library)
            self.add_rxn_insta(id="MDEAH+ = MDEA + H+", library=library)
            #self.add_vapor_pressure_bara(id="MDEA(g) = MDEA(aq)", library=library)

        if np.max(self.info["AMP Mass Fraction"]) > 0:
            self.add_specie(id="AMP", library=library)
            self.add_specie(id="AMPH+", library=library)
            self.add_specie(id="AMPCOO-", library=library)
            self.add_rxn_insta(id="AMPH+ = AMP + H+", library=library)
            self.add_rxn_insta(id="CO2 + AMP = AMPCOO- + H+", library=library)
            #self.add_rxn_reversible(id="CO2 + AMP -> AMPCOO- + H+", library=library)

        if np.max(self.info["PZ Mass Fraction"]) > 0:
            self.add_specie(id="PZ", library=library)
            self.add_specie(id="PZH+", library=library)
            self.add_specie(id="HPZCOO", library=library)
            self.add_specie(id="PZCOO-", library=library)
            self.add_specie(id="PZ(COO)2-2", library=library)
            self.add_specie(id="HOOCPZCOO-", library=library)

            self.add_rxn_insta(id="PZH+ = PZ + H+", library=library)
            self.add_rxn_insta(id="HPZCOO = PZCOO- + H+", library=library)
            self.add_rxn_insta(id="HOOCPZCOO- = PZ(COO)2-2 + H+", library=library)

            self.add_rxn_insta(id="CO2 + PZ = PZCOO- + H+", library=library)
            self.add_rxn_insta(id="CO2 + PZCOO- = PZ(COO)2-2 + H+", library=library)

            #self.add_rxn_reversible(id="CO2 + PZ -> PZCOO- + H+", library=library)
            #self.add_rxn_reversible(id="CO2 + PZCOO- -> PZ(COO)2-2 + H+", library=library)

            #self.add_vapor_pressure_bara(id="PZ(g) = PZ(aq)", library=library)

        if np.max(self.info["KLys Mass Fraction"]) > 0:
            #self.add_specie(id="Lys+2", library=library)
            self.add_specie(id="Lys+", library=library)
            self.add_specie(id="Lys", library=library)
            self.add_specie(id="Lys-", library=library)
            self.add_specie(id="LysCOO-", library=library)
            self.add_specie(id="LysCOO-2", library=library)

            #self.add_rxn_insta(id="Lys+2 = Lys+ + H+", library=library)
            self.add_rxn_insta(id="Lys+ = Lys + H+", library=library)
            self.add_rxn_insta(id="Lys = Lys- + H+", library=library)
            self.add_rxn_insta(id="LysCOO- = LysCOO-2 + H+", library=library)
            self.add_rxn_insta(id="CO2 + Lys- = LysCOO-2 + H+", library=library)
            #self.add_rxn_reversible(id="CO2 + Lys- -> H+ + LysCOO-2", library=library)

        if np.max(self.info["KSar Mass Fraction"]) > 0:
            self.add_specie(id="Sar+", library=library)
            self.add_specie(id="Sar", library=library)
            self.add_specie(id="Sar-", library=library)
            self.add_specie(id="SarCOO-2", library=library)
            if "K+" not in self.specie.keys():
                self.add_specie(id="K+", library=library)
            self.add_rxn_insta(id="Sar+ = Sar + H+", library=library)
            self.add_rxn_insta(id="Sar = Sar- + H+", library=library)
            self.add_rxn_insta(id="CO2 + Sar- = H+ + SarCOO-2", library=library)
            #self.add_rxn_reversible(id="CO2 + Sar- -> H+ + SarCOO-2", library=library)

        if np.max(self.info["KGly Mass Fraction"]) > 0 or np.max(self.info["NGly Mass Fraction"]) > 0:
            self.add_specie(id="Gly+", library=library)
            self.add_specie(id="Gly", library=library)
            self.add_specie(id="Gly-", library=library)
            self.add_specie(id="GlyCOO-2", library=library)
            if "K+" not in self.specie.keys() and np.max(self.info["KGly Mass Fraction"]) > 0:
                self.add_specie(id="K+", library=library)
            if "Na+" not in self.specie.keys() and np.max(self.info["NGly Mass Fraction"]) > 0:
                self.add_specie(id="Na+", library=library)
            self.add_rxn_insta(id="Gly+ = Gly + H+", library=library)
            self.add_rxn_insta(id="Gly = Gly- + H+", library=library)
            self.add_rxn_insta(id="CO2 + Gly- = H+ + GlyCOO-2", library=library)
            #self.add_rxn_reversible(id="CO2 + Gly- -> H+ + GlyCOO-2", library=library)

        if np.max(self.info["K2CO3 Mass Fraction"]) > 0:
            if "K'+" not in self.specie.keys():
                self.add_specie(id="K'+", library=library)

        self.add_specie(id="H2O", library=library)
        self.add_specie(id="CO2", library=library)
        self.add_specie(id="HCO3-", library=library)
        self.add_specie(id="CO3-2", library=library)
        self.add_specie(id="H+", library=library)
        self.add_specie(id="OH-", library=library)

        self.add_rxn_insta(id="H2O = H+ + OH-", library=library)
        self.add_rxn_insta(id="HCO3- = CO3-2 + H+", library=library)
        #self.add_rxn_insta(id="CO2 + H2O = HCO3- + H+", library=library)
        self.add_rxn_reversible(id="CO2 + H2O -> HCO3- + H+", library=library)

        self.add_vapor_pressure_bara(id="CO2(g) = CO2(aq)", library=library)
        self.add_vapor_pressure_bara(id="H2O(g) = H2O(l)", library=library)

        # Temperature and Flow
        self.set_solution_temp_K(value=temp_K)
        self.set_solution_flow_kg_h(value=flow_kg_h)

        # Concentrations
        H2O_mass_fraction = 1 - MDEA_mass_fraction - PZ_mass_fraction - AMP_mass_fraction - MEA_mass_fraction - K2CO3_mass_fraction - KLys_mass_fraction - KSar_mass_fraction - KGly_mass_fraction - NGly_mass_fraction

        # Number of moles in one (1) kg of Unloaded Solution
        n_MDEA = MDEA_mass_fraction / 119
        n_PZ = PZ_mass_fraction / 86
        n_MEA = MEA_mass_fraction / 61
        n_AMP = AMP_mass_fraction / 89
        n_K2CO3 = K2CO3_mass_fraction / (2 * 39 + 60)
        n_KLys = KLys_mass_fraction / (39 + 145)
        n_KSar = KSar_mass_fraction / (39 + 88)
        n_KGly = KGly_mass_fraction / (39 + 74)
        n_NGly = NGly_mass_fraction / (23 + 74)

        # Adding CO2. The total mass of the solution is then above 1 kg. This is handled by normalizing weight fractions in the end.
        n_CO2 = CO2Load * (n_MDEA + n_PZ + n_AMP + n_MEA + n_K2CO3 + n_KLys + n_KSar + n_KGly + n_NGly)
        CO2_mass_fraction = 44 * n_CO2

        # Setting all mass fractions to zero initially
        zero = 0 * np.ones(shape=H2O_mass_fraction.shape)
        for id in self.specie.keys():
            self.set_specie_mass_fraction(id=id, value=zero)

        # Mass Fractions
        self.set_specie_mass_fraction(id="CO2", value=CO2_mass_fraction)
        self.set_specie_mass_fraction(id="H2O", value=H2O_mass_fraction)

        if "MDEA" in self.specie.keys():
            self.set_specie_mass_fraction(id="MDEA", value=MDEA_mass_fraction)

        if "PZ" in self.specie.keys():
            self.set_specie_mass_fraction(id="PZ", value=PZ_mass_fraction)
        self.set_specie_mass_fraction(id="CO3-2", value=K2CO3_mass_fraction * 60 / (2*39 + 60))

        if "AMP" in self.specie.keys():
            self.set_specie_mass_fraction(id="AMP", value=AMP_mass_fraction)

        if "MEA" in self.specie.keys():
            self.set_specie_mass_fraction(id="MEA", value=MEA_mass_fraction)

        if "K'+" in self.specie.keys():
            self.set_specie_mass_fraction(id="K'+", value=K2CO3_mass_fraction * 2*39 / (2*39 + 60))

        if "Lys-" in self.specie.keys():
            self.set_specie_mass_fraction(id="Lys-", value=KLys_mass_fraction * 145 / (39 + 145))

        if "Na+" in self.specie.keys():
            self.set_specie_mass_fraction(id="Na+", value=NGly_mass_fraction * 23 / (23 + 74))

        if "K+" in self.specie.keys():
            self.set_specie_mass_fraction(id="K+", value=KLys_mass_fraction * 39 / (39 + 145) + KSar_mass_fraction * 39 / (39 + 88) + KGly_mass_fraction * 39 / (39 + 74))

        if "Sar-" in self.specie.keys():
            self.set_specie_mass_fraction(id="Sar-", value=KSar_mass_fraction * 88 / (39 + 88))

        if "Gly-" in self.specie.keys():
            self.set_specie_mass_fraction(id="Gly-", value=KGly_mass_fraction * 74 / (39 + 74) + NGly_mass_fraction * 74 / (23 + 74))

        self.normalize_mass_fractions()

    def CO2Load(self, LiquidStream):

        x_CO2 = 0
        x_Amine = 0

        # K'+ are Potassium from K2CO3

        C = {"CO2": 1, "HCO3-": 1, "CO3-2": 1,
             "MEACOO-": 1,
             "PZCOO-": 1, "PZ(COO)2-2": 2, "HOOCPZCOO-": 2, "HPZCOO": 1,
             "AMPCOO-": 1,
             "GlyCOO-2": 1,
             "SarCOO-2": 1,
             "LysCOO-2": 1, "LysCOO-": 1,
             "K'+": -0.5}


        A = {"MEA": 1, "MEACOO-": 1, "MEAH+": 1,
             "MDEA": 1, "MDEAH+": 1,
             "PZ": 1, "PZCOO-": 1, "PZ(COO)2-2": 1, "PZH+": 1, "HOOCPZCOO-": 1, "HPZCOO": 1,
             "AMP": 1, "AMPH+": 1, "AMPCOO-": 1,
             "Gly+": 1, "Gly": 1, "Gly-": 1, "GlyCOO-2": 1, "NGly": 1, "KGly": 1,
             "Sar+": 1, "Sar": 1, "Sar-": 1, "SarCOO-2": 1, "KSar": 1,
             "Lys+": 1, "Lys": 1, "Lys-": 1, "LysCOO-": 1, "LysCOO-2": 1, "KLys": 1,
             "K'+": 0.5}

        for c in C.keys():
            if c in LiquidStream.specie.keys():
                x_CO2 = x_CO2 + C[c] * LiquidStream.get_specie_molar_fraction(id=c)
        for a in A.keys():
            if a in LiquidStream.specie.keys():
                x_Amine = x_Amine + A[a] * LiquidStream.get_specie_molar_fraction(id=a)
        alpha = np.maximum(x_CO2 / np.maximum(x_Amine, 10 ** (-9)), 10**(-9))
        return alpha


class MakeUpWater(lab.LiquidStream):

    def __init__(self, temp_K, flow_kg_h):
        super().__init__(stream_id="Water",solvent_id="H2O")

        self.load_density(id="Water", library=library)
        self.load_heat_capacity(id="Water", library=library)
        self.load_viscosity(id="Water", library=library)
        self.load_activity_coefficient(id="Truesdell Jones", library=library)
        self.load_diffusivity(id="All", library=library)
        self.load_surface_tension_N_m(id="Water", library=library)

        self.add_specie(id="H2O", library=library)
        self.add_vapor_pressure_bara(id="H2O(g) = H2O(l)", library=library)
        shape = np.ones(shape=temp_K.shape)
        self.set_solution_flow_kg_h(value=flow_kg_h)
        self.set_solution_temp_K(value=temp_K)
        self.set_specie_mass_fraction(id="H2O", value=1.0 * shape)


# ---------------------------------------------------------------------------------

class __StructuredPacking__():

    def __init__(self):
        pass

    def Liquid_Holdup_m3_m3(self, Column):
        theta = Column.get_corrugation_angle_degree() * (np.pi / 180)
        ap = Column.get_packing_area_m2_m3()
        mu = Column.LiquidStream.get_solution_viscosity_Pas()
        rho = Column.LiquidStream.get_solution_density_kg_m3()
        v_LS = Column.get_superficial_liquid_velocity_m_s()
        g = 9.81
        enu = 3 * np.abs(v_LS) * mu * ap ** 2
        den = rho * g * np.sin(theta) ** 2
        h_liq = (enu / den) ** (1 / 3)
        return h_liq

    def Mass_Transfer_H2O_kmol_m3s(self, Column):

        p_H2O = Column.GasStream.get_specie_pressure_bara(id="H2O")
        p_H2O_vap = Column.LiquidStream.get_specie_vapor_pressure_bara(gas_id="H2O") #/ Column.LiquidStream.get_specie_molar_fraction(id="H2O")
        T = Column.GasStream.get_gas_temp_K()
        R = 0.08314

        kL, kG, kH, ae = self.Transfer_Coefficients(Column=Column, gas_specie_id="H2O", liq_specie_id="H2O")

        # Overall Mass Transfer Coefficient [kmol/m3.s.bar]
        KGa = kG * ae / (R * T)

        # Condensation Rate [kmol/m3.s]
        r = KGa * (p_H2O - p_H2O_vap)
        return r

    def Mass_Transfer_H2O_kJ_kmol(self, Column):
        T0 = Column.LiquidStream.temp_K
        p0 = Column.LiquidStream.vapor_pressure_bara["H2O(g) = H2O(l)"]["p0"](Column.LiquidStream)
        Column.LiquidStream.temp_K = Column.LiquidStream.temp_K + 0.05
        T1 = Column.LiquidStream.temp_K
        p1 = Column.LiquidStream.vapor_pressure_bara["H2O(g) = H2O(l)"]["p0"](Column.LiquidStream)
        Column.LiquidStream.temp_K = Column.LiquidStream.temp_K - 0.05
        h = 8.314 * (np.log(1/p1) - np.log(1/p0)) / ((1 / T1) - (1 / T0))
        return h

    def Mass_Transfer_CO2_kmol_m3s(self, Column):

        # Packing Geometry
        ap = Column.get_packing_area_m2_m3()

        # Various
        D_CO2 = Column.LiquidStream.get_specie_diffusivity_m2_s(id="CO2")
        T_liq = Column.LiquidStream.get_solution_temp_K()
        rho_L = Column.LiquidStream.get_solution_density_kg_m3() / 1000
        sigma = 72 * 10 ** (-3)
        H_CO2 = Column.LiquidStream.vapor_pressure_bara["CO2(g) = CO2(aq)"]["H"](Column.LiquidStream)

        # Constants
        g = 9.81
        R = 0.08314

        # Driving Force
        p_CO2_vap = Column.LiquidStream.get_specie_vapor_pressure_bara(gas_id="CO2")
        p_CO2 = Column.GasStream.get_specie_pressure_bara(id="CO2")

        # Mass Transfer Coefficients
        kL, kG, kH, ae = self.Transfer_Coefficients(Column=Column, gas_specie_id="H2O", liq_specie_id="H2O")

        # Reaction Rate Coefficients
        k = {}
        k["MEA"] = 5993 * np.exp(-5400 * (1 / T_liq - 1 / 298.15))
        k["PZ"] = 58000 * np.exp(-4209 * (1 / T_liq - 1 / 298.15))
        k["PZCOO-"] = k["PZ"] / 4
        k["OH-"] = 8416 * np.exp(-6667 * (1 / T_liq - 1 / 298.15))
        k["AMP"] = 1500 * np.exp(-5176 * (1 / T_liq - 1 / 298.15))
        k["Sar-"] = 15000 * np.exp(-3127 * (1 / T_liq - 1 / 298.15))  # Need Verification....
        k["Gly-"] = 13900 * np.exp(-5459 * (1 / T_liq - 1 / 298.15))  # Need Verification...
        k["Lys-"] = 13900 * np.exp(-5459 * (1 / T_liq - 1 / 298.15))  # Need Verification...

        # Hatta Number
        Ha = {}
        c = {}
        for id in k.keys():
            if id in Column.LiquidStream.specie.keys():
                c[id] = np.maximum(Column.LiquidStream.get_specie_molarity_kmol_m3(id=id), 10**(-18))
            else:
                c[id] = 10**(-18)
            Ha[id] = (1/kL) * np.sqrt(D_CO2 * k[id] * c[id])

        # Enhancement Factor in Case of Instantaneous Reaction
        E_inf = {}
        for id in k.keys():
            E_inf[id] = 1 + (1/2) * c[id] / (p_CO2 * H_CO2)

        # Enhancement Factor
        E = 1
        for id in k.keys():
            E = E + (np.maximum((1 / E_inf[id] + 1 / Ha[id])**(-1), 1.0) - 1)

        #print(E)

        KG = ((R * T_liq / kG) + 1 / (H_CO2 * E * kL)) ** (-1)
        r = KG * ae * (p_CO2 - p_CO2_vap)
        return r

    def Mass_Transfer_CO2_kJ_kmol(self, Column):
        T0 = Column.LiquidStream.temp_K
        H0 = Column.LiquidStream.vapor_pressure_bara["CO2(g) = CO2(aq)"]["H"](Column.LiquidStream)
        Column.LiquidStream.temp_K = Column.LiquidStream.temp_K + 0.05
        T1 = Column.LiquidStream.temp_K
        H1 = Column.LiquidStream.vapor_pressure_bara["CO2(g) = CO2(aq)"]["H"](Column.LiquidStream)
        Column.LiquidStream.temp_K = Column.LiquidStream.temp_K - 0.05
        q = 8.314 * (np.log(H1) - np.log(H0)) / ((1 / T1) - (1 / T0))
        return q

    def Heat_Transfer_kW_m3(self, Column):
        T_gas = Column.GasStream.get_gas_temp_K()
        T_liq = Column.LiquidStream.get_solution_temp_K()
        kL, kG, kH, ae = self.Transfer_Coefficients(Column=Column, gas_specie_id="H2O", liq_specie_id="H2O")
        q = kH * ae * (T_gas - T_liq)
        return q

    def Transfer_Coefficients(self, Column, gas_specie_id, liq_specie_id):

        # Dimensions
        ap = Column.get_packing_area_m2_m3()
        theta = Column.get_corrugation_angle_degree() * (np.pi / 180)
        Mi = ((3 * ap ** 3 * np.sin(theta) * np.cos(theta)) / (16 * (np.sin(theta) ** 2 + 1) ** (3 / 2))) / ap ** 3

        # Constants of Nature
        g = 9.81
        R = 0.08314

        # GasStream
        n_gas = Column.GasStream.get_gas_flow_kmol_h() / 3600
        rho_gas = Column.GasStream.get_gas_density_kg_m3()
        T_gas = Column.GasStream.get_gas_temp_K()
        p_gas = Column.GasStream.get_gas_pressure_bara()
        mu_gas = Column.GasStream.get_gas_viscosity_Pas()
        D_gas = Column.GasStream.get_specie_diffusivity_m2_s(id=gas_specie_id)
        kappa_gas = Column.GasStream.get_gas_thermal_conductivity_kW_mK()
        c_gas = Column.GasStream.get_gas_molarity_kmol_m3()
        cp_gas = Column.GasStream.get_gas_heat_capacity_kJ_kmolK() * c_gas / rho_gas

        # LiquidStream
        m_liq = Column.LiquidStream.get_solution_flow_kg_h() / 3600
        rho_liq = Column.LiquidStream.get_solution_density_kg_m3()
        T_liq = Column.LiquidStream.get_solution_temp_K()
        mu_liq = Column.LiquidStream.get_solution_viscosity_Pas()
        D_liq = Column.LiquidStream.get_specie_diffusivity_m2_s(id=liq_specie_id)
        sigma = Column.LiquidStream.get_solution_surface_tension_N_m()

        # Superficial Gas Velocity [m/s] and Liquid Load [m/s]
        v_gas = Column.get_superficial_gas_velocity_m_s()
        v_liq = Column.get_superficial_liquid_velocity_m_s()

        # Dimensionless Numbers
        Ar = g * rho_liq * (rho_liq - rho_gas) / (mu_liq ** 2 * ap ** 3)        # Archimedes Number
        Re_gas = (rho_gas * v_gas) / (mu_gas * ap)                              # Reynolds Number
        Sc_gas = mu_gas / (D_gas * rho_gas)                                     # Schmidt Number
        Pr_gas = cp_gas * mu_gas / kappa_gas                                    # Prandtl Number
        Re_liq = (rho_liq * v_liq) / (mu_liq * ap)                              # Reynolds Number
        Sc_liq = mu_liq / (D_liq * rho_liq)                                     # Schmidt Number

        Sh_gas = 0.83 * Re_gas**0.58 * Mi**0.30 * Sc_gas**0.5                   # Sherwood Number
        Sh_liq = 1.79 * Re_liq**0.74 * Mi**0.42 * Sc_liq**0.5                   # Sherwood Number
        Nu_gas = 0.83 * Re_gas**0.58 * Mi**0.3 * Pr_gas ** 0.5                  # Nusselt Number

        # Mass & Heat Transfer Coefficients
        kL = Sh_liq * ap * D_liq            # Liquid Mass Transfer Coefficient [m/s]
        kG = Sh_gas * ap * D_gas            # Gas Mass Transfer Coefficient [m/s]
        kH = Nu_gas * ap * kappa_gas        # Heat Transfer Coefficient [kW/m2.K]

        # Effective Interface Area
        ae = ap * 1.41 * ((rho_liq / sigma) * g ** (1 / 3) * (v_liq / ap) ** (4 / 3)) ** 0.116

        return kL, kG, kH, ae


class Absorber(lab.reactor.Column_StructuredPacking):

    def __init__(self, height_m=5.6, cross_sectional_area_m2=0.5, void_fraction_m3_m3=0.98, packing_area_m2_m3=350, corrugation_angle_degree=60):

        super().__init__(height_m=height_m,
                         cross_sectional_area_m2=cross_sectional_area_m2,
                         void_fraction_m3_m3=void_fraction_m3_m3,
                         packing_area_m2_m3=packing_area_m2_m3,
                         corrugation_angle_degree=corrugation_angle_degree)

        self.add_mass_transfer_kmol_m3s(id="CO2(g) -> CO2(aq)",
                                        stoch_gas={"CO2": -1},
                                        stoch_liq={"CO2": 1},
                                        rate_kmol_m3s=__StructuredPacking__().Mass_Transfer_CO2_kmol_m3s,
                                        exothermic_heat_kJ_kmol=__StructuredPacking__().Mass_Transfer_CO2_kJ_kmol)
        
        self.add_mass_transfer_kmol_m3s(id="H2O(g) -> H2O(aq)",
                                        stoch_gas={"H2O": -1},
                                        stoch_liq={"H2O": 1},
                                        rate_kmol_m3s=__StructuredPacking__().Mass_Transfer_H2O_kmol_m3s,
                                        exothermic_heat_kJ_kmol=__StructuredPacking__().Mass_Transfer_H2O_kJ_kmol)

        self.add_heat_transfer_kW_m3(heat_transfer_kW_m3=__StructuredPacking__().Heat_Transfer_kW_m3)
        self.add_liquid_holdup_m3_m3(liquid_holdup_m3_m3=__StructuredPacking__().Liquid_Holdup_m3_m3)
        self.add_pressure_drop_Pa_m(presure_drop_Pa_m=None)


class WashWaterSection(lab.reactor.Column_StructuredPacking):

    def __init__(self, height_m=1.6, cross_sectional_area_m2=0.5, void_fraction_m3_m3=0.98, packing_area_m2_m3=350, corrugation_angle_degree=60):

        super().__init__(height_m=height_m,
                         cross_sectional_area_m2=cross_sectional_area_m2,
                         void_fraction_m3_m3=void_fraction_m3_m3,
                         packing_area_m2_m3=packing_area_m2_m3,
                         corrugation_angle_degree=corrugation_angle_degree)

        self.add_mass_transfer_kmol_m3s(id="H2O(g) -> H2O(aq)",
                                        stoch_gas={"H2O": -1},
                                        stoch_liq={"H2O": 1},
                                        rate_kmol_m3s=__StructuredPacking__().Mass_Transfer_H2O_kmol_m3s,
                                        exothermic_heat_kJ_kmol=__StructuredPacking__().Mass_Transfer_H2O_kJ_kmol)

        self.add_heat_transfer_kW_m3(heat_transfer_kW_m3=__StructuredPacking__().Heat_Transfer_kW_m3)
        self.add_liquid_holdup_m3_m3(liquid_holdup_m3_m3=__StructuredPacking__().Liquid_Holdup_m3_m3)
        self.add_pressure_drop_Pa_m(presure_drop_Pa_m=None)


class Stripper(lab.reactor.Column_StructuredPacking):

    def __init__(self, height_m=3.4, cross_sectional_area_m2=0.096, void_fraction_m3_m3=0.98, packing_area_m2_m3=250, corrugation_angle_degree=60):
        super().__init__(height_m=height_m,
                         cross_sectional_area_m2=cross_sectional_area_m2,
                         void_fraction_m3_m3=void_fraction_m3_m3,
                         packing_area_m2_m3=packing_area_m2_m3,
                         corrugation_angle_degree=corrugation_angle_degree)

        self.add_mass_transfer_kmol_m3s(id="CO2(g) -> CO2(aq)",
                                        stoch_gas={"CO2": -1},
                                        stoch_liq={"CO2": 1},
                                        rate_kmol_m3s=__StructuredPacking__().Mass_Transfer_CO2_kmol_m3s,
                                        exothermic_heat_kJ_kmol=__StructuredPacking__().Mass_Transfer_CO2_kJ_kmol)

        self.add_mass_transfer_kmol_m3s(id="H2O(g) -> H2O(aq)",
                                        stoch_gas={"H2O": -1},
                                        stoch_liq={"H2O": 1},
                                        rate_kmol_m3s=__StructuredPacking__().Mass_Transfer_H2O_kmol_m3s,
                                        exothermic_heat_kJ_kmol=__StructuredPacking__().Mass_Transfer_H2O_kJ_kmol)

        self.add_heat_transfer_kW_m3(heat_transfer_kW_m3=__StructuredPacking__().Heat_Transfer_kW_m3)
        self.add_liquid_holdup_m3_m3(liquid_holdup_m3_m3=__StructuredPacking__().Liquid_Holdup_m3_m3)
        self.add_pressure_drop_Pa_m(presure_drop_Pa_m=None)


class Rectifier(lab.reactor.Column_StructuredPacking):

    def __init__(self, height_m=0.8, cross_sectional_area_m2=0.096, void_fraction_m3_m3=0.98, packing_area_m2_m3=500, corrugation_angle_degree=60):
        super().__init__(height_m=height_m,
                         cross_sectional_area_m2=cross_sectional_area_m2,
                         void_fraction_m3_m3=void_fraction_m3_m3,
                         packing_area_m2_m3=packing_area_m2_m3,
                         corrugation_angle_degree=corrugation_angle_degree)

        self.add_mass_transfer_kmol_m3s(id="H2O(g) -> H2O(aq)",
                                        stoch_gas={"H2O": -1},
                                        stoch_liq={"H2O": 1},
                                        rate_kmol_m3s=__StructuredPacking__().Mass_Transfer_H2O_kmol_m3s,
                                        exothermic_heat_kJ_kmol=__StructuredPacking__().Mass_Transfer_H2O_kJ_kmol)

        self.add_heat_transfer_kW_m3(heat_transfer_kW_m3=__StructuredPacking__().Heat_Transfer_kW_m3)
        self.add_liquid_holdup_m3_m3(liquid_holdup_m3_m3=__StructuredPacking__().Liquid_Holdup_m3_m3)
        self.add_pressure_drop_Pa_m(presure_drop_Pa_m=None)


# ---------------------------------------------------------------------------------


class HeatExchanger(lab.reactor.PlateHeatExchanger_CounterCurrent):

    def __init__(self, num_of_heights):
        super().__init__(num_of_heights=num_of_heights)
        self.add_heat_transfer_coefficient_kW_m2K(heat_transfer_coefficient_kW_m2K=self.Heat_Transfer_Coefficient_kW_m2K)

    def Heat_Transfer_Coefficient_kW_m2K(self, HX):
        mu1 = HX.LiquidStream1.get_solution_viscosity_Pas()
        mu2 = HX.LiquidStream2.get_solution_viscosity_Pas()
        A = HX.get_interface_area_m2()
        m1 = np.abs(HX.LiquidStream1.get_solution_flow_kg_h())
        m2 = np.abs(HX.LiquidStream2.get_solution_flow_kg_h())
        k1 = 5.3 * (1000 * mu1) ** (-0.5) * (18.52 / A) ** 0.8 * (m1 / 5000) ** 0.8
        k2 = 5.3 * (1000 * mu2) ** (-0.5) * (18.52 / A) ** 0.8 * (m2 / 5000) ** 0.8
        k = (1 / k1 + 1 / k2) ** (-1)
        return k


class Reboiler(lab.reactor.Reboiler_CSTR):

    def __init__(self):
        pass



# ---------------------------------------------------------------------------------


class CCS():

    def __init__(self):
        pass


class CCS_Jarles_WashWater():

    def __init__(self):
        pass






