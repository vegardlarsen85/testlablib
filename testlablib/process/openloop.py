import numpy as np
import testlablib as lab
from testlablib.library.library import library
import time


class ExhaustGas(lab.GasStream):

    def __init__(self, flow_Nm3_h_dry, pressure_bara, temp_K, SO2_ppm_dry, CO2_pct_dry, O2_pct_dry, H2O_pct):
        super().__init__()
        self.load_viscosity_Pas(id="Exhaust Gas", library=library)
        self.load_diffusivity_m2_s(id="Exhaust Gas", library=library)
        self.load_thermal_conductivity_kW_mK(id="Exhaust Gas", library=library)
        self.load_heat_capacity_kJ_kmolK(id="Exhaust Gas", library=library)
        self.add_specie(id="CO2", library=library)
        self.add_specie(id="SO2", library=library)
        self.add_specie(id="O2", library=library)
        self.add_specie(id="H2O", library=library)
        self.add_specie(id="N2", library=library)

        self.set_gas_temp_K(value=temp_K)
        self.set_gas_pressure_bara(value=pressure_bara)
        shape = np.ones(shape=temp_K.shape)
        for id in self.specie.keys():
            self.set_specie_molar_fraction(id=id, value=0 * shape)
        y_SO2 = SO2_ppm_dry * 10 ** (-6)
        y_CO2 = CO2_pct_dry / 100
        y_O2 = O2_pct_dry / 100
        y_N2 = 1 - y_SO2 - y_CO2 - y_O2
        self.set_specie_molar_fraction(id="SO2", value=y_SO2)
        self.set_specie_molar_fraction(id="CO2", value=y_CO2)
        self.set_specie_molar_fraction(id="O2", value=y_O2)
        self.set_specie_molar_fraction(id="N2", value=y_N2)
        if H2O_pct is None:
            p_H2O = self.H2O_vapor_pressure_bara(temp_K=temp_K)
            y_H2O = p_H2O / pressure_bara
        else:
            y_H2O = H2O_pct / 100
        flow_kmol_h_dry = flow_Nm3_h_dry / (0.08314 * 273.15)
        flow_kmol_h_H2O = flow_kmol_h_dry * y_H2O
        flow_kmol_h_wet = flow_kmol_h_dry + flow_kmol_h_H2O
        self.set_gas_flow_kmol_h(value=flow_kmol_h_wet)
        self.set_specie_molar_fraction(id="H2O", value=y_H2O)
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


class Seawater(lab.LiquidStream):

    def __init__(self, temp_K, flow_kg_h, alkalinity_umol_kg, SO2_umol_kg):

        super().__init__(stream_id="Seawater",solvent_id="H2O")

        self.load_density(id="Water", library=library)
        self.load_heat_capacity(id="Water", library=library)
        self.load_viscosity(id="Water", library=library)
        self.load_activity_coefficient(id="Truesdell Jones", library=library)
        self.load_diffusivity(id="All", library=library)

        self.add_specie(id="H2O", library=library)
        self.add_specie(id="H+", library=library)
        self.add_specie(id="OH-", library=library)
        self.add_specie(id="CO2", library=library)
        self.add_specie(id="HCO3-", library=library)
        self.add_specie(id="CO3-2", library=library)
        self.add_specie(id="SO2", library=library)
        self.add_specie(id="HSO3-", library=library)
        self.add_specie(id="O2", library=library)
        self.add_specie(id="SO3-2", library=library)
        self.add_specie(id="Na+", library=library)
        self.add_specie(id="Cl-", library=library)
        self.add_specie(id="SO4-2", library=library)
        self.add_specie(id="HSO4-", library=library)
        self.add_specie(id="Mg+2", library=library)
        self.add_specie(id="Ca+2", library=library)

        self.add_rxn_insta(id="SO2 + H2O = HSO3- + H+", library=library)
        self.add_rxn_insta(id="HSO3- = SO3-2 + H+", library=library)
        self.add_rxn_insta(id="HCO3- = CO3-2 + H+", library=library)
        self.add_rxn_insta(id="H2O = H+ + OH-", library=library)
        self.add_rxn_insta(id="HSO4- = SO4-2 + H+", library=library)
        self.add_rxn_reversible(id="CO2 + H2O -> HCO3- + H+", library=library)
        #self.add_rxn_irreversible(id="2SO3-2 + O2 -> 2SO4-2", library=library)
        self.add_vapor_pressure_bara(id="SO2(g) = SO2(aq)", library=library)
        self.add_vapor_pressure_bara(id="H2O(g) = H2O(l)", library=library)

        shape = np.ones(shape=temp_K.shape)
        solutes_molality_mol_kg = {"CO2": 0 * shape,
                                   "HCO3-": alkalinity_umol_kg * 10**(-6),
                                   "CO3-2": 0 * shape,
                                   "SO2": SO2_umol_kg * 10**(-6),
                                   "HSO3-": 0*shape,
                                   "SO3-2": 0*shape,
                                   "H+": 0 * shape,
                                   "OH-": 0 * shape,
                                   "O2": 0.21 * 0.0013 * np.exp(1500 * (1 / temp_K - 1 / 298.15)),
                                   "Na+": (0.469 / 0.0022) * 10 ** (-6) * alkalinity_umol_kg * shape,
                                   "Cl-": (0.546 / 0.0022) * 10 ** (-6) * alkalinity_umol_kg * shape,
                                   "HSO4-":0*shape,
                                   "SO4-2": (0.0282 / 0.0022) * 10 ** (-6) * alkalinity_umol_kg * shape,
                                   "Mg+2": (0.0528 / 0.0022) * 10 ** (-6) * alkalinity_umol_kg * shape,
                                   "Ca+2": (0.0103 / 0.0022) * 10 ** (-6) * alkalinity_umol_kg * shape}

        self.set_solution_flow_kg_h(value=flow_kg_h)
        self.set_solution_temp_K(value=temp_K)
        self.set_species_molality(solutes_molality_mol_kg=solutes_molality_mol_kg)





