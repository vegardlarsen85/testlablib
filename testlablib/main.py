import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from sklearn import neighbors
from sklearn.neighbors import BallTree
from sklearn.neighbors import KDTree
from sklearn.neighbors import KNeighborsRegressor
import time
import scipy.special
from scipy.integrate import odeint
import pickle
from scipy.integrate import solve_ivp


class Library():

    def __init__(self):

        self.LiquidStream_specie_sol = {}
        self.LiquidStream_specie_liq = {}
        self.GasStream_specie_gas = {}

        self.LiquidStream_rxn_insta = {}
        self.LiquidStream_rxn_reversible = {}
        self.LiquidStream_rxn_irreversible = {}

        self.GasStream_rxn_insta = {}
        self.GasStream_rxn_reversible = {}
        self.GasStream_rxn_irreversible = {}

        self.LiquidStream_density_kg_m3 = {}
        self.LiquidStream_viscosity_Pas = {}
        self.LiquidStream_heat_capacity_kJ_kgK = {}
        self.LiquidStream_thermal_conductivity_kW_mK = {}
        self.LiquidStream_activity_coefficient = {}
        self.LiquidStream_diffusivity_m2_s = {}
        self.LiquidStream_surface_tension_N_m = {}
        self.LiquidStream_vapor_pressure_bara = {}

        self.GasStream_viscosity_Pas = {}
        self.GasStream_thermal_conductivity_kW_mK = {}
        self.GasStream_heat_capacity_kJ_kmolK = {}
        self.GasStream_diffusivity_m2_s = {}

    def append(self, library):

        self.LiquidStream_specie_sol = {**self.LiquidStream_specie_sol, **library.LiquidStream_specie_sol}
        self.LiquidStream_specie_liq = {**self.LiquidStream_specie_liq, **library.LiquidStream_specie_liq}
        self.GasStream_specie_gas = {**self.GasStream_specie_gas, **library.GasStream_specie_gas}

        self.LiquidStream_rxn_insta = {**self.LiquidStream_rxn_insta, **library.LiquidStream_rxn_insta}
        self.LiquidStream_rxn_reversible = {**self.LiquidStream_rxn_reversible, **library.LiquidStream_rxn_reversible}
        self.LiquidStream_rxn_irreversible = {**self.LiquidStream_rxn_irreversible, **library.LiquidStream_rxn_irreversible}

        self.GasStream_rxn_insta = {**self.GasStream_rxn_insta, **library.GasStream_rxn_insta}
        self.GasStream_rxn_reversible = {**self.GasStream_rxn_reversible, **library.GasStream_rxn_reversible}
        self.GasStream_rxn_irreversible = {**self.GasStream_rxn_irreversible, **library.GasStream_rxn_irreversible}

        self.LiquidStream_density_kg_m3 = {**self.LiquidStream_density_kg_m3, **library.LiquidStream_density_kg_m3}
        self.LiquidStream_viscosity_Pas = {**self.LiquidStream_viscosity_Pas, **library.LiquidStream_viscosity_Pas}
        self.LiquidStream_heat_capacity_kJ_kgK = {**self.LiquidStream_heat_capacity_kJ_kgK, **library.LiquidStream_heat_capacity_kJ_kgK}
        self.LiquidStream_thermal_conductivity_kW_mK = {**self.LiquidStream_thermal_conductivity_kW_mK, **library.LiquidStream_thermal_conductivity_kW_mK}
        self.LiquidStream_activity_coefficient = {**self.LiquidStream_activity_coefficient, **library.LiquidStream_activity_coefficient}
        self.LiquidStream_diffusivity_m2_s = {**self.LiquidStream_diffusivity_m2_s, **library.LiquidStream_diffusivity_m2_s}
        self.LiquidStream_surface_tension_N_m = {**self.LiquidStream_surface_tension_N_m, **library.LiquidStream_surface_tension_N_m}
        self.LiquidStream_vapor_pressure_bara = {**self.LiquidStream_vapor_pressure_bara, **library.LiquidStream_vapor_pressure_bara}

        self.GasStream_viscosity_Pas = {**self.GasStream_viscosity_Pas, **library.GasStream_viscosity_Pas}
        self.GasStream_thermal_conductivity_kW_mK = {**self.GasStream_thermal_conductivity_kW_mK, **library.GasStream_thermal_conductivity_kW_mK}
        self.GasStream_heat_capacity_kJ_kmolK = {**self.GasStream_heat_capacity_kJ_kmolK, **library.GasStream_heat_capacity_kJ_kmolK}
        self.GasStream_diffusivity_m2_s = {**self.GasStream_diffusivity_m2_s, **library.GasStream_diffusivity_m2_s}

    def add_LiquidStream_specie_sol(self, id, molar_mass_kg_kmol, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None):
        self.LiquidStream_specie_sol[id] = {}
        self.LiquidStream_specie_sol[id]["Molar Mass [kg/kmol]"] = molar_mass_kg_kmol
        self.LiquidStream_specie_sol[id]["Enthalpy of Formation [kJ/mol]"] = enthalpy_of_formation_kJ_mol
        self.LiquidStream_specie_sol[id]["Entropy [J/mol.K]"] = entropy_J_molK
        self.LiquidStream_specie_sol[id]["Gibbs free Energy of Formation [kJ/mol]"] = gibbs_free_energy_of_formation_kJ_mol

    def add_LiquidStream_specie_liq(self, id, molar_mass_kg_kmol, charge, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None):
        self.LiquidStream_specie_liq[id] = {}
        self.LiquidStream_specie_liq[id]["Molar Mass [kg/kmol]"] = molar_mass_kg_kmol
        self.LiquidStream_specie_liq[id]["Charge"] = charge
        self.LiquidStream_specie_liq[id]["Enthalpy of Formation [kJ/mol]"] = enthalpy_of_formation_kJ_mol
        self.LiquidStream_specie_liq[id]["Entropy [J/mol.K]"] = entropy_J_molK
        self.LiquidStream_specie_liq[id]["Gibbs free Energy of Formation [kJ/mol]"] = gibbs_free_energy_of_formation_kJ_mol

    def add_LiquidStream_rxn_insta(self, id, stoch, unit, equilibrium_constant):
        self.LiquidStream_rxn_insta[id] = {}
        self.LiquidStream_rxn_insta[id]["Stoch"] = stoch
        self.LiquidStream_rxn_insta[id]["Unit"] = unit
        self.LiquidStream_rxn_insta[id]["K"] = equilibrium_constant

    def add_LiquidStream_rxn_reversible(self, id, stoch, rate_kmol_m3s, dependencies, exothermic_heat_kJ_kmol):
        self.LiquidStream_rxn_reversible[id] = {}
        self.LiquidStream_rxn_reversible[id]["Stoch"] = stoch
        self.LiquidStream_rxn_reversible[id]["Rate [kmol/m3.s]"] = rate_kmol_m3s
        self.LiquidStream_rxn_reversible[id]["Dependencies"] = dependencies
        self.LiquidStream_rxn_reversible[id]["Exothermic Heat [kJ/kmol]"] = exothermic_heat_kJ_kmol

    def add_LiquidStream_rxn_irreversible(self, id, stoch, rate_kmol_m3s, dependencies, exothermic_heat_kJ_kmol):
        pass

    def add_LiquidStream_density_kg_m3(self, id, function):
        self.LiquidStream_density_kg_m3[id] = function

    def add_LiquidStream_heat_capacity_kJ_kgK(self, id, function):
        self.LiquidStream_heat_capacity_kJ_kgK[id] = function

    def add_LiquidStream_viscosity_Pas(self, id, function):
        self.LiquidStream_viscosity_Pas[id] = function

    def add_LiquidStream_thermal_conductivity_kW_mK(self, id, function):
        self.thermal_conductivity[function_id] = function

    def add_LiquidStream_activity_coefficient(self, id, function):
        self.LiquidStream_activity_coefficient[id] = function

    def add_LiquidStream_diffusivity_m2_s(self, id, function):
        self.LiquidStream_diffusivity_m2_s[id] = function

    def add_LiquidStream_surface_tension_N_m(self, id, function):
        self.LiquidStream_surface_tension_N_m[id] = function

    def add_LiquidStream_vapor_pressure_bara_henry(self, id, gas_id, liq_id, liq_unit, henrys_coefficient):
        self.LiquidStream_vapor_pressure_bara[id] = {}
        self.LiquidStream_vapor_pressure_bara[id]["Stoch Gas"] = {gas_id: -1}
        self.LiquidStream_vapor_pressure_bara[id]["Stoch Liq"] = {liq_id: 1}
        self.LiquidStream_vapor_pressure_bara[id]["Unit Liq"] = {liq_id: liq_unit}
        self.LiquidStream_vapor_pressure_bara[id]["H"] = henrys_coefficient
        self.LiquidStream_vapor_pressure_bara[id]["p0"] = None

    def add_LiquidStream_vapor_pressure_bara_raoult(self, id, gas_id, liq_id, pure_vapor_pressure_bara):
        self.LiquidStream_vapor_pressure_bara[id] = {}
        self.LiquidStream_vapor_pressure_bara[id]["Stoch Gas"] = {gas_id: -1}
        self.LiquidStream_vapor_pressure_bara[id]["Stoch Liq"] = {liq_id: 1}
        self.LiquidStream_vapor_pressure_bara[id]["Unit Liq"] = {liq_id: "x"}
        self.LiquidStream_vapor_pressure_bara[id]["H"] = None
        self.LiquidStream_vapor_pressure_bara[id]["p0"] = pure_vapor_pressure_bara

    # ------------------------------------------------------

    def add_GasStream_specie_gas(self, id, molar_mass_kg_kmol, charge, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None):
        self.GasStream_specie_gas[id] = {}
        self.GasStream_specie_gas[id]["Molar Mass [kg/kmol]"] = molar_mass_kg_kmol
        self.GasStream_specie_gas[id]["Charge"] = charge
        self.GasStream_specie_gas[id]["Enthalpy of Formation [kJ/mol]"] = enthalpy_of_formation_kJ_mol
        self.GasStream_specie_gas[id]["Entropy [J/mol.K]"] = entropy_J_molK
        self.GasStream_specie_gas[id]["Gibbs free Energy of Formation [kJ/mol]"] = gibbs_free_energy_of_formation_kJ_mol

    def add_GasStream_rxn_insta(self, id, stoch, equilibrium_constant):
        self.GasStream_rxn_insta[id] = {}
        self.GasStream_rxn_insta[id]["Stoch"] = stoch
        self.GasStream_rxn_insta[id]["K"] = equilibrium_constant

    def add_GasStream_rxn_reversible(self, id, stoch, rate_kmol_m3s):
        self.GasStream_rxn_reversible[id] = {}
        self.GasStream_rxn_reversible[id]["Stoch"] = stoch
        self.GasStream_rxn_reversible[id]["Rate [kmol/m3.s]"] = rate_kmol_m3s

    def add_GasStream_rxn_irreversible(self, id, stoch, rate_kmol_m3s, exothermic_heat_kJ_kmol):
        pass

    def add_GasStream_viscosity_Pas(self, id, function):
        self.GasStream_viscosity_Pas[id] = function

    def add_GasStream_thermal_conductivity_kW_mK(self, id, function):
        self.GasStream_thermal_conductivity_kW_mK[id] = function

    def add_GasStream_heat_capacity_kJ_kmolK(self, id, function):
        self.GasStream_heat_capacity_kJ_kmolK[id] = function

    def add_GasStream_diffusivity(self, id, function):
        self.GasStream_diffusivity_m2_s[id] = function

    # ------------------------------------------------------

    def print_info(self):


        print("LiquidStream Species (Solids)")
        print("ID\tMolar Mass\tEnthalpy\tEntropy\tGibbs".expandtabs(20))
        for id in self.LiquidStream_specie_sol.keys():
            s = id
            for el in self.LiquidStream_specie_sol[id].keys():
                s = s + "\t" + str(self.LiquidStream_specie_sol[id][el])
            print(s.expandtabs(20))
        print("")
        print("")

        print("LiquidStream Species (Liquids)")
        print("ID\tMolar Mass\tCharge\tEnthalpy\tEntropy\tGibbs".expandtabs(20))
        for id in self.LiquidStream_specie_liq.keys():
            s = id
            for el in self.LiquidStream_specie_liq[id].keys():
                s = s + "\t" + str(self.LiquidStream_specie_liq[id][el])
            print(s.expandtabs(20))
        print("")
        print("")

        print("GasStream Species (Gas)")
        print("ID\tMolar Mass\tCharge\tEnthalpy\tEntropy\tGibbs".expandtabs(20))
        for id in self.GasStream_specie_gas.keys():
            s = id
            for el in self.GasStream_specie_gas[id].keys():
                s = s + "\t" + str(self.GasStream_specie_gas[id][el])
            print(s.expandtabs(20))
        print("")
        print("")

        print("LiquidStream Instantaneous Reactions")
        print("ID\tStochiometry".expandtabs(60))
        for id in self.LiquidStream_rxn_insta.keys():
            s = id + "\t" + str(self.LiquidStream_rxn_insta[id]["Stoch"])
            print(s.expandtabs(60))
        print("")
        print("")

        print("LiquidStream Reversible Reactions")
        print("ID\tStochiometry".expandtabs(60))
        for id in self.LiquidStream_rxn_reversible.keys():
            s = id + "\t" + str(self.LiquidStream_rxn_reversible[id]["Stoch"])
            print(s.expandtabs(60))
        print("")
        print("")

        print("LiquidStream Irreversible Reactions")
        print("ID\tStochiometry".expandtabs(60))
        for id in self.LiquidStream_rxn_irreversible.keys():
            s = id + "\t" + str(self.LiquidStream_rxn_irreversible[id]["Stoch"])
            print(s.expandtabs(60))
        print("")
        print("")

        print("GasStream Instantaneous Reactions")
        print("ID\tStochiometry".expandtabs(60))
        for id in self.GasStream_rxn_insta.keys():
            s = id + "\t" + str(self.GasStream_rxn_insta[id]["Stoch"])
            print(s.expandtabs(60))
        print("")
        print("")

        print("GasStream Reversible Reactions")
        print("ID\tStochiometry".expandtabs(60))
        for id in self.GasStream_rxn_reversible.keys():
            s = id + "\t" + str(self.GasStream_rxn_reversible[id]["Stoch"])
            print(s.expandtabs(60))
        print("")
        print("")

        print("GasStream Irreversible Reactions")
        print("ID\tStochiometry".expandtabs(60))
        for id in self.GasStream_rxn_irreversible.keys():
            s = id + "\t" + str(self.GasStream_rxn_irreversible[id]["Stoch"])
            print(s.expandtabs(60))
        print("")
        print("")

        print("LiquidStream Density")
        print("ID")
        for id in self.LiquidStream_density_kg_m3.keys():
            print(id)
        print("")
        print("")

        print("LiquidStream Viscosity")
        print("ID")
        for id in self.LiquidStream_viscosity_Pas.keys():
            print(id)
        print("")
        print("")

        print("LiquidStream Heat Capacity")
        print("ID")
        for id in self.LiquidStream_heat_capacity_kJ_kgK.keys():
            print(id)
        print("")
        print("")

        print("LiquidStream Thermal Conductivity")
        print("ID")
        for id in self.LiquidStream_thermal_conductivity_kW_mK.keys():
            print(id)
        print("")
        print("")

        print("LiquidStream Activity Coefficient")
        print("ID")
        for id in self.LiquidStream_activity_coefficient.keys():
            print(id)
        print("")
        print("")

        print("LiquidStream Diffusivity")
        print("ID")
        for id in self.LiquidStream_diffusivity_m2_s.keys():
            print(id)
        print("")
        print("")

        print("LiquidStream Surface Tension")
        print("ID")
        for id in self.LiquidStream_surface_tension_N_m.keys():
            print(id)
        print("")
        print("")

        print("GasStream Viscosity")
        print("ID")
        for id in self.GasStream_viscosity_Pas.keys():
            print(id)
        print("")
        print("")

        print("GasStream Thermal Conductivity")
        print("ID")
        for id in self.GasStream_thermal_conductivity_kW_mK.keys():
            print(id)
        print("")
        print("")

        print("GasStream Heat Capacity")
        print("ID")
        for id in self.GasStream_heat_capacity_kJ_kmolK.keys():
            print(id)
        print("")
        print("")

        print("GasStream Diffusivity")
        print("ID")
        for id in self.GasStream_diffusivity_m2_s.keys():
            print(id)
        print("")
        print("")

        print("LiquidStream Vapor Pressure")
        print("ID\tStochiometry (Liquid)\tStochiometry (Gas)".expandtabs(60))
        for id in self.LiquidStream_vapor_pressure_bara.keys():
            s = id + "\t" + str(self.LiquidStream_vapor_pressure_bara[id]["Stoch Liq"]) + "\t" + str(self.LiquidStream_vapor_pressure_bara[id]["Stoch Gas"])
            print(s.expandtabs(60))
        print("")
        print("")


class _Stochiometry:

    #     | R_liq_insta   R_liq_reversible  R_liq_vap     R_liq_mass_transfer   Z0             Z1
    # R = |
    #     | Z2            Z3                R_gas_vap     R_gas_mass_transfer   R_gas_insta    R_gas_reversible

    def __get_the_matrix__(self, GasStreamIn=None, LiquidStreamIn=None, Column=None, liq_rxn_insta=False, liq_rxn_reversible=False, vapor_pressure=False, mass_transfer=False, gas_rxn_insta=False, gas_rxn_reversible=False):

        if LiquidStreamIn is not None:
            num_of_species = len(LiquidStreamIn.specie.keys())
            R_liq_insta = np.zeros(shape=(num_of_species, 0), dtype=np.float64)
            R_liq_reversible = np.zeros(shape=(num_of_species, 0), dtype=np.float64)
            if liq_rxn_insta:
                R_liq_insta = self.__get_the_matrix_insta__(LiquidStreamIn)
            if liq_rxn_reversible:
                R_liq_reversible = self.__get_the_matrix_reversible__(LiquidStreamIn)

        if GasStreamIn is not None:
            num_of_species = len(GasStreamIn.specie.keys())
            R_gas_insta = np.zeros(shape=(num_of_species, 0), dtype=np.float64)
            R_gas_reversible = np.zeros(shape=(num_of_species, 0), dtype=np.float64)
            if gas_rxn_insta:
                R_gas_insta = self.__get_the_matrix_insta__(GasStreamIn)
            if gas_rxn_reversible:
                R_gas_reversible = self.__get_the_matrix_reversible__(GasStreamIn)

        if LiquidStreamIn is not None and GasStreamIn is not None:
            num_of_species_liq = len(LiquidStreamIn.specie.keys())
            num_of_species_gas = len(GasStreamIn.specie.keys())
            R_gas_vap = np.zeros(shape=(num_of_species_gas, 0), dtype=np.float64)
            R_gas_mass_transfer = np.zeros(shape=(num_of_species_gas, 0), dtype=np.float64)
            R_liq_vap = np.zeros(shape=(num_of_species_liq, 0), dtype=np.float64)
            R_liq_mass_transfer = np.zeros(shape=(num_of_species_liq, 0), dtype=np.float64)
            if vapor_pressure:
                R_gas_vap, R_liq_vap = self.__get_the_matrix_vapor__(GasStreamIn, LiquidStreamIn)
            if mass_transfer:
                R_gas_mass_transfer, R_liq_mass_transfer = self.__get_the_matrix_mass_transfer__(GasStreamIn, LiquidStreamIn, Column)

        # -------------------------------------------------

        if LiquidStreamIn is not None and GasStreamIn is None:
            R = np.concatenate((R_liq_insta, R_liq_reversible), axis=1)

        if LiquidStreamIn is None and GasStreamIn is not None:
            R = np.concatenate((R_gas_insta, R_gas_reversible), axis=1)

        if LiquidStreamIn is not None and GasStreamIn is not None:

            num_of_species_liq = len(LiquidStreamIn.specie.keys())
            num_of_species_gas = len(GasStreamIn.specie.keys())
            Z0 = np.zeros(shape=(num_of_species_liq, R_gas_insta.shape[1]), dtype=np.float64)
            Z1 = np.zeros(shape=(num_of_species_liq, R_gas_reversible.shape[1]), dtype=np.float64)
            Z2 = np.zeros(shape=(num_of_species_gas, R_liq_insta.shape[1]), dtype=np.float64)
            Z3 = np.zeros(shape=(num_of_species_gas, R_liq_reversible.shape[1]), dtype=np.float64)
            R1 = np.concatenate((R_liq_insta, R_liq_reversible, R_liq_vap, R_liq_mass_transfer, Z0, Z1), axis=1)
            R2 = np.concatenate((Z2, Z3, R_gas_vap, R_gas_mass_transfer, R_gas_insta, R_gas_reversible), axis=1)
            R = np.concatenate((R1, R2), axis=0)

        matrix = {}
        matrix["R"] = R
        matrix["R+"] = np.linalg.pinv(matrix["R"]).astype(np.float64)
        matrix["P"] = np.einsum("cq,dq->cd", np.einsum("cr,rq->cq", matrix["R"], np.linalg.pinv(np.einsum("cr,cq->rq", matrix["R"], matrix["R"]))), matrix["R"]).astype(np.float64)
        U, D, _ = np.linalg.svd(matrix["R"], full_matrices=True)
        matrix["A"] = U[:, R.shape[1]::].T.astype(np.float64)
        return matrix

    def __get_the_matrix_insta__(self, StreamIn):
        num_of_species = len(StreamIn.specie.keys())
        num_of_rxn_insta = len(StreamIn.rxn_insta.keys())
        num_of_rxn_reversible = len(StreamIn.rxn_reversible.keys())
        R = np.zeros(shape=(num_of_species, num_of_rxn_insta), dtype=np.float64)
        for j, jd in enumerate(StreamIn.rxn_insta.keys()):
            for id in StreamIn.specie.keys():
                i = StreamIn.specie[id]["Index"]
                if id in StreamIn.rxn_insta[jd]["Stoch"].keys():
                    nu = StreamIn.rxn_insta[jd]["Stoch"][id]
                    R[i, j] = nu * StreamIn.get_specie_molar_mass_kg_kmol(id=id)
        return R

    def __get_the_matrix_reversible__(self, StreamIn):
        num_of_species = len(StreamIn.specie.keys())
        num_of_rxn_insta = len(StreamIn.rxn_insta.keys())
        num_of_rxn_reversible = len(StreamIn.rxn_reversible.keys())
        R = np.zeros(shape=(num_of_species, num_of_rxn_reversible), dtype=np.float64)
        for j, jd in enumerate(StreamIn.rxn_reversible.keys()):
            for id in StreamIn.specie.keys():
                i = StreamIn.specie[id]["Index"]
                if id in StreamIn.rxn_reversible[jd]["Stoch"].keys():
                    nu = StreamIn.rxn_reversible[jd]["Stoch"][id]
                    R[i, j] = nu * StreamIn.get_specie_molar_mass_kg_kmol(id=id)
        return R

    def __get_the_matrix_vapor__(self, GasStreamIn, LiquidStreamIn):

        num_of_species_liq = len(LiquidStreamIn.specie.keys())
        num_of_species_gas = len(GasStreamIn.specie.keys())
        num_of_vapor = len(LiquidStreamIn.vapor_pressure_bara.keys())
        R_liq_vap = np.zeros(shape=(num_of_species_liq, num_of_vapor), dtype=np.float64)
        R_gas_vap = np.zeros(shape=(num_of_species_gas, num_of_vapor), dtype=np.float64)

        for j, jd in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            for id in LiquidStreamIn.specie.keys():
                i = LiquidStreamIn.specie[id]["Index"]
                if id in LiquidStreamIn.vapor_pressure_bara[jd]["Stoch Liq"].keys():
                    nu = LiquidStreamIn.vapor_pressure_bara[jd]["Stoch Liq"][id]
                    R_liq_vap[i, j] = nu * LiquidStreamIn.get_specie_molar_mass_kg_kmol(id=id)

        for j, jd in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            for id in GasStreamIn.specie.keys():
                i = GasStreamIn.specie[id]["Index"]
                if id in LiquidStreamIn.vapor_pressure_bara[jd]["Stoch Gas"].keys():
                    nu = LiquidStreamIn.vapor_pressure_bara[jd]["Stoch Gas"][id]
                    R_gas_vap[i, j] = nu * GasStreamIn.get_specie_molar_mass_kg_kmol(id=id)

        return R_gas_vap, R_liq_vap

    def __get_the_matrix_mass_transfer__(self, GasStreamIn, LiquidStreamIn, Column):

        num_of_species_liq = len(LiquidStreamIn.specie.keys())
        num_of_species_gas = len(GasStreamIn.specie.keys())
        num_of_mass_transfer = len(Column.mass_transfer_kmol_m3s.keys())
        R_liq_mass_transfer = np.zeros(shape=(num_of_species_liq, num_of_mass_transfer), dtype=np.float64)
        R_gas_mass_transfer = np.zeros(shape=(num_of_species_gas, num_of_mass_transfer), dtype=np.float64)

        for j, jd in enumerate(Column.mass_transfer_kmol_m3s.keys()):
            for id in LiquidStreamIn.specie.keys():
                i = LiquidStreamIn.specie[id]["Index"]
                if id in Column.mass_transfer_kmol_m3s[jd]["Stoch Liq"].keys():
                    nu = Column.mass_transfer_kmol_m3s[jd]["Stoch Liq"][id]
                    R_liq_mass_transfer[i, j] = nu * LiquidStreamIn.get_specie_molar_mass_kg_kmol(id=id)

        for j, jd in enumerate(Column.mass_transfer_kmol_m3s.keys()):
            for id in GasStreamIn.specie.keys():
                i = GasStreamIn.specie[id]["Index"]
                if id in Column.mass_transfer_kmol_m3s[jd]["Stoch Gas"].keys():
                    nu = Column.mass_transfer_kmol_m3s[jd]["Stoch Gas"][id]
                    R_gas_mass_transfer[i, j] = nu * GasStreamIn.get_specie_molar_mass_kg_kmol(id=id)

        return R_gas_mass_transfer, R_liq_mass_transfer


class _Serializer:

    def save(self, instance, filename):
        file = open(filename + ".obj","wb")
        pickle.dump(instance, file)
        file.close()

    def load(self, instance, filename):
        file = open(filename + ".obj", 'rb')
        instance = pickle.load(file)
        file.close()
        return instance


class _LiquidEquilibrium:

    def __get_LiquidEquilibrium__(self, LiquidStreamIn, b, matrix, lr):

        w = LiquidStreamIn.__mass_fractions_dic2vec__()
        b = np.einsum("cw,sw->sc", matrix["A"], w) if b is None else b

        LiquidStreamOut = deepcopy(LiquidStreamIn)
        converged = False
        epoch = 0
        iterations = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0],))

        while converged == False:

            Kw = self.__get_LiquidEquilibrium_Kw__(LiquidStreamOut)

            f_rxn_insta = self.__get_LiquidEquilibrium_f_rxn_insta__(LiquidStreamOut, Kw)
            f_mass_balance = self.__get_LiquidEquilibrium_f_mass_balance__(LiquidStreamOut, b, matrix)
            f = np.concatenate((f_rxn_insta, f_mass_balance), axis=1)

            dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta__(LiquidStreamOut, Kw, matrix)
            dfdw_mass_balance = self.__get_LiquidEquilibrium_dfdw_mass_balance__(LiquidStreamOut, matrix)
            dfdw = np.concatenate((dfdw_rxn_insta, dfdw_mass_balance), axis=1)

            # dfdw = self.__get_liquid_equilibrium_dfdw_rxn_insta__(LiquidStreamIn, Kw, matrix)
            #         dfdr = np.einsum("sfw,wr->sfr", dfdw, matrix["R"])

            dw_newton = - np.linalg.solve(dfdw, f)
            dw = lr * dw_newton

            w = w + dw
            w = np.maximum(w, 0)
            LiquidStreamOut.__mass_fractions_vec2dic__(w)

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        return LiquidStreamOut

    def __get_LiquidEquilibrium_reduced__(self, LiquidStreamIn, matrix, lr):

        w0 = LiquidStreamIn.__mass_fractions_dic2vec__()
        w = 1.0 * w0

        LiquidStreamOut = deepcopy(LiquidStreamIn)
        converged = False
        epoch = 0
        iterations = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0],))

        r = 0

        while converged == False:
            Kw = self.__get_LiquidEquilibrium_Kw__(LiquidStreamOut)

            f = self.__get_LiquidEquilibrium_f_rxn_insta__(LiquidStreamIn, Kw)

            dfdw = self.__get_LiquidEquilibrium_dfdw_rxn_insta__(LiquidStreamIn, Kw, matrix)
            dfdr = np.einsum("sfw,wr->sfr", dfdw, matrix["R"])

            dr_newton = - np.linalg.solve(dfdr, f)
            dw_newton = np.einsum("wr,sr->sw", matrix["R"], dr_newton)

            dr = lr * dr_newton
            dw = lr * dw_newton

            # Backtrack to Ensure only Positive Concentrations
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)

            dr = tau * dr
            dw = tau * dw

            r = r + dr

            w = w0 + np.einsum("wr,sr->sw", matrix["R"], dr)
            LiquidStreamOut.__mass_fractions_vec2dic__(w)

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        return LiquidStreamOut

    # ----------------------------------------------------------------------------------

    def __get_LiquidEquilibrium_Ka__(self, LiquidStreamIn):
        Ka = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
        for j, jd in enumerate(LiquidStreamIn.rxn_insta.keys()):
            Ka[:, j] = LiquidStreamIn.get_rxn_insta_equilibrium_constant_wrt_activities(jd)
        return Ka

    def __get_LiquidEquilibrium_Kw__(self, LiquidStreamIn):
        Kw = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
        for j, jd in enumerate(LiquidStreamIn.rxn_insta.keys()):
            Kw[:, j] = LiquidStreamIn.get_rxn_insta_equilibrium_constant_wrt_mass_fractions(jd)
        return Kw

    def __get_LiquidEquilibrium_dKadT__(self, LiquidStreamIn, Ka):
        dT = 0.05
        T0 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + dT
        Ka1 = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.rxn_insta.keys()):
            Ka1[:, i] = LiquidStreamIn.get_rxn_insta_equilibrium_constant_wrt_activities(id)
        T1 = LiquidStreamIn.temp_K
        LiquidStreamIntemp_K = LiquidStreamIn.temp_K - dT
        return (Ka1 - Ka) / dT

    def __get_LiquidEquilibrium_dKwdT__(self, LiquidStreamIn, Kw):
        dT = 0.05
        T0 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + dT
        Kw1 = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.rxn_insta.keys()):
            Kw1[:, i] = LiquidStreamIn.get_rxn_insta_equilibrium_constant_wrt_mass_fractions(id)
        T1 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - dT
        return (Kw1 - Kw) / dT

    def __get_LiquidEquilibrium_exothermic_heat_kJ_kmol_approximation__(self, LiquidStreamIn, Ka, dKadT):
        dT = 0.05
        T = LiquidStreamIn.temp_K
        T1 = T + dT
        Ka1 = Ka + dKadT * dT
        q = 8.314 * (np.log(Ka1) - np.log(Ka)) / ((1 / T1) - (1 / T))
        return q

    def __get_LiquidEquilibrium_exothermic_heat_kJ_kmol__(self, LiquidStreamIn, Kw, dKwdT):
        dT = 0.05
        T = LiquidStreamIn.temp_K[:,None]
        T1 = T + dT
        Kw1 = Kw + dKwdT * dT
        q = 8.314 * (np.log(Kw1) - np.log(Kw)) / ((1 / T1) - (1 / T))
        return q

    # ----------------------------------------------------------------------------------

    def __get_LiquidEquilibrium_f_rxn_insta__(self, LiquidStreamIn, Kw):
        f_rxn_insta = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
        for j, jd in enumerate(LiquidStreamIn.rxn_insta.keys()):
            fr = Kw[:, j]
            fp = np.ones(shape=(LiquidStreamIn.temp_K.shape[0],), dtype=np.float64)
            for id in LiquidStreamIn.rxn_insta[jd]["Stoch"].keys():
                nu = LiquidStreamIn.rxn_insta[jd]["Stoch"][id]
                if nu > 0:
                    fp = fp * LiquidStreamIn.get_specie_mass_fraction(id=id) ** np.abs(nu)
                else:
                    fr = fr * LiquidStreamIn.get_specie_mass_fraction(id=id) ** np.abs(nu)
            f_rxn_insta[:, j] = fr - fp
        return f_rxn_insta

    def __get_LiquidEquilibrium_f_rxn_insta_log__(self, LiquidStreamIn, Kw):
        f_rxn_insta = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
        for j, jd in enumerate(LiquidStreamIn.rxn_insta.keys()):
            f_rxn_insta[:, j] = np.log(Kw[:, j])
            for id in LiquidStreamIn.rxn_insta[jd]["Stoch"].keys():
                nu = LiquidStreamIn.rxn_insta[jd]["Stoch"][id]
                f_rxn_insta[:, j] = f_rxn_insta[:, j] - nu * np.log(LiquidStreamIn.get_specie_mass_fraction(id=id))
        return f_rxn_insta

    def __get_LiquidEquilibrium_f_mass_balance__(self, LiquidStreamIn, b, matrix):
        f_mass_balance = np.einsum("cw,sw->sc", matrix["A"], LiquidStreamIn.__mass_fractions_dic2vec__()) - b
        return f_mass_balance

    def __get_LiquidEquilibrium_f_energy_balance__(self, LiquidStreamInit, LiquidStreamIn, rxn_insta_exothermic_heat_kJ_kmol, matrix):
        cp = (LiquidStreamInit.get_solution_heat_capacity_kJ_kgK() + LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()) / 2
        dw = LiquidStreamIn.__mass_fractions_dic2vec__() - LiquidStreamInit.__mass_fractions_dic2vec__()
        dT = LiquidStreamIn.temp_K - LiquidStreamInit.temp_K
        f_energy_balance = (1 / cp) * np.einsum("sr,sr->s", rxn_insta_exothermic_heat_kJ_kmol, np.einsum("rw,sw->sr", matrix["R+"], dw)) - dT
        f_energy_balance = f_energy_balance[:, None]
        return f_energy_balance

    # ----------------------------------------------------------------------------------

    def __get_LiquidEquilibrium_dfdw_rxn_insta__(self, LiquidStreamIn, Kw, matrix):
        dfdw_rxn_insta = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta, LiquidStreamIn.num_of_species),dtype=np.float64)
        for rxn_insta_i, rxn_insta_id in enumerate(LiquidStreamIn.rxn_insta.keys()):
            for specie_i, specie_id in enumerate(LiquidStreamIn.specie.keys()):
                if specie_id in list(LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys()):
                    nu = LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_id]
                    if nu < 0:
                        dfdw_rxn_insta[:, rxn_insta_i, specie_i] = Kw[:, rxn_insta_i] * np.abs(
                            nu) * LiquidStreamIn.get_specie_mass_fraction(id=specie_id) ** (np.abs(nu) - 1)
                        for specie_jd in LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys():
                            nu = LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_jd]
                            if specie_id != specie_jd and nu < 0:
                                dfdw_rxn_insta[:, rxn_insta_i, specie_i] = dfdw_rxn_insta[:, rxn_insta_i, specie_i] * LiquidStreamIn.get_specie_mass_fraction(id=specie_jd) ** np.abs(nu)
                    else:
                        dfdw_rxn_insta[:, rxn_insta_i, specie_i] = - np.abs(nu) * LiquidStreamIn.get_specie_mass_fraction(id=specie_id) ** (np.abs(nu) - 1)
                        for specie_jd in LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys():
                            nu = LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_jd]
                            if specie_id != specie_jd and nu > 0:
                                dfdw_rxn_insta[:, rxn_insta_i, specie_i] = dfdw_rxn_insta[:, rxn_insta_i, specie_i] * LiquidStreamIn.get_specie_mass_fraction(id=specie_jd) ** np.abs(nu)
        return dfdw_rxn_insta

    def __get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self, LiquidStreamIn):
        dfdw_rxn_insta = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta, LiquidStreamIn.num_of_species),dtype=np.float64)
        for rxn_insta_i, rxn_insta_id in enumerate(LiquidStreamIn.rxn_insta.keys()):
            for _, specie_id in enumerate(LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys()):
                specie_i = LiquidStreamIn.specie[specie_id]["Index"]
                nu = LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_id]
                dfdw_rxn_insta[:, rxn_insta_i, specie_i] = - nu / LiquidStreamIn.get_specie_mass_fraction(id=specie_id)
        return dfdw_rxn_insta

    def __get_LiquidEquilibrium_dfdw_mass_balance__(self, LiquidStreamIn, matrix):
        dfdw_mass_balance = np.broadcast_to(array=matrix["A"][None, :, :], shape=(LiquidStreamIn.temp_K.shape[0], matrix["A"].shape[0], matrix["A"].shape[1])).copy()
        return dfdw_mass_balance

    def __get_LiquidEquilibrium_dfdw_energy_balance__(self, LiquidStreamInit, LiquidStreamIn, rxn_insta_exothermic_heat_kJ_kmol, matrix):
        cp = (LiquidStreamInit.get_solution_heat_capacity_kJ_kgK() + LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()) / 2
        dfdw_energy_balance = np.einsum("s,sw->sw", 1 / cp, np.einsum("sr,rw->sw", rxn_insta_exothermic_heat_kJ_kmol, matrix["R+"]))[:,None, :]
        return dfdw_energy_balance

    # ----------------------------------------------------------------------------------

    def __get_LiquidEquilibrium_dfdT_rxn_insta__(self, LiquidStreamIn, dKwdT):
        dfdT_rxn_insta = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta, 1), dtype=np.float64)
        for j, jd in enumerate(LiquidStreamIn.rxn_insta.keys()):
            fr = dKwdT[:, j]
            for id in LiquidStreamIn.rxn_insta[jd]["Stoch"].keys():
                nu = LiquidStreamIn.rxn_insta[jd]["Stoch"][id]
                if nu < 0:
                    fr = fr * LiquidStreamIn.get_specie_mass_fraction(id=id) ** np.abs(nu)
            dfdT_rxn_insta[:, j, 0] = fr
        return dfdT_rxn_insta

    def __get_LiquidEquilibrium_dfdT_rxn_insta_log__(self, Kw, dKwdT):
        dfdT = dKwdT / Kw
        return dfdT[:,:,None]

    def __get_LiquidEquilibrium_dfdT_mass_balance__(self, LiquidStreamIn, matrix):
        dfdT_mass_balance = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], matrix["A"].shape[0], 1), dtype=np.float64)

        return dfdT_mass_balance

    def __get_LiquidEquilibrium_dfdT_energy_balance__(self, LiquidStreamIn):
        return - np.ones(shape=(LiquidStreamIn.temp_K.shape[0], 1, 1), dtype=np.float64)

    # ----------------------------------------------------------------------------------

    def __get_LiquidEquilibrium_sensitivities__(self, LiquidStreamIn, Kw, matrix):

        # Partial Derivatives
        dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(LiquidStreamIn)
        dfdw_mass_balance = self.__get_LiquidEquilibrium_dfdw_mass_balance__(LiquidStreamIn, matrix)
        dfdw = np.concatenate((dfdw_rxn_insta, dfdw_mass_balance), axis=1)

        # Hessian
        H = np.linalg.inv(dfdw)

        # Sensitivity Matrix (dw/db) with respect to "b"
        dwdb = H[:, :, LiquidStreamIn.num_of_rxn_insta::]

        # Sensitivity Matrix (dw/dT) with respect to "T"
        dKwdT = self.__get_LiquidEquilibrium_dKwdT__(LiquidStreamIn, Kw)
        Hr = H[:, :, :LiquidStreamIn.num_of_rxn_insta:]
        dwdT = - np.einsum("scr,sr->sc", Hr, (1 / Kw) * dKwdT)

        return dwdb, dwdT

    def __get_LiquidEquilibrium_heat_dissipation_kW__(self, LiquidStreamIn, KwIn, dKwdTIn, LiquidStreamOut, KwOut, dKwdTOut, matrix):
        h_in = self.__get_liquid_equilibrium_exothermic_heat_kJ_kmol__(LiquidStreamIn, KwIn, dKwdTin)
        h_out = self.__get_liquid_equilibrium_exothermic_heat_kJ_kmol__(LiquidStreamOut, KwOut, dKwdTOut)
        h = (h_in + h_out) / 2
        cp_in = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        cp_out = LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        cp = (cp_in + cp_out) / 2
        m = LiquidStreamIn.get_solution_flow_kg_h() / 3600
        dw = LiquidStreamOut.__mass_fractions_dic2vec__() - LiquidStreamIn.__mass_fractions_dic2vec__()
        dT = LiquidStreamOut.temp_K - LiquidStreamIn.temp_K
        Q = m * np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", matrix_["R+"], dw)) - m * cp *dT
        return Q


class _GasEquilibrium:

    def __get_GasEquilibrium_Kp__(self, GasStreamIn):
        Kp = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta), dtype=np.float64)
        for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
            Kp[:, j] = GasStreamIn.get_rxn_insta_equilibrium_constant_wrt_partial_pressures(jd)
        return Kp

    def __get_GasEquilibrium_Ky__(self, GasStreamIn):
        Ky = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta), dtype=np.float64)
        for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
            Ky[:, j] = GasStreamIn.get_rxn_insta_equilibrium_constant_wrt_molar_fractions(jd)
        return Ky

    def __get_GasEquilibrium_dKpdT__(self, GasStreamIn, Kp):
        dT = 0.05
        T0 = GasStreamIn.temp_K
        GasStreamIn.temp_K = GasStreamIn.temp_K + dT
        Kp1 = self.__get_GasEquilibrium_Kp__(GasStreamIn)
        T1 = GasStreamIn.temp_K
        GasStreamIn.temp_K = GasStreamIn.temp_K - dT
        return (Kp1 - Kp) / dT

    def __get_GasEquilibrium_dKydT__(self, GasStreamIn, Ky):
        dT = 0.05
        T0 = GasStreamIn.temp_K
        GasStreamIn.temp_K = GasStreamIn.temp_K + dT
        Ky1 = self.__get_GasEquilibrium_Ky__(GasStreamIn)
        T1 = GasStreamIn.temp_K
        GasStreamIn.temp_K = GasStreamIn.temp_K - dT
        return (Ky1 - Ky) / dT

    def __get_GasEquilibrium_exothermic_heat_kJ_kmol__(self, GasStreamIn, Ky, dKydT):
        dT = 0.05
        T = GasStreamIn.temp_K[:,None]
        T1 = T + dT
        Ky1 = Ky + dKydT * dT
        q = 8.314 * (np.log(Ky1) - np.log(Ky)) / ((1 / T1) - (1 / T))
        return q

    # ----------------------------------------------------------------------------------------------------

    def __get_GasEquilibrium_f_rxn_insta__(self, GasStreamIn, Kp):
        f_rxn_insta = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta), dtype=np.float64)
        for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
            fr = Kp[:, j]
            fp = np.ones(shape=(GasStreamIn.temp_K.shape[0],), dtype=np.float64)
            for id in GasStreamIn.rxn_insta[jd]["Stoch"].keys():
                nu = GasStreamIn.rxn_insta[jd]["Stoch"][id]
                if nu > 0:
                    fp = fp * GasStreamIn.get_specie_pressure_bara(id=id) ** np.abs(nu)
                else:
                    fr = fr * GasStreamIn.get_specie_pressure_bara(id=id) ** np.abs(nu)
            f_rxn_insta[:, j] = fr - fp
        return f_rxn_insta

    def __get_GasEquilibrium_f_rxn_insta_log__(self, GasStreamIn, Kp):
        f_rxn_insta = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta), dtype=np.float64)
        for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
            f_rxn_insta[:, j] = np.log(Kp[:, j])
            for id in GasStreamIn.rxn_insta[jd]["Stoch"].keys():
                nu = GasStreamIn.rxn_insta[jd]["Stoch"][id]
                f_rxn_insta[:, j] = f_rxn_insta[:, j] - nu * np.log(GasStreamIn.get_specie_pressure_bara(id=id))
        return f_rxn_insta

    # Ongoing
    def __get_GasEquilibrium_f_mass_balance__(self, GasStreamIn, b, matrix):
        f_mass_balance = np.einsum("cw,sw->sc", matrix["A"], GasStreamIn.__mass_fractions_dic2vec__()) - b
        return f_mass_balance

    # Ongoing
    def __get_GasEquilibrium_f_energy_balance__(self, GasStreamInit, GasStreamIn, rxn_insta_exothermic_heat_kJ_kmol, matrix):
        cp = (LiquidStreamInit.get_solution_heat_capacity_kJ_kgK() + LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()) / 2
        dw = LiquidStreamIn.__mass_fractions_dic2vec__() - LiquidStreamInit.__mass_fractions_dic2vec__()
        dT = LiquidStreamIn.temp_K - LiquidStreamInit.temp_K
        f_energy_balance = (1 / cp) * np.einsum("sr,sr->s", rxn_insta_exothermic_heat_kJ_kmol, np.einsum("rw,sw->sr", matrix["R+"], dw)) - dT
        f_energy_balance = f_energy_balance[:, None]
        return f_energy_balance

    # ----------------------------------------------------------------------------------------------------

    def __get_GasEquilibrium_dfdy_rxn_insta_log__(self, GasStreamIn):
        dfdy_rxn_insta = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta, GasStreamIn.num_of_species), dtype=np.float64)
        for rxn_insta_i, rxn_insta_id in enumerate(GasStreamIn.rxn_insta.keys()):
            for _, specie_id in enumerate(GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys()):
                specie_i = GasStreamIn.specie[specie_id]["Index"]
                nu = GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_id]
                dfdy_rxn_insta[:, rxn_insta_i, specie_i] = - nu / GasStreamIn.get_specie_molar_fraction(id=specie_id)
        return dfdy_rxn_insta

    def __get_GasEquilibrium_dfdT_rxn_insta_log__(self, Ky, dKydT):
        dfdT = dKydT / Ky
        return dfdT[:,:,None]


class _VaporLiquidEquilibrium:

    def __get_VaporLiquidEquilibrium_Kpw__(self, LiquidStreamIn):
        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]

        c = LiquidStreamIn.get_solution_molarity_kmol_m3()
        rho = LiquidStreamIn.get_solution_density_kg_m3()
        den = 0
        for specie_id in LiquidStreamIn.specie.keys():
            den = den + LiquidStreamIn.get_specie_mass_fraction(
                specie_id) / LiquidStreamIn.get_specie_molar_mass_kg_kmol(specie_id)

        Kyw = np.zeros(shape=(num_of_samples, num_of_vap), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            unit = LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][liq_id]
            gamma = LiquidStreamIn.get_specie_activity_coefficient(id=liq_id)
            if LiquidStreamIn.vapor_pressure_bara[id]["H"] == None:
                Kyw[:, i] = 1 / LiquidStreamIn.vapor_pressure_bara[id]["p0"](LiquidStreamIn)
            else:
                Kyw[:, i] = LiquidStreamIn.vapor_pressure_bara[id]["H"](LiquidStreamIn)
            Kyw[:, i] = Kyw[:, i] / gamma

            for r in LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"].keys():
                power = LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"][r]
                if LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][r] == "c":
                    Kyw[:, i] = Kyw[:, i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r) / c) ** power
                elif LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][r] == "x":
                    Kyw[:, i] = Kyw[:, i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r)) ** power
                elif LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][r] == "m":
                    Kyw[:, i] = Kyw[:, i] * ((rho * LiquidStreamIn.get_specie_mass_fraction(
                        id=LiquidStreamIn.solvent_id) * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r) * den) / (1000 * c)) ** power
                elif LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][r] == None:
                    Kyw[:, i] = Kyw[:, i] * (1 / LiquidStreamIn.get_specie_mass_fraction(id=r)) ** power
                elif LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][r] == "w":
                    pass
        return Kyw

    def __get_VaporLiquidEquilibrium_Kyw__(self, GasStreamIn, LiquidStreamIn):

        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]

        c = LiquidStreamIn.get_solution_molarity_kmol_m3()
        rho = LiquidStreamIn.get_solution_density_kg_m3()
        den = 0
        for specie_id in LiquidStreamIn.specie.keys():
            den = den + LiquidStreamIn.get_specie_mass_fraction(
                specie_id) / LiquidStreamIn.get_specie_molar_mass_kg_kmol(specie_id)

        Kyw = np.zeros(shape=(num_of_samples, num_of_vap), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            unit = LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][liq_id]
            gamma = LiquidStreamIn.get_specie_activity_coefficient(id=liq_id)
            if LiquidStreamIn.vapor_pressure_bara[id]["H"] == None:
                Kyw[:,i] = 1 / LiquidStreamIn.vapor_pressure_bara[id]["p0"](LiquidStreamIn)
            else:
                Kyw[:,i] = LiquidStreamIn.vapor_pressure_bara[id]["H"](LiquidStreamIn)
            Kyw[:, i] = Kyw[:, i] / gamma
            Kyw[:, i] = Kyw[:, i] * GasStreamIn.get_gas_pressure_bara()

            for r in LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"].keys():
                power = LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"][r]
                if LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][r] == "c":
                    Kyw[:,i] = Kyw[:,i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r) / c) ** power
                elif LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][r] == "x":
                    Kyw[:,i] = Kyw[:,i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r)) ** power
                elif LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][r] == "m":
                    Kyw[:,i] = Kyw[:,i] * ((rho * LiquidStreamIn.get_specie_mass_fraction(id=LiquidStreamIn.solvent_id) * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r) * den) / (1000 * c)) ** power
                elif LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][r] == None:
                    Kyw[:,i] = Kyw[:,i] * (1 / LiquidStreamIn.get_specie_mass_fraction(id=r)) ** power
                elif LiquidStreamIn.vapor_pressure_bara[id]["Unit Liq"][r] == "w":
                    pass
        return Kyw

    def __get_VaporLiquidEquilibrium_dKywdT__(self, GasStreamIn, LiquidStreamIn, Kyw):
        dT = 0.05
        T0 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + dT
        Kyw1 = self.__get_VaporLiquidEquilibrium_Kyw__(GasStreamIn, LiquidStreamIn)
        T1 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - dT
        return (Kyw1 - Kyw) / dT

    def __get_VaporLiquidEquilibrium_exothermic_heat_kJ_kmol__(self, LiquidStreamIn, Kyw, dKywdT):
        dT = 0.05
        T = LiquidStreamIn.temp_K[:,None]
        T1 = T + dT
        Kyw1 = Kyw + dKywdT * dT
        q = 8.314 * (np.log(Kyw1) - np.log(Kyw)) / ((1 / T1) - (1 / T))
        return q

    def __get_VaporLiquidEquilibrium_f_log__(self, GasStreamIn, LiquidStreamIn, Kyw):
        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        f_vap = np.log(Kyw)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            f_vap[:, i] = f_vap[:, i] + np.log(GasStreamIn.get_specie_molar_fraction(id=gas_id)) - np.log(LiquidStreamIn.get_specie_mass_fraction(id=liq_id))
        return f_vap

    def __get_VaporLiquidEquilibrium_dfdw_log__(self, LiquidStreamIn):
        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_species = LiquidStreamIn.num_of_species
        dfdw = np.zeros(shape=(num_of_samples, num_of_vap, num_of_species), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            liq_i = LiquidStreamIn.specie[liq_id]["Index"]
            dfdw[:, i, liq_i] = - 1 / LiquidStreamIn.get_specie_mass_fraction(id=liq_id)
        return dfdw

    def __get_VaporLiquidEquilibrium_dfdy_log__(self, GasStreamIn, LiquidStreamIn):
        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_species = GasStreamIn.num_of_species
        dfdy = np.zeros(shape=(num_of_samples, num_of_vap, num_of_species), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            gas_i = GasStreamIn.specie[gas_id]["Index"]
            dfdy[:, i, gas_i] = 1 / GasStreamIn.get_specie_molar_fraction(id=gas_id)
        return dfdy

    def __get_VaporLiquidEquilibrium_dfdT_log__(self, Kyw, dKywdT):
        dfdT = dKywdT / Kyw
        return dfdT[:,:,None]


class _LiquidCSTR:

    def __get_LiquidCSTR_rxn_reversible_rates_kmol_m3s__(self, LiquidStreamIn):
        r = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.rxn_reversible.keys()):
            r_forward, r_backward = LiquidStreamIn.rxn_reversible[id]["Rate [kmol/m3.s]"](LiquidStreamIn)
            r[:, i] = r_forward - r_backward
        return r

    # -------------------------------------------------------------------------------------------

    def __get_LiquidCSTR_f_rxn_reversible__(self, LiquidStreamInit, LiquidStreamIn, matrix, rates_kmol_m3s, phi):
        RI = matrix["R+"][LiquidStreamIn.num_of_rxn_insta::, :]
        dw = LiquidStreamIn.__mass_fractions_dic2vec__() - LiquidStreamInit.__mass_fractions_dic2vec__()
        f_rxn_reversible = rates_kmol_m3s - phi[:,None] * np.einsum("rw,sw->sr", RI, dw)
        return f_rxn_reversible

    def __get_LiquidCSTR_f_energy_balance__(self, LiquidStreamInit, LiquidStreamIn, matrix, heat_kW):

        h_init_insta = LiquidStreamInit.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_init_reversible = LiquidStreamInit.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()

        h_in_insta = LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_in_reversible = LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()

        h_init = np.concatenate((h_init_insta, h_init_reversible), axis=1)
        h_in = np.concatenate((h_in_insta, h_in_reversible), axis=1)
        h = (h_init + h_in) / 2

        cp_init = LiquidStreamInit.get_solution_heat_capacity_kJ_kgK()
        cp_in = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        cp = (cp_init + cp_in) / 2

        dw = LiquidStreamIn.__mass_fractions_dic2vec__() - LiquidStreamInit.__mass_fractions_dic2vec__()
        dT = LiquidStreamIn.temp_K - LiquidStreamInit.temp_K

        m = self.LiquidStreamInit.get_solution_flow_kg_h()/3600

        f_energy_balance = np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw)) / cp - dT + heat_kW/(cp * m)
        return f_energy_balance[:,None]

    # -------------------------------------------------------------------------------------------

    def __get_LiquidCSTR_dfdw_rxn_reversible__(self, LiquidStreamIn, matrix, rates_kmol_m3s, phi):
        dfdw_rxn_reversible = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible, LiquidStreamIn.num_of_species), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(LiquidStreamIn.rxn_reversible.keys()):
            for specie_id in LiquidStreamIn.rxn_reversible[rxn_id]["Dependencies"]:
                specie_i = LiquidStreamIn.specie[specie_id]["Index"]
                w_probe = LiquidStreamIn.get_specie_mass_fraction(id=specie_id)
                dw_probe = np.maximum(0.01 * w_probe, 10 ** (-14))
                LiquidStreamIn.set_specie_mass_fraction(id=specie_id, value=w_probe + dw_probe)
                forward_rate_pertubation, backward_rate_pertubation = LiquidStreamIn.rxn_reversible[rxn_id]["Rate [kmol/m3.s]"](LiquidStreamIn)
                rate_pertubation = forward_rate_pertubation - backward_rate_pertubation
                LiquidStreamIn.set_specie_mass_fraction(id=specie_id, value=w_probe)
                dfdw_rxn_reversible[:, rxn_i, specie_i] = (rate_pertubation - rates_kmol_m3s[:, rxn_i]) / dw_probe
        RI = matrix["R+"][LiquidStreamIn.num_of_rxn_insta::, :]
        dfdw_rxn_reversible = dfdw_rxn_reversible - np.einsum("s,rw->srw", phi, RI)
        return dfdw_rxn_reversible

    def __get_LiquidCSTR_dfdT_rxn_reversible__(self, LiquidStreamIn, matrix, rates_kmol_m3s):
        dfdT_rxn_reversible = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible, 1), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(LiquidStreamIn.rxn_reversible.keys()):
            dT = 0.05
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + dT
            forward_rate_pertubation, backward_rate_pertubation = LiquidStreamIn.rxn_reversible[rxn_id]["Rate [kmol/m3.s]"](LiquidStreamIn)
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - dT
            rate_pertubation = forward_rate_pertubation - backward_rate_pertubation
            dfdT_rxn_reversible[:, rxn_i, 0] = (rate_pertubation - rates_kmol_m3s[:, rxn_i]) / dT
        return dfdT_rxn_reversible

    def __get_LiquidCSTR_dfdw_energy_balance__(self, LiquidStreamInit, LiquidStreamIn, matrix):

        h_init_insta = LiquidStreamInit.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_init_reversible = LiquidStreamInit.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()

        h_in_insta = LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_in_reversible = LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()

        h_init = np.concatenate((h_init_insta, h_init_reversible), axis=1)
        h_in = np.concatenate((h_in_insta, h_in_reversible), axis=1)
        h = (h_init + h_in) / 2

        cp_init = LiquidStreamInit.get_solution_heat_capacity_kJ_kgK()
        cp_in = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        cp = (cp_init + cp_in) / 2

        dfdw_energy_balance = np.einsum("sr,rw->sw", h, self.matrix["R+"]) / cp[:,None]
        return dfdw_energy_balance[:,None,:]

    def __get_LiquidCSTR_dfdT_energy_balance__(self, LiquidStreamIn):
        dfdT_energy_balance = - np.ones(shape=(LiquidStreamIn.temp_K.shape[0], 1, 1), dtype=np.float64)
        return dfdT_energy_balance


class _StreamFunctions:

    # ----------------------------------------------------------------------------

    def __broadcast_to_LiquidProfile__(self, LiquidStreamIn, num_of_heights, LiquidProfile=None):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        if LiquidProfile is None:
            LiquidProfile = deepcopy(LiquidStreamIn)
        LiquidProfile.temp_K = np.broadcast_to(array=LiquidStreamIn.temp_K[None, :], shape=(num_of_heights, num_of_samples)).copy()
        LiquidProfile.flow_kg_h = np.broadcast_to(array=LiquidStreamIn.flow_kg_h[None, :], shape=(num_of_heights, num_of_samples)).copy()
        for id in LiquidProfile.specie.keys():
            LiquidProfile.specie[id]["Mass Fraction"] = np.broadcast_to(array=LiquidStreamIn.specie[id]["Mass Fraction"][None, :], shape=(num_of_heights, num_of_samples)).copy()
        for id in LiquidProfile.info.keys():
            LiquidProfile.info[id] = np.broadcast_to(array=LiquidStreamIn.info[id][None, :] ,shape=(num_of_heights, num_of_samples)).copy()
        return LiquidProfile

    def __broadcast_to_GasProfile__(self, GasStreamIn, num_of_heights, GasProfile=None):
        num_of_samples = GasStreamIn.temp_K.shape[0]
        if GasProfile is None:
            GasProfile = deepcopy(GasStreamIn)
        GasProfile.temp_K = np.broadcast_to(array=GasStreamIn.temp_K[None, :], shape=(num_of_heights, num_of_samples)).copy()
        GasProfile.flow_kmol_h = np.broadcast_to(array=GasStreamIn.flow_kmol_h[None, :], shape=(num_of_heights, num_of_samples)).copy()
        GasProfile.pressure_bara = np.broadcast_to(array=GasStreamIn.pressure_bara[None, :], shape=(num_of_heights, num_of_samples)).copy()

        for id in GasProfile.specie.keys():
            GasProfile.specie[id]["Molar Fraction"] = np.broadcast_to(array=GasStreamIn.specie[id]["Molar Fraction"][None, :], shape=(num_of_heights, num_of_samples)).copy()
        for id in GasProfile.info.keys():
            GasProfile.info[id] = np.broadcast_to(array=GasStreamIn.info[id][None, :], shape=(num_of_heights, num_of_samples)).copy()
        return GasProfile

    # ----------------------------------------------------------------------------

    def __get_slice_from_LiquidProfile__(self, LiquidStreamProfile, height_index, LiquidStreamOut=None):
        if LiquidStreamOut is None:
            LiquidStreamOut = deepcopy(LiquidStreamProfile)
        LiquidStreamOut.temp_K = LiquidStreamProfile.temp_K[height_index, :]
        LiquidStreamOut.flow_kg_h = LiquidStreamProfile.flow_kg_h[height_index, :]
        for id in LiquidStreamOut.specie.keys():
            LiquidStreamOut.specie[id]["Mass Fraction"] = LiquidStreamProfile.specie[id]["Mass Fraction"][height_index, :]
        for id in LiquidStreamOut.info.keys():
            LiquidStreamOut.info[id] = LiquidStreamProfile.info[id][height_index, :]
        return LiquidStreamOut

    def __get_slice_from_GasProfile__(self, GasStreamProfile, height_index, GasStreamOut=None):
        if GasStreamOut is None:
            GasStreamOut = deepcopy(GasStreamProfile)
        GasStreamOut.temp_K = GasStreamProfile.temp_K[height_index, :]
        GasStreamOut.flow_kmol_h = GasStreamProfile.flow_kmol_h[height_index, :]
        for id in GasStreamOut.specie.keys():
            GasStreamOut.specie[id]["Molar Fraction"] = GasStreamProfile.specie[id]["Molar Fraction"][height_index, :]
        for id in GasStreamOut.info.keys():
            GasStreamOut.info[id] = GasStreamProfile.info[id][height_index, :]
        return GasStreamOut

    # ----------------------------------------------------------------------------

    def __insert_slice_to_LiquidProfile__(self, LiquidProfile, height_index, LiquidStreamIn):
        LiquidProfile.flow_kg_h[height_index, :] = LiquidStreamIn.flow_kg_h
        LiquidProfile.temp_K[height_index, :] = LiquidStreamIn.temp_K
        for id in LiquidProfile.specie.keys():
            LiquidProfile.specie[id]["Mass Fraction"][height_index, :] = LiquidStreamIn.specie[id]["Mass Fraction"]
        for id in LiquidProfile.info.keys():
            LiquidProfile.info[id][height_index, :] = LiquidStreamIn.info[id]
        return LiquidProfile

    def __insert_slice_to_GasProfile__(self, GasProfile, height_index, GasStreamIn):
        GasProfile.flow_kmol_h[height_index, :] = GasStreamIn.flow_kmol_h
        GasProfile.temp_K[height_index, :] = GasStreamIn.temp_K
        for id in GasProfile.specie.keys():
            GasProfile.specie[id]["Molar Fraction"][height_index, :] = GasStreamIn.specie[id]["Molar Fraction"]
        for id in GasProfile.info.keys():
            GasProfile.info[id][height_index, :] = GasStreamIn.info[id]
        return GasProfile

    # ----------------------------------------------------------------------------

    def __get_compressed_LiquidStream__(self, LiquidStreamIn, condition, LiquidStreamOut=None):
        if LiquidStreamOut is None:
            LiquidStreamOut = deepcopy(LiquidStreamIn)
        LiquidStreamOut.temp_K = np.compress(condition=condition, a=LiquidStreamIn.temp_K, axis=0)
        LiquidStreamOut.flow_kg_h = np.compress(condition=condition, a=LiquidStreamIn.flow_kg_h, axis=0)
        for id in LiquidStreamOut.specie.keys():
            LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.compress(condition=condition, a=LiquidStreamIn.get_specie_mass_fraction(id=id), axis=0))
        for id in LiquidStreamOut.info.keys():
            LiquidStreamOut.info[id] = np.compress(condition=condition, a=LiquidStreamIn.info[id], axis=0)
        return LiquidStreamOut

    def __get_compressed_GasStream__(self, GasStreamIn, condition, GasStreamOut=None):
        if GasStreamOut is None:
            GasStreamOut = deepcopy(GasStreamIn)
        GasStreamOut.temp_K = np.compress(condition=condition, a=GasStreamIn.temp_K, axis=0)
        GasStreamOut.flow_kmol_h = np.compress(condition=condition, a=GasStreamIn.flow_kmol_h, axis=0)
        for id in GasStreamOut.specie.keys():
            GasStreamOut.set_specie_molar_fraction(id=id, value=np.compress(condition=condition, a=GasStreamIn.get_specie_molar_fraction(id=id), axis=0))
        for id in GasStreamOut.info.keys():
            GasStreamOut.info[id] = np.compress(condition=condition, a=GasStreamIn.info[id], axis=0)
        return GasStreamOut

    # ----------------------------------------------------------------------------

    def __insert_into_LiquidStream__(self, LiquidStreamIn, condition, LiquidStreamOut):
        i = np.where(condition)[0]
        LiquidStreamOut.temp_K[i] = LiquidStreamIn.temp_K
        LiquidStreamOut.flow_kg_h[i] = LiquidStreamIn.flow_kg_h
        for id in LiquidStreamOut.specie.keys():
            LiquidStreamOut.specie[id]["Mass Fraction"][i] = LiquidStreamIn.specie[id]["Mass Fraction"]
        for id in LiquidStreamOut.info.keys():
            LiquidStreamOut.info[id][i] = LiquidStreamIn.info[id]
        return LiquidStreamOut

    def __insert_into_GasStream__(self, GasStreamIn, condition, GasStreamOut):
        i = np.where(condition)[0]
        GasStreamOut.temp_K[i] = GasStreamIn.temp_K
        GasStreamOut.flow_kmol_h[i] = GasStreamIn.flow_kmol_h
        for id in GasStreamOut.specie.keys():
            GasStreamOut.specie[id]["Molar Fraction"][i] = GasStreamIn.specie[id]["Molar Fraction"]
        for id in GasStreamOut.info.keys():
            GasStreamOut.info[id][i] = GasStreamIn.info[id]
        return GasStreamOut

    # ----------------------------------------------------------------------------

    def __append_LiquidStream__(self, LiquidStream1, LiquidStream2):
        LiquidStream1.temp_K = np.concatenate((LiquidStream1.temp_K, LiquidStream2.temp_K), axis=0)
        LiquidStream1.flow_kg_h = np.concatenate((LiquidStream1.flow_kg_h, LiquidStream2.flow_kg_h), axis=0)
        for id in LiquidStream1.specie.keys():
            LiquidStream1.set_specie_mass_fraction(id=id, value=np.concatenate((LiquidStream1.get_specie_mass_fraction(id=id), LiquidStream2.get_specie_mass_fraction(id=id)), axis=0))
        for id in LiquidStream1.info.keys():
            LiquidStream1.info[id] = np.concatenate((LiquidStream1.info[id], LiquidStream2.info[id]), axis=0)
        return LiquidStream1

    def __append_GasStream__(self, GasStream1, GasStream2):
        GasStream1.temp_K = np.concatenate((GasStream1.temp_K, GasStream2.temp_K), axis=0)
        GasStream1.flow_kmol_h = np.concatenate((GasStream1.flow_kmol_h, GasStream2.flow_kmol_h), axis=0)
        for id in GasStream1.specie.keys():
            GasStream1.set_specie_molar_fraction(id=id, value=np.concatenate((GasStream1.get_specie_molar_fraction(id=id), GasStream2.get_specie_molar_fraction(id=id)), axis=0))
        for id in GasStream1.info.keys():
            GasStream1.info[id] = np.concatenate((GasStream1.info[id], GasStream2.info[id]), axis=0)
        return GasStream1

    def __interpolate_from_LiquidProfile__(self, LiquidStreamProfile, profile_position_m, position_m, LiquidStream=None):

        # Assume "profile_position_m" is descending

        if LiquidStream is None:
            LiquidStream = deepcopy(LiquidStreamProfile)
        #i = np.argmax(profile_position_m >= position_m) - 1
        i = np.argmax(profile_position_m <= position_m) - 1
        LiquidStream.temp_K = LiquidStreamProfile.temp_K[i,:] + ((LiquidStreamProfile.temp_K[i+1,:] - LiquidStreamProfile.temp_K[i,:]) / (profile_position_m[i+1] - profile_position_m[i])) * (position_m - profile_position_m[i])
        LiquidStream.flow_kg_h = LiquidStreamProfile.flow_kg_h[i,:] + ((LiquidStreamProfile.flow_kg_h[i+1,:] - LiquidStreamProfile.flow_kg_h[i,:]) / (profile_position_m[i+1] - profile_position_m[i])) * (position_m - profile_position_m[i])
        for id in LiquidStreamProfile.specie.keys():
            LiquidStream.specie[id]["Mass Fraction"] = LiquidStreamProfile.specie[id]["Mass Fraction"][i, :] + ((LiquidStreamProfile.specie[id]["Mass Fraction"][i + 1, :] - LiquidStreamProfile.specie[id]["Mass Fraction"][i, :]) / (profile_position_m[i + 1] - profile_position_m[i])) * (position_m - profile_position_m[i])
        for id in LiquidStreamProfile.info.keys():
            LiquidStream.info[id] = LiquidStreamProfile.info[id][i, :] + ((LiquidStreamProfile.info[id][i + 1,:] - LiquidStreamProfile.info[id][i,:]) / (profile_position_m[i + 1] -profile_position_m[i])) * (position_m - profile_position_m[i])
        return LiquidStream

    def __interpolate_from_GasProfile__(self, GasStreamProfile, profile_position_m, position_m, GasStream=None):
        # Assume "profile_position_m" is ascending
        if GasStream is None:
            GasStream = deepcopy(GasStreamProfile)
        i = np.argmax(profile_position_m >= position_m) - 1
        GasStream.temp_K = GasStreamProfile.temp_K[i,:] + ((GasStreamProfile.temp_K[i+1,:] - GasStreamProfile.temp_K[i,:]) / (profile_position_m[i+1] - profile_position_m[i])) * (position_m - profile_position_m[i])
        GasStream.flow_kmol_h = GasStreamProfile.flow_kmol_h[i,:] + ((GasStreamProfile.flow_kmol_h[i+1,:] - GasStreamProfile.flow_kmol_h[i,:]) / (profile_position_m[i+1] - profile_position_m[i])) * (position_m - profile_position_m[i])
        GasStream.pressure_bara = GasStreamProfile.pressure_bara[i,:] + ((GasStreamProfile.pressure_bara[i+1,:] - GasStreamProfile.pressure_bara[i,:]) / (profile_position_m[i+1] - profile_position_m[i])) * (position_m - profile_position_m[i])
        for id in GasStreamProfile.specie.keys():
            GasStream.specie[id]["Molar Fraction"] = GasStreamProfile.specie[id]["Molar Fraction"][i,:] + ((GasStreamProfile.specie[id]["Molar Fraction"][i + 1,:] - GasStreamProfile.specie[id]["Molar Fraction"][i]) / (profile_position_m[i + 1] - profile_position_m[i])) * (position_m - profile_position_m[i])
        for id in GasStreamProfile.info.keys():
            GasStream.info[id] = GasStreamProfile.info[id][i,:] + ((GasStreamProfile.info[id][i + 1,:] - GasStreamProfile.info[id][i,:]) / (profile_position_m[i + 1] -profile_position_m[i])) * (position_m - profile_position_m[i])
        return GasStream


class _Flash:

    def __get_Flash_GasStream__(self, LiquidStreamIn):
        GasStreamOut = GasStream()
        l = Library()
        for id in LiquidStreamIn.vapor_pressure_bara.keys():
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"].keys())[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"].keys())[0]
            molar_mass = LiquidStreamIn.get_specie_molar_mass_kg_kmol(id=liq_id)
            charge = LiquidStreamIn.get_specie_charge(id=liq_id)
            l.add_GasStream_specie_gas(id=gas_id, molar_mass_kg_kmol=molar_mass, charge=charge)
        for id in LiquidStreamIn.vapor_pressure_bara.keys():
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"].keys())[0]
            GasStreamOut.add_specie(id=gas_id, library=l)
        GasStreamOut.set_gas_flow_kmol_h(value=np.zeros(shape=(LiquidStreamIn.temp_K.shape[0],), dtype=np.float64))
        for id in GasStreamOut.specie.keys():
            GasStreamOut.set_specie_molar_fraction(id=id, value=np.ones(shape=(LiquidStreamIn.temp_K.shape[0],), dtype=np.float64) / GasStreamOut.num_of_species)
        return GasStreamOut

    def __get_Flash_f_vap_log__(self, LiquidStreamIn, GasStreamIn, Kpw):
        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        f_vap = np.log(Kpw)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            f_vap[:, i] = f_vap[:, i] + np.log(GasStreamIn.get_specie_molar_fraction(id=gas_id))
            f_vap[:, i] = f_vap[:, i] + np.log(GasStreamIn.get_gas_pressure_bara())
            f_vap[:, i] = f_vap[:, i] - np.log(LiquidStreamIn.get_specie_mass_fraction(id=liq_id))
        return f_vap

    def __get_Flash_dfdw_vap_log__(self, LiquidStreamIn):
        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_species = LiquidStreamIn.num_of_species
        dfdw = np.zeros(shape=(num_of_samples, num_of_vap, num_of_species), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            liq_i = LiquidStreamIn.specie[liq_id]["Index"]
            dfdw[:, i, liq_i] = - 1 / LiquidStreamIn.get_specie_mass_fraction(id=liq_id)
        return dfdw

    def __get_Flash_dfdy_vap_loq__(self, GasStreamIn, LiquidStreamIn):
        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_species = GasStreamIn.num_of_species
        dfdy = np.zeros(shape=(num_of_samples, num_of_vap, num_of_species), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            gas_i = GasStreamIn.specie[gas_id]["Index"]
            dfdy[:, i, gas_i] = 1 / GasStreamIn.get_specie_molar_fraction(id=gas_id)
        return dfdy

    def __get_Flash_dfdT_vap_log__(self, Kyw, dKywdT):
        dfdT = dKywdT / Kyw
        return dfdT[:, :, None]



# ---------------------------------------------------------------------------------------


class GasStream(_Serializer):

    def __init__(self):
        super().__init__()

        # Container for Various Info
        self.info = {}

        # Universal Gas Constant (bar.m3 / kmol.K)
        self.R = 8.314 * 10**(-2)

        # Concentration
        self.specie = {}
        self.rxn_insta = {}
        self.rxn_reversible = {}
        self.rxn_irreversible = {}

        self.num_of_species = 0
        self.num_of_rxn_insta = 0
        self.num_of_rxn_reversible = 0
        self.num_of_rxn_irreversible = 0

        # Values
        self.temp_K = None
        self.flow_kmol_h = None
        self.pressure_bara = None

        # Functions
        self.viscosity_Pas = None
        self.thermal_conductivity_kW_mK = None
        self.heat_capacity_kJ_kmolK = None
        self.diffusivity_m2_s = None

    def add_info(self, key, value):
        self.info[key] = value

    def add_specie(self, id, library):
        self.specie[id] = {}
        self.specie[id]["Index"] = self.num_of_species
        for k in library.GasStream_specie_gas[id].keys():
            self.specie[id][k] = library.GasStream_specie_gas[id][k]
        self.specie[id]["Molar Fraction"] = None
        self.num_of_species = self.num_of_species + 1

    def add_rxn_insta(self, id, library):
        self.rxn_insta[id] = {}
        for k in library.GasStream_rxn_insta[id].keys():
            self.rxn_insta[id][k] = library.GasStream_rxn_insta[id][k]
        self.rxn_insta[id]["Index"] = self.num_of_rxn_insta
        self.num_of_rxn_insta = self.num_of_rxn_insta + 1

    def add_rxn_reversible(self, id, library):
        self.rxn_reversible[id] = {}
        for k in library.GasStream_rxn_reversible[id].keys():
            self.rxn_reversible[id][k] = library.GasStream_rxn_reversible[id][k]
        self.rxn_reversible[id]["Index"] = self.num_of_rxn_reversible
        self.num_of_rxn_reversible = self.num_of_rxn_reversible + 1

    def add_rxn_irreversible(self, id, library):
        self.rxn_irreversible[id] = {}
        for k in library.GasStream_rxn_irreversible[id].keys():
            self.rxn_irreversible[id][k] = library.GasStream_rxn_irreversible[id][k]
        self.rxn_irreversible[id]["Index"] = self.num_of_rxn_irreversible
        self.num_of_rxn_irreversible = self.num_of_rxn_irreversible + 1

    # --------------------------------------------------------------------------------------------

    def load_viscosity_Pas(self, id, library):
        self.viscosity_Pas = library.GasStream_viscosity_Pas[id]

    def load_thermal_conductivity_kW_mK(self, id, library):
        self.thermal_conductivity_kW_mK = library.GasStream_thermal_conductivity_kW_mK[id]

    def load_heat_capacity_kJ_kmolK(self, id, library):
        self.heat_capacity_kJ_kmolK = library.GasStream_heat_capacity_kJ_kmolK[id]

    def load_diffusivity_m2_s(self, id, library):
        self.diffusivity_m2_s = library.GasStream_diffusivity_m2_s[id]

    # --------------------------------------------------------------------------------------------

    def set_gas_temp_K(self, value):
        self.temp_K = value

    def set_gas_flow_kmol_h(self, value):
        self.flow_kmol_h = value

    def set_gas_pressure_bara(self, value):
        self.pressure_bara = value

    def set_specie_molar_fraction(self, id, value):
        self.specie[id]["Molar Fraction"] = value

    # --------------------------------------------------------------------------------------------

    def get_gas_temp_K(self):
        return self.temp_K

    def get_gas_flow_kmol_h(self):
        return self.flow_kmol_h

    def get_gas_flow_kg_h(self):
        m = 0
        for id in self.specie.keys():
            y = self.get_specie_molar_fraction(id=id)
            M = self.get_specie_molar_mass_kg_kmol(id=id)
            m = m + y * M * self.get_gas_flow_kmol_h()
        return m

    def get_gas_pressure_bara(self):
        return self.pressure_bara

    def get_gas_molarity_kmol_m3(self):
        return self.pressure_bara / (0.08314 * self.get_gas_temp_K())

    def get_gas_density_kg_m3(self):
        c = self.get_gas_molarity_kmol_m3()
        rho = 0
        for id in self.specie.keys():
            rho = rho + self.specie[id]["Molar Fraction"] * c * self.specie[id]["Molar Mass [kg/kmol]"]
        return rho

    def get_gas_heat_capacity_kJ_kmolK(self):
        cp = 0
        for id in self.specie.keys():
            cp = cp + self.get_specie_molar_fraction(id=id) * self.get_specie_heat_capacity_kJ_kmolK(id=id)
        return cp

    def get_gas_viscosity_Pas(self):
        return self.viscosity_Pas(self)

    def get_gas_thermal_conductivity_kW_mK(self):
        return  self.thermal_conductivity_kW_mK(self)

    # --------------------------------------------------------------------------------------------

    def get_specie_molar_mass_kg_kmol(self, id):
        return self.specie[id]["Molar Mass [kg/kmol]"]

    def get_specie_flow_kmol_h(self, id):
        return self.get_gas_flow_kmol_h() * self.get_specie_molar_fraction(id=id)

    def get_specie_flow_kg_h(self, id):
        return self.get_specie_flow_kmol_h(id=id) * self.get_specie_molar_mass_kg_kmol(id=id)

    def get_specie_molar_fraction(self, id):
        return self.specie[id]["Molar Fraction"]

    def get_specie_pressure_bara(self, id):
        return self.specie[id]["Molar Fraction"] * self.get_gas_pressure_bara()

    def get_specie_molarity_kmol_m3(self, id):
        return self.get_specie_molar_fraction(id=id) * self.molarity_kmol_m3

    def get_specie_heat_capacity_kJ_kmolK(self, id):
        return self.heat_capacity_kJ_kmolK(self, id)

    def get_specie_ppm_dry(self, id, water_id):
        y_wet = self.get_specie_molar_fraction(id)
        y_water = self.get_specie_molar_fraction(water_id)
        y_dry = y_wet / (1 - y_water)
        return 10**6 * y_dry

    def get_specie_diffusivity_m2_s(self, id):
        return self.diffusivity_m2_s(self, id)

    # --------------------------------------------------------------------------------------------

    def get_rxn_insta_equilibrium_constant_wrt_partial_pressures(self, id):
        Kp = self.rxn_insta[id]["K"](self)
        return Kp

    def get_rxn_insta_equilibrium_constant_wrt_molar_fractions(self, id):
        Kp = self.get_rxn_insta_equilibrium_constant_wrt_partial_pressures(id)
        p = self.get_gas_pressure_bara()
        s = 0
        for k in self.rxn_insta[id]["Stoch"].keys():
            s = s + self.rxn_insta[id]["Stoch"][k]
        Ky = Kp * p ** (-s)
        return Ky

    def get_rxn_insta_exothermic_heat_kJ_kmol(self, id):
        dT = 0.05

        K0 = self.get_rxn_insta_equilibrium_constant_wrt_partial_pressures(id)
        T0 = self.temp_K

        self.temp_K = self.temp_K + dT

        K1 = self.get_rxn_insta_equilibrium_constant_wrt_partial_pressures(id)
        T1 = self.temp_K

        self.temp_K = self.temp_K - dT

        q = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))
        return q

    # --------------------------------------------------------------------------------------------

    def get_rxn_reversible_exothermic_heat_kJ_kmol(self, id):
        return self.rxn_reversible[id]["Exothermic Heat [kJ/kmol]"](self)

    def get_rxn_reversible_rate_kmol_m3s(self, id):
        return self.rxn_reversible[id]["Rate [kmol/m3.s]"](self)

    def get_rxn_reversible_rate_kmol_m3s_as_specie_vector(self):
        r_specie = np.zeros(shape=(self.temp_K.shape[0], self.num_of_species), dtype=np.float64)
        for id in self.rxn_reversible.keys():
            f, b = self.get_rxn_reversible_rate_kmol_m3s(id)
            r_rxn = f - b
            nu = self.rxn_reversible[id]["Stoch"]
            for el in nu.keys():
                i = self.specie[el]["Index"]
                r_specie[:, i] = r_specie[:, i] + r_rxn * nu[el]
        return r_specie

    # ------------------------------------------------------------------------------------

    def normalize_molar_fractions(self):
        y_gas = 0
        y_specie = {}
        for id in self.specie.keys():
            y_specie[id] = self.get_specie_molar_fraction(id=id)
            y_gas = y_gas + y_specie[id]
        for id in self.specie.keys():
            self.set_specie_molar_fraction(id=id, value=y_specie[id] / y_gas)
        return y_gas

    def __molar_fractions_dic2vec__(self):
        y = np.zeros(shape=(len(self.temp_K), self.num_of_species), dtype=np.float64)
        for id in self.specie.keys():
            i = self.specie[id]["Index"]
            y[:, i] = self.get_specie_molar_fraction(id=id)
        return y

    def __molar_fractions_vec2dic__(self, y):
        for i, id in enumerate(self.specie.keys()):
            self.set_specie_molar_fraction(id=id, value=y[:, i])

    def __massfrac2molefrac__(self, w):
        den = 0
        y = 1.0 * w
        for i, id in enumerate(self.specie.keys()):
            den = den + w[:, i] / self.get_specie_molar_mass_kg_kmol(id)
        for i, id in enumerate(self.specie.keys()):
            y[:, i] = (w[:, i] / self.get_specie_molar_mass_kg_kmol(id)) / den
        return y

    def __molefrac2massfrac__(self, y):
        den = 0
        w = 1.0 * y
        for i, id in enumerate(self.specie.keys()):
            den = den + y[:, i] * self.get_specie_molar_mass_kg_kmol(id)
        for i, id in enumerate(self.specie.keys()):
            w[:, i] = (y[:, i] * self.get_specie_molar_mass_kg_kmol(id)) / den
        return w

    def __dydw__(self, w):
        dydw = np.zeros(shape=(w.shape[0], w.shape[1], w.shape[1]), dtype=np.float64)
        M = np.zeros(shape=(w.shape[0], w.shape[1]), dtype=np.float64)
        for i, id in enumerate(self.specie.keys()):
            M[:, i] = self.get_specie_molar_mass_kg_kmol(id=id)
        den = 0
        for i in range(self.num_of_species):
            den = den + w[:, i] / M[:, i]
        for i in range(self.num_of_species):
            for j in range(self.num_of_species):
                if i == j:
                    dydw[:, i, j] = (1 / M[:,i]) / den - (w[:,i] / M[:,i] ** 2) / den**2
                else:
                    dydw[:, i, j] = - (w[:,i] / (M[:,i] * M[:,j])) / den
        return dydw


class LiquidStream(_Serializer):

    def __init__(self, stream_id, solvent_id):

        super().__init__()

        # Unique Solvent Identification
        self.id = stream_id

        # Solvent (E.g. H2O) when calulating molality
        self.solvent_id = solvent_id

        # Container for Various Info
        self.info = {}

        # Concentration
        self.specie = {}
        self.rxn_insta = {}
        self.rxn_reversible = {}
        self.rxn_irreversible = {}
        self.vapor_pressure_bara = {}

        self.num_of_species = 0
        self.num_of_rxn_insta = 0
        self.num_of_rxn_reversible = 0
        self.num_of_rxn_irreversible = 0
        self.num_of_vapor_pressure_bara = 0

        # ------------------------

        # Values
        self.temp_K = None
        self.flow_kg_h = None

        # Functions
        self.density_kg_m3 = None
        self.heat_capacity_kJ_kgK = None
        self.viscosity_Pas = None
        self.thermal_conductivity_kW_mK = None
        self.activity_coefficient = None
        self.diffusivity_m2_s = None
        self.surface_tension_N_m = None

    # --------------------------------------------------------------------

    def load_density(self, id, library):
        self.density_kg_m3 = library.LiquidStream_density_kg_m3[id]

    def load_heat_capacity(self, id, library):
        self.heat_capacity_kJ_kgK = library.LiquidStream_heat_capacity_kJ_kgK[id]

    def load_viscosity(self, id, library):
        self.viscosity_Pas = library.LiquidStream_viscosity_Pas[id]

    def load_diffusivity(self, id, library):
        self.diffusivity_m2_s = library.LiquidStream_diffusivity_m2_s[id]

    def load_activity_coefficient(self, id, library):
        self.activity_coefficient = library.LiquidStream_activity_coefficient[id]

    def load_thermal_conductivity_kW_mK(self, id, library):
        self.thermal_conductivity_kW_mK = library.LiquidStream_thermal_conductivity_kW_mK[id]

    def load_surface_tension_N_m(self, id, library):
        self.surface_tension_N_m = library.LiquidStream_surface_tension_N_m[id]

    # --------------------------------------------------------------------

    def add_info(self, key, value):
        self.info[key] = value

    def add_specie(self, id, library):
        self.specie[id] = {}
        for k in library.LiquidStream_specie_liq[id].keys():
            self.specie[id][k] = library.LiquidStream_specie_liq[id][k]
        self.specie[id]["Mass Fraction"] = None
        self.specie[id]["Index"] = self.num_of_species
        self.num_of_species = self.num_of_species + 1

    def add_rxn_insta(self, id, library):
        self.rxn_insta[id] = {}
        for k in library.LiquidStream_rxn_insta[id].keys():
            self.rxn_insta[id][k] = library.LiquidStream_rxn_insta[id][k]
        self.rxn_insta[id]["Index"] = self.num_of_rxn_insta
        self.num_of_rxn_insta = self.num_of_rxn_insta + 1

    def add_rxn_reversible(self, id, library):
        self.rxn_reversible[id] = {}
        for k in library.LiquidStream_rxn_reversible[id].keys():
            self.rxn_reversible[id][k] = library.LiquidStream_rxn_reversible[id][k]
        self.rxn_reversible[id]["Index"] = self.num_of_rxn_reversible
        self.num_of_rxn_reversible = self.num_of_rxn_reversible + 1

    def add_rxn_irreversible(self, id, library):
        self.rxn_irreversible[id] = {}
        for k in library.LiquidStream_rxn_irreversible[id].keys():
            self.rxn_irreversible[id][k] = library.LiquidStream_rxn_irreversible[id][k]
        self.rxn_irreversible[id]["Index"] = self.num_of_rxn_irreversible
        self.num_of_rxn_irreversible = self.num_of_rxn_irreversible + 1

    def add_vapor_pressure_bara(self, id, library):
        self.vapor_pressure_bara[id] = {}
        for k in library.LiquidStream_vapor_pressure_bara[id].keys():
            self.vapor_pressure_bara[id][k] = library.LiquidStream_vapor_pressure_bara[id][k]
        self.vapor_pressure_bara[id]["Index"] = self.num_of_vapor_pressure_bara
        self.num_of_vapor_pressure_bara = self.num_of_vapor_pressure_bara + 1

    # --------------------------------------------------------------------

    def print_info(self):

        for id in self.specie.keys():
            print("Liquid Specie: \t\t\t" + id)

        for id in self.equilibrium_liq.keys():
            print("Equilibrium Liquid: \t" + id)

        for id in self.reaction_liq.keys():
            print("Reaction Liquid:    \t" + id)

        for id in self.equilibrium_vap.keys():
            print("Equilibrium Vapor:  \t" + id)

    # --------------------------------------------------------------------

    def set_solution_temp_K(self, value):
        self.temp_K = value.astype(np.float64)

    def set_solution_flow_kg_h(self, value):
        self.flow_kg_h = value.astype(np.float64)

    def set_specie_mass_fraction(self, id, value):
        self.specie[id]["Mass Fraction"] = value.astype(np.float64)

    def set_species_molar_fraction(self, molar_fractions):

        m = 0
        for id in molar_fractions.keys():
            self.specie[id]["Mass Fraction"] = molar_fractions[id] * self.specie[id]["Molar Mass [kg/kmol]"]
            m = m + self.specie[id]["Mass Fraction"]

        for id in molar_fractions.keys():
            self.specie[id]["Mass Fraction"] = self.specie[id]["Mass Fraction"] / m

    def set_species_molality(self, solutes_molality_mol_kg):

        # Assuming we start with one kg of solvent
        mass_solvent = 1.0
        mass_tot = mass_solvent
        mass_solute = {}

        # Total Mass and Mass of Solutes
        for id in solutes_molality_mol_kg.keys():
            M = self.get_specie_molar_mass_kg_kmol(id=id)
            mass_tot = mass_tot + solutes_molality_mol_kg[id] * (M * 0.001) * mass_solvent                  # kg = mol/kg H2O * kg/mol * kg H2O
            mass_solute[id] = solutes_molality_mol_kg[id] * (M * 0.001) * mass_solvent                      # kg = mol/kg H2O * kg/mol * kg H2O

        # Mass Fractions
        self.specie[self.solvent_id]["Mass Fraction"] = mass_solvent / mass_tot
        for id in self.specie.keys():
            if id != self.solvent_id:
                if id not in solutes_molality_mol_kg.keys():
                    self.specie[id]["Mass Fraction"] = 0 * self.specie[solvent_id]["Mass Fraction"]
                else:
                    M = self.get_specie_molar_mass_kg_kmol(id=id)
                    self.specie[id]["Mass Fraction"] = mass_solute[id] / mass_tot

    # --------------------------------------------------------------------

    def copy_mass_fractions(self, LiquidStreamIn):
        for id in LiquidStreamIn.specie.keys():
            self.set_specie_mass_fraction(id=id, value=LiquidStreamIn.get_specie_mass_fraction(id))

    # --------------------------------------------------------------------

    def get_solution_temp_K(self):
        T = self.temp_K
        return self.temp_K

    def get_solution_flow_kg_h(self):
        return self.flow_kg_h

    def get_solution_flow_m3_h(self):
        return self.get_solution_flow_kg_h() / self.get_solution_density_kg_m3()

    def get_solution_flow_kmol_h(self):
        return self.get_solution_molarity_kmol_m3() * self.get_solution_flow_m3_h()

    def get_solution_density_kg_m3(self):
        return self.density_kg_m3(self)

    def get_solution_molarity_kmol_m3(self):
        rho = self.get_solution_density_kg_m3()
        c = 0
        for id in self.specie.keys():
            c = c + self.specie[id]["Mass Fraction"] * rho / self.specie[id]["Molar Mass [kg/kmol]"]
        return c

    def get_solution_heat_capacity_kJ_kgK(self):
        return self.heat_capacity_kJ_kgK(self)

    def get_solution_ionic_strength_kmol_m3(self):
        # Ions are applied induvidually and not in pairs... The factor 1/2 are therefore included
        I = 0
        for id in self.specie.keys():
            z = self.specie[id]["Charge"]
            I = I + 0.5 * self.get_specie_molarity_kmol_m3(id=id) * z ** 2
        return I

    def get_solution_ionic_strength_mol_kg(self):
        # Ions are applied induvidually and not in pairs... The factor 1/2 are therefore included
        I = 0
        for id in self.specie.keys():
            z = self.specie[id]["Charge"]
            m = self.get_specie_molality_mol_kg(id=id)
            I = I + 0.5 * m * z ** 2
        return I

    def get_solution_viscosity_Pas(self):
        return self.viscosity_Pas(self)

    def get_solution_surface_tension_N_m(self):
        return self.surface_tension_N_m(self)

    def get_solution_vapor_pressure_bara(self):
        p = 0
        for reaction_id in self.vapor_pressure_bara.keys():
            gas_id = list(self.vapor_pressure_bara[reaction_id]["Stoch Gas"].keys())[0]
            p = p + self.get_specie_vapor_pressure_bara(gas_id)
        return p

    # --------------------------------------------------------------------

    def get_specie_molar_mass_kg_kmol(self, id):
        return self.specie[id]["Molar Mass [kg/kmol]"]

    def get_specie_charge(self, id):
        return self.specie[id]["Charge"]

    def get_specie_molarity_kmol_m3(self, id):
        rho = self.get_solution_density_kg_m3()
        return rho * self.get_specie_mass_fraction(id) / self.get_specie_molar_mass_kg_kmol(id)

    def get_specie_molality_mol_kg(self, id):
        m_tot = 1
        m_solvent = self.get_specie_mass_fraction(id=self.solvent_id)
        m_solute = self.get_specie_mass_fraction(id=id)
        M_solute  = self.get_specie_molar_mass_kg_kmol(id=id) * 0.001
        n_solute = (1 / M_solute) * m_solute
        molality = n_solute / m_solvent
        return molality

    def get_specie_activity_kmol_m3(self, id):
        return self.get_specie_molarity_kmol_m3(id=id) * self.get_specie_activity_coefficient(id=id)

    def get_specie_activity_coefficient(self, id):
        return self.activity_coefficient(self, id)

    def get_specie_molar_fraction(self, id):
        enu = self.specie[id]["Mass Fraction"] / self.specie[id]["Molar Mass [kg/kmol]"]
        den = 0
        for k in self.specie.keys():
            den = den + self.specie[k]["Mass Fraction"] / self.specie[k]["Molar Mass [kg/kmol]"]
        x = enu / den
        return x

    def get_specie_mass_fraction(self, id):
        w = self.specie[id]["Mass Fraction"]
        return w

    def get_specie_diffusivity_m2_s(self, id):
        return self.diffusivity_m2_s(self, id)

    def get_specie_flow_kg_h(self, id):
        return self.get_specie_mass_fraction(id=id) * self.get_solution_flow_kg_h()

    def get_specie_flow_kmol_h(self, id):
        return self.get_specie_molarity_kmol_m3(id=id) * self.get_solution_flow_m3_h()

    def get_specie_vapor_pressure_bara(self, gas_id):

        for id in self.vapor_pressure_bara.keys():
            if gas_id in self.vapor_pressure_bara[id]["Stoch Gas"].keys():
                reaction_id = id
        liq_id = list(self.vapor_pressure_bara[reaction_id]["Stoch Liq"].keys())[0]
        unit = self.vapor_pressure_bara[reaction_id]["Unit Liq"][liq_id]
        gamma = self.get_specie_activity_coefficient(id=liq_id)

        if self.vapor_pressure_bara[reaction_id]["H"] == None:
            K = 1 / self.vapor_pressure_bara[reaction_id]["p0"](self)
        else:
            K = self.vapor_pressure_bara[reaction_id]["H"](self)

        if unit == "c":
            concentration = gamma * self.get_specie_molarity_kmol_m3(id=liq_id)
        elif unit == "m":
            concentration = gamma * self.get_specie_molality_mol_kg(id=liq_id)
        elif unit == "x":
            concentration = gamma * self.get_specie_molar_fraction(id=liq_id)
        elif unit == "w":
            concentration = gamma * self.get_specie_mass_fraction(id=id)
        elif unit == None:
            concentration = gamma * np.ones(shape=(self.temp_K.shape))
        p = concentration/ K
        return p

    # ------------------------------------------------------------------------------------

    def get_rxn_insta_equilibrium_constant_wrt_activities(self, id):
        Ka = self.rxn_insta[id]["K"](self)
        return Ka

    def get_rxn_insta_equilibrium_constant_wrt_concentrations(self, id):
        Kc = self.rxn_insta[id]["K"](self)
        for r in self.rxn_insta[id]["Unit"].keys():
            power = self.rxn_insta[id]["Stoch"][r]
            Kc = Kc * self.get_specie_activity_coefficient(id=r) ** (-power)
        return Kc

    def get_rxn_insta_equilibrium_constant_wrt_mass_fractions(self, id):
        den = 0
        for specie_id in self.specie.keys():
            den = den + self.get_specie_mass_fraction(specie_id) / self.get_specie_molar_mass_kg_kmol(specie_id)
        c = self.get_solution_molarity_kmol_m3()
        rho = self.get_solution_density_kg_m3()
        Kw = self.get_rxn_insta_equilibrium_constant_wrt_concentrations(id=id)
        for r in self.rxn_insta[id]["Unit"].keys():
            power = self.rxn_insta[id]["Stoch"][r]
            if self.rxn_insta[id]["Unit"][r] == "c":
                Kw = Kw * (den * self.get_specie_molar_mass_kg_kmol(r) / c) ** power
            elif self.rxn_insta[id]["Unit"][r] == "x":
                Kw = Kw * (den * self.get_specie_molar_mass_kg_kmol(r)) ** power
            elif self.rxn_insta[id]["Unit"][r] == "m":
                Kw = Kw * ((rho * self.get_specie_mass_fraction(id=self.solvent_id) * self.get_specie_molar_mass_kg_kmol(r) * den) / (1000 * c)) ** power
            elif self.rxn_insta[id]["Unit"][r] == None:
                Kw = Kw * (1 / self.get_specie_mass_fraction(id=r)) ** power
            elif self.rxn_insta[id]["Unit"][r] == "w":
                pass
        return Kw

    def get_rxn_insta_equilibrium_constant_wrt_mass_fractions_temperature_gradient(self, id):
        dT = 0.05
        K0 = self.get_rxn_insta_equilibrium_constant_wrt_mass_fractions(id)
        T0 = self.temp_K
        self.temp_K = self.temp_K + dT
        K1 = self.get_rxn_insta_equilibrium_constant_wrt_mass_fractions(id)
        T1 = self.temp_K
        self.temp_K = self.temp_K - dT
        return (K1 - K0) / dT

    def get_rxn_insta_exothermic_heat_kJ_kmol(self, id):
        dT = 0.05
        K0 = self.get_rxn_insta_equilibrium_constant_wrt_activities(id)
        T0 = 1.0*self.temp_K
        self.temp_K = self.temp_K + dT
        K1 = self.get_rxn_insta_equilibrium_constant_wrt_activities(id)
        T1 = self.temp_K
        self.temp_K = self.temp_K - dT
        q = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))
        return q

    # ------------------------------------------------------------------------------------

    def get_rxn_reversible_exothermic_heat_kJ_kmol(self, id):
        return self.rxn_reversible[id]["Exothermic Heat [kJ/kmol]"](self)

    def get_rxn_reversible_rate_kmol_m3s(self, id):
        return self.rxn_reversible[id]["Rate [kmol/m3.s]"](self)

    def get_rxn_reversible_rate_kmol_m3s_as_specie_vector(self):
        r_specie = np.zeros(shape=(self.temp_K.shape[0], self.num_of_species), dtype=np.float64)
        for id in self.rxn_reversible.keys():
            f, b = self.get_rxn_reversible_rate_kmol_m3s(id)
            r_rxn = f - b
            nu = self.rxn_reversible[id]["Stoch"]
            for el in nu.keys():
                i = self.specie[el]["Index"]
                r_specie[:,i] = r_specie[:,i] + r_rxn * nu[el]
        return r_specie

    # ------------------------------------------------------------------------------------

    def normalize_mass_fractions(self):
        w_sol = 0
        w_specie = {}
        for id in self.specie.keys():
            w_specie[id] = self.get_specie_mass_fraction(id=id)
            w_sol = w_sol + w_specie[id]
        for id in self.specie.keys():
            self.set_specie_mass_fraction(id=id, value=w_specie[id] / w_sol)

    # ------------------------------------------------------------------------------------

    def __mass_fractions_dic2vec__(self):
        w = np.zeros(shape=(len(self.temp_K), self.num_of_species), dtype=np.float64)
        for id in self.specie.keys():
            i = self.specie[id]["Index"]
            w[:, i] = self.get_specie_mass_fraction(id=id)
        return w

    def __mass_fractions_vec2dic__(self, w):
        for i, id in enumerate(self.specie.keys()):
            self.set_specie_mass_fraction(id=id, value=w[:, i])

    def __exothermic_heat_kJ_kmol_as_vector_rxn_insta__(self):
        q = np.zeros(shape=(self.temp_K.shape[0], self.num_of_rxn_insta),dtype=np.float64)
        for i, id in enumerate(self.rxn_insta.keys()):
            q[:, i] = self.get_rxn_insta_exothermic_heat_kJ_kmol(id)
        return q

    def __exothermic_heat_kJ_kmol_as_vector_rxn_reversible__(self):
        q = np.zeros(shape=(self.temp_K.shape[0], self.num_of_rxn_reversible), dtype=np.float64)
        for i, id in enumerate(self.rxn_reversible.keys()):
            q[:, i] = self.get_rxn_reversible_exothermic_heat_kJ_kmol(id)
        return q


def LiquidSum(streams):
    LiquidStreamOut = deepcopy(streams[0])
    LiquidStreamOut.temp_K = 0
    LiquidStreamOut.flow_kg_h = 0

    for stream in streams:
        LiquidStreamOut.flow_kg_h = LiquidStreamOut.flow_kg_h + np.abs(stream.flow_kg_h)

    for id in streams[0].specie.keys():
        m = 0
        for stream in streams:
            if id in stream.specie.keys():
                m = m + stream.get_specie_flow_kg_h(id=id)
        LiquidStreamOut.specie[id]["Mass Fraction"] = m / LiquidStreamOut.flow_kg_h

    for stream in streams:
        T = stream.temp_K
        r = np.abs(stream.flow_kg_h) / LiquidStreamOut.flow_kg_h
        LiquidStreamOut.temp_K = LiquidStreamOut.temp_K + r * T

    LiquidStreamOut.normalize_mass_fractions()

    return LiquidStreamOut


def GasSum(streams):
    GasStreamOut = deepcopy(streams[0])
    GasStreamOut.temp_K = 0
    GasStreamOut.flow_kmol_h = 0

    for stream in streams:
        GasStreamOut.flow_kmol_h = GasStreamOut.flow_kmol_h + np.abs(stream.flow_kmol_h)

    for id in streams[0].specie.keys():
        n = 0
        for stream in streams:
            n = n + stream.get_specie_flow_kmol_h(id=id)
        GasStreamOut.specie[id]["Molar Fraction"] = n / GasStreamOut.flow_kmol_h

    for stream in streams:
        T = stream.temp_K
        r = np.abs(stream.flow_kmol_h) / GasStreamOut.flow_kmol_h
        GasStreamOut.temp_K = GasStreamOut.temp_K + r * T

    GasStreamOut.normalize()
    return GasStreamOut


# ---------------------------------------------------------------------------------------


class LiquidEquilibrium_Isothermal(_Stochiometry, _Serializer, _LiquidEquilibrium):

    def __init__(self):
        self.firstscan = True

    def react(self, LiquidStreamIn, lr=0.75):

        if self.firstscan:
            self.matrix = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True)

        self.LiquidStreamInit = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamInit.specie.keys():
            self.LiquidStreamInit.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamInit.specie[id]["Mass Fraction"], 10**(-18))
        self.LiquidStreamInit.normalize_mass_fractions()
        self.LiquidStreamOut = self.__get_LiquidEquilibrium__(self.LiquidStreamInit, b=None, matrix=self.matrix, lr=lr)
        return self.LiquidStreamOut

    def get_adiabatic_temperature_increase_K(self):
        h_in = self.LiquidStreamInit.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h = (h_in + h_out) / 2
        cp_in = self.LiquidStreamInit.get_solution_heat_capacity_kJ_kgK()
        cp_out = self.LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        cp = (cp_in + cp_out) / 2
        dw = self.LiquidStreamOut.__mass_fractions_dic2vec__() - self.LiquidStreamInit.__mass_fractions_dic2vec__()
        dT = np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw)) / cp
        return dT


class LiquidEquilibrium_Adiabatic(_Stochiometry, _Serializer, _LiquidEquilibrium):

    def __init__(self):
        self.firstscan = True

    def react(self, LiquidStreamIn, lr=0.75):

        # Load Inlet Stream
        self.LiquidStreamInit = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamInit.specie.keys():
            self.LiquidStreamInit.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamInit.specie[id]["Mass Fraction"], 10**(-18))
        self.LiquidStreamInit.normalize_mass_fractions()

        # Reactions Exothermic Heat (h0) at Inlet Condition
        Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamInit)
        dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStreamInit, Kw)
        h0 = self.__get_LiquidEquilibrium_exothermic_heat_kJ_kmol__(self.LiquidStreamInit, Kw, dKwdT)

        if self.firstscan:
            self.matrix = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True)
            self.firstscan = False

        # Conserved Quantity
        b = np.einsum("cw,sw->sc", self.matrix["A"], self.LiquidStreamInit.__mass_fractions_dic2vec__())

        # Initial Guess
        eq_iso = LiquidEquilibrium_Isothermal()
        self.LiquidStreamOut = deepcopy(self.LiquidStreamInit)
        #self.LiquidStreamOut.temp_K = self.LiquidStreamOut.temp_K + 3600 * heat_kW / (self.LiquidStreamOut.get_solution_flow_kg_h() * self.LiquidStreamOut.get_solution_heat_capacity_kJ_kgK())
        self.LiquidStreamOut = eq_iso.react(self.LiquidStreamOut, lr)
        self.LiquidStreamOut.temp_K = self.LiquidStreamOut.temp_K + eq_iso.get_adiabatic_temperature_increase_K()
        w = self.LiquidStreamOut.__mass_fractions_dic2vec__()


        # Iterate Until Convergence
        converged = False
        epoch = 0
        iterations = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0],))

        while converged == False:

            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamOut)
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStreamOut, Kw)
            h = (self.__get_LiquidEquilibrium_exothermic_heat_kJ_kmol__(self.LiquidStreamOut, Kw, dKwdT) + h0) / 2

            f_rxn_insta = self.__get_LiquidEquilibrium_f_rxn_insta__(self.LiquidStreamOut, Kw)
            f_mass_balance = self.__get_LiquidEquilibrium_f_mass_balance__(self.LiquidStreamOut, b, self.matrix)
            f_energy_balance = self.__get_LiquidEquilibrium_f_energy_balance__(self.LiquidStreamInit, self.LiquidStreamOut, h, self.matrix)
            f = np.hstack((f_rxn_insta, f_mass_balance, f_energy_balance))

            dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta__(self.LiquidStreamOut, Kw, self.matrix)
            dfdw_mass_balance = self.__get_LiquidEquilibrium_dfdw_mass_balance__(self.LiquidStreamOut, self.matrix)
            dfdw_energy_balance = self.__get_LiquidEquilibrium_dfdw_energy_balance__(self.LiquidStreamInit, self.LiquidStreamOut, h, self.matrix)
            dfdw = np.concatenate((dfdw_rxn_insta, dfdw_mass_balance, dfdw_energy_balance), axis=1)

            dfdT_rxn_insta = self.__get_LiquidEquilibrium_dfdT_rxn_insta__(self.LiquidStreamOut, dKwdT)
            dfdT_mass_balance = self.__get_LiquidEquilibrium_dfdT_mass_balance__(self.LiquidStreamOut, self.matrix)
            dfdT_energy_balance = self.__get_LiquidEquilibrium_dfdT_energy_balance__(self.LiquidStreamOut)
            dfdT = np.concatenate((dfdT_rxn_insta, dfdT_mass_balance, dfdT_energy_balance), axis=1)

            dfdX = np.concatenate((dfdw, dfdT), axis=2)

            dX_newton = - np.linalg.solve(dfdX, f)
            dw_newton = dX_newton[:, :self.LiquidStreamOut.num_of_species:]
            dT_newton = dX_newton[:, self.LiquidStreamOut.num_of_species]

            dw = lr * dw_newton
            dT = lr * dT_newton

            w = w + dw
            w = np.maximum(w, 0)
            self.LiquidStreamOut.__mass_fractions_vec2dic__(w)
            self.LiquidStreamOut.temp_K = self.LiquidStreamOut.temp_K + dT

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        #print("Epochs: \t" + str(epoch))
        return self.LiquidStreamOut


class LiquidEquilibrium_QPFlash(_Stochiometry, _Serializer, _LiquidEquilibrium, _GasEquilibrium, _VaporLiquidEquilibrium, _Flash, _StreamFunctions):

    def __init__(self):
        self.firstscan = True
        self.eq_adiabatic = LiquidEquilibrium_Adiabatic()

    def react(self, LiquidStreamIn, pressure_bara, heat_kW, lr=0.25):

        # Heat Up Solution and Perform Equilibrium Calculation
        #LiquidStreamHeated = deepcopy(LiquidStreamIn)
        #LiquidStreamHeated.temp_K = LiquidStreamHeated.temp_K + (3600 * heat_kW) / (LiquidStreamHeated.get_solution_heat_capacity_kJ_kgK() * LiquidStreamHeated.get_solution_flow_kg_h())
        #LiquidStreamHeated = self.eq_adiabatic.react(LiquidStreamHeated, lr=lr)

        # Continue with Samples where Flashing Occur
        #vapor_pressure_bara = 0
        #for id in LiquidStreamHeated.vapor_pressure_bara.keys():
        #    gas_id = list(LiquidStreamHeated.vapor_pressure_bara[id]["Stoch Gas"].keys())[0]
        #    vapor_pressure_bara = vapor_pressure_bara + LiquidStreamHeated.get_specie_vapor_pressure_bara(gas_id)

        #LiquidStreamBoiled = self.__get_compressed_LiquidStream__(LiquidStreamIn=LiquidStreamHeated, condition=(vapor_pressure_bara > pressure_bara))
        #LiquidStreamHeated = self.__get_compressed_LiquidStream__(LiquidStreamIn=LiquidStreamHeated, condition=(vapor_pressure_bara <= pressure_bara))
        #GasStreamBoiled, LiquidStreamBoiled = self.__react__(LiquidStreamIn=LiquidStreamIn,
        #                                                     pressure_bara=np.compress(condition=(vapor_pressure_bara > pressure_bara), a=pressure_bara, axis=0),
        #                                                     heat_kW=np.compress(condition=(vapor_pressure_bara > pressure_bara), a=heat_kW, axis=0),
        #                                                     lr=lr)

        #LiquidStreamOut = deepcopy(LiquidStreamIn)
        #LiquidStreamOut = self.__insert_into_LiquidStream__(LiquidStreamIn=LiquidStreamBoiled, condition=(vapor_pressure_bara > pressure_bara), LiquidStreamOut=LiquidStreamOut)
        #LiquidStreamOut = self.__insert_into_LiquidStream__(LiquidStreamIn=LiquidStreamHeated, condition=(vapor_pressure_bara <= pressure_bara), LiquidStreamOut=LiquidStreamOut)

        GasStreamBoiled, LiquidStreamOut = self.__react__(LiquidStreamIn=LiquidStreamIn,
                                                             pressure_bara=pressure_bara,
                                                             heat_kW=heat_kW,
                                                             lr=lr)

        return GasStreamBoiled, LiquidStreamOut

    def __react__(self, LiquidStreamIn, pressure_bara, heat_kW, lr):

        GasStreamIn = self.__get_Flash_GasStream__(LiquidStreamIn)
        GasStreamIn.set_gas_flow_kmol_h(value=10**(-12) * LiquidStreamIn.get_solution_flow_kmol_h())
        GasStreamIn.set_gas_pressure_bara(value=pressure_bara)
        GasStreamIn.temp_K = LiquidStreamIn.temp_K

        self.heat_kW = heat_kW
        num_of_samples = GasStreamIn.temp_K.shape[0]

        if self.firstscan:
            self.matrix = self.__get_the_matrix__(GasStreamIn=GasStreamIn,
                                                  LiquidStreamIn=LiquidStreamIn,
                                                  liq_rxn_insta=True,
                                                  vapor_pressure=True,
                                                  gas_rxn_insta=True)
            self.firstscan = False

        LiquidStreamOut = deepcopy(LiquidStreamIn)
        for id in LiquidStreamOut.specie.keys():
            LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.maximum(LiquidStreamOut.get_specie_mass_fraction(id=id), 10 ** (-18)))

        GasStreamOut = deepcopy(GasStreamIn)
        for id in GasStreamOut.specie.keys():
            GasStreamOut.set_specie_molar_fraction(id=id, value=np.maximum(GasStreamOut.get_specie_molar_fraction(id=id), 10**(-18)))

        # Inlet
        m_in_tot = LiquidStreamOut.get_solution_flow_kg_h() + GasStreamOut.get_gas_flow_kg_h()
        m_in_liq_tot = LiquidStreamOut.get_solution_flow_kg_h()
        m_in_liq = LiquidStreamOut.get_solution_flow_kg_h()[:, None] * LiquidStreamOut.__mass_fractions_dic2vec__()
        m_in_gas = GasStreamOut.get_gas_flow_kg_h()[:, None] * GasStreamOut.__molefrac2massfrac__(y=GasStreamOut.__molar_fractions_dic2vec__())
        n_in_gas_tot = GasStreamOut.get_gas_flow_kmol_h()
        cp_in_liq = LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        #cp_in_gas = GasStreamOut.get_gas_heat_capacity_kJ_kmolK()
        T_in_liq = LiquidStreamOut.get_solution_temp_K()
        T_in_gas = GasStreamOut.get_gas_temp_K()

        # Initial Guess
        T = T_in_liq
        #T = T_in_liq + self.heat_kW/ (cp_in_liq * m_in_liq_tot/3600)

        # Initiate Mass Fracs...
        w_liq = LiquidStreamOut.__mass_fractions_dic2vec__()
        w_gas = GasStreamOut.__molefrac2massfrac__(y=GasStreamOut.__molar_fractions_dic2vec__())

        # Preparations...
        converged = False
        epoch = 0
        error = []
        iterations = np.zeros(shape=(num_of_samples,))

        # Iterate Until Convergence
        while converged == False:

            # Mass Flow
            m_liq_tot = LiquidStreamOut.get_solution_flow_kg_h()
            m_gas_tot = GasStreamOut.get_gas_flow_kg_h()
            m_liq = m_liq_tot[:, None] * w_liq
            m_gas = m_gas_tot[:, None] * w_gas

            # Equilibrium Constants
            Kw = self.__get_LiquidEquilibrium_Kw__(LiquidStreamOut)
            Kpw = self.__get_VaporLiquidEquilibrium_Kpw__(LiquidStreamOut)
            Kyw = self.__get_VaporLiquidEquilibrium_Kyw__(GasStreamOut ,LiquidStreamOut)
            Ky = self.__get_GasEquilibrium_Kp__(GasStreamOut)

            # Equilibrium Constants; Temperature Gradients
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(LiquidStreamOut, Kw)
            dKywdT = self.__get_VaporLiquidEquilibrium_dKywdT__(GasStreamOut, LiquidStreamOut, Kyw)
            dKydT = self.__get_GasEquilibrium_dKydT__(GasStreamOut, Ky)

            # Heat of Reactions
            dh_liq = self.__get_LiquidEquilibrium_exothermic_heat_kJ_kmol__(LiquidStreamOut, Kw, dKwdT)
            dh_vap = self.__get_VaporLiquidEquilibrium_exothermic_heat_kJ_kmol__(LiquidStreamOut, Kyw, dKywdT)
            dh_gas = self.__get_GasEquilibrium_exothermic_heat_kJ_kmol__(GasStreamOut, Ky, dKydT)
            dh = np.concatenate((dh_liq, dh_vap, dh_gas), axis=1)

            # Objective Functions; Equilibrium Constraints
            f_liq = self.__get_LiquidEquilibrium_f_rxn_insta_log__(LiquidStreamOut, Kw)
            f_vap = self.__get_Flash_f_vap_log__(LiquidStreamOut, GasStreamOut, Kpw)
            f_gas = self.__get_GasEquilibrium_f_rxn_insta_log__(GasStreamOut, Ky)

            # Objective Functions; Energy Balance
            dm_liq = m_liq - m_in_liq
            dm_gas = m_gas - m_in_gas
            dm = np.concatenate((dm_liq, dm_gas), axis=1)
            f_energy = (m_in_liq_tot / 3600) * cp_in_liq * (T - T_in_liq)
            #f_energy = f_energy + (n_in_gas_tot / 3600) * cp_in_gas * (T - T_in_gas)
            f_energy = f_energy - np.einsum("sm,sm->s", np.einsum("sr,rm->sm", dh, self.matrix["R+"]), dm / 3600)
            f_energy = f_energy - self.heat_kW

            # Objective Functions
            f = np.concatenate((f_liq, f_vap, f_gas, f_energy[:, None]), axis=1)

            # Partial Derivatives
            dfdw_liq = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(LiquidStreamOut)
            dfdw_vap = self.__get_Flash_dfdw_vap_log__(LiquidStreamOut)
            dfdy_vap = self.__get_Flash_dfdy_vap_loq__(GasStreamOut, LiquidStreamOut)
            dfdy_gas = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(GasStreamOut)

            # Partial derivatives; Mass Flows w.r.t Reaction Extent (z)
            dmdz = np.einsum("s,wr->swr", m_in_tot, self.matrix["R"])
            dmdz_liq = dmdz[:, 0:LiquidStreamOut.num_of_species:1,:]
            dmdz_gas = dmdz[:, LiquidStreamOut.num_of_species::,:]
            dmdz_liq_tot = np.sum(dmdz_liq, axis=1, keepdims=False)
            dmdz_gas_tot = np.sum(dmdz_gas, axis=1, keepdims=False)

            # Partial Derivatives; Quotient Rule
            dwdz_liq = np.einsum("s,smz->smz", m_liq_tot ** (-2), np.einsum("s,smz->smz", m_liq_tot, dmdz_liq) - np.einsum("sm,sz->smz", m_liq, dmdz_liq_tot))
            dwdz_gas = np.einsum("s,smz->smz", m_gas_tot ** (-2), np.einsum("s,smz->smz", m_gas_tot, dmdz_gas) - np.einsum("sm,sz->smz", m_gas, dmdz_gas_tot))

            # Partial Derivative; Gas Phase - Molar Fraction vs. Mass Fraction
            dydw_gas = GasStreamOut.__dydw__(w_gas)

            # Partial Derivatives; Objectives w.r.t. Reaction Extent
            dfdz_liq = np.einsum("slw,swz->slz", dfdw_liq, dwdz_liq)
            dfdz_vap = np.einsum("slw,swz->slz", dfdw_vap, dwdz_liq) + np.einsum("slw,swz->slz", np.einsum("sly,syw->slw", dfdy_vap, dydw_gas), dwdz_gas)
            dfdz_gas = np.einsum("slw,swz->slz", np.einsum("sly,syw->slw", dfdy_gas, dydw_gas), dwdz_gas)
            dfdz_energy = - np.einsum("s,sr->sr", m_in_tot / 3600, dh)[:, None, :]

            # Partial Derivatives; Objectives w.r.t. Temperature
            dfdT_liq = self.__get_LiquidEquilibrium_dfdT_rxn_insta_log__(Kw, dKwdT)
            dfdT_vap = self.__get_Flash_dfdT_vap_log__(Kyw, dKywdT)
            dfdT_gas = self.__get_GasEquilibrium_dfdT_rxn_insta_log__(Ky, dKydT)
            dfdT_energy = (m_in_liq_tot / 3600) * cp_in_liq

            # Partial Derivatives; Finalizing where X = (z,T) = (R*dm,T)
            dfdz = np.concatenate((dfdz_liq, dfdz_vap, dfdz_gas, dfdz_energy), axis=1)
            dfdT = np.concatenate((dfdT_liq, dfdT_vap, dfdT_gas, dfdT_energy[:, None, None]), axis=1)
            dfdX = np.concatenate((dfdz, dfdT), axis=2)

            # Damped Newton's Method
            dX_newton = - np.linalg.solve(dfdX, f)
            dz_newton = dX_newton[:, 0:dX_newton.shape[1] - 1:]
            dT_newton = dX_newton[:, dX_newton.shape[1] - 1]
            dm_newton = m_in_tot[:,None] * np.einsum("ij,sj->si", self.matrix["R"], dz_newton)
            dm = lr * dm_newton
            dT = lr * dT_newton

            # Backtrack to Ensure only Positive Concentrations
            m = np.concatenate((m_liq, m_gas), axis=1)
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * m / dm, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (m > 0) + 1.0 * (m <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            tau = np.broadcast_to(tau, shape=(num_of_samples, GasStreamOut.num_of_species + LiquidStreamOut.num_of_species)).copy()
            dm = dm * tau
            dT = dT * tau[:,0]

            # Mass Flow of individual Species; Split into Liquid and Gas
            dm_liq = dm[:, 0:LiquidStreamOut.num_of_species:1]
            dm_gas = dm[:, LiquidStreamOut.num_of_species::]

            # Mass Flow Total
            dm_liq_tot = np.sum(dm_liq, axis=1, keepdims=False)
            dm_gas_tot = np.sum(dm_gas, axis=1, keepdims=False)

            # Molar Flow Total
            dn_gas_tot = 0
            for i, id in enumerate(GasStreamOut.specie.keys()):
                dn_gas_tot = dn_gas_tot + dm_gas[:, i] / GasStreamOut.get_specie_molar_mass_kg_kmol(id)

            # Update Mass Fractions and Temp
            w_liq = (m_liq + dm_liq) / (m_liq_tot[:, None] + dm_liq_tot[:, None])
            w_gas = (m_gas + dm_gas) / (m_gas_tot[:, None] + dm_gas_tot[:, None])
            T = T + dT

            # Update Gas Stream
            GasStreamOut.__molar_fractions_vec2dic__(y=GasStreamOut.__massfrac2molefrac__(w_gas))
            GasStreamOut.flow_kmol_h = GasStreamOut.flow_kmol_h + dn_gas_tot
            GasStreamOut.temp_K = T
            GasStreamOut.normalize_molar_fractions()

            # Update Liquid Stream
            LiquidStreamOut.__mass_fractions_vec2dic__(w=w_liq)
            LiquidStreamOut.flow_kg_h = LiquidStreamOut.flow_kg_h + dm_liq_tot
            LiquidStreamOut.temp_K = T
            LiquidStreamOut.normalize_mass_fractions()

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dm_newton) < 0.005 * np.abs(m), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            sample_converged = sample_converged * (np.abs(dT_newton) < 10 ** (-2))
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            epoch = epoch + 1

        #print("Epochs: \t" + str(epoch))
        return GasStreamOut, LiquidStreamOut


# ---------------------------------------------------------------------------------------


class VaporLiquidEquilibrium_Isothermal(_Stochiometry, _Serializer, _LiquidEquilibrium, _GasEquilibrium, _VaporLiquidEquilibrium, _StreamFunctions):

    def __init__(self):
        self.firstscan = True

    def react(self, GasStreamIn, LiquidStreamIn, temp_K, lr=0.75):

        GasStreamIn.temp_K = temp_K
        LiquidStreamIn.temp_K = temp_K

        num_of_samples = GasStreamIn.temp_K.shape[0]

        if self.firstscan:
            self.matrix = self.__get_the_matrix__(GasStreamIn=GasStreamIn,
                                                  LiquidStreamIn=LiquidStreamIn,
                                                  liq_rxn_insta=True,
                                                  vapor_pressure=True,
                                                  gas_rxn_insta=True)
            self.firstscan = False

        LiquidStreamOut = deepcopy(LiquidStreamIn)
        for id in LiquidStreamOut.specie.keys():
            LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.maximum(LiquidStreamOut.get_specie_mass_fraction(id=id), 10 ** (-18)))

        GasStreamOut = deepcopy(GasStreamIn)
        for id in GasStreamOut.specie.keys():
            GasStreamOut.set_specie_molar_fraction(id=id, value=np.maximum(GasStreamOut.get_specie_molar_fraction(id=id), 10**(-18)))

        # Inlet
        m_in_tot = LiquidStreamOut.get_solution_flow_kg_h() + GasStreamOut.get_gas_flow_kg_h()
        m_in_liq_tot = LiquidStreamOut.get_solution_flow_kg_h()
        m_in_liq = LiquidStreamOut.get_solution_flow_kg_h()[:, None] * LiquidStreamOut.__mass_fractions_dic2vec__()
        m_in_gas = GasStreamOut.get_gas_flow_kg_h()[:, None] * GasStreamOut.__molefrac2massfrac__(y=GasStreamOut.__molar_fractions_dic2vec__())
        n_in_gas_tot = GasStreamOut.get_gas_flow_kmol_h()
        cp_in_liq = LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        cp_in_gas = GasStreamOut.get_gas_heat_capacity_kJ_kmolK()

        # Initiate Mass Fracs...
        w_liq = LiquidStreamOut.__mass_fractions_dic2vec__()
        w_gas = GasStreamOut.__molefrac2massfrac__(y=GasStreamOut.__molar_fractions_dic2vec__())

        # Preparations...
        converged = False
        epoch = 0
        error = []
        iterations = np.zeros(shape=(num_of_samples,))

        # Iterate Until Convergence
        while converged == False:

            # Mass Flow
            m_liq_tot = LiquidStreamOut.get_solution_flow_kg_h()
            m_gas_tot = GasStreamOut.get_gas_flow_kg_h()
            m_liq = m_liq_tot[:, None] * w_liq
            m_gas = m_gas_tot[:, None] * w_gas

            # Equilibrium Constants
            Kw = self.__get_LiquidEquilibrium_Kw__(LiquidStreamOut)
            Kyw = self.__get_VaporLiquidEquilibrium_Kyw__(GasStreamOut ,LiquidStreamOut)
            Ky = self.__get_GasEquilibrium_Kp__(GasStreamOut)

            # Objective Functions; Equilibrium Constraints
            f_liq = self.__get_LiquidEquilibrium_f_rxn_insta_log__(LiquidStreamOut, Kw)
            f_vap = self.__get_VaporLiquidEquilibrium_f_log__(GasStreamOut, LiquidStreamOut, Kyw)
            f_gas = self.__get_GasEquilibrium_f_rxn_insta_log__(GasStreamOut, Ky)

            # Objective Functions
            f = np.concatenate((f_liq, f_vap, f_gas), axis=1)

            # Partial Derivatives
            dfdw_liq = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(LiquidStreamOut)
            dfdw_vap = self.__get_VaporLiquidEquilibrium_dfdw_log__(LiquidStreamOut)
            dfdy_vap = self.__get_VaporLiquidEquilibrium_dfdy_log__(GasStreamOut, LiquidStreamOut)
            dfdy_gas = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(GasStreamOut)

            # Partial derivatives; Mass Flows w.r.t Reaction Extent (z)
            dmdz = np.einsum("s,wr->swr", m_in_tot, self.matrix["R"])
            dmdz_liq = dmdz[:, 0:LiquidStreamOut.num_of_species:1,:]
            dmdz_gas = dmdz[:, LiquidStreamOut.num_of_species::,:]
            dmdz_liq_tot = np.sum(dmdz_liq, axis=1, keepdims=False)
            dmdz_gas_tot = np.sum(dmdz_gas, axis=1, keepdims=False)

            # Partial Derivatives; Quotient Rule
            dwdz_liq = np.einsum("s,smz->smz", m_liq_tot ** (-2), np.einsum("s,smz->smz", m_liq_tot, dmdz_liq) - np.einsum("sm,sz->smz", m_liq, dmdz_liq_tot))
            dwdz_gas = np.einsum("s,smz->smz", m_gas_tot ** (-2), np.einsum("s,smz->smz", m_gas_tot, dmdz_gas) - np.einsum("sm,sz->smz", m_gas, dmdz_gas_tot))

            # Partial Derivative; Gas Phase - Molar Fraction vs. Mass Fraction
            dydw_gas = GasStreamOut.__dydw__(w_gas)

            # Partial Derivatives; Objectives w.r.t. Reaction Extent
            dfdz_liq = np.einsum("slw,swz->slz", dfdw_liq, dwdz_liq)
            dfdz_vap = np.einsum("slw,swz->slz", dfdw_vap, dwdz_liq) + np.einsum("slw,swz->slz", np.einsum("sly,syw->slw", dfdy_vap, dydw_gas), dwdz_gas)
            dfdz_gas = np.einsum("slw,swz->slz", np.einsum("sly,syw->slw", dfdy_gas, dydw_gas), dwdz_gas)
            dfdz = np.concatenate((dfdz_liq, dfdz_vap, dfdz_gas), axis=1)

            # Damped Newton's Method
            dz_newton = - np.linalg.solve(dfdz, f)
            dm_newton = m_in_tot[:,None] * np.einsum("ij,sj->si", self.matrix["R"], dz_newton)
            dm = lr * dm_newton

            # Backtrack to Ensure only Positive Concentrations
            m = np.concatenate((m_liq, m_gas), axis=1)
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * m / dm, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (m > 0) + 1.0 * (m <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            tau = np.broadcast_to(tau, shape=(num_of_samples, GasStreamOut.num_of_species + LiquidStreamOut.num_of_species)).copy()
            dm = dm * tau

            # Mass Flow of individual Species; Split into Liquid and Gas
            dm_liq = dm[:, 0:LiquidStreamOut.num_of_species:1]
            dm_gas = dm[:, LiquidStreamOut.num_of_species::]

            # Mass Flow Total
            dm_liq_tot = np.sum(dm_liq, axis=1, keepdims=False)
            dm_gas_tot = np.sum(dm_gas, axis=1, keepdims=False)

            # Molar Flow Total
            dn_gas_tot = 0
            for i, id in enumerate(GasStreamOut.specie.keys()):
                dn_gas_tot = dn_gas_tot + dm_gas[:, i] / GasStreamOut.get_specie_molar_mass_kg_kmol(id)

            # Update Mass Fractions and Temp
            w_liq = (m_liq + dm_liq) / (m_liq_tot[:, None] + dm_liq_tot[:, None])
            w_gas = (m_gas + dm_gas) / (m_gas_tot[:, None] + dm_gas_tot[:, None])

            # Update Gas Stream
            GasStreamOut.__molar_fractions_vec2dic__(y=GasStreamOut.__massfrac2molefrac__(w_gas))
            GasStreamOut.flow_kmol_h = GasStreamOut.flow_kmol_h + dn_gas_tot
            GasStreamOut.normalize_molar_fractions()

            # Update Liquid Stream
            LiquidStreamOut.__mass_fractions_vec2dic__(w=w_liq)
            LiquidStreamOut.flow_kg_h = LiquidStreamOut.flow_kg_h + dm_liq_tot
            LiquidStreamOut.normalize_mass_fractions()

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dm_newton) < 0.005 * np.abs(m), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            epoch = epoch + 1

        #print("Epochs: \t" + str(epoch))
        return GasStreamOut, LiquidStreamOut


class VaporLiquidEquilibrium_Adiabatic(_Stochiometry, _Serializer, _LiquidEquilibrium, _GasEquilibrium, _VaporLiquidEquilibrium, _StreamFunctions):

    def __init__(self):
        self.firstscan = True

    def react(self, GasStreamIn, LiquidStreamIn, lr=0.75):

        num_of_samples = GasStreamIn.temp_K.shape[0]

        if self.firstscan:
            self.matrix = self.__get_the_matrix__(GasStreamIn=GasStreamIn,
                                                  LiquidStreamIn=LiquidStreamIn,
                                                  liq_rxn_insta=True,
                                                  vapor_pressure=True,
                                                  gas_rxn_insta=True)
            self.firstscan = False

        LiquidStreamOut = deepcopy(LiquidStreamIn)
        for id in LiquidStreamOut.specie.keys():
            LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.maximum(LiquidStreamOut.get_specie_mass_fraction(id=id), 10 ** (-18)))

        GasStreamOut = deepcopy(GasStreamIn)
        for id in GasStreamOut.specie.keys():
            GasStreamOut.set_specie_molar_fraction(id=id, value=np.maximum(GasStreamOut.get_specie_molar_fraction(id=id), 10**(-18)))

        # Inlet
        m_in_tot = LiquidStreamOut.get_solution_flow_kg_h() + GasStreamOut.get_gas_flow_kg_h()
        m_in_liq_tot = LiquidStreamOut.get_solution_flow_kg_h()
        m_in_liq = LiquidStreamOut.get_solution_flow_kg_h()[:, None] * LiquidStreamOut.__mass_fractions_dic2vec__()
        m_in_gas = GasStreamOut.get_gas_flow_kg_h()[:, None] * GasStreamOut.__molefrac2massfrac__(y=GasStreamOut.__molar_fractions_dic2vec__())
        n_in_gas_tot = GasStreamOut.get_gas_flow_kmol_h()
        cp_in_liq = LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        cp_in_gas = GasStreamOut.get_gas_heat_capacity_kJ_kmolK()
        T_in_liq = LiquidStreamOut.get_solution_temp_K()
        T_in_gas = GasStreamOut.get_gas_temp_K()

        # Equilibrium Temp before taking reactions into account
        T = (n_in_gas_tot * cp_in_gas * T_in_gas + m_in_liq_tot * cp_in_liq * T_in_liq) / (n_in_gas_tot * cp_in_gas + m_in_liq_tot * cp_in_liq)
        LiquidStreamOut.temp_K = T
        GasStreamOut.temp_K = T
        T_in_liq = T
        T_in_gas = T

        # Initiate Mass Fracs...
        w_liq = LiquidStreamOut.__mass_fractions_dic2vec__()
        w_gas = GasStreamOut.__molefrac2massfrac__(y=GasStreamOut.__molar_fractions_dic2vec__())

        # Preparations...
        converged = False
        epoch = 0
        error = []
        iterations = np.zeros(shape=(num_of_samples,))

        # Iterate Until Convergence
        while converged == False:

            # Mass Flow
            m_liq_tot = LiquidStreamOut.get_solution_flow_kg_h()
            m_gas_tot = GasStreamOut.get_gas_flow_kg_h()
            m_liq = m_liq_tot[:, None] * w_liq
            m_gas = m_gas_tot[:, None] * w_gas

            # Equilibrium Constants
            Kw = self.__get_LiquidEquilibrium_Kw__(LiquidStreamOut)
            Kyw = self.__get_VaporLiquidEquilibrium_Kyw__(GasStreamOut ,LiquidStreamOut)
            Ky = self.__get_GasEquilibrium_Kp__(GasStreamOut)

            # Equilibrium Constants; Temperature Gradients
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(LiquidStreamOut, Kw)
            dKywdT = self.__get_VaporLiquidEquilibrium_dKywdT__(GasStreamOut, LiquidStreamOut, Kyw)
            dKydT = self.__get_GasEquilibrium_dKydT__(GasStreamOut, Ky)

            # Heat of Reactions
            dh_liq = self.__get_LiquidEquilibrium_exothermic_heat_kJ_kmol__(LiquidStreamOut, Kw, dKwdT)
            dh_vap = self.__get_VaporLiquidEquilibrium_exothermic_heat_kJ_kmol__(LiquidStreamOut, Kyw, dKywdT)
            dh_gas = self.__get_GasEquilibrium_exothermic_heat_kJ_kmol__(GasStreamOut, Ky, dKydT)
            dh = np.concatenate((dh_liq, dh_vap, dh_gas), axis=1)

            # Objective Functions; Equilibrium Constraints
            f_liq = self.__get_LiquidEquilibrium_f_rxn_insta_log__(LiquidStreamOut, Kw)
            f_vap = self.__get_VaporLiquidEquilibrium_f_log__(GasStreamOut, LiquidStreamOut, Kyw)
            f_gas = self.__get_GasEquilibrium_f_rxn_insta_log__(GasStreamOut, Ky)

            # Objective Functions; Energy Balance
            dm_liq = m_liq - m_in_liq
            dm_gas = m_gas - m_in_gas
            dm = np.concatenate((dm_liq, dm_gas), axis=1)
            f_energy = (m_in_liq_tot / 3600) * cp_in_liq * (T - T_in_liq)
            f_energy = f_energy + (n_in_gas_tot / 3600) * cp_in_gas * (T - T_in_gas)
            f_energy = f_energy - np.einsum("sm,sm->s", np.einsum("sr,rm->sm", dh, self.matrix["R+"]), dm / 3600)

            # Objective Functions
            f = np.concatenate((f_liq, f_vap, f_gas, f_energy[:, None]), axis=1)

            # Partial Derivatives
            dfdw_liq = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(LiquidStreamOut)
            dfdw_vap = self.__get_VaporLiquidEquilibrium_dfdw_log__(LiquidStreamOut)
            dfdy_vap = self.__get_VaporLiquidEquilibrium_dfdy_log__(GasStreamOut, LiquidStreamOut)
            dfdy_gas = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(GasStreamOut)

            # Partial derivatives; Mass Flows w.r.t Reaction Extent (z)
            dmdz = np.einsum("s,wr->swr", m_in_tot, self.matrix["R"])
            dmdz_liq = dmdz[:, 0:LiquidStreamOut.num_of_species:1,:]
            dmdz_gas = dmdz[:, LiquidStreamOut.num_of_species::,:]
            dmdz_liq_tot = np.sum(dmdz_liq, axis=1, keepdims=False)
            dmdz_gas_tot = np.sum(dmdz_gas, axis=1, keepdims=False)

            # Partial Derivatives; Quotient Rule
            dwdz_liq = np.einsum("s,smz->smz", m_liq_tot ** (-2), np.einsum("s,smz->smz", m_liq_tot, dmdz_liq) - np.einsum("sm,sz->smz", m_liq, dmdz_liq_tot))
            dwdz_gas = np.einsum("s,smz->smz", m_gas_tot ** (-2), np.einsum("s,smz->smz", m_gas_tot, dmdz_gas) - np.einsum("sm,sz->smz", m_gas, dmdz_gas_tot))

            # Partial Derivative; Gas Phase - Molar Fraction vs. Mass Fraction
            dydw_gas = GasStreamOut.__dydw__(w_gas)

            # Partial Derivatives; Objectives w.r.t. Reaction Extent
            dfdz_liq = np.einsum("slw,swz->slz", dfdw_liq, dwdz_liq)
            dfdz_vap = np.einsum("slw,swz->slz", dfdw_vap, dwdz_liq) + np.einsum("slw,swz->slz", np.einsum("sly,syw->slw", dfdy_vap, dydw_gas), dwdz_gas)
            dfdz_gas = np.einsum("slw,swz->slz", np.einsum("sly,syw->slw", dfdy_gas, dydw_gas), dwdz_gas)
            dfdz_energy = - np.einsum("s,sr->sr", m_in_tot / 3600, dh)[:, None, :]

            # Partial Derivatives; Objectives w.r.t. Temperature
            dfdT_liq = self.__get_LiquidEquilibrium_dfdT_rxn_insta_log__(Kw, dKwdT)
            dfdT_vap = self.__get_VaporLiquidEquilibrium_dfdT_log__(Kyw, dKywdT)
            dfdT_gas = self.__get_GasEquilibrium_dfdT_rxn_insta_log__(Ky, dKydT)
            dfdT_energy = (m_in_liq_tot / 3600) * cp_in_liq + (n_in_gas_tot / 3600) * cp_in_gas

            # Partial Derivatives; Finalizing where X = (z,T) = (R*dm,T)
            dfdz = np.concatenate((dfdz_liq, dfdz_vap, dfdz_gas, dfdz_energy), axis=1)
            dfdT = np.concatenate((dfdT_liq, dfdT_vap, dfdT_gas, dfdT_energy[:, None, None]), axis=1)
            dfdX = np.concatenate((dfdz, dfdT), axis=2)

            # Damped Newton's Method
            dX_newton = - np.linalg.solve(dfdX, f)
            dz_newton = dX_newton[:, 0:dX_newton.shape[1] - 1:]
            dT_newton = dX_newton[:, dX_newton.shape[1] - 1]
            dm_newton = m_in_tot[:,None] * np.einsum("ij,sj->si", self.matrix["R"], dz_newton)
            dm = lr * dm_newton
            dT = lr * dT_newton

            # Backtrack to Ensure only Positive Concentrations
            m = np.concatenate((m_liq, m_gas), axis=1)
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * m / dm, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (m > 0) + 1.0 * (m <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            tau = np.broadcast_to(tau, shape=(num_of_samples, GasStreamOut.num_of_species + LiquidStreamOut.num_of_species)).copy()
            dm = dm * tau
            dT = dT * tau[:,0]

            # Mass Flow of individual Species; Split into Liquid and Gas
            dm_liq = dm[:, 0:LiquidStreamOut.num_of_species:1]
            dm_gas = dm[:, LiquidStreamOut.num_of_species::]

            # Mass Flow Total
            dm_liq_tot = np.sum(dm_liq, axis=1, keepdims=False)
            dm_gas_tot = np.sum(dm_gas, axis=1, keepdims=False)

            # Molar Flow Total
            dn_gas_tot = 0
            for i, id in enumerate(GasStreamOut.specie.keys()):
                dn_gas_tot = dn_gas_tot + dm_gas[:, i] / GasStreamOut.get_specie_molar_mass_kg_kmol(id)

            # Update Mass Fractions and Temp
            w_liq = (m_liq + dm_liq) / (m_liq_tot[:, None] + dm_liq_tot[:, None])
            w_gas = (m_gas + dm_gas) / (m_gas_tot[:, None] + dm_gas_tot[:, None])
            T = T + dT

            # Update Gas Stream
            GasStreamOut.__molar_fractions_vec2dic__(y=GasStreamOut.__massfrac2molefrac__(w_gas))
            GasStreamOut.flow_kmol_h = GasStreamOut.flow_kmol_h + dn_gas_tot
            GasStreamOut.temp_K = T
            GasStreamOut.normalize_molar_fractions()

            # Update Liquid Stream
            LiquidStreamOut.__mass_fractions_vec2dic__(w=w_liq)
            LiquidStreamOut.flow_kg_h = LiquidStreamOut.flow_kg_h + dm_liq_tot
            LiquidStreamOut.temp_K = T
            LiquidStreamOut.normalize_mass_fractions()

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dm_newton) < 0.005 * np.abs(m), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            sample_converged = sample_converged * (np.abs(dT_newton) < 10 ** (-2))
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            epoch = epoch + 1

        #print("Epochs: \t" + str(epoch))
        return GasStreamOut, LiquidStreamOut


class VaporLiquidEquilibrium_EquilibriumStages(_Serializer):

    def __init__(self, num_of_stages):
        self.firstscan = True
        self.num_of_stages = num_of_stages
        self.GasStreams = [None] * (num_of_stages + 1)
        self.LiquidStreams = [None] * (num_of_stages + 1)
        self.Stage = VaporLiquidEquilibrium_Adiabatic()

    def react(self, GasStreamIn, LiquidStreamIn):

        if self.firstscan or True:
            for i in range(self.num_of_stages + 1):
                self.GasStreams[i] = deepcopy(GasStreamIn)
                self.LiquidStreams[i] = deepcopy(LiquidStreamIn)
        else:
            self.GasStreams[0] = deepcopy(GasStreamIn)
            self.LiquidStreams[self.num_of_stages] = deepcopy(LiquidStreamIn)

        w_old = self.LiquidStreams[0].__mass_fractions_dic2vec__()
        y_old = self.GasStreams[self.num_of_stages].__molar_fractions_dic2vec__()

        converged = False
        epochs = 0
        while converged == False:

            # Iterate LiquidStreams Downward
            for i in range(self.num_of_stages, 0, -1):
                _, self.LiquidStreams[i-1] = self.Stage.react(self.GasStreams[i-1], self.LiquidStreams[i])

            # Iterate GasStreams Upward
            for i in range(0, self.num_of_stages, 1):
                self.GasStreams[i+1], _ = self.Stage.react(self.GasStreams[i], self.LiquidStreams[i+1])

            # Check if Algorithm have Converged
            w = self.LiquidStreams[0].__mass_fractions_dic2vec__()
            y = self.GasStreams[self.num_of_stages].__molar_fractions_dic2vec__()
            dw = w - w_old
            dy = y - y_old
            liq_converged = np.array((np.abs(dw) < 0.01 * np.abs(w)) | (w < 10**(-17)), dtype=np.float32)
            liq_converged = np.min(liq_converged, axis=1, keepdims=False)
            gas_converged = np.array((np.abs(dy) < 0.01 * np.abs(y)) | (y < 10**(-17)), dtype=np.float32)
            gas_converged = np.min(gas_converged, axis=1, keepdims=False)
            sample_converged = np.minimum(liq_converged, gas_converged)
            converged = np.min(sample_converged, axis=0)
            w_old = 1.0 * w
            y_old = 1.0 * y
            epochs = epochs + 1

        #print(epochs)

        return self.GasStreams[self.num_of_stages], self.LiquidStreams[0]


# ---------------------------------------------------------------------------------------


class LiquidCSTR_Isothermal(_Stochiometry, _Serializer, _LiquidEquilibrium, _LiquidCSTR):

    def __init__(self):
        self.firstscan = True

    def react(self, LiquidStreamIn, volume_m3, lr=0.75):

        # Load Inlet Stream, and Initiate Outlet Stream
        self.LiquidStreamInit = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamInit.specie.keys():
            self.LiquidStreamInit.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamInit.specie[id]["Mass Fraction"], 10**(-18))
        self.LiquidStreamInit.normalize_mass_fractions()

        if self.firstscan:
            self.matrix = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True, liq_rxn_reversible=True)
            self.LiquidStreamOut = deepcopy(self.LiquidStreamInit)
            self.firstscan = False

        if self.LiquidStreamOut.temp_K.shape[0] != self.LiquidStreamInit.temp_K.shape[0]:
            self.LiquidStreamOut = deepcopy(self.LiquidStreamInit)

        if self.LiquidStreamOut.id != self.LiquidStreamInit.id:
            self.matrix = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True, liq_rxn_reversible=True)
            self.LiquidStreamOut = deepcopy(self.LiquidStreamInit)

        # Conserved Quantity
        b = np.einsum("cw,sw->sc", self.matrix["A"], self.LiquidStreamInit.__mass_fractions_dic2vec__())

        # Initial Guess or Mass Fractions
        w = self.LiquidStreamOut.__mass_fractions_dic2vec__()

        # Reactor Dimension
        self.volume_m3 = volume_m3
        phi = self.LiquidStreamInit.get_solution_flow_kg_h() / (self.volume_m3 * 3600)

        # Iterate Until Convergence
        converged = False
        epoch = 0
        iterations = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0],))

        while converged == False:

            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamOut)
            rates_kmol_m3s = self.__get_LiquidCSTR_rxn_reversible_rates_kmol_m3s__(self.LiquidStreamOut)

            f_rxn_insta = self.__get_LiquidEquilibrium_f_rxn_insta__(self.LiquidStreamOut, Kw)
            f_rxn_reversible = self.__get_LiquidCSTR_f_rxn_reversible__(self.LiquidStreamInit, self.LiquidStreamOut, self.matrix, rates_kmol_m3s, phi)
            f_mass_balance = self.__get_LiquidEquilibrium_f_mass_balance__(self.LiquidStreamOut, b, self.matrix)

            f = np.hstack((f_rxn_insta, f_rxn_reversible, f_mass_balance))

            dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta__(self.LiquidStreamOut, Kw, self.matrix)
            dfdw_rxn_reversible = self.__get_LiquidCSTR_dfdw_rxn_reversible__(self.LiquidStreamOut, self.matrix, rates_kmol_m3s, phi)
            dfdw_mass_balance = self.__get_LiquidEquilibrium_dfdw_mass_balance__(self.LiquidStreamOut, self.matrix)

            dfdw = np.concatenate((dfdw_rxn_insta, dfdw_rxn_reversible, dfdw_mass_balance), axis=1)

            dw_newton = - np.linalg.solve(dfdw, f)
            dw = lr * dw_newton

            w = w + dw
            w = np.maximum(w, 0)
            self.LiquidStreamOut.__mass_fractions_vec2dic__(w)

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        #print("Epochs: \t" + str(epoch))
        return self.LiquidStreamOut

    def get_heat_dissipation_kW(self):

        h_in_insta = self.LiquidStreamInit.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_in_reversible = self.LiquidStreamInit.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()

        h_out_insta = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out_reversible = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()

        h_in = np.concatenate((h_in_insta, h_in_reversible), axis=1)
        h_out = np.concatenate((h_out_insta, h_out_reversible), axis=1)
        h = (h_in + h_out) / 2

        cp_in = self.LiquidStreamInit.get_solution_heat_capacity_kJ_kgK()
        cp_out = self.LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        cp = (cp_in + cp_out) / 2

        dw = self.LiquidStreamOut.__mass_fractions_dic2vec__() - self.LiquidStreamInit.__mass_fractions_dic2vec__()
        dT = self.LiquidStreamOut.temp_K - self.LiquidStreamInit.temp_K

        m = self.LiquidStreamInit.get_solution_flow_kg_h() / 3600

        heat_dissipation_kW = m * np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw)) - m * cp * dT
        return heat_dissipation_kW

    def get_adiabatic_temperature_increase_K(self):

        h_in_insta = self.LiquidStreamInit.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_in_reversible = self.LiquidStreamInit.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()

        h_out_insta = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out_reversible = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()

        h_in = np.concatenate((h_in_insta, h_in_reversible), axis=1)
        h_out = np.concatenate((h_out_insta, h_out_reversible), axis=1)
        h = (h_in + h_out) / 2

        cp_in = self.LiquidStreamInit.get_solution_heat_capacity_kJ_kgK()
        cp_out = self.LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        cp = (cp_in + cp_out) / 2

        dw = self.LiquidStreamOut.__mass_fractions_dic2vec__() - self.LiquidStreamInit.__mass_fractions_dic2vec__()
        dT = np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw)) / cp

        return dT


class LiquidCSTR_Adiabatic(_Stochiometry, _Serializer, _LiquidEquilibrium, _LiquidCSTR):

    def __init__(self):
        self.firstscan = True

    def react(self, LiquidStreamIn, volume_m3, heat_kW, lr=0.75):

        # Load Inlet Stream
        self.LiquidStreamInit = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamInit.specie.keys():
            self.LiquidStreamInit.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamInit.specie[id]["Mass Fraction"], 10**(-18))
        self.LiquidStreamInit.normalize_mass_fractions()

        if self.firstscan:
            self.matrix = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True, liq_rxn_reversible=True)
            self.firstscan = False

        # Conserved Quantity
        b = np.einsum("cw,sw->sc", self.matrix["A"], self.LiquidStreamInit.__mass_fractions_dic2vec__())

        # Initial Guess
        eq_iso = LiquidCSTR_Isothermal()
        self.LiquidStreamOut = deepcopy(self.LiquidStreamInit)
        self.LiquidStreamOut.temp_K = self.LiquidStreamOut.temp_K + 3600 * heat_kW / (self.LiquidStreamOut.get_solution_flow_kg_h() * self.LiquidStreamOut.get_solution_heat_capacity_kJ_kgK())
        self.LiquidStreamOut = eq_iso.react(self.LiquidStreamOut, volume_m3=volume_m3)
        self.LiquidStreamOut.temp_K = self.LiquidStreamOut.temp_K + eq_iso.get_adiabatic_temperature_increase_K()
        w = self.LiquidStreamOut.__mass_fractions_dic2vec__()

        # Reactor Dimension
        self.volume_m3 = volume_m3
        self.heat_kW = heat_kW
        phi = self.LiquidStreamInit.get_solution_flow_kg_h() / (self.volume_m3 * 3600)

        # Iterate Until Convergence
        converged = False
        epoch = 0
        iterations = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0],))

        while converged == False:

            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamOut)
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStreamOut, Kw)

            rates_kmol_m3s = self.__get_LiquidCSTR_rxn_reversible_rates_kmol_m3s__(self.LiquidStreamOut)

            f_rxn_insta = self.__get_LiquidEquilibrium_f_rxn_insta__(self.LiquidStreamOut, Kw)
            f_rxn_reversible = self.__get_LiquidCSTR_f_rxn_reversible__(self.LiquidStreamInit, self.LiquidStreamOut, self.matrix, rates_kmol_m3s, phi)
            f_mass_balance = self.__get_LiquidEquilibrium_f_mass_balance__(self.LiquidStreamOut, b, self.matrix)
            f_energy_balance = self.__get_LiquidCSTR_f_energy_balance__(self.LiquidStreamInit, self.LiquidStreamOut, self.matrix, self.heat_kW)
            f = np.hstack((f_rxn_insta, f_rxn_reversible, f_mass_balance, f_energy_balance))

            dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta__(self.LiquidStreamOut, Kw, self.matrix)
            dfdw_rxn_reversible = self.__get_LiquidCSTR_dfdw_rxn_reversible__(self.LiquidStreamOut, self.matrix, rates_kmol_m3s, phi)
            dfdw_mass_balance = self.__get_LiquidEquilibrium_dfdw_mass_balance__(self.LiquidStreamOut, self.matrix)
            dfdw_energy_balance = self.__get_LiquidCSTR_dfdw_energy_balance__(self.LiquidStreamInit, self.LiquidStreamOut, self.matrix)
            dfdw = np.concatenate((dfdw_rxn_insta, dfdw_rxn_reversible, dfdw_mass_balance, dfdw_energy_balance), axis=1)

            dfdT_rxn_insta = self.__get_LiquidEquilibrium_dfdT_rxn_insta__(self.LiquidStreamOut, dKwdT)
            dfdT_rxn_reversible = self.__get_LiquidCSTR_dfdT_rxn_reversible__(self.LiquidStreamOut, self.matrix, rates_kmol_m3s)
            dfdT_mass_balance = self.__get_LiquidEquilibrium_dfdT_mass_balance__(self.LiquidStreamOut, self.matrix)
            dfdT_energy_balance = self.__get_LiquidCSTR_dfdT_energy_balance__(self.LiquidStreamOut)
            dfdT = np.concatenate((dfdT_rxn_insta, dfdT_rxn_reversible, dfdT_mass_balance, dfdT_energy_balance), axis=1)

            dfdX = np.concatenate((dfdw, dfdT), axis=2)

            dX_newton = - np.linalg.solve(dfdX, f)
            dw_newton = dX_newton[:, :self.LiquidStreamOut.num_of_species:]
            dT_newton = dX_newton[:, self.LiquidStreamOut.num_of_species]

            dw = lr * dw_newton
            dT = lr * dT_newton

            w = w + dw
            w = np.maximum(w, 0)
            self.LiquidStreamOut.__mass_fractions_vec2dic__(w)
            self.LiquidStreamOut.temp_K = self.LiquidStreamOut.temp_K + dT

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        print("Epochs: \t" + str(epoch))
        return self.LiquidStreamOut


class LiquidCSTR_QPFlash(_Stochiometry, _Serializer, _LiquidEquilibrium, _LiquidCSTR, _Flash):

    def __init__(self):
        self.firstscan = True

    def react(self, LiquidStreamIn, volume_m3, heat_kW, pressure_bara, lr=0.75):

        # Load Inlet Stream
        self.LiquidStreamInit = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamInit.specie.keys():
            self.LiquidStreamInit.specie[id]["Mass Fraction"] = np.maximum(
                self.LiquidStreamInit.specie[id]["Mass Fraction"], 10 ** (-18))
        self.LiquidStreamInit.normalize_mass_fractions()

        # Generate "GasStream"
        self.GasStreamOut = self.__get_Flash_GasStream__(self.LiquidStreamInit)

        if self.firstscan:
            self.matrix = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, GasStreamIn=self.GasStreamOut, liq_rxn_insta=True, liq_rxn_reversible=True, vapor_pressure=True)
            self.firstscan = False

        # Heat Up Liquid without Flash
        self.LiquidStreamHeated = LiquidCSTR_Adiabatic().react(LiquidStreamIn=self.LiquidStreamInit, volume_m3=volume_m3, heat_kW=heat_kW, lr=lr)

        # Continue with Samples where Flashing Occur
        vapor_pressure_bara = 0
        for id in self.LiquidStreamHeated.vapor_pressure_bara.keys():
            gas_id = list(self.LiquidStreamHeated.vapor_pressure_bara[id]["Stoch Gas"].keys())[0]
            vapor_pressure_bara = vapor_pressure_bara + self.LiquidStreamHeated.get_specie_vapor_pressure_bara(gas_id)
        self.LiquidStreamBoiled = deepcopy(self.LiquidStreamHeated)
        self.LiquidStreamBoiled.__compress__(condition=(vapor_pressure_bara > pressure_bara))
        self.LiquidStreamHeated.__compress__(condition=(vapor_pressure_bara <= pressure_bara))


        # -------------------------------------------------------------------------------------------

        # Iterate Until Convergence
        converged = False
        epoch = 0
        iterations = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0],))


        b = None

        while converged == False:

            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamBoiled)
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStreamBoiled, Kw)
            rates_kmol_m3s = self.__get_LiquidCSTR_rxn_reversible_rates_kmol_m3s__(self.LiquidStreamBoiled)

            f_rxn_insta = self.__get_LiquidEquilibrium_f_rxn_insta__(self.LiquidStreamBoiled, Kw)
            f_rxn_reversible = self.__get_LiquidCSTR_f_rxn_reversible__(self.LiquidStreamInit, self.LiquidStreamBoiled, self.matrix, rates_kmol_m3s, phi)


            f_mass_balance = self.__get_FlashEvaporation_f_mass_balance__(self.LiquidStreamBoiled, self.GasStreamOut, self.matrix, b)



            f_energy_balance = self.__get_LiquidCSTR_f_energy_balance__(self.LiquidStreamInit, self.LiquidStreamOut, self.matrix, self.heat_kW)




            f = np.hstack((f_rxn_insta, f_rxn_reversible, f_mass_balance, f_energy_balance))

            dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta__(self.LiquidStreamOut, Kw, self.matrix)
            dfdw_rxn_reversible = self.__get_LiquidCSTR_dfdw_rxn_reversible__(self.LiquidStreamOut, self.matrix,
                                                                              rates_kmol_m3s, phi)
            dfdw_mass_balance = self.__get_LiquidEquilibrium_dfdw_mass_balance__(self.LiquidStreamOut, self.matrix)
            dfdw_energy_balance = self.__get_LiquidCSTR_dfdw_energy_balance__(self.LiquidStreamInit,
                                                                              self.LiquidStreamOut, self.matrix)
            dfdw = np.concatenate((dfdw_rxn_insta, dfdw_rxn_reversible, dfdw_mass_balance, dfdw_energy_balance), axis=1)

            dfdT_rxn_insta = self.__get_LiquidEquilibrium_dfdT_rxn_insta__(self.LiquidStreamOut, dKwdT)
            dfdT_rxn_reversible = self.__get_LiquidCSTR_dfdT_rxn_reversible__(self.LiquidStreamOut, self.matrix,
                                                                              rates_kmol_m3s)
            dfdT_mass_balance = self.__get_LiquidEquilibrium_dfdT_mass_balance__(self.LiquidStreamOut, self.matrix)
            dfdT_energy_balance = self.__get_LiquidCSTR_dfdT_energy_balance__(self.LiquidStreamOut)
            dfdT = np.concatenate((dfdT_rxn_insta, dfdT_rxn_reversible, dfdT_mass_balance, dfdT_energy_balance), axis=1)

            dfdX = np.concatenate((dfdw, dfdT), axis=2)

            dX_newton = - np.linalg.solve(dfdX, f)
            dw_newton = dX_newton[:, :self.LiquidStreamOut.num_of_species:]
            dT_newton = dX_newton[:, self.LiquidStreamOut.num_of_species]

            dw = lr * dw_newton
            dT = lr * dT_newton

            w = w + dw
            w = np.maximum(w, 0)
            self.LiquidStreamOut.__mass_fractions_vec2dic__(w)
            self.LiquidStreamOut.temp_K = self.LiquidStreamOut.temp_K + dT

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1




        # -------------------------------------------------------------------------------------------




        # Combine
        self.LiquidStreamOut = deepcopy(self.LiquidStreamInit)
        self.LiquidStreamOut.__where__(condition=(vapor_pressure_bara > pressure_bara), LiquidStreamIn=self.LiquidStreamBoiled)
        self.LiquidStreamOut.__where__(condition=(vapor_pressure_bara <= pressure_bara), LiquidStreamIn=self.LiquidStreamHeated)

        return self.LiquidStreamOut, vapor_pressure_bara

    def get_vapor_flow_kmol_h(self):
        return self.GasStreamOut.get_gas_flow_kmol_h()

    def get_vapor_temp_K(self):
        return self.GasStreamOut.get_gas_temp_K()

    def get_vapor_molar_fractions(self):
        y = {}
        for id in self.GasStreamOut.specie.keys():
            y[id] = self.GasStreamOut.get_specie_molar_fraction(id)
        return y


# ---------------------------------------------------------------------------------------


class LiquidHeatExchanger_CounterCurrent(_Stochiometry, _Serializer, _LiquidEquilibrium, _StreamFunctions):

    def __init__(self, num_of_heights):
        self.num_of_heights = num_of_heights
        self.firstscan = True

    def add_heat_transfer_coefficient_kW_m2K(self, heat_transfer_coefficient_kW_m2K):
        self._heat_transfer_coefficient_kW_m2K_ = heat_transfer_coefficient_kW_m2K

    def react(self, LiquidStreamIn1, LiquidStreamIn2, interface_area_m2, volume1_m3=None, volume2_m3=None, lr=0.75):

        _, _  = self.approx(LiquidStreamIn1, LiquidStreamIn2, interface_area_m2)

        # Geometry
        self.volume1_m3 = 1.8 * 10**(-3) * interface_area_m2 if volume1_m3 is None else volume1_m3
        self.volume2_m3 = 1.8 * 10**(-3) * interface_area_m2 if volume2_m3 is None else volume2_m3
        self.interface_area_m2 = interface_area_m2
        self.length1_m = np.ones(shape=self.volume1_m3.shape)       # Assumption
        self.length2_m = np.ones(shape=self.volume2_m3.shape)       # Assumption
        self.cross_sectional_area1_m2 = self.volume1_m3 / self.length1_m
        self.cross_sectional_area2_m2 = self.volume2_m3 / self.length2_m
        self.dz = 1 / self.num_of_heights

        # Strictly Positive Mass Fractions
        for id in LiquidStreamIn1.specie.keys():
            LiquidStreamIn1.specie[id]["Mass Fraction"] = np.maximum(LiquidStreamIn1.specie[id]["Mass Fraction"], 10 ** (-18))
        LiquidStreamIn1.normalize_mass_fractions()

        for id in LiquidStreamIn2.specie.keys():
            LiquidStreamIn2.specie[id]["Mass Fraction"] = np.maximum(LiquidStreamIn2.specie[id]["Mass Fraction"], 10 ** (-18))
        LiquidStreamIn2.normalize_mass_fractions()

        if self.firstscan:
            self.matrix1 = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn1, liq_rxn_insta=True)
            self.matrix2 = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn2, liq_rxn_insta=True)

        self.LiquidStreamIn1 = self.__get_LiquidEquilibrium_reduced__(LiquidStreamIn1, matrix=self.matrix1, lr=lr)
        self.LiquidStreamIn2 = self.__get_LiquidEquilibrium_reduced__(LiquidStreamIn2, matrix=self.matrix2, lr=lr)

        self.LiquidStream1 = deepcopy(self.LiquidStreamIn1)
        self.LiquidStream2 = deepcopy(self.LiquidStreamIn2)

        for epoch in range(10):

            # Integrate LiquidStream1 "Upward"
            self.LiquidStream1 = deepcopy(self.LiquidStreamIn1)
            for i in range(1, self.num_of_heights):

                # Load LiquidStream2 at height = z[i-1]
                self.LiquidStream2 = self.__get_slice_from_LiquidProfile__(self.LP2, i-1, self.LiquidStream2)

                # Calculate Gradients
                dwdz, dTdz = self.__get_gradients_liq1__()

                # Integrate
                #self.LiquidStream1.__mass_fractions_vec2dic__(self.LiquidStream1.__mass_fractions_dic2vec__() + dwdz * dz[:, None])
                self.LiquidStream1.temp_K = self.LiquidStream1.temp_K + dTdz * self.dz
                #self.LiquidStream1.normalize_mass_fractions()

                # Update Liquid Profile
                self.LP1 = self.__insert_slice_to_LiquidProfile__(self.LP1, i, self.LiquidStream1)


            # Integrate LiquidStream2 "Downward"
            self.LiquidStream2 = deepcopy(self.LiquidStreamIn2)
            for i in range(self.num_of_heights - 2, -1, -1):


                # Load LiquidStream1 at height = z[i+1]
                self.LiquidStream1 = self.__get_slice_from_LiquidProfile__(self.LP1, i + 1, self.LiquidStream1)

                # Calculate Gradients
                dwdz, dTdz = self.__get_gradients_liq2__()

                # Integrate
                #self.LiquidStream2.__mass_fractions_vec2dic__(self.LiquidStream2.__mass_fractions_dic2vec__() + dwdz * dz[:, None])
                self.LiquidStream2.temp_K = self.LiquidStream2.temp_K + dTdz * self.dz
                #self.LiquidStream2.normalize_mass_fractions()

                # Update Liquid Profile
                self.LP2 = self.__insert_slice_to_LiquidProfile__(self.LP2, i, self.LiquidStream2)

        return self.LiquidStream1, self.LiquidStream2

    def approx(self, LiquidStreamIn1, LiquidStreamIn2, interface_area_m2):

        self.interface_area_m2 = interface_area_m2
        self.LP1 = self.__broadcast_to_LiquidProfile__(LiquidStreamIn1, self.num_of_heights)
        self.LP2 = self.__broadcast_to_LiquidProfile__(LiquidStreamIn2, self.num_of_heights)


        LiquidStreamOut1 = deepcopy(LiquidStreamIn1)
        LiquidStreamOut2 = deepcopy(LiquidStreamIn2)

        self.LiquidStream1 = deepcopy(LiquidStreamIn1)
        self.LiquidStream2 = deepcopy(LiquidStreamIn2)

        for epoch in range(5):

            cp1 = LiquidStreamIn1.get_solution_heat_capacity_kJ_kgK()
            cp2 = LiquidStreamIn2.get_solution_heat_capacity_kJ_kgK()

            m1 = LiquidStreamIn1.get_solution_flow_kg_h() / 3600
            m2 = LiquidStreamIn2.get_solution_flow_kg_h() / 3600

            # "Capacity"
            C1 = cp1 * m1
            C2 = cp2 * m2
            Cmin = np.minimum(C1, C2)
            Cmax = np.maximum(C1, C2)
            Cr = Cmin / Cmax

            # Maximum Possible Heat Transferred
            Qmax = Cmin * (LiquidStreamIn1.temp_K - LiquidStreamIn2.temp_K)

            # Heat Transfer Coefficient
            self.LiquidStream1.temp_K = 1.0 * LiquidStreamIn1.temp_K
            self.LiquidStream2.temp_K = 1.0 * LiquidStreamOut2.temp_K
            kH_left = self._heat_transfer_coefficient_kW_m2K_(self)

            self.LiquidStream1.temp_K = 1.0 * LiquidStreamOut1.temp_K
            self.LiquidStream2.temp_K = 1.0 * LiquidStreamIn2.temp_K
            kH_right = self._heat_transfer_coefficient_kW_m2K_(self)

            # Average
            kH = (kH_left + kH_right) / 2

            # Number of Transfer Units
            NTU = kH * self.interface_area_m2 / Cmin

            # Actual Heat Transferred
            effectiveness = (1 - np.exp(-NTU * (1 - Cr))) / (1 - Cr * np.exp(- NTU * (1 - Cr)))
            Q = effectiveness * Qmax

            # Calculate Outlet Temperatures
            LiquidStreamOut1.temp_K = LiquidStreamIn1.temp_K - Q / C1
            LiquidStreamOut2.temp_K = LiquidStreamIn2.temp_K + Q / C2
            dT_left = np.abs(LiquidStreamIn1.temp_K - LiquidStreamOut2.temp_K)
            dT_right = np.abs(LiquidStreamOut1.temp_K - LiquidStreamIn2.temp_K)
            #LMTD = (dT_left - dT_right) / (np.log(dT_left) - np.log(dT_right))

        for i in range(LiquidStreamIn1.temp_K.shape[0]):
            self.LP1.temp_K[:,i] = np.linspace(LiquidStreamIn1.temp_K[i], LiquidStreamOut1.temp_K[i], self.num_of_heights)
            self.LP2.temp_K[:, i] = np.linspace(LiquidStreamOut2.temp_K[i], LiquidStreamIn2.temp_K[i], self.num_of_heights)

        return LiquidStreamOut1, LiquidStreamOut2

    def __get_gradients_liq1__(self):
        m1 = self.LiquidStream1.get_solution_flow_kg_h() / 3600
        cp1 = self.LiquidStream1.get_solution_heat_capacity_kJ_kgK()
        Ae = self.interface_area_m2
        dT = self.LiquidStream1.get_solution_temp_K() - self.LiquidStream2.get_solution_temp_K()
        kH = self._heat_transfer_coefficient_kW_m2K_(self)
        q = -kH * dT
        Z = 1.0
        dTdz = (Ae * q) / (m1 * cp1 * Z)
        dwdz = None
        return dwdz, dTdz

    def __get_gradients_liq2__(self):
        m2 = self.LiquidStream2.get_solution_flow_kg_h() / 3600
        cp2 = self.LiquidStream2.get_solution_heat_capacity_kJ_kgK()
        Ae = self.interface_area_m2
        dT = self.LiquidStream1.get_solution_temp_K() - self.LiquidStream2.get_solution_temp_K()
        kH = self._heat_transfer_coefficient_kW_m2K_(self)
        q = kH * dT
        Z = 1.0
        dTdz = (Ae * q) / (m2 * cp2 * Z)
        dwdz = None
        return dwdz, dTdz

    def get_profile_stream1_temp_K(self, sample_index):
        T = self.LP1.get_solution_temp_K()[:, sample_index]
        z = np.linspace(0, 1, self.num_of_heights)
        return T, z

    def get_profile_stream2_temp_K(self, sample_index):
        T = self.LP2.get_solution_temp_K()[:, sample_index]
        z = np.linspace(0, 1, self.num_of_heights)
        return T, z

    def get_profile_stream1_total_vapor_pressure_bara(self):
        pass

    def get_profile_stream2_total_vapor_pressure_bara(self):
        pass

    def get_interface_area_m2(self):
        return self.interface_area_m2


class LiquidHeatExchanger_CounterCurrent_NTU(_Serializer):

    def __init__(self, id, area_m2, heat_transfer_coefficient_kW_m2K):
        super().__init__()
        self.id = id
        self.area_m2 = area_m2
        self.heat_transfer_coefficient_kW_m2K = heat_transfer_coefficient_kW_m2K

    def react(self, LiquidStreamIn1, LiquidStreamIn2):
        self.LiquidStream1 = deepcopy(LiquidStreamIn1)
        self.LiquidStream2 = deepcopy(LiquidStreamIn2)
        LiquidStreamOut1 = deepcopy(LiquidStreamIn1)
        LiquidStreamOut2 = deepcopy(LiquidStreamIn2)
        T1_in = LiquidStreamIn1.get_solution_temp_K()
        T2_in = LiquidStreamIn2.get_solution_temp_K()
        cp1 = LiquidStreamIn1.get_solution_heat_capacity_kJ_kgK()
        cp2 = LiquidStreamIn2.get_solution_heat_capacity_kJ_kgK()
        m1 = LiquidStreamIn1.get_solution_flow_kg_h() / 3600
        m2 = LiquidStreamIn2.get_solution_flow_kg_h() / 3600
        C1 = cp1 * m1
        C2 = cp2 * m2
        Cmin = np.minimum(C1, C2)
        Cmax = np.maximum(C1, C2)
        Cr = Cmin / Cmax
        Qmax = Cmin * (T1_in - T2_in)
        self.kH = self.heat_transfer_coefficient_kW_m2K(self)
        self.NTU = self.kH * self.area_m2 / Cmin
        effectiveness = (1 - np.exp(-self.NTU *(1 - Cr))) / (1 - Cr * np.exp(- self.NTU * (1 - Cr)))
        Q = effectiveness * Qmax
        T1_out = T1_in - Q / C1
        T2_out = T2_in + Q / C2
        dT_left = np.abs(T1_in - T2_out)
        dT_right = np.abs(T1_out - T2_in)
        self.LMTD = (dT_left - dT_right) / (np.log(dT_left) - np.log(dT_right))
        LiquidStreamOut1.temp_K = T1_out
        LiquidStreamOut2.temp_K = T2_out
        return LiquidStreamOut1, LiquidStreamOut2


# ---------------------------------------------------------------------------------------


class GasCompressor_Isentropic(_Serializer):

    def __init__(self):
        super().__init__()
        self.compressor_power_kW = None

    def react(self, GasStreamIn, pressure_out_bara):
        cp = GasStreamIn.get_gas_heat_capacity_kJ_kmolK()
        cv = cp - 8.314
        y = cp / cv
        shape = np.ones(shape=GasStreamIn.temp_K.shape)

        # Outlet Stream
        GasStreamOut = deepcopy(GasStreamIn)
        GasStreamOut.pressure_bara = pressure_out_bara
        GasStreamOut.temp_K = GasStreamIn.temp_K * (pressure_out_bara / GasStreamIn.pressure_bara) ** ((y - 1) / y)

        # Compressor Power
        R = 8.314           # kJ/kmol.K
        dPdn = (y/(y-1)) * R * GasStreamIn.temp_K * ((GasStreamOut.pressure_bara/GasStreamIn.pressure_bara)**((y - 1) / y) - 1)
        self.compressor_power_kW = dPdn * GasStreamIn.flow_kmol_h / 3600

        return GasStreamOut

    def get_compressor_power_kW(self):
        return self.compressor_power_kW


# ---------------------------------------------------------------------------------------


class GasLiquidContactor_PFR(_Stochiometry, _Serializer, _LiquidEquilibrium, _StreamFunctions):

    def __init__(self, position_m, cross_sectional_area_m2, void_fraction_m3_m3):
        self.info = {}
        self.mass_transfer_kmol_m3s = {}
        self.num_of_mass_transfer = 0
        self.firstscan = True
        self.position_m = position_m
        self.cross_sectional_area_m2 = cross_sectional_area_m2
        self.void_fraction_m3_m3 = void_fraction_m3_m3
        self.num_of_heights = len(position_m)
        self.flow = ""

    # --------------------------------------------------------------------

    def add_info(self, key, value):
        self.info[key] = value

    def add_mass_transfer_kmol_m3s(self, id, stoch_gas, stoch_liq, rate_kmol_m3s, exothermic_heat_kJ_kmol=None):
        self.mass_transfer_kmol_m3s[id] = {}
        self.mass_transfer_kmol_m3s[id]["Stoch Gas"] = stoch_gas
        self.mass_transfer_kmol_m3s[id]["Stoch Liq"] = stoch_liq
        self.mass_transfer_kmol_m3s[id]["Rate [kmol/m3.s]"] = rate_kmol_m3s
        self.mass_transfer_kmol_m3s[id]["Exothermic Heat [kJ/kmol]"] = exothermic_heat_kJ_kmol
        self.mass_transfer_kmol_m3s[id]["Index"] = self.num_of_mass_transfer
        self.num_of_mass_transfer = self.num_of_mass_transfer + 1

    def add_heat_transfer_kW_m3(self, heat_transfer_kW_m3):
        self.heat_transfer_kW_m3 = heat_transfer_kW_m3

    def add_liquid_holdup_m3_m3(self, liquid_holdup_m3_m3):
        self.liquid_holdup_m3_m3 = liquid_holdup_m3_m3

    def add_pressure_drop_Pa_m(self, presure_drop_Pa_m):
        pass

    # --------------------------------------------------------------------

    def get_superficial_gas_velocity_m_s(self):
        T = self.GasStream.get_gas_temp_K()
        n = self.GasStream.get_gas_flow_kmol_h() / 3600
        p = self.GasStream.get_gas_pressure_bara()
        R = 0.08314
        V = n * R * T / p
        v = V / self._cross_sectional_area_m2
        return v

    def get_superficial_liquid_velocity_m_s(self):
        rho = self.LiquidStream.get_solution_density_kg_m3()
        m = self.LiquidStream.get_solution_flow_kg_h() / 3600
        V = m / rho
        v = V / self._cross_sectional_area_m2
        return v

    def get_liquid_holdup_m3_m3(self):
        return self.liquid_holdup_m3_m3(self)

    def get_pressure_drop_Pa_m(self):
        return self.pressure_drop_Pa_m(self)

    def get_mass_transfer_kmol_m3s(self, id):
        return self.mass_transfer_kmol_m3s[id]["Rate [kmol/m3.s]"](self)

    def get_heat_transfer_kW_m3(self):
        return self.heat_transfer_kW_m3(self)

    # --------------------------------------------------------------------

    def __react_countercurrent__(self, GasStreamIn, LiquidStreamIn, step_size_m, max_epochs=np.inf, lr=0.75):

        self.flow = "countercurrent"
        self.__prepocessing__(GasStreamIn, LiquidStreamIn, step_size_m, max_epochs, lr)

        # Convergence Stuff...
        w_old = self.LiquidStreamIn.__mass_fractions_dic2vec__()
        y_old = self.GasStreamIn.__molar_fractions_dic2vec__()
        T_old_liq = self.LiquidStreamIn.temp_K
        T_old_gas = self.GasStreamIn.temp_K

        converged = False
        self.epochs = 0

        while converged == False:

            # -----------------------------------------------------------------
            # Integrate LiquidStream Downward
            # -----------------------------------------------------------------

            z = self.position_m[-1]

            self.LiquidStream = deepcopy(self.LiquidStreamIn)
            self.LP = deepcopy(self.LiquidStreamIn)
            self.LP.position_m = [z]

            while z > self.position_m[0]:

                m = self.LiquidStream.get_solution_flow_kg_h()/3600
                T = self.LiquidStream.get_solution_temp_K()
                w = self.LiquidStream.__mass_fractions_dic2vec__()
                X = np.concatenate((m[:,None],T[:,None],w),axis=1)

                X_new, z_new = self.__runge_kutta_liq__(self.__system_dynamics_liq__, X, z, order=1, direction="backward", step_size_fixed=False)
                m_new = X_new[:,0]
                T_new = X_new[:,1]
                w_new = X_new[:,2::]

                self.LiquidStream.__mass_fractions_vec2dic__(w_new)
                self.LiquidStream.temp_K = T_new
                self.LiquidStream.flow_kg_h = 3600*m_new
                self.LiquidStream.normalize_mass_fractions()
                z = z_new

                # Update Liquid Profile
                self.LP.temp_K = np.vstack((self.LP.temp_K, self.LiquidStream.temp_K))
                self.LP.flow_kg_h = np.vstack((self.LP.flow_kg_h, self.LiquidStream.flow_kg_h))
                for id in self.LiquidStream.specie.keys():
                    self.LP.specie[id]["Mass Fraction"] = np.vstack((self.LP.specie[id]["Mass Fraction"], self.LiquidStream.specie[id]["Mass Fraction"]))
                for id in self.LiquidStream.info.keys():
                    self.LP.info[id] = np.vstack((self.LP.info[id], self.LiquidStream.info[id]))
                self.LP.position_m.append(z)

            # Arrify
            self.LP.position_m = np.array(self.LP.position_m)

            # -----------------------------------------------------------------
            # Integrate GasStream Upward
            # -----------------------------------------------------------------
            z = self.position_m[0]

            self.GasStream = deepcopy(self.GasStreamIn)
            self.GP = deepcopy(self.GasStreamIn)
            self.GP.position_m = [z]

            while z < self.position_m[-1]:

                n = self.GasStream.get_gas_flow_kmol_h() / 3600
                T = self.GasStream.get_gas_temp_K()
                y = self.GasStream.__molar_fractions_dic2vec__()

                X = np.concatenate((n[:, None], T[:, None], y), axis=1)
                X_new, z_new = self.__runge_kutta_gas__(self.__system_dynamics_gas__, X, z, order=1, direction="forward", step_size_fixed=False)
                n_new = X_new[:, 0]
                T_new = X_new[:, 1]
                y_new = X_new[:, 2::]

                self.GasStream.__molar_fractions_vec2dic__(y_new)
                self.GasStream.temp_K = T_new
                self.GasStream.flow_kmol_h = 3600 * n_new
                self.GasStream.normalize_molar_fractions()
                z = z_new

                # Update Gas Profile
                self.GP.temp_K = np.vstack((self.GP.temp_K, self.GasStream.temp_K))
                self.GP.flow_kmol_h = np.vstack((self.GP.flow_kmol_h, self.GasStream.flow_kmol_h))
                self.GP.pressure_bara = np.vstack((self.GP.pressure_bara, self.GasStream.pressure_bara))
                for id in self.GasStream.specie.keys():
                    self.GP.specie[id]["Molar Fraction"] = np.vstack((self.GP.specie[id]["Molar Fraction"], self.GasStream.specie[id]["Molar Fraction"]))
                for id in self.GasStream.info.keys():
                    self.GP.info[id] = np.vstack((self.GP.info[id], self.GasStream.info[id]))
                self.GP.position_m.append(z)

            # Arrify
            self.GP.position_m = np.array(self.GP.position_m)

            # Check if Algorithm have Converged
            w = self.LiquidStream.__mass_fractions_dic2vec__()
            y = self.GasStream.__molar_fractions_dic2vec__()

            T_liq = self.LiquidStream.temp_K
            T_gas = self.GasStream.temp_K
            dw = w - w_old
            dy = y - y_old
            dT_gas = T_gas - T_old_gas
            dT_liq = T_liq - T_old_liq
            liq_converged = np.array((np.abs(dw) < 0.01 * np.abs(w)) | (w <= 10 ** (-18)), dtype=np.float32)
            liq_converged = np.min(liq_converged, axis=1, keepdims=False)
            gas_converged = np.array((np.abs(dy) < 0.01 * np.abs(y)) | (y <= 10 ** (-18)), dtype=np.float32)
            gas_converged = np.min(gas_converged, axis=1, keepdims=False)
            sample_converged = np.minimum(liq_converged, gas_converged)
            sample_converged = sample_converged * (np.abs(dT_gas) < 0.05)
            sample_converged = sample_converged * (np.abs(dT_liq) < 0.05)
            converged = np.min(sample_converged, axis=0)
            w_old = 1.0 * w
            y_old = 1.0 * y
            T_old_liq = 1.0 * T_liq
            T_old_gas = 1.0 * T_gas
            self.epochs = self.epochs + 1
            converged = converged or (self.epochs>=max_epochs)

        print(self.epochs)
        return self.GasStream, self.equilibrium_liquid.react(self.LiquidStream, lr=lr)

    def __react_cocurrent__(self, GasStreamIn, LiquidStreamIn, step_size_m, max_epochs=np.inf, lr=0.75):

        self.flow = "cocurrent"
        self.__prepocessing__(GasStreamIn, LiquidStreamIn, step_size_m, max_epochs, lr)

        # Integrate both Gas and LiquidStream Downward
        z = self.position_m[-1]

        self.LiquidStream = deepcopy(self.LiquidStreamIn)
        self.LP = deepcopy(self.LiquidStreamIn)
        self.LP.position_m = [z]

        self.GasStream = deepcopy(self.GasStreamIn)
        self.GP = deepcopy(self.GasStreamIn)
        self.GP.position_m = [z]

        while z > self.position_m[0]:

            m = self.LiquidStream.get_solution_flow_kg_h() / 3600
            T_liq = self.LiquidStream.get_solution_temp_K()
            w = self.LiquidStream.__mass_fractions_dic2vec__()
            X = np.concatenate((m[:, None], T_liq[:, None], w), axis=1)

            X_new, z_new = self.__runge_kutta_liq__(self.__system_dynamics_liq__, X, z, order=1, direction="backward", step_size_fixed=True)
            m_new = X_new[:, 0]
            T_liq_new = X_new[:, 1]
            w_new = X_new[:, 2::]

            # --------------------------------------------------------------

            n = self.GasStream.get_gas_flow_kmol_h() / 3600
            T_gas = self.GasStream.get_gas_temp_K()
            y = self.GasStream.__molar_fractions_dic2vec__()

            X = np.concatenate((n[:, None], T_gas[:, None], y), axis=1)
            X_new, z_new = self.__runge_kutta_gas__(self.__system_dynamics_gas__, X, z, order=1, direction="backward", step_size_fixed=True)
            n_new = X_new[:, 0]
            T_gas_new = X_new[:, 1]
            y_new = X_new[:, 2::]

            # --------------------------------------------------------------

            self.LiquidStream.__mass_fractions_vec2dic__(w_new)
            self.LiquidStream.temp_K = T_liq_new
            self.LiquidStream.flow_kg_h = 3600 * m_new
            self.LiquidStream.normalize_mass_fractions()
            z = z_new

            self.GasStream.__molar_fractions_vec2dic__(y_new)
            self.GasStream.temp_K = T_gas_new
            self.GasStream.flow_kmol_h = 3600 * n_new
            self.GasStream.normalize_molar_fractions()
            z = z_new

            # --------------------------------------------------------------

            self.LP.temp_K = np.vstack((self.LP.temp_K, self.LiquidStream.temp_K))
            self.LP.flow_kg_h = np.vstack((self.LP.flow_kg_h, self.LiquidStream.flow_kg_h))
            for id in self.LiquidStream.specie.keys():
                self.LP.specie[id]["Mass Fraction"] = np.vstack(
                    (self.LP.specie[id]["Mass Fraction"], self.LiquidStream.specie[id]["Mass Fraction"]))
            for id in self.LiquidStream.info.keys():
                self.LP.info[id] = np.vstack((self.LP.info[id], self.LiquidStream.info[id]))
            self.LP.position_m.append(z)

            self.GP.temp_K = np.vstack((self.GP.temp_K, self.GasStream.temp_K))
            self.GP.flow_kmol_h = np.vstack((self.GP.flow_kmol_h, self.GasStream.flow_kmol_h))
            self.GP.pressure_bara = np.vstack((self.GP.pressure_bara, self.GasStream.pressure_bara))
            for id in self.GasStream.specie.keys():
                self.GP.specie[id]["Molar Fraction"] = np.vstack(
                    (self.GP.specie[id]["Molar Fraction"], self.GasStream.specie[id]["Molar Fraction"]))
            for id in self.GasStream.info.keys():
                self.GP.info[id] = np.vstack((self.GP.info[id], self.GasStream.info[id]))
            self.GP.position_m.append(z)

        # Arrify
        self.LP.position_m = np.array(self.LP.position_m)
        self.GP.position_m = np.array(self.GP.position_m)

        return self.GasStream, self.equilibrium_liquid.react(self.LiquidStream, lr=lr)

    def __react_constant_gas__(self, GasStreamIn, LiquidStreamIn, step_size_m, max_epochs=np.inf, lr=0.75):
        pass

    def __react_constant_liq__(self, GasStreamIn, LiquidStreamIn, step_size_m, max_epochs=np.inf, lr=0.75):
        pass

    def __react_stirred_gas__(self, GasStreamIn, LiquidStreamIn, step_size_m, max_epochs=np.inf, lr=0.75):
        pass

    def __react_stirred_liq__(self, GasStreamIn, LiquidStreamIn, step_size_m, max_epochs=np.inf, lr=0.75):
        pass

    # --------------------------------------------------------------------

    def __prepocessing__(self, GasStreamIn, LiquidStreamIn, step_size_m, max_epochs=np.inf, lr=0.75):

        # Step Size of Integrator
        self.step_size_m = step_size_m

        # Only Strictly Positive Mass Fractions Allowed
        for id in LiquidStreamIn.specie.keys():
            LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(LiquidStreamIn.specie[id]["Mass Fraction"],
                                                                    10 ** (-18))
        LiquidStreamIn.normalize_mass_fractions()
        for id in GasStreamIn.specie.keys():
            GasStreamIn.specie[id]["Molar Fraction"] = np.maximum(GasStreamIn.specie[id]["Molar Fraction"], 10 ** (-18))
        GasStreamIn.normalize_molar_fractions()

        if self.firstscan:

            # Extract Molar Masses
            self.M_liq = np.zeros(shape=(LiquidStreamIn.num_of_species,))
            for i, id in enumerate(LiquidStreamIn.specie.keys()):
                self.M_liq[i] = LiquidStreamIn.get_specie_molar_mass_kg_kmol(id)

            self.M_gas = np.zeros(shape=(GasStreamIn.num_of_species))
            for i, id in enumerate(GasStreamIn.specie.keys()):
                self.M_gas[i] = GasStreamIn.get_specie_molar_mass_kg_kmol(id)

            # Matrix or Reaction Stochiometry: Used to Calculate Sensitivity Matrices
            self.matrix_liq = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True)
            self.matrix_gas = self.__get_the_matrix__(GasStreamIn=GasStreamIn, gas_rxn_insta=True)

            # Matrix of Reaction Stochiometry: Used to Calculate Heat of Reactions
            self.matrix_liq_heat = self.__get_the_matrix__(GasStreamIn=GasStreamIn,
                                                           LiquidStreamIn=LiquidStreamIn,
                                                           Column=self,
                                                           liq_rxn_insta=True,
                                                           liq_rxn_reversible=True,
                                                           vapor_pressure=False,
                                                           mass_transfer=True,
                                                           gas_rxn_insta=False,
                                                           gas_rxn_reversible=False)

            self.matrix_gas_heat = self.__get_the_matrix__(GasStreamIn=GasStreamIn,
                                                           gas_rxn_insta=True,
                                                           gas_rxn_reversible=True)

            # Instance for Calculating Chemical Equilibrium
            self.equilibrium_liquid = LiquidEquilibrium_Adiabatic()
            self.equilibrium_gas = None

        # Process Liquid and Gas Inlet Streams
        self.LiquidStreamIn = self.equilibrium_liquid.react(LiquidStreamIn, lr=lr)
        self.GasStreamIn = deepcopy(GasStreamIn)

        # Initiate Streams used for Propagation up and Down
        self.LiquidStream = deepcopy(self.LiquidStreamIn)
        self.GasStream = deepcopy(self.GasStreamIn)

        # Initiate Profiles
        if self.firstscan:
            self.LP = self.__broadcast_to_LiquidProfile__(LiquidStreamIn=self.LiquidStream,
                                                          num_of_heights=self.num_of_heights)
            self.GP = self.__broadcast_to_GasProfile__(GasStreamIn=self.GasStream, num_of_heights=self.num_of_heights)
            self.LP.position_m = self.position_m
            self.GP.position_m = self.position_m
            self.firstscan = False

    def __system_dynamics_liq__(self, X, z, direction):

        m = X[:,0]
        T = X[:,1]
        w = X[:,2::]

        # Geometry at Position = z
        self._cross_sectional_area_m2 = np.interp(z, self.position_m, self.cross_sectional_area_m2)
        self._void_fraction_m3_m3 = np.interp(z, self.position_m, self.void_fraction_m3_m3)

        # LiquidStream at Position = z
        self.LiquidStream.__mass_fractions_vec2dic__(w)
        self.LiquidStream.temp_K = T
        self.LiquidStream.flow_kg_h = 3600 * m

        # GasStream at Position = z
        if self.flow == "countercurrent":
            self.GasStream = self.__interpolate_from_GasProfile__(self.GP, self.GP.position_m, z, self.GasStream)

        # Exothermic Heat [kJ/kmol rxn]
        h = np.zeros(shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.num_of_rxn_insta + self.LiquidStream.num_of_rxn_reversible + self.num_of_mass_transfer), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(self.LiquidStream.rxn_insta.keys()):
            h[:, rxn_i] = self.LiquidStream.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
        for rxn_i, rxn_id in enumerate(self.LiquidStream.rxn_reversible.keys()):
            h[:,rxn_i + self.LiquidStream.num_of_rxn_insta] = self.LiquidStream.get_rxn_reversible_exothermic_heat_kJ_kmol(rxn_id)
        for rxn_i, rxn_id in enumerate(self.mass_transfer_kmol_m3s.keys()):
            h[:, rxn_i + self.LiquidStream.num_of_rxn_insta + self.LiquidStream.num_of_rxn_reversible] = self.mass_transfer_kmol_m3s[rxn_id]["Exothermic Heat [kJ/kmol]"](self)

        # Liquid and Gas Mass Flow
        m_liq_tot = self.LiquidStream.get_solution_flow_kg_h() / 3600
        m_liq = m_liq_tot[:, None] * w
        m_gas_tot = self.GasStream.get_gas_flow_kg_h() / 3600
        m_gas = m_gas_tot[:,None] * self.GasStream.__molefrac2massfrac__(y=self.GasStream.__molar_fractions_dic2vec__())

        # Calculate Absorption Rates and some Related Heat Dissipation
        aJ_rxn, aJ_specie_gas, aJ_specie_liq, q_abs, q_des, q_rxn = self.__get_mass_transfer_vectors__()

        # Calculate Rate-Reactions in Liquid Phase
        r = self.LiquidStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

        # Packing Hydrodynamics
        liquid_holdup_m3_m3 = self.get_liquid_holdup_m3_m3()

        # Differential Equations (dm/dz)
        dmdz_liq = self._cross_sectional_area_m2 * (self.M_liq[None,:] * aJ_specie_liq + self._void_fraction_m3_m3 * liquid_holdup_m3_m3[:,None] * self.M_liq[None,:] * r)
        dmdz_liq_tot = self._cross_sectional_area_m2 * np.sum(self.M_liq[None,:] * aJ_specie_liq, axis=1, keepdims=False)
        dmdz_gas = self._cross_sectional_area_m2 * self.M_gas[None,:] * aJ_specie_gas
        dmdz_gas_tot = - dmdz_liq_tot

        # Sensitivity in Liquid Phase (dmdB, dmdT)
        m_liq_broadcast = np.broadcast_to(array=m_liq[:,:,None], shape=(m_liq.shape[0], m_liq.shape[1], m_liq.shape[1]))
        dwdm_liq = - np.einsum("smn,s->smn", m_liq_broadcast, 1/m_liq_tot**2) + np.einsum("vw,s->svw", np.eye(N=self.LiquidStream.num_of_species), 1 / m_liq_tot)
        dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStream)
        dfdm_rxn_insta = np.einsum("sfw,swv->sfv", dfdw_rxn_insta, dwdm_liq)
        dfdm_mass_balance = np.broadcast_to(array=self.matrix_liq["A"][None, :, :], shape=(self.LiquidStream.temp_K.shape[0], self.matrix_liq["A"].shape[0], self.matrix_liq["A"].shape[1])).copy()
        dfdm = np.concatenate((dfdm_rxn_insta, dfdm_mass_balance), axis=1)
        H = np.linalg.inv(dfdm)
        dmdB = H[:, :, self.LiquidStream.num_of_rxn_insta::]
        Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStream)
        dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStream, Kw)
        Hr = H[:, :, :self.LiquidStream.num_of_rxn_insta:]
        dmdT = - np.einsum("scr,sr->sc", Hr, (1 / Kw) * dKwdT)

        # Change in "B-Vector" for Liquid Phase
        dBdz = np.einsum("bw,sw->sb", self.matrix_liq["A"], dmdz_liq)

        # Direct Heat Transfer (kW/m3)
        q_dir = self.get_heat_transfer_kW_m3()

        # Heat Capacity
        cp = self.LiquidStream.get_solution_heat_capacity_kJ_kgK()

        # Initial Guess
        dmdz_liq_at_eq = dmdz_liq

        # Iterate Until "Convergence"
        for _ in range(10):
            dmdz_at_eq = np.concatenate((dmdz_liq_at_eq, dmdz_gas), axis=1)
            q_lat = np.einsum("sr,sr->s", h, np.einsum("rm,sm->sr", self.matrix_liq_heat["R+"], dmdz_at_eq)) / self._cross_sectional_area_m2
            dTdz = self._cross_sectional_area_m2 * (1 / (cp * m_liq_tot)) * (q_dir + q_abs + q_lat)
            dmdz_liq_at_eq = np.einsum("swb,sb->sw", dmdB, dBdz) + np.einsum("sm,s->sm", dmdT, dTdz)

        dwdz_liq = np.einsum("swm,sm->sw", dwdm_liq, dmdz_liq_at_eq)
        dXdz = np.concatenate((dmdz_liq_tot[:,None], dTdz[:,None], dwdz_liq), axis=1)

        if direction == "backward":
            dXdz = - dXdz
        return dXdz

    def __system_dynamics_gas__(self, X, z, direction):

        n = X[:, 0]
        T = X[:, 1]
        y = X[:, 2::]

        # Geometry at Position = z
        self._cross_sectional_area_m2 = np.interp(z, self.position_m, self.cross_sectional_area_m2)
        self._void_fraction_m3_m3 = np.interp(z, self.position_m, self.void_fraction_m3_m3)

        # GasStream at Position = z
        self.GasStream.__molar_fractions_vec2dic__(y)
        self.GasStream.temp_K = T
        self.GasStream.flow_kmol_h = 3600 * n

        # LiquidStream at Position = z
        if self.flow == "countercurrent":
            self.LiquidStream = self.__interpolate_from_LiquidProfile__(self.LP, self.LP.position_m, z, self.LiquidStream)

        # Other Stuff...
        cp = self.GasStream.get_gas_heat_capacity_kJ_kmolK()
        p_gas = self.GasStream.get_gas_pressure_bara()
        v_GS = self.get_superficial_gas_velocity_m_s()
        R = 0.08314

        # Exothermic Heat
        h = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.num_of_rxn_insta + self.GasStream.num_of_rxn_reversible), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(self.GasStream.rxn_insta.keys()):
            h[:, rxn_i] = self.GasStream.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
        for rxn_i, rxn_id in enumerate(self.GasStream.rxn_reversible.keys()):
            h[:,rxn_i + self.GasStream.num_of_rxn_insta] = self.GasStream.get_rxn_reversible_exothermic_heat_kJ_kmol(rxn_id)

        # Molar Flow
        n_tot = self.GasStream.get_gas_flow_kmol_h() / 3600
        n = n_tot[:,None] * self.GasStream.__molar_fractions_dic2vec__()

        # Calculate Absorption Rates and some Related Heat Dissipation
        aJ_rxn, aJ_specie_gas, aJ_specie_liq, q_abs, q_des, _ = self.__get_mass_transfer_vectors__()

        # Calculate Rate-Reactions in Gas Phase
        r = self.GasStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

        # Differential Equations (dn/dz)
        liquid_holdup_m3_m3 = self.get_liquid_holdup_m3_m3()

        dndz = self._cross_sectional_area_m2 * (aJ_specie_gas + self._void_fraction_m3_m3 * (1 - liquid_holdup_m3_m3[:,None]) * r)
        dndz_tot = np.sum(dndz, axis=1, keepdims=False)

        # Quotient Rule
        dydz = (dndz * n_tot[:,None] - n * dndz_tot[:,None]) / n_tot[:,None] ** 2

        # Direct Heat Transfer (kW/m3)
        q_dir = self.get_heat_transfer_kW_m3()
        q_gas = - q_dir + q_des

        # Temperature Gradient
        dTdz = ((R * T) / (p_gas * cp * v_GS)) * q_gas

        # Output
        dXdz = np.concatenate((dndz_tot[:, None], dTdz[:, None], dydz), axis=1)
        if direction == "backward":
            dXdz = - dXdz

        return dXdz

    def __runge_kutta_liq__(self, f, X, z, order, direction, step_size_fixed):

        if order == 1:

            dXdz = f(X, z, direction)
            h = self.__step_size_liq__(X, dXdz, z, direction, step_size_fixed)

            X_new = X + dXdz * h
            z_new = z + h

        elif order == 2:

            dXdz1 = f(X, z, direction)
            h = self.__step_size_liq__(X, dXdz1, z, direction, step_size_fixed)
            dX1 = h * dXdz1

            dXdz2 = f(X + dX1, z + h, direction)

            X_new = X + 0.5 * (dXdz1 + dXdz2) * h
            z_new = z + h

        elif order == 4:
            dXdz1 = f(X, z, direction)
            h = self.__step_size_liq__(X, dXdz1, z, direction, step_size_fixed)
            dX1 = h * dXdz1

            dXdz2 = f(X + 0.5 * dX1, z + 0.5 * h, direction)
            dX2 = h * dXdz2

            dXdz3 = f(X + 0.5 * dX2, z + 0.5 * h, direction)
            dX3 = h * dXdz3

            dXdz4 = f(X + dX3, z + h, direction)
            dX4 = h * dXdz4

            X_new = X + (dX1 + 2 * dX2 + 2 * dX3 + dX4) / 6
            z_new = z + h

        return X_new, z_new

    def __runge_kutta_gas__(self, f, X, z, order, direction, step_size_fixed):

        if order == 1:
            dXdz = f(X, z, direction)
            h = self.__step_size_gas__(X, dXdz, z, direction, step_size_fixed)
            X_new = X + dXdz * h
            z_new = z + h

        elif order == 2:
            dXdz1 = f(X, z, direction)
            h = self.__step_size_gas__(X, dXdz1, z, direction, step_size_fixed)
            dX1 = h * dXdz1
            dXdz2 = f(X + dX1, z + h, direction)
            X_new = X + 0.5 * (dXdz1 + dXdz2) * h
            z_new = z + h

        elif order == 4:
            dXdz1 = f(X, z, direction)
            h = self.__step_size_gas__(X, dXdz1, z, direction, step_size_fixed)
            dX1 = h * dXdz1

            dXdz2 = f(X + 0.5 * dX1, z + 0.5 * h, direction)
            dX2 = h * dXdz2

            dXdz3 = f(X + 0.5 * dX2, z + 0.5 * h, direction)
            dX3 = h * dXdz3

            dXdz4 = f(X + dX3, z + h, direction)
            dX4 = h * dXdz4

            X_new = X + (dX1 + 2 * dX2 + 2 * dX3 + dX4) / 6
            z_new = z + h

        return X_new, z_new

    def __step_size_liq__(self, X, dXdz, z, direction, fixed):

        if fixed:
            h = self.step_size_m
        else:
            m = X[:, 0]
            T = X[:, 1]
            w = X[:, 2::]
            dmdz = dXdz[:, 0]
            dTdz = dXdz[:, 1]
            dwdz = dXdz[:, 2::]
            with np.errstate(divide='ignore', invalid='ignore'):
                h1 = 0.1 * np.min(np.nan_to_num(w / np.abs(dwdz), nan=np.inf, posinf=np.inf, neginf=np.inf), axis=1, keepdims=False)
            with np.errstate(divide='ignore', invalid='ignore'):
                h2 = 5 / np.abs(dTdz)
            h1 = np.min(h1)
            h2 = np.min(h2)
            h = np.minimum(h1, h2)
            h = np.minimum(h, self.step_size_m)

        if direction == "forward":
            h = np.minimum(h, self.position_m[-1] - z)
        if direction == "backward":
            h = np.minimum(h, z - self.position_m[0])
            h = - h

        return h

    def __step_size_gas__(self, X, dXdz, z, direction, fixed):

        if fixed:
            h = self.step_size_m
        else:
            m = X[:, 0]
            T = X[:, 1]
            w = X[:, 2::]
            dmdz = dXdz[:, 0]
            dTdz = dXdz[:, 1]
            dwdz = dXdz[:, 2::]
            with np.errstate(divide='ignore', invalid='ignore'):
                h1 = 0.1 * np.min(np.nan_to_num(w / np.abs(dwdz), nan=np.inf, posinf=np.inf, neginf=np.inf), axis=1, keepdims=False)
            with np.errstate(divide='ignore', invalid='ignore'):
                h2 = 5 / np.abs(dTdz)
            h1 = np.min(h1)
            h2 = np.min(h2)
            h = np.minimum(h1, h2)
            h = np.minimum(h, self.step_size_m)

        if direction == "forward":
            h = np.minimum(h, self.position_m[-1] - z)
        if direction == "backward":
            h = np.minimum(h, z - self.position_m[0])
            h = - h
        return h

    def __get_mass_transfer_vectors__(self):

        dT = self.GasStream.get_gas_temp_K() - self.LiquidStream.get_solution_temp_K()
        r_specie_liq = np.zeros(shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.num_of_species), dtype=np.float64)
        r_specie_gas = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.num_of_species), dtype=np.float64)
        r_rxn = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.num_of_mass_transfer), dtype=np.float64)
        q_rxn = np.zeros(shape=(self.GasStream.temp_K.shape[0],), dtype=np.float64)

        for k, id in enumerate(self.mass_transfer_kmol_m3s.keys()):
            r_rxn[:, k] = self.mass_transfer_kmol_m3s[id]["Rate [kmol/m3.s]"](self)

        for k, id in enumerate(self.mass_transfer_kmol_m3s.keys()):
            nu = self.mass_transfer_kmol_m3s[id]["Stoch Liq"]
            for el in nu.keys():
                i = self.LiquidStream.specie[el]["Index"]
                r_specie_liq[:, i] = r_specie_liq[:, i] + r_rxn[:, k] * nu[el]
            nu = self.mass_transfer_kmol_m3s[id]["Stoch Gas"]
            for el in nu.keys():
                i = self.GasStream.specie[el]["Index"]
                r_specie_gas[:, i] = r_specie_gas[:, i] + r_rxn[:, k] * nu[el]

        for k, id in enumerate(self.mass_transfer_kmol_m3s.keys()):
            # Heat of Absorption
            q_rxn = q_rxn + r_rxn[:, k] * self.mass_transfer_kmol_m3s[id]["Exothermic Heat [kJ/kmol]"](self)

        # Heat due to Heating/Cooling of Gas during Absorption
        q_abs = np.zeros(shape=(self.GasStream.temp_K.shape[0],), dtype=np.float64)
        q_des = np.zeros(shape=(self.GasStream.temp_K.shape[0],), dtype=np.float64)
        for i, id in enumerate(self.GasStream.specie.keys()):
            cp = self.GasStream.get_specie_heat_capacity_kJ_kmolK(id)
            q_abs = q_abs + dT * cp * np.maximum(r_specie_gas[:, i], 0)
            q_des = q_des + dT * cp * np.minimum(r_specie_gas[:, i], 0)

        return r_rxn, r_specie_gas, r_specie_liq, q_abs, q_des, q_rxn


class GasLiquidContactor_CSTR():

    def __init__(self):
        pass

    def __react__(self, GasStreamIn, LiquidStreamIn):
        pass


# ---------------------------------------------------------------------------------------


class Dynamic_LiquidCSTR_Adiabatic():

    def __init__(self):
        pass


class Dynamic_GasLiquidContactor_CSTR_Adiabatic():

    def __init__(self):
        pass


class Dynamic_GasLiquidContactor_CSTR_Isothermal():

    def __init__(self):
        pass



# ---------------------------------------------------------------------------------------



class OLD:


    class _GasLiquidContactor_PFR_Scipy():

        def __init__(self):
            pass

        # --------------------------------------------------------------------

        def add_info(self, key, value):
            self.info[key] = value

        def __f_euler_liq__(self, z, X):

            # Mass Flow, Temperature and Mass Fractions
            m = np.array([X[0]], dtype=np.float64)
            T = np.array([X[1]], dtype=np.float64)
            w = X[2::][None,:]

            # Geometry at Position = z
            self._cross_sectional_area_m2 = np.interp(z, self._profile_position_m, self._profile_cross_sectional_area_m2)
            self._void_fraction_m3_m3 = np.interp(z, self._profile_position_m, self._profile_void_fraction_m3_m3)

            # LiquidStream at Position = z
            self.LiquidStream.__mass_fractions_vec2dic__(w)
            self.LiquidStream.temp_K = T
            self.LiquidStream.flow_kg_h = 3600 * m
            for id in self.LiquidStream.info.keys():
                self.LiquidStream.info[id] = np.array([self.LiquidStreamIn.info[id][self.s]])       # Ugly...

            # GasStream at Position = z
            self.GasStream = self.__interpolate_from_GasProfile__(self.GP[self.s], self.GP[self.s].position_m, z, self.GasStream)

            cp = self.LiquidStream.get_solution_heat_capacity_kJ_kgK()
            v_LS = self.get_superficial_liquid_velocity_m_s()
            rho = self.LiquidStream.get_solution_density_kg_m3()

            # Exothermic Heat
            h = np.zeros(shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.num_of_rxn_insta + self.LiquidStream.num_of_rxn_reversible + self.num_of_mass_transfer), dtype=np.float64)
            for rxn_i, rxn_id in enumerate(self.LiquidStream.rxn_insta.keys()):
                h[:, rxn_i] = self.LiquidStream.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.LiquidStream.rxn_reversible.keys()):
                h[:,rxn_i + self.LiquidStream.num_of_rxn_insta] = self.LiquidStream.get_rxn_reversible_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.mass_transfer_kmol_m3s.keys()):
                h[:, rxn_i + self.LiquidStream.num_of_rxn_insta + self.LiquidStream.num_of_rxn_reversible] = self.mass_transfer_kmol_m3s[rxn_id]["Exothermic Heat [kJ/kmol]"](self)

            # Liquid and Gas Mass Flow
            m_liq_tot = self.LiquidStream.get_solution_flow_kg_h() / 3600
            m_liq = m_liq_tot[:, None] * w
            m_gas_tot = self.GasStream.get_gas_flow_kg_h() / 3600
            m_gas = m_gas_tot * self.GasStream.__molefrac2massfrac__(y=self.GasStream.__molar_fractions_dic2vec__())

            # Calculate Hessian and Sensitivity Matrix (dw/db)
            Kw = self.__get_LiquidEquilibrium_Kw__(LiquidStreamIn = self.LiquidStream)
            dwdb, dwdT = self.__get_LiquidEquilibrium_sensitivities__(self.LiquidStream, Kw, self.matrix_liq)

            # Calculate Absorption Rates and some Related Heat Dissipation
            aJ_rxn, aJ_specie_gas, aJ_specie_liq, q_abs, q_des = self.__get_mass_transfer_vectors__()

            # Calculate Rate-Reactions in Liquid Phase
            r = self.LiquidStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

            # Packing Hydrodynamics
            liquid_holdup_m3_m3 = 0.1

            # Differential Equations (dm/dz)
            dmdz_liq = self._cross_sectional_area_m2 * (self.M_liq[self.s,:] * aJ_specie_liq + self._void_fraction_m3_m3 * 0.1 * self.M_liq[self.s,:] * r)
            dmdz_liq_tot = self._cross_sectional_area_m2 * np.sum(self.M_liq[self.s,:] * aJ_specie_liq, axis=1, keepdims=False)
            dmdz_gas = self._cross_sectional_area_m2 * self.M_gas[self.s,:] * aJ_specie_gas
            dmdz_gas_tot = - dmdz_liq_tot

            # Quotient Rule
            dwdz_liq = (dmdz_liq * m_liq_tot[:,None] - m_liq * dmdz_liq_tot[:,None]) / m_liq_tot[:,None] ** 2
            dwdz_gas = (dmdz_gas * m_gas_tot[:,None] - m_gas * dmdz_gas_tot[:,None]) / m_gas_tot[:,None] ** 2

            # Derivatives (dw/dz) and (dm/dz) also including equilibrium constraints (but only Isothermal Case)
            dbdz = np.einsum("bw,sw->sb", self.matrix_liq["A"], dwdz_liq)
            dwdz_liq = np.einsum("swb,sb->sw", dwdb, dbdz)

            # Direct Heat Transfer (kW/m3)
            q_dir = self.get_heat_transfer_kW_m3()

            dwdz_liq_isothermal = 1.0 * dwdz_liq

            # Heat and Temperature
            for _ in range(4):

                dmdz_liq = (dwdz_liq * m_liq_tot[:,None] ** 2 + m_liq * dmdz_liq_tot[:,None]) / m_liq_tot[:,None]

                # Heat into Liquid Phase due to Chemical Reactions (kW/m3)
                dmdz = np.concatenate((dmdz_liq, dmdz_gas), axis=1)
                drdz = np.einsum("rm,sm->sr", self.matrix_liq_heat["R+"], dmdz)
                q_lat = (m_gas_tot + m_liq_tot) * np.einsum("sr,sr->s", h, drdz) / 0.5  # [kW/m3] = [kg/s] [kJ/kmol] [kmol/kg.m] * [1/m2]

                # Overall Heat to Liquid Phase
                q_liq = q_dir + q_abs + q_lat

                # Temperature Gradient
                dTdz = (1 / (rho * cp * v_LS)) * q_liq

                # Compensate for temperature effects
                dwdz_liq = dwdz_liq_isothermal + np.einsum("sw,s->sw", dwdT, dTdz)

            dXdz = np.hstack((-dmdz_liq_tot, -dTdz, -dwdz_liq.flatten()))
            return dXdz

        def __f_euler_gas__(self, n, T, y, z):

            # Geometry at Position = z
            self._cross_sectional_area_m2 = np.interp(z, self._profile_position_m, self._profile_cross_sectional_area_m2)
            self._void_fraction_m3_m3 = np.interp(z, self._profile_position_m, self._profile_void_fraction_m3_m3)

            # GasStream at Position = z
            self.GasStream.__molar_fractions_vec2dic__(y)
            self.GasStream.temp_K = T
            self.GasStream.flow_kmol_h = 3600 * n

            # LiquidStream at Position = z
            self.LiquidStream = self.__interpolate_from_LiquidProfile__(self.LP, self.zp_liq, z, self.LiquidStream)

            # Other Stuff...
            cp = self.GasStream.get_gas_heat_capacity_kJ_kmolK()
            p_gas = self.GasStream.get_gas_pressure_bara()
            v_GS = self.get_superficial_gas_velocity_m_s()
            R = 0.08314

            # Exothermic Heat
            h = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.num_of_rxn_insta + self.GasStream.num_of_rxn_reversible), dtype=np.float64)
            for rxn_i, rxn_id in enumerate(self.GasStream.rxn_insta.keys()):
                h[:, rxn_i] = self.GasStream.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.GasStream.rxn_reversible.keys()):
                h[:,rxn_i + self.GasStream.num_of_rxn_insta] = self.GasStream.get_rxn_reversible_exothermic_heat_kJ_kmol(rxn_id)

            # Molar Flow
            n_tot = self.GasStream.get_gas_flow_kmol_h() / 3600
            n = n_tot * self.GasStream.__molar_fractions_dic2vec__()

            # Calculate Absorption Rates and some Related Heat Dissipation
            aJ_rxn, aJ_specie_gas, aJ_specie_liq, q_abs, q_des = self.__get_mass_transfer_vectors__()

            # Calculate Rate-Reactions in Gas Phase
            r = self.GasStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

            # Packing Hydrodynamics
            liquid_holdup_m3_m3 = 0.1

            # Differential Equations (dn/dz)
            dndz = 0.5 * (aJ_specie_gas + 0.98 * (1 - 0.1) * r)
            dndz_tot = np.sum(dndz, axis=1, keepdims=False)

            # Quotient Rule
            dydz = (dndz * n_tot[:,None] - n * dndz_tot[:,None]) / n_tot[:,None] ** 2

            # Direct Heat Transfer (kW/m3)
            q_dir = self.get_heat_transfer_kW_m3()
            q_gas = - q_dir + q_des

            # Temperature Gradient
            dTdz = ((R * T) / (p_gas * cp * v_GS)) * q_gas

            return dndz_tot, dTdz, dydz

        def __RK1_liq__(self, m, T, w, z):

            dmdz, dTdz, dwdz = self.__f_euler_liq__(m, T, w, z)

            h = self.__step_size_liq__(w, dwdz, z)

            m_new = m + dmdz * h
            T_new = T + dTdz * h
            w_new = w + dwdz * h
            z_new = z + h
            return m_new, T_new, w_new, z_new

        def __RK2_liq__(self, m, T, w, z):

            dmdz1, dTdz1, dwdz1 = self.__f_euler_liq__(m, T, w, z)

            h = self.__step_size_liq__(w, dwdz1, z)

            dm1 = h * dmdz1
            dT1 = h * dTdz1
            dw1 = h * dwdz1
            dmdz2, dTdz2, dwdz2 = self.__f_euler_liq__(m + dm1, T + dT1, w + dw1, z + h)
            m_new = m + 0.5 * (dmdz1 + dmdz2) * h
            T_new = T + 0.5 * (dTdz1 + dTdz2) * h
            w_new = w + 0.5 * (dwdz1 + dwdz2) * h
            z_new = z + h
            return m_new, T_new, w_new, z_new

        def __RK4_liq__(self, m, T, w, z):

            dmdz1, dTdz1, dwdz1 = self.__f_euler_liq__(m, T, w, z)

            h = self.__step_size_liq__(w, dwdz1, z)

            dm1 = h * dmdz1
            dT1 = h * dTdz1
            dw1 = h * dwdz1

            dmdz2, dTdz2, dwdz2 = self.__f_euler_liq__(m + 0.5 * dm1, T + 0.5 * dT1, w + 0.5 * dw1, z + 0.5 * h)
            dm2 = h * dmdz2
            dT2 = h * dTdz2
            dw2 = h * dwdz2

            dmdz3, dTdz3, dwdz3 = self.__f_euler_liq__(m + 0.5 * dm2, T + 0.5 * dT2, w + 0.5 * dw2, z + 0.5 * h)
            dm3 = h * dmdz3
            dT3 = h * dTdz3
            dw3 = h * dwdz3

            dmdz4, dTdz4, dwdz4 = self.__f_euler_liq__(m + dm3, T + dT3, w + dw3, z + h)
            dm4 = h * dmdz4
            dT4 = h * dTdz4
            dw4 = h * dwdz4

            m_new = m + (dm1 + 2 * dm2 + 2*dm3 + dm4) / 6
            T_new = T + (dT1 + 2 * dT2 + 2 * dT3 + dT4) / 6
            w_new = w + (dw1 + 2 * dw2 + 2 * dw3 + dw4) / 6
            z_new = z + h
            return m_new, T_new, w_new, z_new

        def __RK1_gas__(self, n, T, y, z):

            dndz, dTdz, dydz = self.__f_euler_gas__(n, T, y, z)

            h = self.__step_size_gas__(y, dydz, z)

            n_new = n + dndz * h
            T_new = T + dTdz * h
            y_new = y + dydz * h
            z_new = z + h
            return n_new, T_new, y_new, z_new

        def __RK2_gas__(self, n, T, y, z):

            dndz1, dTdz1, dydz1 = self.__f_euler_gas__(n, T, y, z)

            h = self.__step_size_gas__(y, dydz1, z)

            dn1 = h * dndz1
            dT1 = h * dTdz1
            dy1 = h * dydz1
            dndz2, dTdz2, dydz2 = self.__f_euler_gas__(n + dn1, T + dT1, y + dy1, z + h)
            n_new = n + 0.5 * (dndz1 + dndz2) * h
            T_new = T + 0.5 * (dTdz1 + dTdz2) * h
            y_new = y + 0.5 * (dydz1 + dydz2) * h
            z_new = z + h
            return n_new, T_new, y_new, z_new

        def __RK4_gas__(self, n, T, y, z):

            dndz1, dTdz1, dydz1 = self.__f_euler_gas__(n, T, y, z)

            h = self.__step_size_gas__(y, dydz1, z)

            dn1 = h * dndz1
            dT1 = h * dTdz1
            dy1 = h * dydz1

            dndz2, dTdz2, dydz2 = self.__f_euler_gas__(n + 0.5 * dn1, T + 0.5 * dT1, y + 0.5 * dy1, z + 0.5 * h)
            dn2 = h * dndz2
            dT2 = h * dTdz2
            dy2 = h * dydz2

            dndz3, dTdz3, dydz3 = self.__f_euler_gas__(n + 0.5 * dn2, T + 0.5 * dT2, y + 0.5 * dy2, z + 0.5 * h)
            dn3 = h * dndz3
            dT3 = h * dTdz3
            dy3 = h * dydz3

            dndz4, dTdz4, dydz4 = self.__f_euler_gas__(n + dn3, T + dT3, y + dy3, z + h)
            dn4 = h * dndz4
            dT4 = h * dTdz4
            dy4 = h * dydz4

            n_new = n + (dn1 + 2 * dn2 + 2 * dn3 + dn4) / 6
            T_new = T + (dT1 + 2 * dT2 + 2 * dT3 + dT4) / 6
            y_new = y + (dy1 + 2 * dy2 + 2 * dy3 + dy4) / 6
            z_new = z + h
            return n_new, T_new, y_new, z_new

        def __step_size_gas__(self, y, dydz, z):
            with np.errstate(divide='ignore', invalid='ignore'):
                h = 0.02 * np.min(np.nan_to_num(y / np.abs(dydz), nan=np.inf, posinf=np.inf, neginf=np.inf), axis=1, keepdims=False)
            h = np.min(h)
            h = np.minimum(h, 0.075)
            h = np.minimum(h, self._profile_position_m[-1] - z)
            return h

        def __step_size_liq__(self, w, dwdz, z):
            with np.errstate(divide='ignore', invalid='ignore'):
                h = 0.05 * np.min(np.nan_to_num(w / np.abs(dwdz), nan=np.inf, posinf=np.inf, neginf=np.inf), axis=1, keepdims=False)
            h = np.min(h)
            h = np.minimum(h, 0.075)
            h = np.minimum(h, z - self._profile_position_m[0])
            h = - h
            return h

        # --------------------------------------------------------------------

        def add_info(self, key, value):
            self.info[key] = value

        def add_mass_transfer_kmol_m3s(self, id, stoch_gas, stoch_liq, rate_kmol_m3s, exothermic_heat_kJ_kmol=None):
            self.mass_transfer_kmol_m3s[id] = {}
            self.mass_transfer_kmol_m3s[id]["Stoch Gas"] = stoch_gas
            self.mass_transfer_kmol_m3s[id]["Stoch Liq"] = stoch_liq
            self.mass_transfer_kmol_m3s[id]["Rate [kmol/m3.s]"] = rate_kmol_m3s
            self.mass_transfer_kmol_m3s[id]["Exothermic Heat [kJ/kmol]"] = exothermic_heat_kJ_kmol
            self.mass_transfer_kmol_m3s[id]["Index"] = self.num_of_mass_transfer
            self.num_of_mass_transfer = self.num_of_mass_transfer + 1

        def add_heat_transfer_kW_m3(self, heat_transfer_kW_m3):
            self.heat_transfer_kW_m3 = heat_transfer_kW_m3

        def add_liquid_holdup_m3_m3(self, liquid_holdup_m3_m3):
            self.liquid_holdup_m3_m3 = liquid_holdup_m3_m3

        def add_pressure_drop_Pa_m(self, presure_drop_Pa_m):
            pass

        # --------------------------------------------------------------------

        def get_superficial_gas_velocity_m_s(self):
            T = self.GasStream.get_gas_temp_K()
            n = self.GasStream.get_gas_flow_kmol_h() / 3600
            p = self.GasStream.get_gas_pressure_bara()
            R = 0.08314
            V = n * R * T / p
            v = V / self._cross_sectional_area_m2
            return v

        def get_superficial_liquid_velocity_m_s(self):
            rho = self.LiquidStream.get_solution_density_kg_m3()
            m = self.LiquidStream.get_solution_flow_kg_h() / 3600
            V = m / rho
            v = V / self._cross_sectional_area_m2
            return v

        def get_liquid_holdup_m3_m3(self):
            return self.liquid_holdup_m3_m3(self)

        def get_pressure_drop_Pa_m(self):
            return self.pressure_drop_Pa_m(self)

        def get_mass_transfer_kmol_m3s(self, id):
            return self.mass_transfer_kmol_m3s[id]["Rate [kmol/m3.s]"](self)

        def get_heat_transfer_kW_m3(self):
            return self.heat_transfer_kW_m3(self)

        # --------------------------------------------------------------------

        def __get_mass_transfer_vectors__(self):

            dT = self.GasStream.get_gas_temp_K() - self.LiquidStream.get_solution_temp_K()
            r_specie_liq = np.zeros(shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.num_of_species), dtype=np.float64)
            r_specie_gas = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.num_of_species), dtype=np.float64)
            r_rxn = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.num_of_mass_transfer), dtype=np.float64)

            #q_lat = np.zeros(shape=(self.GasStream.temp_K.shape[0],), dtype=np.float64)

            #print(self.GasStream.temp_K.shape)
            #print(self.GasStream.flow_kmol_h.shape)
            #print(self.GasStream.specie["CO2"]["Molar Fraction"].shape)
            #print(self.GasStream.get_gas_pressure_bara().shape)
            #print(self.GasStream.info)
            #time.sleep(10)

            for k, id in enumerate(self.mass_transfer_kmol_m3s.keys()):
                r_rxn[:, k] = self.mass_transfer_kmol_m3s[id]["Rate [kmol/m3.s]"](self)
                nu = self.mass_transfer_kmol_m3s[id]["Stoch Liq"]
                for el in nu.keys():
                    i = self.LiquidStream.specie[el]["Index"]
                    r_specie_liq[:, i] = r_specie_liq[:, i] + r_rxn[:, k] * nu[el]
                nu = self.mass_transfer_kmol_m3s[id]["Stoch Gas"]
                for el in nu.keys():
                    i = self.GasStream.specie[el]["Index"]
                    r_specie_gas[:, i] = r_specie_gas[:, i] + r_rxn[:, k] * nu[el]

                # Heat of Absorption
                #q_lat = q_lat + r_rxn[:, k] * self.mass_transfer_kmol_m3s[id]["Exothermic Heat [kJ/kmol]"](self)


            # Heat due to Heating/Cooling of Gas during Absorption
            q_abs = np.zeros(shape=(self.GasStream.temp_K.shape[0],), dtype=np.float64)
            q_des = np.zeros(shape=(self.GasStream.temp_K.shape[0],), dtype=np.float64)
            for i, id in enumerate(self.GasStream.specie.keys()):
                cp = self.GasStream.get_specie_heat_capacity_kJ_kmolK(id)
                q_abs = q_abs + dT * cp * np.maximum(r_specie_gas[:, i], 0)
                q_des = q_des + dT * cp * np.minimum(r_specie_gas[:, i], 0)

            return r_rxn, r_specie_gas, r_specie_liq, q_abs, q_des


    class GasLiquidContactor_PFR_CounterCurrent_Scipy(_GasLiquidContactor_PFR_Scipy ,_Stochiometry, _Serializer, _LiquidEquilibrium, _StreamFunctions):

        def __init__(self, position_m, cross_sectional_area_m2, void_fraction_m3_m3):
            self.info = {}
            self.mass_transfer_kmol_m3s = {}
            self.num_of_mass_transfer = 0
            self.firstscan = True
            self._profile_position_m = position_m
            self._profile_cross_sectional_area_m2 = cross_sectional_area_m2
            self._profile_void_fraction_m3_m3 = void_fraction_m3_m3
            self.num_of_heights = len(position_m)

        def __react__(self, GasStreamIn, LiquidStreamIn, max_epochs=np.inf):

            # Various
            self.num_of_samples = GasStreamIn.temp_K.shape[0]

            # Only Strictly Positive Mass Fractions Allowed
            for id in LiquidStreamIn.specie.keys():
                LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(LiquidStreamIn.specie[id]["Mass Fraction"], 10 ** (-18))
            LiquidStreamIn.normalize_mass_fractions()
            for id in GasStreamIn.specie.keys():
                GasStreamIn.specie[id]["Molar Fraction"] = np.maximum(GasStreamIn.specie[id]["Molar Fraction"], 10 ** (-18))
            GasStreamIn.normalize_molar_fractions()

            if self.firstscan:

                # Extract Molar Masses
                self.M_liq = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_species))
                for i, id in enumerate(LiquidStreamIn.specie.keys()):
                    self.M_liq[:, i] = LiquidStreamIn.get_specie_molar_mass_kg_kmol(id)
                self.M_gas = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_species))
                for i, id in enumerate(GasStreamIn.specie.keys()):
                    self.M_gas[:, i] = GasStreamIn.get_specie_molar_mass_kg_kmol(id)

                # Matrix or Reaction Stochiometry: Used to Calculate Sensitivity Matrices
                self.matrix_liq = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True)
                self.matrix_gas = self.__get_the_matrix__(GasStreamIn=GasStreamIn, gas_rxn_insta=True)

                # Matrix of Reaction Stochiometry: Used to Calculate Heat of Reactions
                self.matrix_liq_heat = self.__get_the_matrix__(GasStreamIn=GasStreamIn,
                                                               LiquidStreamIn=LiquidStreamIn,
                                                               Column=self,
                                                               liq_rxn_insta=True,
                                                               liq_rxn_reversible=True,
                                                               vapor_pressure=False,
                                                               mass_transfer=True,
                                                               gas_rxn_insta=False,
                                                               gas_rxn_reversible=False)

                self.matrix_gas_heat = self.__get_the_matrix__(GasStreamIn=GasStreamIn,
                                                               gas_rxn_insta=True,
                                                               gas_rxn_reversible=True)

                self.equilibrium_liquid = LiquidEquilibrium_Adiabatic()
                self.equilibrium_gas = None

            # Process Liquid and Gas Inlet Streams
            self.LiquidStreamIn = self.equilibrium_liquid.react(LiquidStreamIn)
            self.GasStreamIn = deepcopy(GasStreamIn)

            # Initiate Streams used for Propagation up and Down
            self.LiquidStream = deepcopy(self.LiquidStreamIn)
            self.GasStream = deepcopy(self.GasStreamIn)



            # Memory
            self.LP = [deepcopy(self.LiquidStreamIn)] * self.num_of_samples
            self.GP = [deepcopy(self.GasStreamIn)] * self.num_of_samples

            # Initiate Profiles
            for s in range(self.num_of_samples):
                self.LP[s].position_m = np.flip(self._profile_position_m)
                self.LP[s].flow_kg_h = self.LiquidStreamIn.flow_kg_h[s] * np.ones(shape=(self.num_of_heights,1), dtype=np.float64)
                self.LP[s].temp_K = self.LiquidStreamIn.temp_K[s] * np.ones(shape=(self.num_of_heights,1), dtype=np.float64)
                for id in self.LiquidStreamIn.specie.keys():
                    self.LP[s].specie[id]["Mass Fraction"] = self.LiquidStreamIn.specie[id]["Mass Fraction"][s] * np.ones(shape=(self.num_of_heights,1), dtype=np.float64)
                for id in self.LiquidStreamIn.info.keys():
                    self.LP[s].info[id] = self.LiquidStreamIn.info[id][s] * np.ones(shape=(self.num_of_heights,1), dtype=np.float64)

                self.GP[s].position_m = self._profile_position_m
                self.GP[s].flow_kmol_h = self.GasStreamIn.flow_kmol_h[s] * np.ones(shape=(self.num_of_heights,1), dtype=np.float64)
                self.GP[s].temp_K = self.GasStreamIn.temp_K[s] * np.ones(shape=(self.num_of_heights,1), dtype=np.float64)
                self.GP[s].pressure_bara = self.GasStreamIn.pressure_bara[s] * np.ones(shape=(self.num_of_heights, 1), dtype=np.float64)
                for id in self.GasStreamIn.specie.keys():
                    self.GP[s].specie[id]["Molar Fraction"] = self.GasStreamIn.specie[id]["Molar Fraction"][s] * np.ones(shape=(self.num_of_heights,1), dtype=np.float64)
                for id in self.GasStreamIn.info.keys():
                    self.GP[s].info[id] = self.GasStreamIn.info[id][s] * np.ones(shape=(self.num_of_heights,1), dtype=np.float64)


            for self.s in range(self.num_of_samples):

                # Convergence Stuff...
                converged = False

                while converged == False:

                    # ------------------------------------------------------------------
                    # Integrate LiquidStream Downward
                    # ------------------------------------------------------------------
                    z_start = self._profile_position_m[-1]
                    z_stop = self._profile_position_m[0]
                    m0 = np.array([self.LiquidStreamIn.get_solution_flow_kg_h()[self.s] / 3600])         # shape = (1,)
                    T0 = np.array([self.LiquidStreamIn.get_solution_temp_K()[self.s]])                   # shape = (1,)
                    w0 = self.LiquidStreamIn.__mass_fractions_dic2vec__()[self.s,:]                      # shape = (num_of_species,)
                    Y0 = np.hstack((m0, T0, w0))

                    # Radau
                    # RK23
                    # RK45
                    # LSODA
                    # BDF
                    soln = solve_ivp(fun=self.__f_euler_liq__, t_span=(z_start, z_stop), y0=Y0, method="RK45", dense_output=False)
                    z = soln.t
                    Y = soln.y.T

                    print("Number of Funtion Calls")
                    print(soln.nfev)

                    # Save Liquid Profile
                    self.LP[self.s].position_m = z
                    self.LP[self.s].flow_kg_h = 3600 * Y[:,0][:,None]
                    self.LP[self.s].temp_K = Y[:,1][:,None]
                    for i, id in enumerate(self.LP[self.s].specie.keys()):
                        self.LP[self.s].specie[id]["Mass Fraction"] = Y[:, i + 2][:,None]
                    for i, id in enumerate(self.LP[self.s].info.keys()):
                        self.LP[self.s].info[id] = self.LiquidStreamIn.info[id] * np.ones(shape=z.shape)[:,None]


                    """""""""
                    # ------------------------------------------------------------------
                    # Integrate GasStream Upward
                    # ------------------------------------------------------------------
                    z_start = self._profile_position_m[0]
                    z_stop = self._profile_position_m[-1]
                    n0 = np.array([self.GasStreamIn.get_solution_flow_kmol_h()[s] / 3600])  # shape = (1,)
                    T0 = np.array([self.GasStreamIn.get_solution_temp_K()[s]])  # shape = (1,)
                    y0 = sef.GasStreamIn.__molar_fractions_dic2vec__()[s,:]
                    Y0 = np.hstack((n0, T0, y0))
                    soln = solve_ivp(fun=self.__f_euler_gas__, t_span=(z_start, z_stop), y0=Y0, method="Radau", dense_output=False)
                    """""""""


                    #plt.plot(self.LP[self.s].temp_K[:,0] - 273.15, self.LP[self.s].position_m, marker="o")
                    #plt.grid(True)
                    #plt.show()



                    # Check if Algorithm have Converged
                    converged = True


            return self.GasStream, self.LiquidStream




