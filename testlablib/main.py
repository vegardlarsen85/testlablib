import time
import numpy as np
from copy import deepcopy
import pickle


class Serializer:

    def save(self, instance, filename):
        file = open(filename + ".obj","wb")
        pickle.dump(instance, file)
        file.close()

    def load(self, instance, filename):
        file = open(filename + ".obj", 'rb')
        instance = pickle.load(file)
        file.close()
        return instance


# ---------------------------------------------------------------------------------------


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


class _LiquidEquilibrium:


    def __get_LiquidEquilibrium_Ka__(self, LiquidStreamIn):
        Ka = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
        for j, jd in enumerate(LiquidStreamIn.rxn_insta.keys()):
            Ka[:, j] = LiquidStreamIn.get_rxn_insta_equilibrium_constant_wrt_activities(jd)
        return Ka

    def __get_LiquidEquilibrium_Kw__(self, LiquidStreamIn):

        if LiquidStreamIn.temp_K.ndim == 1:
            Kw = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
            for j, jd in enumerate(LiquidStreamIn.rxn_insta.keys()):
                Kw[:, j] = LiquidStreamIn.get_rxn_insta_equilibrium_constant_wrt_mass_fractions(jd)

        elif LiquidStreamIn.temp_K.ndim == 2:
            Kw = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.temp_K.shape[1], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
            for j, jd in enumerate(LiquidStreamIn.rxn_insta.keys()):
                Kw[:,:, j] = LiquidStreamIn.get_rxn_insta_equilibrium_constant_wrt_mass_fractions(jd)
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

        if LiquidStreamIn.temp_K.ndim == 1:
            dT = 0.05
            T0 = LiquidStreamIn.temp_K
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + dT
            Kw1 = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
            for i, id in enumerate(LiquidStreamIn.rxn_insta.keys()):
                Kw1[:, i] = LiquidStreamIn.get_rxn_insta_equilibrium_constant_wrt_mass_fractions(id)
            T1 = LiquidStreamIn.temp_K
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - dT

        if LiquidStreamIn.temp_K.ndim == 2:
            dT = 0.05
            T0 = LiquidStreamIn.temp_K
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + dT
            Kw1 = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.temp_K.shape[1], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
            for i, id in enumerate(LiquidStreamIn.rxn_insta.keys()):
                Kw1[:, :, i] = LiquidStreamIn.get_rxn_insta_equilibrium_constant_wrt_mass_fractions(id)
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

        if LiquidStreamIn.temp_K.ndim == 1:
            f_rxn_insta = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
            for j, jd in enumerate(LiquidStreamIn.rxn_insta.keys()):
                f_rxn_insta[:, j] = np.log(Kw[:, j])
                for id in LiquidStreamIn.rxn_insta[jd]["Stoch"].keys():
                    nu = LiquidStreamIn.rxn_insta[jd]["Stoch"][id]
                    f_rxn_insta[:, j] = f_rxn_insta[:, j] - nu * np.log(LiquidStreamIn.get_specie_mass_fraction(id=id))

        if LiquidStreamIn.temp_K.ndim == 2:
            f_rxn_insta = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.temp_K.shape[1], LiquidStreamIn.num_of_rxn_insta), dtype=np.float64)
            for j, jd in enumerate(LiquidStreamIn.rxn_insta.keys()):
                f_rxn_insta[:,:, j] = np.log(Kw[:,:, j])
                for id in LiquidStreamIn.rxn_insta[jd]["Stoch"].keys():
                    nu = LiquidStreamIn.rxn_insta[jd]["Stoch"][id]
                    f_rxn_insta[:, :, j] = f_rxn_insta[:, :, j] - nu * np.log(LiquidStreamIn.get_specie_mass_fraction(id=id))

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

    def __get_LiquidEquilibrium_dfdw_rxn_insta__(self, LiquidStreamIn, Kw):
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

        if LiquidStreamIn.temp_K.ndim == 1:
            dfdw_rxn_insta = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_insta, LiquidStreamIn.num_of_species),dtype=np.float64)
            for rxn_insta_i, rxn_insta_id in enumerate(LiquidStreamIn.rxn_insta.keys()):
                for _, specie_id in enumerate(LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys()):
                    specie_i = LiquidStreamIn.specie[specie_id]["Index"]
                    nu = LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_id]
                    dfdw_rxn_insta[:, rxn_insta_i, specie_i] = - nu / LiquidStreamIn.get_specie_mass_fraction(id=specie_id)

        if LiquidStreamIn.temp_K.ndim == 2:
            dfdw_rxn_insta = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.temp_K.shape[1], LiquidStreamIn.num_of_rxn_insta, LiquidStreamIn.num_of_species),dtype=np.float64)
            for rxn_insta_i, rxn_insta_id in enumerate(LiquidStreamIn.rxn_insta.keys()):
                for _, specie_id in enumerate(LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys()):
                    specie_i = LiquidStreamIn.specie[specie_id]["Index"]
                    nu = LiquidStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_id]
                    dfdw_rxn_insta[:,:, rxn_insta_i, specie_i] = - nu / LiquidStreamIn.get_specie_mass_fraction(id=specie_id)
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
        if GasStreamIn.temp_K.ndim == 1:
            Kp = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta), dtype=np.float64)
            for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
                Kp[:, j] = GasStreamIn.get_rxn_insta_equilibrium_constant_wrt_partial_pressures(jd)
        if GasStreamIn.temp_K.ndim == 2:
            Kp = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.temp_K.shape[1], GasStreamIn.num_of_rxn_insta), dtype=np.float64)
            for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
                Kp[:, :, j] = GasStreamIn.get_rxn_insta_equilibrium_constant_wrt_partial_pressures(jd)
        return Kp

    def __get_GasEquilibrium_Ky__(self, GasStreamIn):
        if GasStreamIn.temp_K.ndim == 1:
            Ky = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta), dtype=np.float64)
            for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
                Ky[:, j] = GasStreamIn.get_rxn_insta_equilibrium_constant_wrt_molar_fractions(jd)
        if GasStreamIn.temp_K.ndim == 2:
            Ky = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.temp_K.shape[1], GasStreamIn.num_of_rxn_insta), dtype=np.float64)
            for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
                Ky[:, :, j] = GasStreamIn.get_rxn_insta_equilibrium_constant_wrt_molar_fractions(jd)
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

        if GasStreamIn.temp_K.ndim == 1:
            f_rxn_insta = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta), dtype=np.float64)
            for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
                f_rxn_insta[:, j] = np.log(Kp[:, j])
                for id in GasStreamIn.rxn_insta[jd]["Stoch"].keys():
                    nu = GasStreamIn.rxn_insta[jd]["Stoch"][id]
                    f_rxn_insta[:, j] = f_rxn_insta[:, j] - nu * np.log(GasStreamIn.get_specie_pressure_bara(id=id))

        elif GasStreamIn.temp_K.ndim == 2:
            f_rxn_insta = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.temp_K.shape[1], GasStreamIn.num_of_rxn_insta), dtype=np.float64)
            for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
                f_rxn_insta[:,:, j] = np.log(Kp[:,:, j])
                for id in GasStreamIn.rxn_insta[jd]["Stoch"].keys():
                    nu = GasStreamIn.rxn_insta[jd]["Stoch"][id]
                    f_rxn_insta[:,:, j] = f_rxn_insta[:,:, j] - nu * np.log(GasStreamIn.get_specie_pressure_bara(id=id))

        return f_rxn_insta

    def __get_GasEquilibrium_f_mass_balance__(self, GasStreamIn, b, matrix):
        w = GasStreamIn.__molefrac2massfrac__(GasStreamIn.__molar_fractions_dic2vec__())
        f_mass_balance = np.einsum("cw,sw->sc", matrix["A"], w) - b
        return f_mass_balance

    # Need Verification... dw -> dy and kJ/kg.K -> kJ/kmol.K OK?
    def __get_GasEquilibrium_f_energy_balance__(self, GasStreamInit, GasStreamIn, rxn_insta_exothermic_heat_kJ_kmol, matrix):
        cp = (GasStreamInit.get_gas_heat_capacity_kJ_kmolK() + GasStreamIn.get_gas_heat_capacity_kJ_kmolK()) / 2
        dy = GasStreamIn.__molar_fractions_dic2vec__() - GasStreamInit.__molar_fractions_dic2vec__()
        dT = GasStreamIn.temp_K - GasStreamInit.temp_K
        f_energy_balance = (1 / cp) * np.einsum("sr,sr->s", rxn_insta_exothermic_heat_kJ_kmol, np.einsum("rw,sw->sr", matrix["R+"], dy)) - dT
        f_energy_balance = f_energy_balance[:, None]
        return f_energy_balance

    # ----------------------------------------------------------------------------------------------------

    def __get_GasEquilibrium_dfdy_rxn_insta__(self, GasStreamIn, Ky):

        dfdy_rxn_insta = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta, GasStreamIn.num_of_species),dtype=np.float64)

        for rxn_insta_i, rxn_insta_id in enumerate(GasStreamIn.rxn_insta.keys()):
            for specie_i, specie_id in enumerate(GasStreamIn.specie.keys()):
                if specie_id in list(GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys()):
                    nu = GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_id]
                    if nu < 0:
                        dfdy_rxn_insta[:, rxn_insta_i, specie_i] = Ky[:, rxn_insta_i] * np.abs(nu) * GasStreamIn.get_specie_molar_fraction(id=specie_id) ** (np.abs(nu) - 1)
                        for specie_jd in GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys():
                            nu = GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_jd]
                            if specie_id != specie_jd and nu < 0:
                                dfdy_rxn_insta[:, rxn_insta_i, specie_i] = dfdy_rxn_insta[:, rxn_insta_i, specie_i] * GasStreamIn.get_specie_molar_fraction(id=specie_jd) ** np.abs(nu)
                    else:
                        dfdy_rxn_insta[:, rxn_insta_i, specie_i] = - np.abs(nu) * GasStreamIn.get_specie_molar_fraction(id=specie_id) ** (np.abs(nu) - 1)
                        for specie_jd in GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys():
                            nu = GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_jd]
                            if specie_id != specie_jd and nu > 0:
                                dfdy_rxn_insta[:, rxn_insta_i, specie_i] = dfdy_rxn_insta[:, rxn_insta_i, specie_i] * GasStreamIn.get_specie_molar_fraction(id=specie_jd) ** np.abs(nu)
        return dfdy_rxn_insta

    def __get_GasEquilibrium_dfdT_rxn_insta__(self, GasStreamIn, dKydT):
        dfdT_rxn_insta = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta, 1), dtype=np.float64)
        for j, jd in enumerate(GasStreamIn.rxn_insta.keys()):
            fr = dKydT[:, j]
            for id in GasStreamIn.rxn_insta[jd]["Stoch"].keys():
                nu = GasStreamIn.rxn_insta[jd]["Stoch"][id]
                if nu < 0:
                    fr = fr * GasStreamIn.get_specie_molar_fraction(id=id) ** np.abs(nu)
            dfdT_rxn_insta[:, j, 0] = fr
        return dfdT_rxn_insta

    def __get_GasEquilibrium_dfdy_rxn_insta_log__(self, GasStreamIn):

        if GasStreamIn.temp_K.ndim == 1:
            dfdy_rxn_insta = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.num_of_rxn_insta, GasStreamIn.num_of_species), dtype=np.float64)
            for rxn_insta_i, rxn_insta_id in enumerate(GasStreamIn.rxn_insta.keys()):
                for _, specie_id in enumerate(GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys()):
                    specie_i = GasStreamIn.specie[specie_id]["Index"]
                    nu = GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_id]
                    dfdy_rxn_insta[:, rxn_insta_i, specie_i] = - nu / GasStreamIn.get_specie_molar_fraction(id=specie_id)

        if GasStreamIn.temp_K.ndim == 2:
            dfdy_rxn_insta = np.zeros(shape=(GasStreamIn.temp_K.shape[0], GasStreamIn.temp_K.shape[1], GasStreamIn.num_of_rxn_insta, GasStreamIn.num_of_species), dtype=np.float64)
            for rxn_insta_i, rxn_insta_id in enumerate(GasStreamIn.rxn_insta.keys()):
                for _, specie_id in enumerate(GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"].keys()):
                    specie_i = GasStreamIn.specie[specie_id]["Index"]
                    nu = GasStreamIn.rxn_insta[rxn_insta_id]["Stoch"][specie_id]
                    dfdy_rxn_insta[:,:, rxn_insta_i, specie_i] = - nu / GasStreamIn.get_specie_molar_fraction(id=specie_id)

        return dfdy_rxn_insta

    def __get_GasEquilibrium_dfdT_rxn_insta_log__(self, Ky, dKydT):
        dfdT = dKydT / Ky
        return dfdT[:,:,None]

    # Need Verification... Same as Above...
    def __get_GasEquilibrium_dfdy_energy_balance__(self, GasStreamInit, GasStreamIn, rxn_insta_exothermic_heat_kJ_kmol, matrix):
        cp = (GasStreamInit.get_gas_heat_capacity_kJ_kmolK() + GasStreamIn.get_gas_heat_capacity_kJ_kmolK()) / 2
        dfdy_energy_balance = np.einsum("s,sw->sw", 1 / cp, np.einsum("sr,rw->sw", rxn_insta_exothermic_heat_kJ_kmol, matrix["R+"]))[:,None, :]
        return dfdy_energy_balance

    def __get_GasEquilibrium_dfdT_energy_balance__(self, GasStreamIn):
        return - np.ones(shape=(GasStreamIn.temp_K.shape[0], 1, 1), dtype=np.float64)


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

    def __get_VaporLiquidEquilibrium_dKpwdT__(self, LiquidStreamIn, Kpw):
        dT = 0.05
        T0 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + dT
        Kpw1 = self.__get_VaporLiquidEquilibrium_Kpw__(LiquidStreamIn)
        T1 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - dT
        return (Kpw1 - Kpw) / dT

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

    def __get_LiquidCSTR_rxn_reversible_rates_II_kmol_m3s__(self, LiquidStreamIn):
        r_forward = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible), dtype=np.float64)
        r_backward = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.rxn_reversible.keys()):
            rf, rb = LiquidStreamIn.rxn_reversible[id]["Rate [kmol/m3.s]"](LiquidStreamIn)
            r_forward[:, i] = rf
            r_backward[:, i] = rb
        return r_forward, r_backward

    def __get_LiquidCSTR_drdw__(self, LiquidStreamIn, rates_kmol_m3s):
        drdw = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible, LiquidStreamIn.num_of_species),dtype=np.float64)
        for rxn_i, rxn_id in enumerate(LiquidStreamIn.rxn_reversible.keys()):
            for specie_id in LiquidStreamIn.rxn_reversible[rxn_id]["Dependencies"]:
                specie_i = LiquidStreamIn.specie[specie_id]["Index"]
                w_probe = LiquidStreamIn.get_specie_mass_fraction(id=specie_id)
                dw_probe = np.maximum(0.01 * w_probe, 10 ** (-14))
                LiquidStreamIn.set_specie_mass_fraction(id=specie_id, value=w_probe + dw_probe)
                forward_rate_pertubation, backward_rate_pertubation = LiquidStreamIn.rxn_reversible[rxn_id]["Rate [kmol/m3.s]"](LiquidStreamIn)
                rate_pertubation = forward_rate_pertubation - backward_rate_pertubation
                LiquidStreamIn.set_specie_mass_fraction(id=specie_id, value=w_probe)
                drdw[:, rxn_i, specie_i] = (rate_pertubation - rates_kmol_m3s[:, rxn_i]) / dw_probe
        return drdw

    def __get_LiquidCSTR_drdT__(self, LiquidStreamIn, rates_kmol_m3s):
        drdT = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible, 1), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(LiquidStreamIn.rxn_reversible.keys()):
            dT = 0.05
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + dT
            forward_rate_pertubation, backward_rate_pertubation = LiquidStreamIn.rxn_reversible[rxn_id]["Rate [kmol/m3.s]"](LiquidStreamIn)
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - dT
            rate_pertubation = forward_rate_pertubation - backward_rate_pertubation
            drdT[:, rxn_i, 0] = (rate_pertubation - rates_kmol_m3s[:, rxn_i]) / dT
        return drdT

    # -------------------------------------------------------------------------------------------

    def __get_LiquidCSTR_f_rxn_reversible__(self, LiquidStreamInit, LiquidStreamIn, matrix, rates_kmol_m3s, phi):
        RI = matrix["R+"][LiquidStreamIn.num_of_rxn_insta::, :]
        dw = LiquidStreamIn.__mass_fractions_dic2vec__() - LiquidStreamInit.__mass_fractions_dic2vec__()
        f_rxn_reversible = rates_kmol_m3s - phi[:,None] * np.einsum("rw,sw->sr", RI, dw)
        return f_rxn_reversible

    def __get_LiquidCSTR_f_rxn_reversible_log__(self, forward_rates_kmol_m3s, backward_rates_kmol_m3s):
        f_rxn_reversible_log = np.log(forward_rates_kmol_m3s) - np.log(backward_rates_kmol_m3s)
        return f_rxn_reversible_log

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

        m = LiquidStreamIn.get_solution_flow_kg_h()/3600

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

    def __get_LiquidCSTR_dfdw_rxn_reversible_log__(self, LiquidStreamIn, forward_rates_kmol_m3s, backward_rates_kmol_m3s):
        drfdw = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible, LiquidStreamIn.num_of_species), dtype=np.float64)
        drbdw = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible, LiquidStreamIn.num_of_species), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(LiquidStreamIn.rxn_reversible.keys()):
            for specie_id in LiquidStreamIn.rxn_reversible[rxn_id]["Dependencies"]:
                specie_i = LiquidStreamIn.specie[specie_id]["Index"]
                w_probe = LiquidStreamIn.get_specie_mass_fraction(id=specie_id)
                dw_probe = np.maximum(0.01 * w_probe, 10 ** (-14))
                LiquidStreamIn.set_specie_mass_fraction(id=specie_id, value=w_probe + dw_probe)
                forward_rate_pertubation, backward_rate_pertubation = LiquidStreamIn.rxn_reversible[rxn_id]["Rate [kmol/m3.s]"](LiquidStreamIn)
                LiquidStreamIn.set_specie_mass_fraction(id=specie_id, value=w_probe)
                drfdw[:, rxn_i, specie_i] = (forward_rate_pertubation - forward_rates_kmol_m3s[:, rxn_i]) / dw_probe
                drbdw[:, rxn_i, specie_i] = (backward_rate_pertubation - backward_rates_kmol_m3s[:, rxn_i]) / dw_probe
        dfdw = (1/forward_rates_kmol_m3s[:,:,None]) * drfdw - (1/backward_rates_kmol_m3s[:,:,None]) * drbdw
        # f = ln(rf) - ln(rb)
        # df/dw = (1/rf) drf/dw - (1/rb) drb/dw
        return dfdw

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

    def __get_LiquidCSTR_dfdT_rxn_reversible_log__(self, LiquidStreamIn, forward_rates_kmol_m3s, backward_rates_kmol_m3s):
        drfdT = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible, 1), dtype=np.float64)
        drbdT = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_rxn_reversible, 1), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(LiquidStreamIn.rxn_reversible.keys()):
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + 0.1
            forward_rate_pertubation, backward_rate_pertubation = LiquidStreamIn.rxn_reversible[rxn_id]["Rate [kmol/m3.s]"](LiquidStreamIn)
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - 0.1
            drfdT[:, rxn_i, 0] = (forward_rate_pertubation - forward_rates_kmol_m3s[:, rxn_i]) / 0.1
            drbdT[:, rxn_i, 0] = (backward_rate_pertubation - backward_rates_kmol_m3s[:, rxn_i]) / 0.1
        dfdT = (1/forward_rates_kmol_m3s[:,:,None]) * drfdT - (1/backward_rates_kmol_m3s[:,:,None]) * drbdT
        # f = ln(rf) - ln(rb)
        # df/dT = (1/rf) drf/dT - (1/rb) drb/dT
        return dfdT

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
        GasStreamOut.pressure_bara = np.compress(condition=condition, a=GasStreamIn.pressure_bara, axis=0)
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
        GasStreamOut.pressure_bara[i] = GasStreamIn.pressure_bara
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
        for id in LiquidStreamIn.vapor_pressure_bara.keys():
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"].keys())[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"].keys())[0]
            molar_mass = LiquidStreamIn.get_specie_molar_mass_kg_kmol(id=liq_id)
            charge = LiquidStreamIn.get_specie_charge(id=liq_id)
            GasStreamOut.add_specie(id=gas_id, molar_mass_kg_kmol=molar_mass, charge=charge)
        GasStreamOut.set_gas_flow_kmol_h(value=np.zeros(shape=(LiquidStreamIn.temp_K.shape[0],), dtype=np.float64))
        for id in GasStreamOut.specie.keys():
            GasStreamOut.set_specie_molar_fraction(id=id, value=np.ones(shape=(LiquidStreamIn.temp_K.shape[0],), dtype=np.float64) / GasStreamOut.num_of_species)
        return GasStreamOut

    def __get_Flash_f_vap__(self, LiquidStreamIn, GasStreamIn, Kpw):
        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        f_vap = np.zeros(shape=Kpw.shape, dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            f_vap[:, i] = GasStreamIn.get_specie_molar_fraction(id=gas_id) * GasStreamIn.get_gas_pressure_bara() - LiquidStreamIn.get_specie_mass_fraction(id=liq_id) / Kpw[:,i]
        return f_vap

    def __get_Flash_dfdw_vap__(self, LiquidStreamIn, Kpw):
        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_species = LiquidStreamIn.num_of_species
        dfdw = np.zeros(shape=(num_of_samples, num_of_vap, num_of_species), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            liq_i = LiquidStreamIn.specie[liq_id]["Index"]
            dfdw[:, i, liq_i] = - 1 / Kpw[:, i]
        return dfdw

    def __get_Flash_dfdy_vap__(self, GasStreamIn, LiquidStreamIn):
        num_of_vap = LiquidStreamIn.num_of_vapor_pressure_bara
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_species = GasStreamIn.num_of_species
        dfdy = np.zeros(shape=(num_of_samples, num_of_vap, num_of_species), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            gas_i = GasStreamIn.specie[gas_id]["Index"]
            dfdy[:, i, gas_i] = GasStreamIn.get_gas_pressure_bara()
        return dfdy

    def __get_Flash_dfdT_vap__(self, LiquidStreamIn, Kpw, dKpwdT):
        dfdT = np.zeros(shape=Kpw.shape, dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.vapor_pressure_bara.keys()):
            gas_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Gas"])[0]
            liq_id = list(LiquidStreamIn.vapor_pressure_bara[id]["Stoch Liq"])[0]
            dfdT[:, i] = LiquidStreamIn.get_specie_mass_fraction(id=liq_id) * dKpwdT[:,i] / Kpw[:,i]**2
        return dfdT[:, :, None]

    # -----------------------------------------------------------

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


class _Gas(Serializer):

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

        self.num_of_species = 0
        self.num_of_rxn_insta = 0
        self.num_of_rxn_reversible = 0

        # Values
        self.temp_K = None
        self.pressure_bara = None

        # Functions
        self.viscosity_Pas = None
        self.thermal_conductivity_kW_mK = None
        self.heat_capacity_kJ_kmolK = None
        self.diffusivity_m2_s = None

    def add_info(self, key, value):
        self.info[key] = value

    def add_specie(self, id, molar_mass_kg_kmol, charge):
        self.specie[id] = {}
        self.specie[id]["Molar Mass [kg/kmol]"] = molar_mass_kg_kmol
        self.specie[id]["Charge"] = charge
        self.specie[id]["Molar Fraction"] = None
        self.specie[id]["Index"] = self.num_of_species
        self.num_of_species = self.num_of_species + 1

    def add_rxn_insta(self, id, stoch, equilibrium_constant):
        self.rxn_insta[id] = {}
        self.rxn_insta[id]["Stoch"] = stoch
        self.rxn_insta[id]["K"] = equilibrium_constant
        self.rxn_insta[id]["Index"] = self.num_of_rxn_insta
        self.num_of_rxn_insta = self.num_of_rxn_insta + 1

    def add_rxn_reversible(self, id, stoch, rate_kmol_m3s):
        self.rxn_reversible[id] = {}
        self.rxn_reversible[id]["Stoch"] = stoch
        self.rxn_reversible[id]["Rate [kmol/m3.s]"] = rate_kmol_m3s
        self.rxn_reversible[id]["Index"] = self.num_of_rxn_reversible
        self.num_of_rxn_reversible = self.num_of_rxn_reversible + 1

    # --------------------------------------------------------------------------------------------

    def load_viscosity_Pas(self, function):
        self.viscosity_Pas = function

    def load_thermal_conductivity_kW_mK(self, function):
        self.thermal_conductivity_kW_mK = function

    def load_heat_capacity_kJ_kmolK(self, function):
        self.heat_capacity_kJ_kmolK = function

    def load_diffusivity_m2_s(self, function):
        self.diffusivity_m2_s = function

    # --------------------------------------------------------------------------------------------

    def set_gas_temp_K(self, value):
        self.temp_K = value

    def set_gas_pressure_bara(self, value):
        self.pressure_bara = value

    def set_specie_molar_fraction(self, id, value):
        self.specie[id]["Molar Fraction"] = value

    # --------------------------------------------------------------------------------------------

    def get_gas_temp_K(self):
        return self.temp_K

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

    def get_specie_pct_vol_dry(self, id, water_id):
        y_wet = self.get_specie_molar_fraction(id)
        y_water = self.get_specie_molar_fraction(water_id)
        y_dry = y_wet / (1 - y_water)
        return 100 * y_dry

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

        if self.temp_K.ndim == 1:

            r_specie = np.zeros(shape=(self.temp_K.shape[0], self.num_of_species), dtype=np.float64)
            for id in self.rxn_reversible.keys():
                f, b = self.get_rxn_reversible_rate_kmol_m3s(id)
                r_rxn = f - b
                nu = self.rxn_reversible[id]["Stoch"]
                for el in nu.keys():
                    i = self.specie[el]["Index"]
                    r_specie[:, i] = r_specie[:, i] + r_rxn * nu[el]

        elif self.temp_K.ndim == 2:

            r_specie = np.zeros(shape=(self.temp_K.shape[0], self.temp_K.shape[1], self.num_of_species), dtype=np.float64)
            for id in self.rxn_reversible.keys():
                f, b = self.get_rxn_reversible_rate_kmol_m3s(id)
                r_rxn = f - b
                nu = self.rxn_reversible[id]["Stoch"]
                for el in nu.keys():
                    i = self.specie[el]["Index"]
                    r_specie[:,:, i] = r_specie[:, i] + r_rxn * nu[el]



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

        if self.temp_K.ndim == 1:
            y = np.zeros(shape=(len(self.temp_K), self.num_of_species), dtype=np.float64)
            for id in self.specie.keys():
                i = self.specie[id]["Index"]
                y[:, i] = self.get_specie_molar_fraction(id=id)

        elif self.temp_K.ndim == 2:
            y = np.zeros(shape=(self.temp_K.shape[0], self.temp_K.shape[1], self.num_of_species), dtype=np.float64)
            for id in self.specie.keys():
                i = self.specie[id]["Index"]
                y[:,:, i] = self.get_specie_molar_fraction(id=id)

        return y

    def __molar_fractions_vec2dic__(self, y):

        if y.ndim == 2:
            for i, id in enumerate(self.specie.keys()):
                self.set_specie_molar_fraction(id=id, value=y[:, i])
        if y.ndim == 3:
            for i, id in enumerate(self.specie.keys()):
                self.set_specie_molar_fraction(id=id, value=y[:, :, i])

    def __massfrac2molefrac__(self, w):

        if w.ndim == 2:
            den = 0
            y = 1.0 * w
            for i, id in enumerate(self.specie.keys()):
                den = den + w[:, i] / self.get_specie_molar_mass_kg_kmol(id)
            for i, id in enumerate(self.specie.keys()):
                y[:, i] = (w[:, i] / self.get_specie_molar_mass_kg_kmol(id)) / den

        if w.ndim == 3:
            den = 0
            y = 1.0 * w
            for i, id in enumerate(self.specie.keys()):
                den = den + w[:,:, i] / self.get_specie_molar_mass_kg_kmol(id)
            for i, id in enumerate(self.specie.keys()):
                y[:,:, i] = (w[:,:, i] / self.get_specie_molar_mass_kg_kmol(id)) / den

        return y

    def __molefrac2massfrac__(self, y):

        if y.ndim == 2:
            den = 0
            w = 1.0 * y
            for i, id in enumerate(self.specie.keys()):
                den = den + y[:, i] * self.get_specie_molar_mass_kg_kmol(id)
            for i, id in enumerate(self.specie.keys()):
                w[:, i] = (y[:, i] * self.get_specie_molar_mass_kg_kmol(id)) / den
        if y.ndim == 3:
            den = 0
            w = 1.0 * y
            for i, id in enumerate(self.specie.keys()):
                den = den + y[:,:, i] * self.get_specie_molar_mass_kg_kmol(id)
            for i, id in enumerate(self.specie.keys()):
                w[:,:, i] = (y[:,:, i] * self.get_specie_molar_mass_kg_kmol(id)) / den

        return w

    def __dydw__(self, w):

        if self.temp_K.ndim == 1:
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

        if self.temp_K.ndim == 2:
            dydw = np.zeros(shape=(w.shape[0], w.shape[1], w.shape[2], w.shape[2]), dtype=np.float64)
            M = np.zeros(shape=(w.shape[0], w.shape[1], w.shape[2]), dtype=np.float64)
            for i, id in enumerate(self.specie.keys()):
                M[:,:, i] = self.get_specie_molar_mass_kg_kmol(id=id)
            den = 0
            for i in range(self.num_of_species):
                den = den + w[:,:, i] / M[:,:, i]
            for i in range(self.num_of_species):
                for j in range(self.num_of_species):
                    if i == j:
                        dydw[:,:, i, j] = (1 / M[:,:,i]) / den - (w[:,:,i] / M[:,:,i] ** 2) / den**2
                    else:
                        dydw[:,:, i, j] = - (w[:,:,i] / (M[:,:,i] * M[:,:,j])) / den

        return dydw

    # ------------------------------------------------------------------------------------

    def __exothermic_heat_kJ_kmol_as_vector_rxn_insta__(self):
        q = np.zeros(shape=(self.temp_K.shape[0], self.num_of_rxn_insta),dtype=np.float64)
        for i, id in enumerate(self.rxn_insta.keys()):
            q[:, i] = self.get_rxn_insta_exothermic_heat_kJ_kmol(id)
        return q


class GasStream(_Gas):

    def __init__(self):
        super().__init__()
        self.flow_kmol_h = None

    def set_gas_flow_kmol_h(self, value):
        self.flow_kmol_h = value

    def get_gas_flow_kmol_h(self):
        return self.flow_kmol_h

    def get_gas_flow_kg_h(self):
        m = 0
        for id in self.specie.keys():
            y = self.get_specie_molar_fraction(id=id)
            M = self.get_specie_molar_mass_kg_kmol(id=id)
            m = m + y * M * self.get_gas_flow_kmol_h()
        return m

    def get_specie_flow_kmol_h(self, id):
        return self.get_gas_flow_kmol_h() * self.get_specie_molar_fraction(id=id)

    def get_specie_flow_kg_h(self, id):
        return self.get_specie_flow_kmol_h(id=id) * self.get_specie_molar_mass_kg_kmol(id=id)


def GasStreamSum(streams):
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

    GasStreamOut.normalize_molar_fractions()
    return GasStreamOut


class GasEquilibrium_Isothermal_Isobaric(_Stochiometry, Serializer, _GasEquilibrium):

    def __init__(self):
        self.firstscan = True

    def react(self, GasStreamIn, lr=0.75):

        # Load Matrices defined by Stochiometry of the Reactions
        if self.firstscan:
            self.matrix = self.__get_the_matrix__(GasStreamIn=GasStreamIn, gas_rxn_insta=True)

        # Initiate Outlet Stream
        self.GasStreamIn = deepcopy(GasStreamIn)
        for id in self.GasStreamIn.specie.keys():
            self.GasStreamIn.specie[id]["Molar Fraction"] = np.maximum(self.GasStreamIn.specie[id]["Molar Fraction"], 10**(-18))
        self.GasStreamIn.normalize_molar_fractions()
        self.GasStreamOut = deepcopy(self.GasStreamIn)

        # Initiate
        y0 = self.GasStreamOut.__molar_fractions_dic2vec__()
        w0 = self.GasStreamOut.__molefrac2massfrac__(y0)

        y = 1.0 * y0
        w = 1.0 * w0

        b = np.einsum("cw,sw->sc", self.matrix["A"], w0)
        r = 0

        converged = False
        epoch = 0
        iterations = np.zeros(shape=(self.GasStreamOut.temp_K.shape[0],))

        while converged == False:

            # Equilibrium Constants
            Kp = self.__get_GasEquilibrium_Kp__(self.GasStreamOut)

            # Objective Functions
            f = self.__get_GasEquilibrium_f_rxn_insta_log__(self.GasStreamOut, Kp)

            # Partial Derivatives
            dfdy = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(self.GasStreamOut)
            dydw = self.GasStreamOut.__dydw__(w)
            dfdw = np.einsum("src,scd->srd", dfdy, dydw)
            dfdr = np.einsum("src,cq->srq", dfdw, self.matrix["R"])

            # Newton's Method
            #dr_newton = - np.linalg.solve(dfdr, f)
            dr_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdr), f)

            dw_newton = np.einsum("ij,sj->si", self.matrix["R"], dr_newton)

            # Reducing Step Size (Slightly)
            dr = lr * dr_newton
            dw = np.einsum("ij,sj->si", self.matrix["R"], dr)

            # Backtrack to Ensure only Positive Concentrations
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)

            # Update Mass Fractions...
            r = r + tau * dr
            w = w0 + np.einsum("ij,sj->si", self.matrix["R"], r)
            y = self.GasStreamOut.__massfrac2molefrac__(w)

            # Outlet
            self.GasStreamOut.__molar_fractions_vec2dic__(np.maximum(w, 10**(-18)))

            # Check if Algorithm have Converged
            specie_converged_1 = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
            specie_converged_2 = np.array(np.abs(dw_newton) < 10**(-18), dtype=np.float32)
            specie_converged = specie_converged_1 + specie_converged_2
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = bool(converged) and (epoch > 0)
            converged = converged or (epoch > 1000)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        self.firstscan = False
        return self.GasStreamOut

    def get_heat_dissipation_kW(self):
        h_in = self.GasStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out = self.GasStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h = (h_in + h_out) / 2
        w_in = self.GasStreamIn.__molefrac2massfrac__(self.GasStreamIn.__molar_fractions_dic2vec__())
        w_out = self.GasStreamOut.__molefrac2massfrac__(self.GasStreamOut.__molar_fractions_dic2vec__())
        dw = w_out - w_in
        m = self.GasStreamIn.get_gas_flow_kg_h() / 3600
        Q = m * np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw))
        return Q

    def get_adiabatic_temperature_increase_K(self):
        h_in = self.GasStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out = self.GasStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h = (h_in + h_out) / 2
        cp_in = self.GasStreamIn.get_gas_heat_capacity_kJ_kmolK() * self.GasStreamIn.get_gas_molarity_kmol_m3() / self.GasStreamIn.get_gas_density_kg_m3()
        cp_out = self.GasStreamOut.get_gas_heat_capacity_kJ_kmolK() * self.GasStreamOut.get_gas_molarity_kmol_m3() / self.GasStreamOut.get_gas_density_kg_m3()
        cp = (cp_in + cp_out) / 2
        w_in = self.GasStreamIn.__molefrac2massfrac__(self.GasStreamIn.__molar_fractions_dic2vec__())
        w_out = self.GasStreamOut.__molefrac2massfrac__(self.GasStreamOut.__molar_fractions_dic2vec__())
        dw = w_out - w_in
        dT = np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw)) / cp
        return dT


class GasEquilibrium_Adiabatic_Isobaric(_Stochiometry, Serializer, _GasEquilibrium):

    def __init__(self):
        self.firstscan = True

    def react(self, GasStreamIn, lr=0.75):

        if self.firstscan:
            self.matrix = self.__get_the_matrix__(GasStreamIn=GasStreamIn, gas_rxn_insta=True)
            self.gas_equilibrium_isothermal = GasEquilibrium_Isothermal_Isobaric()

        # Load Inlet Stream
        self.GasStreamIn = deepcopy(GasStreamIn)
        for id in self.GasStreamIn.specie.keys():
            self.GasStreamIn.specie[id]["Molar Fraction"] = np.maximum(self.GasStreamIn.specie[id]["Molar Fraction"], 10**(-18))
        self.GasStreamIn.normalize_molar_fractions()

        # Reactions Exothermic Heat (h0) at Inlet Condition
        Ky = self.__get_GasEquilibrium_Ky__(self.GasStreamIn)
        dKydT = self.__get_GasEquilibrium_dKydT__(self.GasStreamIn, Ky)
        h0 = self.__get_GasEquilibrium_exothermic_heat_kJ_kmol__(self.GasStreamIn, Ky, dKydT)

        # Init
        y0 = self.GasStreamIn.__molar_fractions_dic2vec__()
        w0 = self.GasStreamIn.__molefrac2massfrac__(y0)
        b0 = np.einsum("cw,sw->sc", self.matrix["A"], w0)

        # Initial Guess
        self.GasStreamOut = deepcopy(self.GasStreamIn)
        self.GasStreamOut = self.gas_equilibrium_isothermal.react(self.GasStreamOut, lr)
        self.GasStreamOut.temp_K = self.GasStreamOut.temp_K + self.gas_equilibrium_isothermal.get_adiabatic_temperature_increase_K()
        w = self.GasStreamOut.__molefrac2massfrac__(self.GasStreamOut.__molar_fractions_dic2vec__())
        T = self.GasStreamOut.temp_K
        r = np.einsum("rw,sw->sr", self.matrix["R+"], w - w0)

        # Iterate Until Convergence
        converged = False
        epoch = 0
        iterations = np.zeros(shape=(self.GasStreamIn.temp_K.shape[0],))

        while converged == False:

            Kp = self.__get_GasEquilibrium_Kp__(self.GasStreamOut)
            Ky = self.__get_GasEquilibrium_Ky__(self.GasStreamOut)
            dKydT = self.__get_GasEquilibrium_dKydT__(self.GasStreamOut, Ky)
            h = (self.__get_GasEquilibrium_exothermic_heat_kJ_kmol__(self.GasStreamOut, Ky, dKydT) + h0) / 2

            f_rxn_insta = self.__get_GasEquilibrium_f_rxn_insta_log__(self.GasStreamOut, Kp)
            f_energy_balance = self.__get_GasEquilibrium_f_energy_balance__(self.GasStreamIn, self.GasStreamOut, h, self.matrix)
            f = np.hstack((f_rxn_insta, f_energy_balance))

            dydw = self.GasStreamOut.__dydw__(w)
            dfdy_rxn_insta = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(self.GasStreamOut)
            dfdw_rxn_insta = np.einsum("src,scd->srd", dfdy_rxn_insta, dydw)
            dfdy_energy_balance = self.__get_GasEquilibrium_dfdy_energy_balance__(self.GasStreamIn, self.GasStreamOut, h, self.matrix)
            dfdw_energy_balance = np.einsum("src,scd->srd", dfdy_energy_balance, dydw)

            dfdw = np.concatenate((dfdw_rxn_insta, dfdw_energy_balance), axis=1)
            dfdr = np.einsum("src,cq->srq", dfdw, self.matrix["R"])

            dfdT_rxn_insta = self.__get_GasEquilibrium_dfdT_rxn_insta_log__(Ky, dKydT)
            dfdT_energy_balance = self.__get_GasEquilibrium_dfdT_energy_balance__(self.GasStreamOut)
            dfdT = np.concatenate((dfdT_rxn_insta, dfdT_energy_balance), axis=1)

            dfdX = np.concatenate((dfdr, dfdT), axis=2)
            #dX_newton = - np.linalg.solve(dfdX, f)
            dX_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdX), f)

            dr_newton = dX_newton[:, :self.GasStreamOut.num_of_rxn_insta:]
            dT_newton = dX_newton[:, self.GasStreamOut.num_of_rxn_insta]
            dw_newton = np.einsum("ij,sj->si", self.matrix["R"], dr_newton)

            dr = lr * dr_newton
            dw = lr * dw_newton
            dT = lr * dT_newton

            # Backtrack to Ensure only Positive Concentrations
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)

            # Update Mass Fractions and Temperature...
            r = r + tau * dr
            T = T + tau[:,0] * dT
            w = w0 + np.einsum("ij,sj->si", self.matrix["R"], r)

            # Newton's Method
            self.GasStreamOut.__molar_fractions_vec2dic__(self.GasStreamOut.__massfrac2molefrac__(np.maximum(w, 10**(-18))))
            self.GasStreamOut.temp_K = T

            # Check if Algorithm have Converged
            specie_converged_1 = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
            specie_converged_2 = np.array(np.abs(dw_newton) < 10 ** (-18), dtype=np.float32)
            specie_converged = specie_converged_1 + specie_converged_2
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = bool(converged) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        #print("Epochs: \t" + str(epoch))
        self.firstscan = False
        return self.GasStreamOut


class GasPFR_Isothermal_Isobaric(_Stochiometry, Serializer, _GasEquilibrium, _StreamFunctions):

    def __init__(self, position_m, cross_sectional_area_m2):
        self.info = {}
        self.firstscan = True
        self.position_m = position_m
        self.cross_sectional_area_m2 = cross_sectional_area_m2
        self.equilibrium_gas_inlet = GasEquilibrium_Adiabatic_Isobaric()
        self.equilibrium_gas_outlet = GasEquilibrium_Adiabatic_Isobaric()

    def get_superficial_gas_velocity_m_s(self):
        A = np.interp(x=self.z, xp=self.position_m, fp=self.cross_sectional_area_m2)
        T = self.GasStream.get_gas_temp_K()
        n = self.GasStream.get_gas_flow_kmol_h() / 3600
        p = self.GasStream.get_gas_pressure_bara()
        R = 0.08314
        V = n * R * T / p
        v = V / A
        return v

    def __react__(self, GasStreamIn, step_size_m, lr=0.75):

        # Step Size of Integrator
        self.step_size_m = step_size_m

        # Only Strictly Positive Molar Fractions Allowed
        for id in GasStreamIn.specie.keys():
            GasStreamIn.specie[id]["Molar Fraction"] = np.maximum(GasStreamIn.specie[id]["Molar Fraction"], 10 ** (-18))
        GasStreamIn.normalize_molar_fractions()

        # Process Gas Inlet Stream
        self.GasStreamIn = self.equilibrium_gas_inlet.react(GasStreamIn, lr=lr)

        if self.firstscan:
            self.M_gas = np.zeros(shape=(GasStreamIn.num_of_species))
            for i, id in enumerate(GasStreamIn.specie.keys()):
                self.M_gas[i] = GasStreamIn.get_specie_molar_mass_kg_kmol(id)
            self.matrix_gas = self.__get_the_matrix__(GasStreamIn=GasStreamIn, gas_rxn_insta=True)
            self.firstscan = False

        # Initiate as Inlet
        self.GasStream = self.__broadcast_to_GasProfile__(GasStreamIn=self.GasStreamIn, num_of_heights=1)

        self.GP = deepcopy(self.GasStream)
        self.z = 0
        self.GP.position_m = [self.z]

        while self.z < self.position_m[-1]:

            A = np.interp(x=self.z, xp=self.position_m, fp=self.cross_sectional_area_m2)

            n_gas_tot = self.GasStream.get_gas_flow_kmol_h() / 3600
            y_gas = self.GasStream.__molar_fractions_dic2vec__()
            n_gas = y_gas * n_gas_tot[:, :, None]
            cp_gas = self.GasStream.get_gas_heat_capacity_kJ_kmolK()

            # Calculate Rate-Reactions in Gas Phase
            r_gas = self.GasStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

            # Differential Equations, Equilibrium Reactions not taken into Account
            dn_gas__dz = A * r_gas
            dn_gas_tot__dz = np.sum(dn_gas__dz, axis=2, keepdims=False)

            # Sensitivity in Gas Phase (dndB)
            n_gas_broadcast = np.broadcast_to(array=n_gas[:, :, :, None], shape=(n_gas.shape[0], n_gas.shape[1], n_gas.shape[2], n_gas.shape[2]))
            dy_gas__dn_gas = - np.einsum("zsmn,zs->zsmn", n_gas_broadcast, 1 / n_gas_tot ** 2) + np.einsum("vw,zs->zsvw", np.eye(N=self.GasStream.num_of_species), 1 / n_gas_tot)
            df1__dy_gas = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(self.GasStream)
            df1__dn_gas = np.einsum("zsfw,zswv->zsfv", df1__dy_gas, dy_gas__dn_gas)
            df2__dn_gas = np.broadcast_to(array=self.matrix_gas["A"][None, None, :, :], shape=(
            self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1], self.matrix_gas["A"].shape[0],
            self.matrix_gas["A"].shape[1])).copy()
            df__dn_gas = np.concatenate((df1__dn_gas, df2__dn_gas), axis=2)
            H = np.linalg.inv(df__dn_gas)
            dn_gas__dB = H[:, :, :, self.GasStream.num_of_rxn_insta::]

            # Change in "B-Vector" for Gas Phase
            dB__dz = np.einsum("bw,zsw->zsb", self.matrix_gas["A"], dn_gas__dz)

            # Differential Equations, Equilibrium Reactions taken into Account
            dn_gas__dz = np.einsum("zswb,zsb->zsw", dn_gas__dB, dB__dz)

            """""""""
            Now we have obtained following derivatives.
            1) dn_gas__dz
            """""""""

            # Calculate Step Size
            with np.errstate(divide='ignore', invalid='ignore'):
                dz_max_03 = 0.1 * np.min(np.nan_to_num(n_gas / np.abs(dn_gas__dz), nan=np.inf, posinf=np.inf, neginf=np.inf))
            dz_max_05 = self.position_m[-1] - self.z
            dz = np.min(a=[dz_max_03, dz_max_05, self.step_size_m])

            # Integrate
            dn_gas = dn_gas__dz * dz
            dn_gas_tot = np.sum(dn_gas, axis=2, keepdims=False)
            y_gas_new = (n_gas + dn_gas) / (n_gas_tot + dn_gas_tot)[:, :, None]
            self.GasStream.flow_kmol_h = self.GasStream.flow_kmol_h + 3600 * dn_gas_tot
            self.GasStream.__molar_fractions_vec2dic__(y=y_gas_new)

            self.z = self.z + dz

            # Perform One Single (Isothermal) Equilibrium Step for the sake of Numerical Stability........
            y = self.GasStream.__molar_fractions_dic2vec__()
            w = self.GasStream.__molefrac2massfrac__(y)
            Kp = self.__get_GasEquilibrium_Kp__(self.GasStream)
            f = self.__get_GasEquilibrium_f_rxn_insta_log__(self.GasStream, Kp)
            dfdy = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(self.GasStream)
            dydw = self.GasStream.__dydw__(w)
            dfdw = np.einsum("zsfy,zsyw->zsfw", dfdy, dydw)
            dfdr = np.einsum("zsrc,cq->zsrq", dfdw, self.equilibrium_gas_inlet.matrix["R"])
            dr_newton = - np.linalg.solve(dfdr, f)
            dw_newton = np.einsum("ij,zsj->zsi", self.equilibrium_gas_inlet.matrix["R"], dr_newton)
            dw = lr * dw_newton
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            w = w + dw * tau
            self.GasStream.__molar_fractions_vec2dic__(y=self.GasStream.__massfrac2molefrac__(w))

            # Append to Profile
            self.GP.temp_K = np.vstack((self.GP.temp_K, self.GasStream.temp_K))
            self.GP.flow_kmol_h = np.vstack((self.GP.flow_kmol_h, self.GasStream.flow_kmol_h))
            self.GP.pressure_bara = np.vstack((self.GP.pressure_bara, self.GasStream.pressure_bara))
            for id in self.GasStream.specie.keys():
                self.GP.specie[id]["Molar Fraction"] = np.vstack((self.GP.specie[id]["Molar Fraction"], self.GasStream.specie[id]["Molar Fraction"]))
            for id in self.GasStream.info.keys():
                self.GP.info[id] = np.vstack((self.GP.info[id], self.GasStream.info[id]))
            self.GP.position_m.append(self.z)

        # From List to Numpy Array
        self.GP.position_m = np.array(self.GP.position_m)

        # Return Result
        self.GasStreamOut = self.__get_slice_from_GasProfile__(GasStreamProfile=self.GasStream, height_index=0, GasStreamOut=self.GasStream)
        self.GasStream = deepcopy(self.GP)
        del self.GP
        return self.GasStreamOut


class GasPFR_Adiabatic_Isobaric(_Stochiometry, Serializer, _GasEquilibrium, _StreamFunctions):

    def __init__(self, position_m, cross_sectional_area_m2):
        self.info = {}
        self.firstscan = True
        self.position_m = position_m
        self.cross_sectional_area_m2 = cross_sectional_area_m2
        self.equilibrium_gas_inlet = GasEquilibrium_Adiabatic_Isobaric()
        self.equilibrium_gas_outlet = GasEquilibrium_Adiabatic_Isobaric()

    def get_superficial_gas_velocity_m_s(self):
        A = np.interp(x=self.z, xp=self.position_m, fp=self.cross_sectional_area_m2)
        T = self.GasStream.get_gas_temp_K()
        n = self.GasStream.get_gas_flow_kmol_h() / 3600
        p = self.GasStream.get_gas_pressure_bara()
        R = 0.08314
        V = n * R * T / p
        v = V / A
        return v

    def __react__(self, GasStreamIn, step_size_m, lr=0.75):

        # Step Size of Integrator
        self.step_size_m = step_size_m

        # Only Strictly Positive Molar Fractions Allowed
        for id in GasStreamIn.specie.keys():
            GasStreamIn.specie[id]["Molar Fraction"] = np.maximum(GasStreamIn.specie[id]["Molar Fraction"], 10 ** (-18))
        GasStreamIn.normalize_molar_fractions()

        # Process Gas Inlet Stream
        self.GasStreamIn = self.equilibrium_gas_inlet.react(GasStreamIn, lr=lr)

        if self.firstscan:

            # Extract Molar Masses
            self.M_gas = np.zeros(shape=(GasStreamIn.num_of_species))
            for i, id in enumerate(GasStreamIn.specie.keys()):
                self.M_gas[i] = GasStreamIn.get_specie_molar_mass_kg_kmol(id)

            # Matrix or Reaction Stochiometry: Used to Calculate Sensitivity Matrices
            self.matrix_gas = self.__get_the_matrix__(GasStreamIn=GasStreamIn, gas_rxn_insta=True)

            # Matrix of Reaction Stochiometry: Used to Calculate Heat of Reactions
            self.matrix_gas_heat = self.__get_the_matrix__(GasStreamIn=GasStreamIn,
                                                           gas_rxn_insta=True,
                                                           gas_rxn_reversible=True)
            self.firstscan = False

        # Initiate as Inlet
        self.GasStream = self.__broadcast_to_GasProfile__(GasStreamIn=self.GasStreamIn, num_of_heights=1)

        self.GP = deepcopy(self.GasStream)
        self.z = 0
        self.GP.position_m = [self.z]


        while self.z < self.position_m[-1]:

            A = np.interp(x=self.z, xp=self.position_m, fp=self.cross_sectional_area_m2)

            n_gas_tot = self.GasStream.get_gas_flow_kmol_h() / 3600
            y_gas = self.GasStream.__molar_fractions_dic2vec__()
            n_gas = y_gas * n_gas_tot[:, :, None]
            cp_gas = self.GasStream.get_gas_heat_capacity_kJ_kmolK()
            T_gas = self.GasStream.get_gas_temp_K()

            # Exothermic Heat [kJ/kmol rxn] of the following reactions
            # - Gas Insta Reactions
            # - Gas Reversibe Reactions
            h_gas = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1], self.GasStream.num_of_rxn_insta + self.GasStream.num_of_rxn_reversible),dtype=np.float64)
            for rxn_i, rxn_id in enumerate(self.GasStream.rxn_insta.keys()):
                h_gas[:,:, rxn_i] = self.GasStream.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.GasStream.rxn_reversible.keys()):
                h_gas[:,:,rxn_i + self.GasStream.num_of_rxn_insta] = self.GasStream.get_rxn_reversible_exothermic_heat_kJ_kmol(rxn_id)

            # Calculate Rate-Reactions in Gas Phase
            r_gas = self.GasStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

            # Differential Equations, Equilibrium Reactions not taken into Account
            dn_gas__dz = A * r_gas
            dn_gas_tot__dz = np.sum(dn_gas__dz, axis=2, keepdims=False)

            # Sensitivity in Gas Phase (dndB, dndT)
            n_gas_broadcast = np.broadcast_to(array=n_gas[:, :, :, None], shape=(n_gas.shape[0], n_gas.shape[1], n_gas.shape[2], n_gas.shape[2]))
            dy_gas__dn_gas = - np.einsum("zsmn,zs->zsmn", n_gas_broadcast, 1 / n_gas_tot ** 2) + np.einsum("vw,zs->zsvw", np.eye(N=self.GasStream.num_of_species), 1 / n_gas_tot)
            df1__dy_gas = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(self.GasStream)
            df1__dn_gas = np.einsum("zsfw,zswv->zsfv", df1__dy_gas, dy_gas__dn_gas)
            df2__dn_gas = np.broadcast_to(array=self.matrix_gas["A"][None, None, :, :], shape=(
            self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1], self.matrix_gas["A"].shape[0],
            self.matrix_gas["A"].shape[1])).copy()
            df__dn_gas = np.concatenate((df1__dn_gas, df2__dn_gas), axis=2)
            H = np.linalg.inv(df__dn_gas)
            dn_gas__dB = H[:, :, :, self.GasStream.num_of_rxn_insta::]
            Ky = self.__get_GasEquilibrium_Ky__(self.GasStream)
            dKydT = self.__get_GasEquilibrium_dKydT__(self.GasStream, Ky)
            Hr = H[:, :, :, :self.GasStream.num_of_rxn_insta:]
            dn_gas__dT = - np.einsum("zscr,zsr->zsc", Hr, (1 / Ky) * dKydT)

            # Change in "B-Vector" for Gas Phase
            dB__dz = np.einsum("bw,zsw->zsb", self.matrix_gas["A"], dn_gas__dz)

            # Temperature Gradient
            dn_gas_at_eq__dz = dn_gas__dz
            for _ in range(10):
                q_lat = np.einsum("zsr,zsr->zs", h_gas, np.einsum("rm,zsm->zsr", self.matrix_gas_heat["R+"], dn_gas_at_eq__dz)) / A
                dT_gas__dz = A * (1 / (cp_gas * n_gas_tot)) * (-q_dir + q_des + q_lat)
                dn_gas_at_eq__dz = np.einsum("zswb,zsb->zsw", dn_gas__dB, dB__dz) + np.einsum("zsm,zs->zsm", dn_gas__dT, dT_gas__dz)
            dn_gas__dz = dn_gas_at_eq__dz

            """""""""
            Now we have obtained following derivatives.
            1) dn_gas__dz
            2) dT_gas__dz
            """""""""

            # Calculate Step Size
            with np.errstate(divide='ignore', invalid='ignore'):
                dz_max_03 = 0.1 * np.min(np.nan_to_num(n_gas / np.abs(dn_gas__dz), nan=np.inf, posinf=np.inf, neginf=np.inf))
            with np.errstate(divide='ignore', invalid='ignore'):
                dz_max_04 = np.nan_to_num(5 / np.abs(dT_gas__dz), nan=np.inf, posinf=np.inf, neginf=np.inf)[0,0]
            dz_max_05 = self.position_m[-1] - self.z
            dz = np.min(a=[dz_max_03, dz_max_04, dz_max_05, self.step_size_m])

            # Integrate
            dn_gas = dn_gas__dz * dz
            dT_gas = dT_gas__dz * dz
            dn_gas_tot = np.sum(dn_gas, axis=2, keepdims=False)

            y_gas_new = (n_gas + dn_gas) / (n_gas_tot + dn_gas_tot)[:, :, None]

            self.GasStream.temp_K = self.GasStream.temp_K + dT_gas
            self.GasStream.flow_kmol_h = self.GasStream.flow_kmol_h + 3600 * dn_gas_tot
            self.GasStream.__molar_fractions_vec2dic__(y=y_gas_new)

            self.z = self.z + dz

            # Perform One Single (Isothermal) Equilibrium Step for the sake of Numerical Stability........
            y = self.GasStream.__molar_fractions_dic2vec__()
            w = self.GasStream.__molefrac2massfrac__(y)
            Kp = self.__get_GasEquilibrium_Kp__(self.GasStream)
            f = self.__get_GasEquilibrium_f_rxn_insta_log__(self.GasStream, Kp)
            dfdy = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(self.GasStream)
            dydw = self.GasStream.__dydw__(w)
            dfdw = np.einsum("zsfy,zsyw->zsfw", dfdy, dydw)
            dfdr = np.einsum("zsrc,cq->zsrq", dfdw, self.equilibrium_gas_inlet.matrix["R"])
            dr_newton = - np.linalg.solve(dfdr, f)
            dw_newton = np.einsum("ij,zsj->zsi", self.equilibrium_gas_inlet.matrix["R"], dr_newton)
            dw = lr * dw_newton
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            w = w + dw * tau
            self.GasStream.__molar_fractions_vec2dic__(y=self.GasStream.__massfrac2molefrac__(w))

            # Append to Profile
            self.GP.temp_K = np.vstack((self.GP.temp_K, self.GasStream.temp_K))
            self.GP.flow_kmol_h = np.vstack((self.GP.flow_kmol_h, self.GasStream.flow_kmol_h))
            self.GP.pressure_bara = np.vstack((self.GP.pressure_bara, self.GasStream.pressure_bara))
            for id in self.GasStream.specie.keys():
                self.GP.specie[id]["Molar Fraction"] = np.vstack((self.GP.specie[id]["Molar Fraction"], self.GasStream.specie[id]["Molar Fraction"]))
            for id in self.GasStream.info.keys():
                self.GP.info[id] = np.vstack((self.GP.info[id], self.GasStream.info[id]))
            self.GP.position_m.append(self.z)

        # From List to Numpy Array
        self.GP.position_m = np.array(self.GP.position_m)

        # Return Result
        self.GasStreamOut = self.__get_slice_from_GasProfile__(GasStreamProfile=self.GasStream, height_index=0, GasStreamOut=self.GasStream)
        self.GasStream = deepcopy(self.GP)
        del self.GP
        return self.GasStreamOut


# ---------------------------------------------------------------------------------------


class _Liquid(Serializer):

    def __init__(self, solvent_id):

        super().__init__()

        # Solvent (E.g. H2O) when calulating molality
        self.solvent_id = solvent_id

        # Container for Various Info
        self.info = {}
        self.function = {}

        # Concentration
        self.specie = {}
        self.rxn_insta = {}
        self.rxn_reversible = {}
        self.vapor_pressure_bara = {}

        self.num_of_species = 0
        self.num_of_rxn_insta = 0
        self.num_of_rxn_reversible = 0
        self.num_of_vapor_pressure_bara = 0

        # ------------------------

        # Values
        self.temp_K = None

        # Functions
        self.density_kg_m3 = None
        self.heat_capacity_kJ_kgK = None
        self.viscosity_Pas = None
        self.thermal_conductivity_kW_mK = None
        self.activity_coefficient = None
        self.diffusivity_m2_s = None
        self.surface_tension_N_m = None

    # --------------------------------------------------------------------

    def load_density_kg_m3(self, function):
        self.density_kg_m3 = function

    def load_heat_capacity_kJ_kgK(self, function):
        self.heat_capacity_kJ_kgK = function

    def load_viscosity_Pas(self, function):
        self.viscosity_Pas = function

    def load_diffusivity(self, function):
        self.diffusivity_m2_s = function

    def load_activity_coefficient(self, function):
        self.activity_coefficient = function

    def load_thermal_conductivity_kW_mK(self, function):
        self.thermal_conductivity_kW_mK = function

    def load_surface_tension_N_m(self, function):
        self.surface_tension_N_m = function

    # --------------------------------------------------------------------

    def add_info(self, key, value):
        self.info[key] = value

    def add_function(self, key, function):
        self.function[key] = function

    def add_specie(self, id, molar_mass_kg_kmol, charge):
        self.specie[id] = {}
        self.specie[id]["Molar Mass [kg/kmol]"] = molar_mass_kg_kmol
        self.specie[id]["Charge"] = charge
        self.specie[id]["Mass Fraction"] = None
        self.specie[id]["Index"] = self.num_of_species
        self.num_of_species = self.num_of_species + 1

    def add_rxn_insta(self, id, stoch, unit, equilibrium_constant):
        self.rxn_insta[id] = {}
        self.rxn_insta[id]["Stoch"] = stoch
        self.rxn_insta[id]["Unit"] = unit
        self.rxn_insta[id]["K"] = equilibrium_constant
        self.rxn_insta[id]["Index"] = self.num_of_rxn_insta
        self.num_of_rxn_insta = self.num_of_rxn_insta + 1

    def add_rxn_reversible(self, id, stoch, rate_kmol_m3s, dependencies, exothermic_heat_kJ_kmol):
        self.rxn_reversible[id] = {}
        self.rxn_reversible[id]["Stoch"] = stoch
        self.rxn_reversible[id]["Rate [kmol/m3.s]"] = rate_kmol_m3s
        self.rxn_reversible[id]["Dependencies"] = dependencies
        self.rxn_reversible[id]["Exothermic Heat [kJ/kmol]"] = exothermic_heat_kJ_kmol
        self.rxn_reversible[id]["Index"] = self.num_of_rxn_reversible
        self.num_of_rxn_reversible = self.num_of_rxn_reversible + 1

    def add_vapor_pressure_bara_henry(self, id, gas_id, liq_id, liq_unit, henrys_coefficient):
        self.vapor_pressure_bara[id] = {}
        self.vapor_pressure_bara[id]["Stoch Gas"] = {gas_id: -1}
        self.vapor_pressure_bara[id]["Stoch Liq"] = {liq_id: 1}
        self.vapor_pressure_bara[id]["Unit Liq"] = {liq_id: liq_unit}
        self.vapor_pressure_bara[id]["H"] = henrys_coefficient
        self.vapor_pressure_bara[id]["p0"] = None
        self.vapor_pressure_bara[id]["Index"] = self.num_of_vapor_pressure_bara
        self.num_of_vapor_pressure_bara = self.num_of_vapor_pressure_bara + 1

    def add_vapor_pressure_bara_raoult(self, id, gas_id, liq_id, pure_vapor_pressure_bara):
        self.vapor_pressure_bara[id] = {}
        self.vapor_pressure_bara[id]["Stoch Gas"] = {gas_id: -1}
        self.vapor_pressure_bara[id]["Stoch Liq"] = {liq_id: 1}
        self.vapor_pressure_bara[id]["Unit Liq"] = {liq_id: "x"}
        self.vapor_pressure_bara[id]["H"] = None
        self.vapor_pressure_bara[id]["p0"] = pure_vapor_pressure_bara
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
                    self.specie[id]["Mass Fraction"] = 0 * self.specie[self.solvent_id]["Mass Fraction"]
                else:
                    M = self.get_specie_molar_mass_kg_kmol(id=id)
                    self.specie[id]["Mass Fraction"] = mass_solute[id] / mass_tot

    # --------------------------------------------------------------------

    def copy_mass_fractions(self, LiquidStreamIn):
        for id in LiquidStreamIn.specie.keys():
            self.set_specie_mass_fraction(id=id, value=LiquidStreamIn.get_specie_mass_fraction(id))

    # --------------------------------------------------------------------

    def get_info(self, id):
        return self.info[id]

    def get_function(self, id):
        return self.function[id]

    # --------------------------------------------------------------------

    def get_function_value(self, id):
        return self.function[id](self)

    def get_solution_temp_K(self):
        T = self.temp_K
        return self.temp_K

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
        w_solvent = self.get_specie_mass_fraction(id=self.solvent_id)
        w = self.get_specie_mass_fraction(id=id)
        M = self.get_specie_molar_mass_kg_kmol(id=id)
        m = (1000 / w_solvent) * (w / M)
        return m

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

        c_present = 0
        for r in self.rxn_insta[id]["Unit"].keys():
            c_present = c_present + (self.rxn_insta[id]["Unit"][r] == "c")
        if c_present > 0:
            c = self.get_solution_molarity_kmol_m3()
            rho = self.get_solution_density_kg_m3()

        Kw = self.get_rxn_insta_equilibrium_constant_wrt_concentrations(id=id)
        for r in self.rxn_insta[id]["Unit"].keys():
            power = self.rxn_insta[id]["Stoch"][r]
            if self.rxn_insta[id]["Unit"][r] == "c":
                Kw = Kw * (self.get_specie_molar_mass_kg_kmol(r) / rho) ** power
            elif self.rxn_insta[id]["Unit"][r] == "x":
                Kw = Kw * (den * self.get_specie_molar_mass_kg_kmol(r)) ** power
            elif self.rxn_insta[id]["Unit"][r] == "m":
                Kw = Kw * (self.get_specie_mass_fraction(id=self.solvent_id) * self.get_specie_molar_mass_kg_kmol(r) / 1000) ** power
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

        if self.temp_K.ndim == 1:
            r_specie = np.zeros(shape=(self.temp_K.shape[0], self.num_of_species), dtype=np.float64)
            for id in self.rxn_reversible.keys():
                f, b = self.get_rxn_reversible_rate_kmol_m3s(id)
                r_rxn = f - b
                nu = self.rxn_reversible[id]["Stoch"]
                for el in nu.keys():
                    i = self.specie[el]["Index"]
                    r_specie[:,i] = r_specie[:,i] + r_rxn * nu[el]

        if self.temp_K.ndim == 2:
            r_specie = np.zeros(shape=(self.temp_K.shape[0], self.temp_K.shape[1], self.num_of_species), dtype=np.float64)
            for id in self.rxn_reversible.keys():
                f, b = self.get_rxn_reversible_rate_kmol_m3s(id)
                r_rxn = f - b
                nu = self.rxn_reversible[id]["Stoch"]
                for el in nu.keys():
                    i = self.specie[el]["Index"]
                    r_specie[:,:,i] = r_specie[:,:,i] + r_rxn * nu[el]

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

        if self.temp_K.ndim == 1:
            w = np.zeros(shape=(len(self.temp_K), self.num_of_species), dtype=np.float64)
            for id in self.specie.keys():
                i = self.specie[id]["Index"]
                w[:, i] = self.get_specie_mass_fraction(id=id)

        if self.temp_K.ndim == 2:
            w = np.zeros(shape=(self.temp_K.shape[0], self.temp_K.shape[1], self.num_of_species), dtype=np.float64)
            for id in self.specie.keys():
                i = self.specie[id]["Index"]
                w[:,:, i] = self.get_specie_mass_fraction(id=id)
        return w

    def __mass_fractions_vec2dic__(self, w):

        if w.ndim == 2:
            for i, id in enumerate(self.specie.keys()):
                self.set_specie_mass_fraction(id=id, value=w[:, i])
        if w.ndim == 3:
            for i, id in enumerate(self.specie.keys()):
                self.set_specie_mass_fraction(id=id, value=w[:, :, i])

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


class LiquidStream(_Liquid):

    def __init__(self, solvent_id):
        super().__init__(solvent_id)
        self.flow_kg_h = None

    def set_solution_flow_kg_h(self, value):
        self.flow_kg_h = value.astype(np.float64)

    def get_solution_flow_kg_h(self):
        return self.flow_kg_h

    def get_solution_flow_m3_h(self):
        return self.get_solution_flow_kg_h() / self.get_solution_density_kg_m3()

    def get_solution_flow_kmol_h(self):
        return self.get_solution_molarity_kmol_m3() * self.get_solution_flow_m3_h()

    def get_specie_flow_kg_h(self, id):
        return self.get_specie_mass_fraction(id=id) * self.get_solution_flow_kg_h()

    def get_specie_flow_kmol_h(self, id):
        return self.get_specie_molarity_kmol_m3(id=id) * self.get_solution_flow_m3_h()


def LiquidStreamSum(streams):
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


class LiquidEquilibrium_Isothermal(_Stochiometry, Serializer, _LiquidEquilibrium):

    def __init__(self):
        self.firstscan = True

    def react(self, LiquidStreamIn, lr=0.75):

        # Load Matrices defined by Stochiometry of the Reactions
        if self.firstscan:
            self.matrix = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True)

        # Initiate Outlet Stream
        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamIn.specie.keys():
            self.LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamIn.specie[id]["Mass Fraction"], 10**(-18))
        self.LiquidStreamIn.normalize_mass_fractions()
        self.LiquidStreamOut = deepcopy(self.LiquidStreamIn)

        # Initiate
        w0 = self.LiquidStreamOut.__mass_fractions_dic2vec__()
        w = 1.0 * w0
        b = np.einsum("cw,sw->sc", self.matrix["A"], w0)
        r = 0

        converged = False
        epoch = 0
        iterations = np.zeros(shape=(self.LiquidStreamOut.temp_K.shape[0],))

        while converged == False:

            # Equilibrium Constants
            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamOut)

            # Objective Functions
            f = self.__get_LiquidEquilibrium_f_rxn_insta_log__(self.LiquidStreamOut, Kw)

            # Partial Derivatives
            dfdw = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStreamOut)
            dfdr = np.einsum("src,cq->srq", dfdw, self.matrix["R"])

            # Newton's Method
            dr_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdr), f)

            dw_newton = np.einsum("ij,sj->si", self.matrix["R"], dr_newton)

            # Reducing Step Size (Slightly)
            dr = lr * dr_newton

            dw = np.einsum("ij,sj->si", self.matrix["R"], dr)
            # Backtrack to Ensure only Positive Concentrations
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)

            # Update Mass Fractions...
            r = r + tau * dr
            w = w0 + np.einsum("ij,sj->si", self.matrix["R"], r)

            # Outlet
            self.LiquidStreamOut.__mass_fractions_vec2dic__(np.maximum(w, 10**(-18)))

            # Check if Algorithm have Converged
            specie_converged_1 = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
            specie_converged_2 = np.array(np.abs(dw_newton) < 10**(-18), dtype=np.float32)
            specie_converged = specie_converged_1 + specie_converged_2
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = bool(converged) and (epoch > 0)
            converged = converged or (epoch > 1000)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        self.firstscan = False
        return self.LiquidStreamOut

    def get_heat_dissipation_kW(self):
        h_in = self.LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h = (h_in + h_out) / 2
        dw = self.LiquidStreamOut.__mass_fractions_dic2vec__() - self.LiquidStreamIn.__mass_fractions_dic2vec__()
        m = self.LiquidStreamIn.get_solution_flow_kg_h() / 3600
        Q = m * np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw))
        return Q

    def get_adiabatic_temperature_increase_K(self):
        h_in = self.LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h = (h_in + h_out) / 2
        cp_in = self.LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        cp_out = self.LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        cp = (cp_in + cp_out) / 2
        dw = self.LiquidStreamOut.__mass_fractions_dic2vec__() - self.LiquidStreamIn.__mass_fractions_dic2vec__()
        dT = np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw)) / cp
        return dT

    def get_heat_of_vaporization_kJ_kmol(self, LiquidStreamIn, gas_id, lr=0.75):
        LiquidStreamOut = self.react(LiquidStreamIn, lr=lr)
        LiquidStreamOut2 = deepcopy(LiquidStreamOut)
        LiquidStreamOut2.temp_K = LiquidStreamOut.temp_K + 0.1
        LiquidStreamOut2 = self.react(LiquidStreamOut2, lr=lr)
        p = LiquidStreamOut.get_specie_vapor_pressure_bara(gas_id=gas_id)
        p2 = LiquidStreamOut2.get_specie_vapor_pressure_bara(gas_id=gas_id)
        T = LiquidStreamOut.temp_K
        T2 = LiquidStreamOut2.temp_K
        H = 8.314 * np.log(p / p2) * (1 / T2 - 1 / T) ** (-1)
        return H


class LiquidEquilibrium_Adiabatic(_Stochiometry, Serializer, _LiquidEquilibrium):

    def __init__(self):
        self.firstscan = True

    def react(self, LiquidStreamIn, lr=0.75):

        if self.firstscan:
            self.matrix = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True)
            self.liquid_equilibrium_isothermal = LiquidEquilibrium_Isothermal()

        # Load Inlet Stream
        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamIn.specie.keys():
            self.LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamIn.specie[id]["Mass Fraction"], 10**(-18))
        self.LiquidStreamIn.normalize_mass_fractions()

        # Reactions Exothermic Heat (h0) at Inlet Condition
        Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamIn)
        dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStreamIn, Kw)
        h0 = self.__get_LiquidEquilibrium_exothermic_heat_kJ_kmol__(self.LiquidStreamIn, Kw, dKwdT)

        # Init
        b0 = np.einsum("cw,sw->sc", self.matrix["A"], self.LiquidStreamIn.__mass_fractions_dic2vec__())
        w0 = self.LiquidStreamIn.__mass_fractions_dic2vec__()

        # Initial Guess
        self.LiquidStreamOut = deepcopy(self.LiquidStreamIn)
        self.LiquidStreamOut = self.liquid_equilibrium_isothermal.react(self.LiquidStreamOut, lr)
        self.LiquidStreamOut.temp_K = self.LiquidStreamOut.temp_K + self.liquid_equilibrium_isothermal.get_adiabatic_temperature_increase_K()
        w = self.LiquidStreamOut.__mass_fractions_dic2vec__()
        T = self.LiquidStreamOut.temp_K
        r = np.einsum("rw,sw->sr", self.matrix["R+"], w - w0)

        # Iterate Until Convergence
        converged = False
        epoch = 0
        iterations = np.zeros(shape=(self.LiquidStreamIn.temp_K.shape[0],))

        while converged == False:

            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamOut)
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStreamOut, Kw)
            h = (self.__get_LiquidEquilibrium_exothermic_heat_kJ_kmol__(self.LiquidStreamOut, Kw, dKwdT) + h0) / 2

            f_rxn_insta = self.__get_LiquidEquilibrium_f_rxn_insta_log__(self.LiquidStreamOut, Kw)
            f_energy_balance = self.__get_LiquidEquilibrium_f_energy_balance__(self.LiquidStreamIn, self.LiquidStreamOut, h, self.matrix)
            f = np.hstack((f_rxn_insta, f_energy_balance))

            dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStreamOut)
            dfdw_energy_balance = self.__get_LiquidEquilibrium_dfdw_energy_balance__(self.LiquidStreamIn, self.LiquidStreamOut, h, self.matrix)
            dfdw = np.concatenate((dfdw_rxn_insta, dfdw_energy_balance), axis=1)
            dfdr = np.einsum("src,cq->srq", dfdw, self.matrix["R"])

            dfdT_rxn_insta = self.__get_LiquidEquilibrium_dfdT_rxn_insta__(self.LiquidStreamOut, dKwdT)
            dfdT_energy_balance = self.__get_LiquidEquilibrium_dfdT_energy_balance__(self.LiquidStreamOut)
            dfdT = np.concatenate((dfdT_rxn_insta, dfdT_energy_balance), axis=1)

            dfdX = np.concatenate((dfdr, dfdT), axis=2)


            dX_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdX), f)


            dr_newton = dX_newton[:, :self.LiquidStreamOut.num_of_rxn_insta:]
            dT_newton = dX_newton[:, self.LiquidStreamOut.num_of_rxn_insta]
            dw_newton = np.einsum("ij,sj->si", self.matrix["R"], dr_newton)

            dr = lr * dr_newton
            dw = lr * dw_newton
            dT = lr * dT_newton

            # Backtrack to Ensure only Positive Concentrations
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)

            # Update Mass Fractions and Temperature...
            r = r + tau * dr
            T = T + tau[:,0] * dT
            w = w0 + np.einsum("ij,sj->si", self.matrix["R"], r)

            # Newton's Method
            self.LiquidStreamOut.__mass_fractions_vec2dic__(np.maximum(w, 10 ** (-18)))
            self.LiquidStreamOut.temp_K = T

            # Check if Algorithm have Converged
            specie_converged_1 = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
            specie_converged_2 = np.array(np.abs(dw_newton) < 10 ** (-18), dtype=np.float32)
            specie_converged = specie_converged_1 + specie_converged_2
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = bool(converged) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        #print("Epochs: \t" + str(epoch))
        self.firstscan = False
        return self.LiquidStreamOut


class LiquidEquilibrium_QPFlash(_Stochiometry, Serializer, _LiquidEquilibrium, _GasEquilibrium, _VaporLiquidEquilibrium, _Flash, _StreamFunctions):

    def __init__(self, GasStreamRef=None):
        self.firstscan = True
        self.liquid_equilibrium_adiabatic = LiquidEquilibrium_Adiabatic()
        self.GasStreamRef = GasStreamRef

    def react(self, LiquidStreamIn, heat_kW, pressure_bara, lr, flash_always_occur):

        # Load Inlet Stream
        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamIn.specie.keys():
            self.LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamIn.specie[id]["Mass Fraction"], 10 ** (-18))
        self.LiquidStreamIn.normalize_mass_fractions()

        # Load GasStream
        if self.GasStreamRef is None:
            self.GasStreamIn = self.__get_Flash_GasStream__(self.LiquidStreamIn)
        else:
            self.GasStreamIn = deepcopy(self.GasStreamRef)
            n = self.GasStreamIn.num_of_species
            for id in self.GasStreamIn.specie.keys():
                self.GasStreamIn.set_specie_molar_fraction(id=id, value = np.ones(shape=heat_kW.shape) / n)
        self.GasStreamIn.set_gas_flow_kmol_h(value=10 ** (-12) * self.LiquidStreamIn.get_solution_flow_kmol_h())
        self.GasStreamIn.set_gas_pressure_bara(value=pressure_bara)
        self.GasStreamIn.temp_K = self.LiquidStreamIn.temp_K

        # Initial Guess...
        if self.firstscan:
            self.matrix = self.__get_the_matrix__(GasStreamIn=self.GasStreamIn,
                                                  LiquidStreamIn=self.LiquidStreamIn,
                                                  liq_rxn_insta=True,
                                                  vapor_pressure=True,
                                                  gas_rxn_insta=True)
        
        # Solution is Heated Regardless...
        LiquidStreamHeated = deepcopy(self.LiquidStreamIn)
        GasStreamHeated = deepcopy(self.GasStreamIn)
        
        # Check for Flash....
        if flash_always_occur:
            flashing = (self.LiquidStreamIn.temp_K < np.inf)
        else:
            for i in range(4):
                LiquidStreamHeated.temp_K = LiquidStreamHeated.temp_K + (1 / 4) * (3600 * heat_kW) / (LiquidStreamHeated.get_solution_heat_capacity_kJ_kgK() * LiquidStreamHeated.get_solution_flow_kg_h())
            LiquidStreamHeated = self.liquid_equilibrium_adiabatic.react(LiquidStreamHeated, lr=lr)
            GasStreamHeated.temp_K = LiquidStreamHeated.temp_K
            vapor_pressure_bara = 0
            for id in LiquidStreamHeated.vapor_pressure_bara.keys():
                gas_id = list(LiquidStreamHeated.vapor_pressure_bara[id]["Stoch Gas"].keys())[0]
                vapor_pressure_bara = vapor_pressure_bara + LiquidStreamHeated.get_specie_vapor_pressure_bara(gas_id)
            flashing = (vapor_pressure_bara > pressure_bara)

        # Calculate Equilibrium w/QP Flash
        if np.sum(flashing) > 0:
            GasStreamFlashed, LiquidStreamFlashed = self.__flash__(GasStreamIn=self.__get_compressed_GasStream__(GasStreamIn=self.GasStreamIn, condition=flashing),
                                                                   LiquidStreamIn=self.__get_compressed_LiquidStream__(LiquidStreamIn=self.LiquidStreamIn,condition=flashing),
                                                                   heat_kW=np.compress(condition=flashing, a=heat_kW, axis=0),
                                                                   pressure_bara=np.compress(condition=flashing, a=pressure_bara, axis=0),
                                                                   lr=lr)
            self.LiquidStreamOut = self.__insert_into_LiquidStream__(LiquidStreamIn=LiquidStreamFlashed, condition=flashing, LiquidStreamOut=LiquidStreamHeated)
            self.GasStreamOut = self.__insert_into_GasStream__(GasStreamIn=GasStreamFlashed, condition=flashing, GasStreamOut=GasStreamHeated)
        else:
            self.LiquidStreamOut = deepcopy(LiquidStreamHeated)
            self.GasStreamOut = deepcopy(GasStreamHeated)
        

        self.firstscan = False
        return self.GasStreamOut, self.LiquidStreamOut

    def __flash__(self, GasStreamIn, LiquidStreamIn, heat_kW, pressure_bara, lr):

        num_of_samples = GasStreamIn.temp_K.shape[0]

        LiquidStreamOut = deepcopy(LiquidStreamIn)
        GasStreamOut = deepcopy(GasStreamIn)

        # Inlet
        m_in_tot = LiquidStreamIn.get_solution_flow_kg_h() + GasStreamIn.get_gas_flow_kg_h()
        m_in_liq_tot = LiquidStreamIn.get_solution_flow_kg_h()
        m_in_liq = LiquidStreamIn.get_solution_flow_kg_h()[:, None] * LiquidStreamIn.__mass_fractions_dic2vec__()
        m_in_gas = GasStreamIn.get_gas_flow_kg_h()[:, None] * GasStreamIn.__molefrac2massfrac__(y=GasStreamIn.__molar_fractions_dic2vec__())
        m_in = np.concatenate((m_in_liq, m_in_gas), axis=1)
        n_in_gas_tot = GasStreamIn.get_gas_flow_kmol_h()
        cp_in_liq = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        # cp_in_gas = GasStreamIn.get_gas_heat_capacity_kJ_kmolK()
        T_in_liq = LiquidStreamIn.get_solution_temp_K()
        T_in_gas = GasStreamIn.get_gas_temp_K()

        # Initial Guess
        T = 1.0 * T_in_liq

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
            f_energy = f_energy - heat_kW

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
            dX_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdX), f)

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
            converged = bool(converged) and (epoch > 0)
            epoch = epoch + 1

        #print(epoch)
        return GasStreamOut, LiquidStreamOut


class LiquidEquilibrium_TPFlash(_Stochiometry, Serializer, _LiquidEquilibrium, _GasEquilibrium, _VaporLiquidEquilibrium, _Flash, _StreamFunctions):

    def __init__(self, GasStreamRef=None):
        self.firstscan = True
        self.liquid_equilibrium_isothermal = LiquidEquilibrium_Isothermal()
        self.GasStreamRef = GasStreamRef

    def react(self, LiquidStreamIn, temp_K, pressure_bara, lr, flash_always_occur):

        # Load Inlet Stream
        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamIn.specie.keys():
            self.LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamIn.specie[id]["Mass Fraction"], 10 ** (-18))
        self.LiquidStreamIn.normalize_mass_fractions()

        # Load GasStream
        if self.GasStreamRef is None:
            self.GasStreamIn = self.__get_Flash_GasStream__(self.LiquidStreamIn)
        else:
            self.GasStreamIn = deepcopy(self.GasStreamRef)
            n = self.GasStreamIn.num_of_species
            for id in self.GasStreamIn.specie.keys():
                self.GasStreamIn.set_specie_molar_fraction(id=id, value=np.ones(shape=temp_K.shape) / n)
        self.GasStreamIn.set_gas_flow_kmol_h(value=10 ** (-12) * self.LiquidStreamIn.get_solution_flow_kmol_h())
        self.GasStreamIn.set_gas_pressure_bara(value=pressure_bara)
        self.GasStreamIn.temp_K = self.LiquidStreamIn.temp_K

        # Initial Guess...
        if self.firstscan:
            self.matrix = self.__get_the_matrix__(GasStreamIn=self.GasStreamIn,
                                                  LiquidStreamIn=self.LiquidStreamIn,
                                                  liq_rxn_insta=True,
                                                  vapor_pressure=True,
                                                  gas_rxn_insta=True)

        # Solution is Heated Regardless...
        LiquidStreamHeated = deepcopy(self.LiquidStreamIn)
        GasStreamHeated = deepcopy(self.GasStreamIn)

        # Check for Flash....
        if flash_always_occur:
            flashing = (self.LiquidStreamIn.temp_K < np.inf)
        else:
            LiquidStreamHeated.temp_K = temp_K
            LiquidStreamHeated = self.liquid_equilibrium_isothermal.react(LiquidStreamHeated, lr=lr)
            GasStreamHeated.temp_K = LiquidStreamHeated.temp_K
            vapor_pressure_bara = 0
            for id in LiquidStreamHeated.vapor_pressure_bara.keys():
                gas_id = list(LiquidStreamHeated.vapor_pressure_bara[id]["Stoch Gas"].keys())[0]
                vapor_pressure_bara = vapor_pressure_bara + LiquidStreamHeated.get_specie_vapor_pressure_bara(gas_id)
            flashing = (vapor_pressure_bara > pressure_bara)

        # Calculate Equilibrium w/QP Flash
        if np.sum(flashing) > 0:

            GasStreamFlashed, LiquidStreamFlashed = self.__flash__(GasStreamIn=self.__get_compressed_GasStream__(GasStreamIn=self.GasStreamIn, condition=flashing),
                                                                   LiquidStreamIn=self.__get_compressed_LiquidStream__(LiquidStreamIn=self.LiquidStreamIn, condition=flashing),
                                                                   temp_K=np.compress(condition=flashing, a=temp_K, axis=0),
                                                                   pressure_bara=np.compress(condition=flashing, a=pressure_bara, axis=0),
                                                                   lr=lr)

            self.LiquidStreamOut = self.__insert_into_LiquidStream__(LiquidStreamIn=LiquidStreamFlashed,
                                                                     condition=flashing,
                                                                     LiquidStreamOut=LiquidStreamHeated)

            self.GasStreamOut = self.__insert_into_GasStream__(GasStreamIn=GasStreamFlashed, condition=flashing,
                                                               GasStreamOut=GasStreamHeated)
        else:
            self.LiquidStreamOut = deepcopy(LiquidStreamHeated)
            self.GasStreamOut = deepcopy(GasStreamHeated)

        # Heat Input
        Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamOut)
        Kpw = self.__get_VaporLiquidEquilibrium_Kpw__(self.LiquidStreamOut)
        Kyw = self.__get_VaporLiquidEquilibrium_Kyw__(self.GasStreamOut, self.LiquidStreamOut)
        Ky = self.__get_GasEquilibrium_Kp__(self.GasStreamOut)

        dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStreamOut, Kw)
        dKywdT = self.__get_VaporLiquidEquilibrium_dKywdT__(self.GasStreamOut, self.LiquidStreamOut, Kyw)
        dKydT = self.__get_GasEquilibrium_dKydT__(self.GasStreamOut, Ky)

        w_liq = self.LiquidStreamOut.__mass_fractions_dic2vec__()
        w_gas = self.GasStreamOut.__molefrac2massfrac__(y=self.GasStreamOut.__molar_fractions_dic2vec__())

        m_in_liq = self.LiquidStreamIn.get_solution_flow_kg_h()[:, None] * self.LiquidStreamIn.__mass_fractions_dic2vec__()
        m_in_liq_tot = self.LiquidStreamIn.get_solution_flow_kg_h()

        cp_in_liq = self.LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        cp_in_gas = self.GasStreamIn.get_gas_heat_capacity_kJ_kmolK()

        m_in_gas = self.GasStreamIn.get_gas_flow_kg_h()[:, None] * self.GasStreamIn.__molefrac2massfrac__(y=self.GasStreamIn.__molar_fractions_dic2vec__())
        n_in_gas_tot = self.GasStreamIn.get_gas_flow_kmol_h()

        m_liq_tot = self.LiquidStreamOut.get_solution_flow_kg_h()
        m_gas_tot = self.GasStreamOut.get_gas_flow_kg_h()
        m_liq = m_liq_tot[:, None] * w_liq
        m_gas = m_gas_tot[:, None] * w_gas

        dh_liq = self.__get_LiquidEquilibrium_exothermic_heat_kJ_kmol__(self.LiquidStreamOut, Kw, dKwdT)
        dh_vap = self.__get_VaporLiquidEquilibrium_exothermic_heat_kJ_kmol__(self.LiquidStreamOut, Kyw, dKywdT)
        dh_gas = self.__get_GasEquilibrium_exothermic_heat_kJ_kmol__(self.GasStreamOut, Ky, dKydT)
        dh = np.concatenate((dh_liq, dh_vap, dh_gas), axis=1)
        dm_liq = m_liq - m_in_liq
        dm_gas = m_gas - m_in_gas
        dm = np.concatenate((dm_liq, dm_gas), axis=1)
        self.heat_kW = (m_in_liq_tot / 3600) * cp_in_liq * (temp_K - self.LiquidStreamIn.get_solution_temp_K())
        self.heat_kW = self.heat_kW + (n_in_gas_tot / 3600) * cp_in_gas * (temp_K - self.GasStreamIn.get_gas_temp_K())
        self.heat_kW = self.heat_kW - np.einsum("sm,sm->s", np.einsum("sr,rm->sm", dh, self.matrix["R+"]), dm / 3600)


        self.firstscan = False
        return self.GasStreamOut, self.LiquidStreamOut

    def __flash__(self, GasStreamIn, LiquidStreamIn, temp_K, pressure_bara, lr):

        num_of_samples = GasStreamIn.temp_K.shape[0]

        LiquidStreamOut = deepcopy(LiquidStreamIn)
        GasStreamOut = deepcopy(GasStreamIn)

        # Inlet
        m_in_tot = LiquidStreamIn.get_solution_flow_kg_h() + GasStreamIn.get_gas_flow_kg_h()
        m_in_liq_tot = LiquidStreamIn.get_solution_flow_kg_h()
        m_in_liq = LiquidStreamIn.get_solution_flow_kg_h()[:, None] * LiquidStreamIn.__mass_fractions_dic2vec__()
        m_in_gas = GasStreamIn.get_gas_flow_kg_h()[:, None] * GasStreamIn.__molefrac2massfrac__(y=GasStreamIn.__molar_fractions_dic2vec__())
        m_in = np.concatenate((m_in_liq, m_in_gas), axis=1)
        n_in_gas_tot = GasStreamIn.get_gas_flow_kmol_h()
        cp_in_liq = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        cp_in_gas = GasStreamIn.get_gas_heat_capacity_kJ_kmolK()
        T_in_liq = LiquidStreamIn.get_solution_temp_K()
        T_in_gas = GasStreamIn.get_gas_temp_K()

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
            Kyw = self.__get_VaporLiquidEquilibrium_Kyw__(GasStreamOut, LiquidStreamOut)
            Ky = self.__get_GasEquilibrium_Kp__(GasStreamOut)

            # Equilibrium Constants; Temperature Gradients
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(LiquidStreamOut, Kw)
            dKywdT = self.__get_VaporLiquidEquilibrium_dKywdT__(GasStreamOut, LiquidStreamOut, Kyw)
            dKydT = self.__get_GasEquilibrium_dKydT__(GasStreamOut, Ky)

            # Objective Functions
            f_liq = self.__get_LiquidEquilibrium_f_rxn_insta_log__(LiquidStreamOut, Kw)
            f_vap = self.__get_Flash_f_vap_log__(LiquidStreamOut, GasStreamOut, Kpw)
            f_gas = self.__get_GasEquilibrium_f_rxn_insta_log__(GasStreamOut, Ky)
            f = np.concatenate((f_liq, f_vap, f_gas), axis=1)

            # Partial Derivatives
            dfdw_liq = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(LiquidStreamOut)
            dfdw_vap = self.__get_Flash_dfdw_vap_log__(LiquidStreamOut)
            dfdy_vap = self.__get_Flash_dfdy_vap_loq__(GasStreamOut, LiquidStreamOut)
            dfdy_gas = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(GasStreamOut)

            # Partial derivatives; Mass Flows w.r.t Reaction Extent (z)
            dmdz = np.einsum("s,wr->swr", m_in_tot, self.matrix["R"])
            dmdz_liq = dmdz[:, 0:LiquidStreamOut.num_of_species:1, :]
            dmdz_gas = dmdz[:, LiquidStreamOut.num_of_species::, :]
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

            # Partial Derivatives; Objectives w.r.t. Temperature
            dfdT_liq = self.__get_LiquidEquilibrium_dfdT_rxn_insta_log__(Kw, dKwdT)
            dfdT_vap = self.__get_Flash_dfdT_vap_log__(Kyw, dKywdT)
            dfdT_gas = self.__get_GasEquilibrium_dfdT_rxn_insta_log__(Ky, dKydT)

            # Partial Derivatives; Finalizing where X = (z,T) = (R*dm,T)
            dfdz = np.concatenate((dfdz_liq, dfdz_vap, dfdz_gas), axis=1)

            # Damped Newton's Method
            #dz_newton = - np.linalg.solve(dfdz, f)
            dz_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdz), f)

            dm_newton = m_in_tot[:, None] * np.einsum("ij,sj->si", self.matrix["R"], dz_newton)
            dm = lr * dm_newton

            # Backtrack to Ensure only Positive Concentrations
            m = np.concatenate((m_liq, m_gas), axis=1)
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * m / dm, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (m > 0) + 1.0 * (m <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            tau = np.broadcast_to(tau, shape=(
            num_of_samples, GasStreamOut.num_of_species + LiquidStreamOut.num_of_species)).copy()
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
            GasStreamOut.temp_K = temp_K
            GasStreamOut.normalize_molar_fractions()

            # Update Liquid Stream
            LiquidStreamOut.__mass_fractions_vec2dic__(w=w_liq)
            LiquidStreamOut.flow_kg_h = LiquidStreamOut.flow_kg_h + dm_liq_tot
            LiquidStreamOut.temp_K = temp_K
            LiquidStreamOut.normalize_mass_fractions()

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dm_newton) < 0.005 * np.abs(m), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = bool(converged) and (epoch > 0)
            epoch = epoch + 1


        return GasStreamOut, LiquidStreamOut

    def get_heat_kW(self):
        return self.heat_kW


class LiquidCSTR_Isothermal(_Stochiometry, Serializer, _LiquidEquilibrium, _LiquidCSTR):

    def __init__(self, volume_m3):
        self.firstscan = True
        self.volume_m3 = volume_m3

    def react(self, LiquidStreamIn, lr=0.25):

        # Initiate Outlet Stream
        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamIn.specie.keys():
            self.LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamIn.specie[id]["Mass Fraction"], 10 ** (-18))
        self.LiquidStreamIn.normalize_mass_fractions()

        # Load Matrices defined by Stochiometry of the Reactions
        if self.firstscan:
            self.matrix = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn,
                                                  liq_rxn_insta=True,
                                                  liq_rxn_reversible=True)

        # Quantities from Inlet Stream
        b = np.einsum("cw,sw->sc", self.matrix["A"], self.LiquidStreamIn.__mass_fractions_dic2vec__())
        w0 = self.LiquidStreamIn.__mass_fractions_dic2vec__()
        w = 1.0 * w0
        r = 0

        # Initial Guess
        self.LiquidStreamOut = deepcopy(self.LiquidStreamIn)

        # Reactor Dimension
        phi = self.LiquidStreamIn.get_solution_flow_kg_h() / (self.volume_m3 * 3600)

        #reactors = [np.inf, 1000, 100, 10, 5, 2, 1]
        #reactors = [np.inf, 2**16, 2**15, 2**14,2**13, 2**12, 2**11, 2**10, 2**9, 2**8, 2**7, 2**6, 2**5, 2**4, 2**3, 2**2, 2**1, 1]
        reactors = [np.inf, 1.0]

        for reactor in reactors:

            converged = False
            epoch = 0
            iterations = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0],))

            while converged == False:

                Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamOut)

                forward_rates_kmol_m3s, bacward_rates_kmol_m3s = self.__get_LiquidCSTR_rxn_reversible_rates_II_kmol_m3s__(self.LiquidStreamOut)
                rates_kmol_m3s = forward_rates_kmol_m3s - bacward_rates_kmol_m3s

                if reactor == np.inf:
                    f_rxn_insta = self.__get_LiquidEquilibrium_f_rxn_insta_log__(self.LiquidStreamOut, Kw)
                    f_rxn_reversible = self.__get_LiquidCSTR_f_rxn_reversible_log__(forward_rates_kmol_m3s, bacward_rates_kmol_m3s)
                    f = np.hstack((f_rxn_insta, f_rxn_reversible))
                    dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStreamOut)
                    dfdw_rxn_reversible = self.__get_LiquidCSTR_dfdw_rxn_reversible_log__(self.LiquidStreamOut, forward_rates_kmol_m3s, bacward_rates_kmol_m3s)
                    dfdw = np.concatenate((dfdw_rxn_insta, dfdw_rxn_reversible), axis=1)
                    dfdr = np.einsum("src,cq->srq", dfdw, self.matrix["R"])
                else:
                    f_rxn_insta = self.__get_LiquidEquilibrium_f_rxn_insta__(self.LiquidStreamOut, Kw)
                    f_rxn_reversible = self.__get_LiquidCSTR_f_rxn_reversible__(self.LiquidStreamIn, self.LiquidStreamOut, self.matrix, rates_kmol_m3s, phi / reactor)
                    f = np.hstack((f_rxn_insta, f_rxn_reversible))
                    dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta__(self.LiquidStreamOut, Kw)
                    dfdw_rxn_reversible = self.__get_LiquidCSTR_dfdw_rxn_reversible__(self.LiquidStreamOut, self.matrix, rates_kmol_m3s, phi / reactor)
                    dfdw = np.concatenate((dfdw_rxn_insta, dfdw_rxn_reversible), axis=1)
                    dfdr = np.einsum("src,cq->srq", dfdw, self.matrix["R"])

                #dr_newton = - np.linalg.solve(dfdr, f)
                dr_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdr), f)

                dw_newton = np.einsum("ij,sj->si", self.matrix["R"], dr_newton)

                # Reducing Step Size (Slightly)
                dr = lr * dr_newton
                dw = np.einsum("ij,sj->si", self.matrix["R"], dr)

                # Backtrack to Ensure only Positive Concentrations
                with np.errstate(divide='ignore', invalid='ignore'):
                    tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
                tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
                tau = tau * (w > 0) + 1.0 * (w <= 0)
                tau = np.min(tau, axis=1, keepdims=True)

                # Update Mass Fractions...
                r = r + tau * dr
                w = w0 + np.einsum("ij,sj->si", self.matrix["R"], r)

                # Outlet
                self.LiquidStreamOut.__mass_fractions_vec2dic__(np.maximum(w, 10**(-18)))

                # Check if Algorithm have Converged
                specie_converged_1 = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
                specie_converged_2 = np.array(np.abs(dw_newton) <= 10 ** (-18), dtype=np.float32)
                specie_converged = specie_converged_1 + specie_converged_2
                sample_converged = np.min(specie_converged, axis=1)
                converged = np.min(sample_converged, axis=0)
                converged = bool(converged) and (epoch > 0)
                converged = converged or (epoch > 500)
                iterations = iterations + (1 - sample_converged)
                epoch = epoch + 1

            #print(epoch)


        self.firstscan = False
        return self.LiquidStreamOut

    def get_heat_dissipation_kW(self):
        h_in_1 = self.LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_in_2 = self.LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()
        h_in = np.concatenate((h_in_1, h_in_2), axis=1)

        h_out_1 = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out_2 = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()
        h_out = np.concatenate((h_out_1, h_out_2), axis=1)

        h = (h_in + h_out) / 2
        dw = self.LiquidStreamOut.__mass_fractions_dic2vec__() - self.LiquidStreamIn.__mass_fractions_dic2vec__()
        m = self.LiquidStreamIn.get_solution_flow_kg_h() / 3600
        Q = m * np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw))
        return Q

    def get_adiabatic_temperature_increase_K(self):
        h_in_1 = self.LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_in_2 = self.LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()
        h_in = np.concatenate((h_in_1, h_in_2), axis=1)

        h_out_1 = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out_2 = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()
        h_out = np.concatenate((h_out_1, h_out_2), axis=1)

        h = (h_in + h_out) / 2
        cp_in = self.LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        cp_out = self.LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        cp = (cp_in + cp_out) / 2
        dw = self.LiquidStreamOut.__mass_fractions_dic2vec__() - self.LiquidStreamIn.__mass_fractions_dic2vec__()
        dT = np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw)) / cp
        return dT


class LiquidCSTR_Adiabatic(_Stochiometry, Serializer, _LiquidEquilibrium, _LiquidCSTR):

    def __init__(self, volume_m3):
        self.firstscan = True
        self.volume_m3 = volume_m3

    def react(self, LiquidStreamIn, heat_kW, lr=0.25):

        self.heat_kW = heat_kW

        # Initiate Outlet Stream
        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamIn.specie.keys():
            self.LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamIn.specie[id]["Mass Fraction"], 10 ** (-18))
        self.LiquidStreamIn.normalize_mass_fractions()

        # Load Matrices defined by Stochiometry of the Reactions
        if self.firstscan:
            self.matrix = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn,
                                                  liq_rxn_insta=True,
                                                  liq_rxn_reversible=True)

        # Quantities from Inlet Stream
        b = np.einsum("cw,sw->sc", self.matrix["A"], self.LiquidStreamIn.__mass_fractions_dic2vec__())
        w0 = self.LiquidStreamIn.__mass_fractions_dic2vec__()
        w = 1.0 * w0
        r = 0
        T = self.LiquidStreamIn.temp_K

        # Initial Guess
        self.LiquidStreamOut = deepcopy(self.LiquidStreamIn)

        # Reactor Dimension
        phi = self.LiquidStreamIn.get_solution_flow_kg_h() / (self.volume_m3 * 3600)

        reactors = ["Infinite", "Correct Size"]

        for reactor in reactors:

            converged = False
            epoch = 0
            iterations = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0],))

            while converged == False:

                Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStreamOut)
                dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStreamOut, Kw)

                forward_rates_kmol_m3s, bacward_rates_kmol_m3s = self.__get_LiquidCSTR_rxn_reversible_rates_II_kmol_m3s__(self.LiquidStreamOut)
                rates_kmol_m3s = forward_rates_kmol_m3s - bacward_rates_kmol_m3s

                if reactor == "Infinite":
                    f_rxn_insta = self.__get_LiquidEquilibrium_f_rxn_insta_log__(self.LiquidStreamOut, Kw)
                    f_rxn_reversible = self.__get_LiquidCSTR_f_rxn_reversible_log__(forward_rates_kmol_m3s, bacward_rates_kmol_m3s)
                    f_energy_balance = self.__get_LiquidCSTR_f_energy_balance__(self.LiquidStreamIn, self.LiquidStreamOut, self.matrix, self.heat_kW)
                    f = np.hstack((f_rxn_insta, f_rxn_reversible, f_energy_balance))
                    dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStreamOut)
                    dfdw_rxn_reversible = self.__get_LiquidCSTR_dfdw_rxn_reversible_log__(self.LiquidStreamOut, forward_rates_kmol_m3s, bacward_rates_kmol_m3s)
                    dfdw_energy_balance = self.__get_LiquidCSTR_dfdw_energy_balance__(self.LiquidStreamIn, self.LiquidStreamOut, self.matrix)
                    dfdw = np.concatenate((dfdw_rxn_insta, dfdw_rxn_reversible, dfdw_energy_balance), axis=1)
                    dfdr = np.einsum("src,cq->srq", dfdw, self.matrix["R"])
                    dfdT_rxn_insta = self.__get_LiquidEquilibrium_dfdT_rxn_insta_log__(Kw, dKwdT)
                    dfdT_rxn_reversible = self.__get_LiquidCSTR_dfdT_rxn_reversible_log__(self.LiquidStreamOut, forward_rates_kmol_m3s, bacward_rates_kmol_m3s)
                    dfdT_energy_balance = self.__get_LiquidCSTR_dfdT_energy_balance__(self.LiquidStreamOut)
                    dfdT = np.concatenate((dfdT_rxn_insta, dfdT_rxn_reversible, dfdT_energy_balance), axis=1)
                    dfdX = np.concatenate((dfdr, dfdT), axis=2)

                if reactor == "Correct Size":
                    f_rxn_insta = self.__get_LiquidEquilibrium_f_rxn_insta__(self.LiquidStreamOut, Kw)
                    f_rxn_reversible = self.__get_LiquidCSTR_f_rxn_reversible__(self.LiquidStreamIn, self.LiquidStreamOut, self.matrix, rates_kmol_m3s, phi)
                    f_energy_balance = self.__get_LiquidCSTR_f_energy_balance__(self.LiquidStreamIn, self.LiquidStreamOut, self.matrix, self.heat_kW)
                    f = np.hstack((f_rxn_insta, f_rxn_reversible, f_energy_balance))
                    dfdw_rxn_insta = self.__get_LiquidEquilibrium_dfdw_rxn_insta__(self.LiquidStreamOut, Kw)
                    dfdw_rxn_reversible = self.__get_LiquidCSTR_dfdw_rxn_reversible__(self.LiquidStreamOut, self.matrix, rates_kmol_m3s, phi)
                    dfdw_energy_balance = self.__get_LiquidCSTR_dfdw_energy_balance__(self.LiquidStreamIn, self.LiquidStreamOut, self.matrix)
                    dfdw = np.concatenate((dfdw_rxn_insta, dfdw_rxn_reversible, dfdw_energy_balance), axis=1)
                    dfdr = np.einsum("src,cq->srq", dfdw, self.matrix["R"])
                    dfdT_rxn_insta = self.__get_LiquidEquilibrium_dfdT_rxn_insta__(self.LiquidStreamOut, dKwdT)
                    dfdT_rxn_reversible = self.__get_LiquidCSTR_dfdT_rxn_reversible__(self.LiquidStreamOut, self.matrix, rates_kmol_m3s)
                    dfdT_energy_balance = self.__get_LiquidCSTR_dfdT_energy_balance__(self.LiquidStreamOut)
                    dfdT = np.concatenate((dfdT_rxn_insta, dfdT_rxn_reversible, dfdT_energy_balance), axis=1)
                    dfdX = np.concatenate((dfdr, dfdT), axis=2)

                # Newton's Method
                #dX_newton = - np.linalg.solve(dfdX, f)
                dX_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdX), f)

                dr_newton = dX_newton[:, :self.LiquidStreamOut.num_of_rxn_insta + self.LiquidStreamOut.num_of_rxn_reversible:]
                dT_newton = dX_newton[:, self.LiquidStreamOut.num_of_rxn_insta + self.LiquidStreamOut.num_of_rxn_reversible]
                dw_newton = np.einsum("ij,sj->si", self.matrix["R"], dr_newton)

                dr = lr * dr_newton
                dw = lr * dw_newton
                dT = lr * dT_newton

                # Backtrack to Ensure only Positive Concentrations
                with np.errstate(divide='ignore', invalid='ignore'):
                    tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
                tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
                tau = tau * (w > 0) + 1.0 * (w <= 0)
                tau = np.min(tau, axis=1, keepdims=True)

                # Update Mass Fractions and Temperature...
                r = r + tau * dr
                T = T + tau[:, 0] * dT
                w = w0 + np.einsum("ij,sj->si", self.matrix["R"], r)

                # Newton's Method
                self.LiquidStreamOut.__mass_fractions_vec2dic__(np.maximum(w, 10 ** (-18)))
                self.LiquidStreamOut.temp_K = T

                # Check if Algorithm have Converged
                specie_converged_1 = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)
                specie_converged_2 = np.array(np.abs(dw_newton) <= 10 ** (-18), dtype=np.float32)
                specie_converged = specie_converged_1 + specie_converged_2
                sample_converged = np.min(specie_converged, axis=1)
                converged = np.min(sample_converged, axis=0)
                converged = bool(converged) and (epoch > 0)
                iterations = iterations + (1 - sample_converged)
                epoch = epoch + 1

        self.firstscan = False
        return self.LiquidStreamOut

    def get_heat_dissipation_kW(self):
        h_in_1 = self.LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_in_2 = self.LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()
        h_in = np.concatenate((h_in_1, h_in_2), axis=1)

        h_out_1 = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out_2 = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()
        h_out = np.concatenate((h_out_1, h_out_2), axis=1)

        h = (h_in + h_out) / 2
        dw = self.LiquidStreamOut.__mass_fractions_dic2vec__() - self.LiquidStreamIn.__mass_fractions_dic2vec__()
        m = self.LiquidStreamIn.get_solution_flow_kg_h() / 3600
        Q = m * np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw))
        return Q

    def get_adiabatic_temperature_increase_K(self):
        h_in_1 = self.LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_in_2 = self.LiquidStreamIn.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()
        h_in = np.concatenate((h_in_1, h_in_2), axis=1)

        h_out_1 = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_insta__()
        h_out_2 = self.LiquidStreamOut.__exothermic_heat_kJ_kmol_as_vector_rxn_reversible__()
        h_out = np.concatenate((h_out_1, h_out_2), axis=1)

        h = (h_in + h_out) / 2
        cp_in = self.LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        cp_out = self.LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        cp = (cp_in + cp_out) / 2
        dw = self.LiquidStreamOut.__mass_fractions_dic2vec__() - self.LiquidStreamIn.__mass_fractions_dic2vec__()
        dT = np.einsum("sr,sr->s", h, np.einsum("rw,sw->sr", self.matrix["R+"], dw)) / cp
        return dT


class LiquidCSTR_QPFlash(_Stochiometry, Serializer, _LiquidEquilibrium, _LiquidCSTR, _Flash, _StreamFunctions, _VaporLiquidEquilibrium, _GasEquilibrium):

    def __init__(self, volume_m3):
        self.firstscan = True
        self.volume_m3 = volume_m3
        self.liquid_cstr_adiabatic = LiquidCSTR_Adiabatic(volume_m3=volume_m3)

    def react(self, LiquidStreamIn, heat_kW, pressure_bara, lr, flash_always_occur):

        # Load Inlet Stream
        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamIn.specie.keys():
            self.LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStreamIn.specie[id]["Mass Fraction"], 10 ** (-18))
        self.LiquidStreamIn.normalize_mass_fractions()

        # Load GasStream
        self.GasStreamIn = self.__get_Flash_GasStream__(self.LiquidStreamIn)
        self.GasStreamIn.set_gas_flow_kmol_h(value=10 ** (-12) * self.LiquidStreamIn.get_solution_flow_kmol_h())
        self.GasStreamIn.set_gas_pressure_bara(value=pressure_bara)
        self.GasStreamIn.temp_K = self.LiquidStreamIn.temp_K

        # Initial Guess...
        if self.firstscan:
            self.matrix = self.__get_the_matrix__(GasStreamIn=self.GasStreamIn,
                                                  LiquidStreamIn=self.LiquidStreamIn,
                                                  liq_rxn_insta=True,
                                                  liq_rxn_reversible=True,
                                                  vapor_pressure=True,
                                                  gas_rxn_insta=True)

        # Solution is Heated Regardless...
        LiquidStreamHeated = deepcopy(self.LiquidStreamIn)
        GasStreamHeated = deepcopy(self.GasStreamIn)

        # Check for Flash....
        if flash_always_occur:
            flashing = (self.LiquidStreamIn.temp_K < np.inf)
        else:
            LiquidStreamHeated = self.liquid_cstr_adiabatic.react(self.LiquidStreamIn, heat_kW=heat_kW, lr=lr)
            GasStreamHeated.temp_K = LiquidStreamHeated.temp_K
            vapor_pressure_bara = 0
            for id in LiquidStreamHeated.vapor_pressure_bara.keys():
                gas_id = list(LiquidStreamHeated.vapor_pressure_bara[id]["Stoch Gas"].keys())[0]
                vapor_pressure_bara = vapor_pressure_bara + LiquidStreamHeated.get_specie_vapor_pressure_bara(gas_id)
            flashing = (vapor_pressure_bara > pressure_bara)

        # Calculate Equilibrium w/QP Flash
        if np.sum(flashing) > 0:
            GasStreamFlashed, LiquidStreamFlashed = self.__flash__(GasStreamIn=self.__get_compressed_GasStream__(GasStreamIn=self.GasStreamIn, condition=flashing),
                                                                   LiquidStreamIn=self.__get_compressed_LiquidStream__(LiquidStreamIn=self.LiquidStreamIn, condition=flashing),
                                                                   heat_kW=np.compress(condition=flashing, a=heat_kW, axis=0),
                                                                   pressure_bara=np.compress(condition=flashing, a=pressure_bara, axis=0),
                                                                   lr=lr)

            self.LiquidStreamOut = self.__insert_into_LiquidStream__(LiquidStreamIn=LiquidStreamFlashed,
                                                                     condition=flashing,
                                                                     LiquidStreamOut=LiquidStreamHeated)

            self.GasStreamOut = self.__insert_into_GasStream__(GasStreamIn=GasStreamFlashed, condition=flashing,
                                                               GasStreamOut=GasStreamHeated)

        else:
            self.LiquidStreamOut = deepcopy(LiquidStreamHeated)
            self.GasStreamOut = deepcopy(GasStreamHeated)

        self.firstscan = False
        return self.GasStreamOut, self.LiquidStreamOut

    def __flash__(self, GasStreamIn, LiquidStreamIn, heat_kW, pressure_bara, lr):

        num_of_samples = GasStreamIn.temp_K.shape[0]

        LiquidStreamOut = deepcopy(LiquidStreamIn)
        GasStreamOut = deepcopy(GasStreamIn)

        # Inlet
        m_in_tot = LiquidStreamIn.get_solution_flow_kg_h() + GasStreamIn.get_gas_flow_kg_h()
        m_in_liq_tot = LiquidStreamIn.get_solution_flow_kg_h()
        m_in_liq = LiquidStreamIn.get_solution_flow_kg_h()[:, None] * LiquidStreamIn.__mass_fractions_dic2vec__()
        m_in_gas = GasStreamIn.get_gas_flow_kg_h()[:, None] * GasStreamIn.__molefrac2massfrac__(y=GasStreamIn.__molar_fractions_dic2vec__())
        m_in = np.concatenate((m_in_liq, m_in_gas), axis=1)
        n_in_gas_tot = GasStreamIn.get_gas_flow_kmol_h()
        cp_in_liq = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        # cp_in_gas = GasStreamIn.get_gas_heat_capacity_kJ_kmolK()
        T_in_liq = LiquidStreamIn.get_solution_temp_K()
        T_in_gas = GasStreamIn.get_gas_temp_K()

        # Initial Guess
        T = 1.0 * T_in_liq

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

            # Reaction Rates, Liquid Phase
            forward_rates_kmol_m3s, backward_rates_kmol_m3s = self.__get_LiquidCSTR_rxn_reversible_rates_II_kmol_m3s__(LiquidStreamOut)

            # Equilibrium Constants
            Kw = self.__get_LiquidEquilibrium_Kw__(LiquidStreamOut)
            Kpw = self.__get_VaporLiquidEquilibrium_Kpw__(LiquidStreamOut)
            Kyw = self.__get_VaporLiquidEquilibrium_Kyw__(GasStreamOut, LiquidStreamOut)
            Ky = self.__get_GasEquilibrium_Kp__(GasStreamOut)

            # Equilibrium Constants; Temperature Gradients
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(LiquidStreamOut, Kw)
            dKywdT = self.__get_VaporLiquidEquilibrium_dKywdT__(GasStreamOut, LiquidStreamOut, Kyw)
            dKydT = self.__get_GasEquilibrium_dKydT__(GasStreamOut, Ky)

            # Heat of Reactions
            dh_liq = self.__get_LiquidEquilibrium_exothermic_heat_kJ_kmol__(LiquidStreamOut, Kw, dKwdT)
            dh_liq_2 = np.zeros(shape=(num_of_samples, LiquidStreamOut.num_of_rxn_reversible), dtype=np.float64)
            for i, id in enumerate(LiquidStreamOut.rxn_reversible.keys()):
                dh_liq_2[:, i] = LiquidStreamOut.rxn_reversible[id]["Exothermic Heat [kJ/kmol]"](LiquidStreamOut)
            dh_vap = self.__get_VaporLiquidEquilibrium_exothermic_heat_kJ_kmol__(LiquidStreamOut, Kyw, dKywdT)
            dh_gas = self.__get_GasEquilibrium_exothermic_heat_kJ_kmol__(GasStreamOut, Ky, dKydT)
            dh = np.concatenate((dh_liq, dh_liq_2, dh_vap, dh_gas), axis=1)

            # Objective Functions; Equilibrium Constraints
            f_liq = self.__get_LiquidEquilibrium_f_rxn_insta_log__(LiquidStreamOut, Kw)
            f_liq_2 = self.__get_LiquidCSTR_f_rxn_reversible_log__(forward_rates_kmol_m3s, backward_rates_kmol_m3s)
            f_vap = self.__get_Flash_f_vap_log__(LiquidStreamOut, GasStreamOut, Kpw)
            f_gas = self.__get_GasEquilibrium_f_rxn_insta_log__(GasStreamOut, Ky)

            # Objective Functions; Energy Balance
            dm_liq = m_liq - m_in_liq
            dm_gas = m_gas - m_in_gas
            dm = np.concatenate((dm_liq, dm_gas), axis=1)
            f_energy = (m_in_liq_tot / 3600) * cp_in_liq * (T - T_in_liq) - np.einsum("sm,sm->s", np.einsum("sr,rm->sm", dh, self.matrix["R+"]), dm / 3600) - heat_kW

            # Objective Functions
            f = np.concatenate((f_liq, f_liq_2, f_vap, f_gas, f_energy[:, None]), axis=1)

            # Partial Derivatives
            dfdw_liq = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(LiquidStreamOut)
            dfdw_liq_2 = self.__get_LiquidCSTR_dfdw_rxn_reversible_log__(LiquidStreamOut, forward_rates_kmol_m3s, backward_rates_kmol_m3s)
            dfdw_vap = self.__get_Flash_dfdw_vap_log__(LiquidStreamOut)
            dfdy_vap = self.__get_Flash_dfdy_vap_loq__(GasStreamOut, LiquidStreamOut)
            dfdy_gas = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(GasStreamOut)

            # Partial derivatives; Mass Flows w.r.t Reaction Extent (z)
            dmdz = np.einsum("s,wr->swr", m_in_tot, self.matrix["R"])
            dmdz_liq = dmdz[:, 0:LiquidStreamOut.num_of_species:1, :]
            dmdz_gas = dmdz[:, LiquidStreamOut.num_of_species::, :]
            dmdz_liq_tot = np.sum(dmdz_liq, axis=1, keepdims=False)
            dmdz_gas_tot = np.sum(dmdz_gas, axis=1, keepdims=False)

            # Partial Derivatives; Quotient Rule
            dwdz_liq = np.einsum("s,smz->smz", m_liq_tot ** (-2), np.einsum("s,smz->smz", m_liq_tot, dmdz_liq) - np.einsum("sm,sz->smz", m_liq, dmdz_liq_tot))
            dwdz_gas = np.einsum("s,smz->smz", m_gas_tot ** (-2), np.einsum("s,smz->smz", m_gas_tot, dmdz_gas) - np.einsum("sm,sz->smz", m_gas, dmdz_gas_tot))

            # Partial Derivative; Gas Phase - Molar Fraction vs. Mass Fraction
            dydw_gas = GasStreamOut.__dydw__(w_gas)

            # Partial Derivatives; Objectives w.r.t. Reaction Extent
            dfdz_liq = np.einsum("slw,swz->slz", dfdw_liq, dwdz_liq)
            dfdz_liq_2 = np.einsum("slw,swz->slz", dfdw_liq_2, dwdz_liq)
            dfdz_vap = np.einsum("slw,swz->slz", dfdw_vap, dwdz_liq) + np.einsum("slw,swz->slz", np.einsum("sly,syw->slw", dfdy_vap, dydw_gas), dwdz_gas)
            dfdz_gas = np.einsum("slw,swz->slz", np.einsum("sly,syw->slw", dfdy_gas, dydw_gas), dwdz_gas)
            dfdz_energy = - np.einsum("s,sr->sr", m_in_tot / 3600, dh)[:, None, :]

            # Partial Derivatives; Objectives w.r.t. Temperature
            dfdT_liq = self.__get_LiquidEquilibrium_dfdT_rxn_insta_log__(Kw, dKwdT)
            dfdT_liq_2 = self.__get_LiquidCSTR_dfdT_rxn_reversible_log__(LiquidStreamOut, forward_rates_kmol_m3s, backward_rates_kmol_m3s)
            dfdT_vap = self.__get_Flash_dfdT_vap_log__(Kyw, dKywdT)
            dfdT_gas = self.__get_GasEquilibrium_dfdT_rxn_insta_log__(Ky, dKydT)
            dfdT_energy = (m_in_liq_tot / 3600) * cp_in_liq

            # Partial Derivatives; Finalizing where X = (z,T) = (R*dm,T)
            dfdz = np.concatenate((dfdz_liq, dfdz_liq_2, dfdz_vap, dfdz_gas, dfdz_energy), axis=1)
            dfdT = np.concatenate((dfdT_liq, dfdT_liq_2, dfdT_vap, dfdT_gas, dfdT_energy[:, None, None]), axis=1)
            dfdX = np.concatenate((dfdz, dfdT), axis=2)

            # Damped Newton's Method
            #dX_newton = - np.linalg.solve(dfdX, f)
            dX_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdX), f)

            dz_newton = dX_newton[:, 0:dX_newton.shape[1] - 1:]
            dT_newton = dX_newton[:, dX_newton.shape[1] - 1]
            dm_newton = m_in_tot[:, None] * np.einsum("ij,sj->si", self.matrix["R"], dz_newton)
            dm = lr * dm_newton
            dT = lr * dT_newton

            # Backtrack to Ensure only Positive Concentrations
            m = np.concatenate((m_liq, m_gas), axis=1)
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * m / dm, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (m > 0) + 1.0 * (m <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            tau = np.broadcast_to(tau, shape=(
            num_of_samples, GasStreamOut.num_of_species + LiquidStreamOut.num_of_species)).copy()
            dm = dm * tau
            dT = dT * tau[:, 0]

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
            converged = bool(converged) and (epoch > 0)
            epoch = epoch + 1

        return GasStreamOut, LiquidStreamOut


class LiquidPFR_Isothermal(_Stochiometry, Serializer, _LiquidEquilibrium, _StreamFunctions):

    def __init__(self, position_m, cross_sectional_area_m2):

        self.info = {}

        self.firstscan = True
        self.position_m = position_m
        self.cross_sectional_area_m2 = cross_sectional_area_m2

        self.equilibrium_liquid_inlet = LiquidEquilibrium_Isothermal()
        self.equilibrium_liquid_outlet = LiquidEquilibrium_Isothermal()

    def react(self, LiquidStreamIn, step_size_m, lr=0.75):

        # Step Size of Integrator
        self.step_size_m = step_size_m

        # Only Strictly Positive Mass Fractions Allowed
        for id in LiquidStreamIn.specie.keys():
            LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(LiquidStreamIn.specie[id]["Mass Fraction"], 10 ** (-18))
        LiquidStreamIn.normalize_mass_fractions()

        # Process Liquid Inlet Stream
        self.LiquidStreamIn = self.equilibrium_liquid_inlet.react(LiquidStreamIn, lr=lr)

        if self.firstscan:
            self.M_liq = np.zeros(shape=(LiquidStreamIn.num_of_species,))
            for i, id in enumerate(LiquidStreamIn.specie.keys()):
                self.M_liq[i] = LiquidStreamIn.get_specie_molar_mass_kg_kmol(id)
            self.matrix_liq = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True)
            self.firstscan = False

        # Initiate as Inlet
        self.LiquidStream = self.__broadcast_to_LiquidProfile__(LiquidStreamIn=self.LiquidStreamIn, num_of_heights=1)
        self.LP = deepcopy(self.LiquidStream)
        self.z = 0
        self.LP.position_m = [self.z]

        while self.z < self.position_m[-1]:

            A = np.interp(x=self.z, xp=self.position_m, fp=self.cross_sectional_area_m2)

            m_liq_tot = self.LiquidStream.get_solution_flow_kg_h() / 3600
            w_liq = self.LiquidStream.__mass_fractions_dic2vec__()
            m_liq = w_liq * m_liq_tot[:, :, None]
            cp_liq = self.LiquidStream.get_solution_heat_capacity_kJ_kgK()

            # Calculate Rate-Reactions in Liquid Phase
            # [r] = [kmol/m3.s]
            r_liq = self.LiquidStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

            # Differential Equations (dm/dz)
            dm_liq__dz = A * self.M_liq[None, None, :] * r_liq
            dm_liq_tot__dz = A * np.sum(dm_liq__dz, axis=2, keepdims=False)

            # Sensitivity in Liquid Phase (dmdB, dmdT)
            m_liq_broadcast = np.broadcast_to(array=m_liq[:, :, :, None], shape=(m_liq.shape[0], m_liq.shape[1], m_liq.shape[2], m_liq.shape[2]))
            dw_liq__dm_liq = - np.einsum("zsmn,zs->zsmn", m_liq_broadcast, 1 / m_liq_tot ** 2) + np.einsum("vw,zs->zsvw", np.eye(N=self.LiquidStream.num_of_species), 1 / m_liq_tot)
            df1__dw_liq = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStream)
            df1__dm_liq = np.einsum("zsfw,zswv->zsfv", df1__dw_liq, dw_liq__dm_liq)
            df2__dm_liq = np.broadcast_to(array=self.matrix_liq["A"][None, None, :, :], shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.temp_K.shape[1], self.matrix_liq["A"].shape[0], self.matrix_liq["A"].shape[1])).copy()
            df__dm_liq = np.concatenate((df1__dm_liq, df2__dm_liq), axis=2)
            H = np.linalg.inv(df__dm_liq)
            dm_liq__dB = H[:, :, :, self.LiquidStream.num_of_rxn_insta::]

            # Change in "B-Vector" for Liquid Phase
            dB__dz = np.einsum("bw,zsw->zsb", self.matrix_liq["A"], dm_liq__dz)

            # Change in Flow wrt Position after Equilibrium Constraints is taken into Consideration
            dm_liq__dz = np.einsum("zswb,zsb->zsw", dm_liq__dB, dB__dz)

            # Calculate Step Size
            with np.errstate(divide='ignore', invalid='ignore'):
                dz_max_01 = 0.1 * np.min(np.nan_to_num(m_liq / np.abs(dm_liq__dz), nan=np.inf, posinf=np.inf, neginf=np.inf))
            dz_max_03 = self.position_m[-1] - self.z
            dz = np.min(a=[dz_max_01, dz_max_03, self.step_size_m])

            # Integrate
            dm_liq = dm_liq__dz * dz
            dm_liq_tot = np.sum(dm_liq, axis=2, keepdims=False)
            w_liq_new = (m_liq + dm_liq) / (m_liq_tot + dm_liq_tot)[:, :, None]

            self.LiquidStream.flow_kg_h = self.LiquidStream.flow_kg_h + 3600 * dm_liq_tot
            self.LiquidStream.__mass_fractions_vec2dic__(w=w_liq_new)

            self.z = self.z + dz

            # Append to Profile
            self.LP.temp_K = np.vstack((self.LP.temp_K, self.LiquidStream.temp_K))
            self.LP.flow_kg_h = np.vstack((self.LP.flow_kg_h, self.LiquidStream.flow_kg_h))
            for id in self.LiquidStream.specie.keys():
                self.LP.specie[id]["Mass Fraction"] = np.vstack((self.LP.specie[id]["Mass Fraction"], self.LiquidStream.specie[id]["Mass Fraction"]))
            for id in self.LiquidStream.info.keys():
                self.LP.info[id] = np.vstack((self.LP.info[id], self.LiquidStream.info[id]))
            self.LP.position_m.append(self.z)

        # From List to Numpy Array
        self.LP.position_m = np.array(self.LP.position_m)

        # Return Result
        self.LiquidStreamOut = self.__get_slice_from_LiquidProfile__(LiquidStreamProfile=self.LiquidStream, height_index=0, LiquidStreamOut=self.LiquidStream)
        self.LiquidStream = deepcopy(self.LP)
        del self.LP

        return self.LiquidStreamOut


class LiquidPFR_Adiabatic(_Stochiometry, Serializer, _LiquidEquilibrium, _StreamFunctions):

    def __init__(self, position_m, cross_sectional_area_m2):

        self.info = {}

        self.firstscan = True
        self.position_m = position_m
        self.cross_sectional_area_m2 = cross_sectional_area_m2

        self.equilibrium_liquid_inlet = LiquidEquilibrium_Adiabatic()
        self.equilibrium_liquid_outlet = LiquidEquilibrium_Adiabatic()
        self.equilibrium_gas = None

    def react(self, LiquidStreamIn, step_size_m, lr=0.75):

        # Step Size of Integrator
        self.step_size_m = step_size_m

        # Only Strictly Positive Mass Fractions Allowed
        for id in LiquidStreamIn.specie.keys():
            LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(LiquidStreamIn.specie[id]["Mass Fraction"], 10 ** (-18))
        LiquidStreamIn.normalize_mass_fractions()

        # Process Liquid and Gas Inlet Streams
        self.LiquidStreamIn = self.equilibrium_liquid_inlet.react(LiquidStreamIn, lr=lr)

        if self.firstscan:

            # Extract Molar Masses
            self.M_liq = np.zeros(shape=(LiquidStreamIn.num_of_species,))
            for i, id in enumerate(LiquidStreamIn.specie.keys()):
                self.M_liq[i] = LiquidStreamIn.get_specie_molar_mass_kg_kmol(id)

            # Matrix or Reaction Stochiometry: Used to Calculate Sensitivity Matrices
            self.matrix_liq = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn, liq_rxn_insta=True)

            # Matrix of Reaction Stochiometry: Used to Calculate Heat of Reactions
            self.matrix_liq_heat = self.__get_the_matrix__(LiquidStreamIn=LiquidStreamIn,
                                                           liq_rxn_insta=True,
                                                           liq_rxn_reversible=True)
            self.firstscan = False

        # Initiate as Inlet
        self.LiquidStream = self.__broadcast_to_LiquidProfile__(LiquidStreamIn=self.LiquidStreamIn, num_of_heights=1)
        self.LP = deepcopy(self.LiquidStream)
        self.z = 0
        self.LP.position_m = [self.z]

        while self.z < self.position_m[-1]:

            A = np.interp(x=self.z, xp=self.position_m, fp=self.cross_sectional_area_m2)

            m_liq_tot = self.LiquidStream.get_solution_flow_kg_h() / 3600
            w_liq = self.LiquidStream.__mass_fractions_dic2vec__()
            m_liq = w_liq * m_liq_tot[:, :, None]
            cp_liq = self.LiquidStream.get_solution_heat_capacity_kJ_kgK()
            T_liq = self.LiquidStream.get_solution_temp_K()

            # Exothermic Heat [kJ/kmol rxn] of the following reactions
            # - Liquid Insta Reactions
            # - Liquid Reversibe Reactions
            # - Mass Transfer
            h = np.zeros(shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.temp_K.shape[1], self.LiquidStream.num_of_rxn_insta + self.LiquidStream.num_of_rxn_reversible), dtype=np.float64)
            for rxn_i, rxn_id in enumerate(self.LiquidStream.rxn_insta.keys()):
                h[:, :, rxn_i] = self.LiquidStream.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.LiquidStream.rxn_reversible.keys()):
                h[:, :, rxn_i + self.LiquidStream.num_of_rxn_insta] = self.LiquidStream.get_rxn_reversible_exothermic_heat_kJ_kmol(rxn_id)

            # Calculate Rate-Reactions in Liquid Phase
            # [r] = [kmol/m3.s]
            r_liq = self.LiquidStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

            # Differential Equations (dm/dz)
            dm_liq__dz = A * self.M_liq[None, None, :] * r_liq
            dm_liq_tot__dz = A * np.sum(dm_liq__dz, axis=2, keepdims=False)

            # Sensitivity in Liquid Phase (dmdB, dmdT)
            m_liq_broadcast = np.broadcast_to(array=m_liq[:, :, :, None], shape=(m_liq.shape[0], m_liq.shape[1], m_liq.shape[2], m_liq.shape[2]))
            dw_liq__dm_liq = - np.einsum("zsmn,zs->zsmn", m_liq_broadcast, 1 / m_liq_tot ** 2) + np.einsum("vw,zs->zsvw", np.eye(N=self.LiquidStream.num_of_species), 1 / m_liq_tot)
            df1__dw_liq = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStream)
            df1__dm_liq = np.einsum("zsfw,zswv->zsfv", df1__dw_liq, dw_liq__dm_liq)
            df2__dm_liq = np.broadcast_to(array=self.matrix_liq["A"][None, None, :, :], shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.temp_K.shape[1], self.matrix_liq["A"].shape[0], self.matrix_liq["A"].shape[1])).copy()
            df__dm_liq = np.concatenate((df1__dm_liq, df2__dm_liq), axis=2)
            H = np.linalg.inv(df__dm_liq)
            dm_liq__dB = H[:, :, :, self.LiquidStream.num_of_rxn_insta::]
            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStream)
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStream, Kw)
            Hr = H[:, :, :, :self.LiquidStream.num_of_rxn_insta:]
            dm_liq__dT = - np.einsum("zscr,zsr->zsc", Hr, (1 / Kw) * dKwdT)

            # Change in "B-Vector" for Liquid Phase
            dB__dz = np.einsum("bw,zsw->zsb", self.matrix_liq["A"], dm_liq__dz)


            # Initial Guess
            dm_liq_at_eq__dz = dm_liq__dz

            # Iterate Until "Convergence"
            for _ in range(10):
                q_lat = np.einsum("zsr,zsr->zs", h, np.einsum("rm,zsm->zsr", self.matrix_liq_heat["R+"], dm_liq_at_eq__dz)) / A
                dT_liq__dz = A * (1 / (cp_liq * m_liq_tot)) * q_lat
                dm_liq_at_eq__dz = np.einsum("zswb,zsb->zsw", dm_liq__dB, dB__dz) + np.einsum("zsm,zs->zsm", dm_liq__dT, dT_liq__dz)
            dm_liq__dz = dm_liq_at_eq__dz

            # Calculate Step Size
            with np.errstate(divide='ignore', invalid='ignore'):
                dz_max_01 = 0.1 * np.min(np.nan_to_num(m_liq / np.abs(dm_liq__dz), nan=np.inf, posinf=np.inf, neginf=np.inf))
            with np.errstate(divide='ignore', invalid='ignore'):
                dz_max_02 = np.nan_to_num(5 / np.abs(dT_liq__dz), nan=np.inf, posinf=np.inf, neginf=np.inf)[0,0]
            dz_max_03 = self.position_m[-1] - self.z
            dz = np.min(a=[dz_max_01, dz_max_02, dz_max_03, self.step_size_m])

            # Integrate
            dm_liq = dm_liq__dz * dz
            dT_liq = dT_liq__dz * dz
            dm_liq_tot = np.sum(dm_liq, axis=2, keepdims=False)
            w_liq_new = (m_liq + dm_liq) / (m_liq_tot + dm_liq_tot)[:, :, None]

            self.LiquidStream.temp_K = self.LiquidStream.temp_K + dT_liq
            self.LiquidStream.flow_kg_h = self.LiquidStream.flow_kg_h + 3600 * dm_liq_tot
            self.LiquidStream.__mass_fractions_vec2dic__(w=w_liq_new)

            self.z = self.z + dz

            # Append to Profile
            self.LP.temp_K = np.vstack((self.LP.temp_K, self.LiquidStream.temp_K))
            self.LP.flow_kg_h = np.vstack((self.LP.flow_kg_h, self.LiquidStream.flow_kg_h))
            for id in self.LiquidStream.specie.keys():
                self.LP.specie[id]["Mass Fraction"] = np.vstack((self.LP.specie[id]["Mass Fraction"], self.LiquidStream.specie[id]["Mass Fraction"]))
            for id in self.LiquidStream.info.keys():
                self.LP.info[id] = np.vstack((self.LP.info[id], self.LiquidStream.info[id]))
            self.LP.position_m.append(self.z)

        # From List to Numpy Array
        self.LP.position_m = np.array(self.LP.position_m)

        # Return Result
        self.LiquidStreamOut = self.__get_slice_from_LiquidProfile__(LiquidStreamProfile=self.LiquidStream, height_index=0, LiquidStreamOut=self.LiquidStream)
        self.LiquidStream = deepcopy(self.LP)
        del self.LP

        return self.LiquidStreamOut


# ---------------------------------------------------------------------------------------


class VaporLiquidEquilibrium_Isothermal(_Stochiometry, Serializer, _LiquidEquilibrium, _GasEquilibrium, _VaporLiquidEquilibrium, _StreamFunctions):

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
            #dz_newton = - np.linalg.solve(dfdz, f)
            dz_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdz), f)

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


class VaporLiquidEquilibrium_Adiabatic(_Stochiometry, Serializer, _LiquidEquilibrium, _GasEquilibrium, _VaporLiquidEquilibrium, _StreamFunctions):

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
            #dX_newton = - np.linalg.solve(dfdX, f)
            dX_newton = - np.einsum("scd,sd->sc", np.linalg.pinv(dfdX), f)

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


class VaporLiquidEquilibrium_EquilibriumStages(Serializer):

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


class LiquidHeatExchanger_CounterCurrent(Serializer):

    def __init__(self, area_m2):
        super().__init__()
        self.area_m2 = area_m2

    def load_heat_transfer_coefficient_kW_m2K(self, function):
        self.heat_transfer_coefficient_kW_m2K = function

    def react(self, LiquidStreamIn1, LiquidStreamIn2):

        self.LiquidStream1 = deepcopy(LiquidStreamIn1)
        self.LiquidStream2 = deepcopy(LiquidStreamIn2)
        LiquidStreamOut1 = deepcopy(LiquidStreamIn1)
        LiquidStreamOut2 = deepcopy(LiquidStreamIn2)
        T1_in = LiquidStreamIn1.get_solution_temp_K()
        T2_in = LiquidStreamIn2.get_solution_temp_K()
        m1 = LiquidStreamIn1.get_solution_flow_kg_h() / 3600
        m2 = LiquidStreamIn2.get_solution_flow_kg_h() / 3600

        for _ in range(5):
            cp1 = self.LiquidStream1.get_solution_heat_capacity_kJ_kgK()
            cp2 = self.LiquidStream2.get_solution_heat_capacity_kJ_kgK()
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
            LiquidStreamOut1.temp_K = T1_out
            LiquidStreamOut2.temp_K = T2_out
            self.LiquidStream1.temp_K = 0.5 * (LiquidStreamOut1.temp_K + LiquidStreamIn1.temp_K)
            self.LiquidStream2.temp_K = 0.5 * (LiquidStreamOut2.temp_K + LiquidStreamIn2.temp_K)

        # Additional info; LMTD
        dT_left = np.abs(T1_in - T2_out)
        dT_right = np.abs(T1_out - T2_in)
        self.LMTD = (dT_left - dT_right) / (np.log(dT_left) - np.log(dT_right))

        return LiquidStreamOut1, LiquidStreamOut2

    def get_interface_area_m2(self):
        return self.area_m2

    def get_logarithmic_mean_temperature_difference_K(self):
        return self.LMTD

    def get_heat_transfer_coefficient_kW_m2K(self):
        return self.kH

    def get_number_of_transfer_units(self):
        return self.NTU


class LiquidHeatExchanger_CoCurrent(Serializer):

    def __init__(self, area_m2):
        super().__init__()
        self.area_m2 = area_m2

    def load_heat_transfer_coefficient_kW_m2K(self, function):
        self.heat_transfer_coefficient_kW_m2K = function

    def react(self, LiquidStreamIn1, LiquidStreamIn2):

        self.LiquidStream1 = deepcopy(LiquidStreamIn1)
        self.LiquidStream2 = deepcopy(LiquidStreamIn2)
        LiquidStreamOut1 = deepcopy(LiquidStreamIn1)
        LiquidStreamOut2 = deepcopy(LiquidStreamIn2)
        T1_in = LiquidStreamIn1.get_solution_temp_K()
        T2_in = LiquidStreamIn2.get_solution_temp_K()
        m1 = LiquidStreamIn1.get_solution_flow_kg_h() / 3600
        m2 = LiquidStreamIn2.get_solution_flow_kg_h() / 3600

        for _ in range(5):
            cp1 = self.LiquidStream1.get_solution_heat_capacity_kJ_kgK()
            cp2 = self.LiquidStream2.get_solution_heat_capacity_kJ_kgK()
            C1 = cp1 * m1
            C2 = cp2 * m2
            Cmin = np.minimum(C1, C2)
            Cmax = np.maximum(C1, C2)
            Cr = Cmin / Cmax
            Qmax = Cmin * (T1_in - T2_in)
            self.kH = self.heat_transfer_coefficient_kW_m2K(self)
            self.NTU = self.kH * self.area_m2 / Cmin
            effectiveness = (1 - np.exp(-self.NTU * (1 + Cr))) / (1 + Cr)
            Q = effectiveness * Qmax
            T1_out = T1_in - Q / C1
            T2_out = T2_in + Q / C2
            LiquidStreamOut1.temp_K = T1_out
            LiquidStreamOut2.temp_K = T2_out
            self.LiquidStream1.temp_K = 0.5 * (LiquidStreamOut1.temp_K + LiquidStreamIn1.temp_K)
            self.LiquidStream2.temp_K = 0.5 * (LiquidStreamOut2.temp_K + LiquidStreamIn2.temp_K)

        # Additional info; LMTD
        dT_left = np.abs(T1_in - T2_out)
        dT_right = np.abs(T1_out - T2_in)
        self.LMTD = (dT_left - dT_right) / (np.log(dT_left) - np.log(dT_right))

        return LiquidStreamOut1, LiquidStreamOut2

    def get_interface_area_m2(self):
        return self.area_m2

    def get_logarithmic_mean_temperature_difference_K(self):
        return self.LMTD

    def get_heat_transfer_coefficient_kW_m2K(self):
        return self.kH

    def get_number_of_transfer_units(self):
        return self.NTU


# ---------------------------------------------------------------------------------------


class GasCompressor_Isentropic(Serializer):

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


class GasCompressor_Isothermal(Serializer):

    def __init__(self):
        super().__init__()
        self.compressor_power_kW = None

    def react(self, GasStreamIn, pressure_out_bara):
        GasStreamOut = deepcopy(GasStreamIn)
        GasStreamOut.pressure_bara = pressure_out_bara
        self.compressor_power_kW = 8.314 * GasStreamIn.temp_K * np.log(pressure_out_bara / GasStreamIn.get_gas_pressure_bara()) * GasStreamIn.flow_kmol_h / 3600
        return GasStreamOut

    def get_compressor_power_kW(self):
        return self.compressor_power_kW


# ---------------------------------------------------------------------------------------


class _GasLiquidContactor:

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

    def get_mass_transfer_kmol_m3s(self, id):
        return self.mass_transfer_kmol_m3s[id]["Rate [kmol/m3.s]"](self)

    def get_heat_transfer_kW_m3(self):
        return self.heat_transfer_kW_m3(self)

    def get_liquid_holdup_m3_m3(self):
        return self.liquid_holdup_m3_m3(self)

    def __get_mass_transfer_vectors__(self):

        dT = self.GasStream.get_gas_temp_K() - self.LiquidStream.get_solution_temp_K()
        aJ_specie_liq = np.zeros(shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.temp_K.shape[1], self.LiquidStream.num_of_species), dtype=np.float64)
        aJ_specie_gas = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1], self.GasStream.num_of_species), dtype=np.float64)
        aJ_rxn = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1], self.num_of_mass_transfer), dtype=np.float64)
        q_rxn = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1]), dtype=np.float64)

        for k, id in enumerate(self.mass_transfer_kmol_m3s.keys()):
            aJ_rxn[:,:, k] = self.mass_transfer_kmol_m3s[id]["Rate [kmol/m3.s]"](self)

        for k, id in enumerate(self.mass_transfer_kmol_m3s.keys()):
            nu = self.mass_transfer_kmol_m3s[id]["Stoch Liq"]
            for el in nu.keys():
                i = self.LiquidStream.specie[el]["Index"]
                aJ_specie_liq[:,:, i] = aJ_specie_liq[:,:, i] + aJ_rxn[:,:, k] * nu[el]
            nu = self.mass_transfer_kmol_m3s[id]["Stoch Gas"]
            for el in nu.keys():
                i = self.GasStream.specie[el]["Index"]
                aJ_specie_gas[:,:, i] = aJ_specie_gas[:,:, i] + aJ_rxn[:,:, k] * nu[el]

        # Heat due to Heating/Cooling of Gas during Absorption
        q_abs = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1]), dtype=np.float64)
        q_des = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1]), dtype=np.float64)
        for i, id in enumerate(self.GasStream.specie.keys()):
            cp = self.GasStream.get_specie_heat_capacity_kJ_kmolK(id)
            q_abs = q_abs + dT * cp * np.maximum(aJ_specie_gas[:,:, i], 0)
            q_des = q_des + dT * cp * np.minimum(aJ_specie_gas[:,:, i], 0)

        return aJ_rxn, aJ_specie_gas, aJ_specie_liq, q_abs, q_des

    def __get_mass_transfer_vectors_bulk__(self):

        dT = self.GasBulk.get_gas_temp_K() - self.LiquidBulk.get_solution_temp_K()
        aJ_specie_liq = np.zeros(shape=(self.LiquidBulk.temp_K.shape[0], self.LiquidBulk.temp_K.shape[1], self.LiquidBulk.num_of_species), dtype=np.float64)
        aJ_specie_gas = np.zeros(shape=(self.GasBulk.temp_K.shape[0], self.GasBulk.temp_K.shape[1], self.GasBulk.num_of_species), dtype=np.float64)
        aJ_rxn = np.zeros(shape=(self.GasBulk.temp_K.shape[0], self.GasBulk.temp_K.shape[1], self.num_of_mass_transfer), dtype=np.float64)
        q_rxn = np.zeros(shape=(self.GasBulk.temp_K.shape[0], self.GasBulk.temp_K.shape[1]), dtype=np.float64)

        for k, id in enumerate(self.mass_transfer_kmol_m3s.keys()):
            aJ_rxn[:,:, k] = self.mass_transfer_kmol_m3s[id]["Rate [kmol/m3.s]"](self)

        for k, id in enumerate(self.mass_transfer_kmol_m3s.keys()):
            nu = self.mass_transfer_kmol_m3s[id]["Stoch Liq"]
            for el in nu.keys():
                i = self.LiquidBulk.specie[el]["Index"]
                aJ_specie_liq[:,:, i] = aJ_specie_liq[:,:, i] + aJ_rxn[:,:, k] * nu[el]
            nu = self.mass_transfer_kmol_m3s[id]["Stoch Gas"]
            for el in nu.keys():
                i = self.GasBulk.specie[el]["Index"]
                aJ_specie_gas[:,:, i] = aJ_specie_gas[:,:, i] + aJ_rxn[:,:, k] * nu[el]

        # Heat due to Heating/Cooling of Gas during Absorption
        q_abs = np.zeros(shape=(self.GasBulk.temp_K.shape[0], self.GasBulk.temp_K.shape[1]), dtype=np.float64)
        q_des = np.zeros(shape=(self.GasBulk.temp_K.shape[0], self.GasBulk.temp_K.shape[1]), dtype=np.float64)
        for i, id in enumerate(self.GasBulk.specie.keys()):
            cp = self.GasBulk.get_specie_heat_capacity_kJ_kmolK(id)
            q_abs = q_abs + dT * cp * np.maximum(aJ_specie_gas[:,:, i], 0)
            q_des = q_des + dT * cp * np.minimum(aJ_specie_gas[:,:, i], 0)

        return aJ_rxn, aJ_specie_gas, aJ_specie_liq, q_abs, q_des


class _GasLiquidContactor_PFR:

    def add_pressure_drop_Pa_m(self, pressure_drop_Pa_m):
        self.pressure_drop_Pa_m = pressure_drop_Pa_m

    def get_pressure_drop_Pa_m(self):
        return self.pressure_drop_Pa_m(self)


class GasLiquidContactor_PFR_CounterCurrent(_Stochiometry, Serializer, _LiquidEquilibrium, _GasEquilibrium, _StreamFunctions, _GasLiquidContactor, _GasLiquidContactor_PFR):

    def __init__(self, position_m, cross_sectional_area_m2, void_fraction_m3_m3):
        self.info = {}
        self.mass_transfer_kmol_m3s = {}
        self.num_of_mass_transfer = 0
        self.firstscan = True
        self.position_m = position_m
        self.cross_sectional_area_m2 = cross_sectional_area_m2
        self.void_fraction_m3_m3 = void_fraction_m3_m3
        self.num_of_heights = len(position_m)

        # Instance for Calculating Chemical Equilibrium
        self.equilibrium_liquid_inlet = LiquidEquilibrium_Adiabatic()
        self.equilibrium_liquid_outlet = LiquidEquilibrium_Adiabatic()
        self.equilibrium_gas_inlet = GasEquilibrium_Adiabatic_Isobaric()
        self.equilibrium_gas_outlet = GasEquilibrium_Adiabatic_Isobaric()

    def get_superficial_gas_velocity_m_s(self):
        T = self.GasStream.get_gas_temp_K()
        n = self.GasStream.get_gas_flow_kmol_h() / 3600
        p = self.GasStream.get_gas_pressure_bara()
        R = 0.08314
        V = n * R * T / p
        v = V / self.cross_sectional_area_m2[:,None]
        return v

    def get_superficial_liquid_velocity_m_s(self):
        rho = self.LiquidStream.get_solution_density_kg_m3()
        m = self.LiquidStream.get_solution_flow_kg_h() / 3600
        V = m / rho
        v = V / self.cross_sectional_area_m2[:,None]
        return v

    def get_cross_sectional_area_m2(self):
        return self.cross_sectional_area_m2[:,None]

    def __react__(self, GasStreamIn, LiquidStreamIn, epochs, lr=0.75):

        # Only Strictly Positive Mass Fractions Allowed
        for id in LiquidStreamIn.specie.keys():
            LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(LiquidStreamIn.specie[id]["Mass Fraction"], 10 ** (-18))
        LiquidStreamIn.normalize_mass_fractions()
        for id in GasStreamIn.specie.keys():
            GasStreamIn.specie[id]["Molar Fraction"] = np.maximum(GasStreamIn.specie[id]["Molar Fraction"], 10 ** (-18))
        GasStreamIn.normalize_molar_fractions()

        # Process Liquid and Gas Inlet Streams
        self.LiquidStreamIn = self.equilibrium_liquid_inlet.react(LiquidStreamIn, lr=lr)
        self.GasStreamIn = self.equilibrium_gas_inlet.react(GasStreamIn, lr=lr)

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

            self.LiquidStream = self.__broadcast_to_LiquidProfile__(LiquidStreamIn=self.LiquidStreamIn, num_of_heights=self.num_of_heights)
            self.GasStream = self.__broadcast_to_GasProfile__(GasStreamIn=self.GasStreamIn, num_of_heights=self.num_of_heights)
            self.firstscan = False
            self.GasStreamOut = deepcopy(self.GasStreamIn)
            self.LiquidStreamOut = deepcopy(self.LiquidStreamIn)

        # Update Top and Bottom as Inlet Streams
        self.LiquidStream.temp_K[-1,:] = self.LiquidStreamIn.temp_K
        self.LiquidStream.flow_kg_h[-1,:] = self.LiquidStreamIn.flow_kg_h
        for id in self.LiquidStream.specie.keys():
            self.LiquidStream.specie[id]["Mass Fraction"][-1,:] = self.LiquidStreamIn.specie[id]["Mass Fraction"]
        for id in self.LiquidStream.info.keys():
            self.LiquidStream.info[id] = np.broadcast_to(array=self.LiquidStreamIn.info[id][None,:], shape=(self.num_of_heights, self.LiquidStreamIn.temp_K.shape[0])).copy()

        self.GasStream.temp_K[0, :] = self.GasStreamIn.temp_K
        self.GasStream.flow_kmol_h[0, :] = self.GasStreamIn.flow_kmol_h
        for id in self.GasStream.specie.keys():
            self.GasStream.specie[id]["Molar Fraction"][0, :] = self.GasStreamIn.specie[id]["Molar Fraction"]
        for id in self.GasStream.info.keys():
            self.GasStream.info[id] = np.broadcast_to(array=self.GasStreamIn.info[id][None, :], shape=(self.num_of_heights, self.GasStreamIn.temp_K.shape[0])).copy()


        for self.epochs in range(epochs):

            """""""""
            LIQUID STREAM
            """""""""

            m_liq_tot = self.LiquidStream.get_solution_flow_kg_h() / 3600
            w_liq = self.LiquidStream.__mass_fractions_dic2vec__()
            m_liq = w_liq * m_liq_tot[:,:,None]
            cp_liq = self.LiquidStream.get_solution_heat_capacity_kJ_kgK()

            # Exothermic Heat [kJ/kmol rxn] of the following reactions
            # - Liquid Insta Reactions
            # - Liquid Reversibe Reactions
            # - Mass Transfer
            h = np.zeros(shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.temp_K.shape[1], self.LiquidStream.num_of_rxn_insta + self.LiquidStream.num_of_rxn_reversible + self.num_of_mass_transfer), dtype=np.float64)
            for rxn_i, rxn_id in enumerate(self.LiquidStream.rxn_insta.keys()):
                h[:,:, rxn_i] = self.LiquidStream.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.LiquidStream.rxn_reversible.keys()):
                h[:,:,rxn_i + self.LiquidStream.num_of_rxn_insta] = self.LiquidStream.get_rxn_reversible_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.mass_transfer_kmol_m3s.keys()):
                h[:,:, rxn_i + self.LiquidStream.num_of_rxn_insta + self.LiquidStream.num_of_rxn_reversible] = self.mass_transfer_kmol_m3s[rxn_id]["Exothermic Heat [kJ/kmol]"](self)

            # Calculate Absorption Rates and some Related Heat Dissipation.
            # [aJ] = [kmol/m3.s]
            # [q] = [kW/m3]
            aJ_rxn, aJ_specie_gas, aJ_specie_liq, q_abs, q_des = self.__get_mass_transfer_vectors__()

            # Calculate Rate-Reactions in Liquid Phase
            # [r] = [kmol/m3.s]
            r_liq = self.LiquidStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

            # Packing Hydrodynamics
            liquid_holdup_m3_m3 = self.get_liquid_holdup_m3_m3()

            # Differential Equations (dm/dz)
            dm_liq__dz = self.cross_sectional_area_m2[:,None,None] * (self.M_liq[None, None, :] * aJ_specie_liq + self.void_fraction_m3_m3[:,None,None] * liquid_holdup_m3_m3[:,:,None] * self.M_liq[None,None,:] * r_liq)
            dm_liq_tot__dz = self.cross_sectional_area_m2[:,None] * np.sum(dm_liq__dz , axis=2, keepdims=False)
            dm_gas__dz = self.cross_sectional_area_m2[:,None,None] * self.M_gas[None, None, :] * aJ_specie_gas
            dm_gas_tot__dz = - dm_liq_tot__dz

            # Sensitivity in Liquid Phase (dmdB, dmdT)
            m_liq_broadcast = np.broadcast_to(array=m_liq[:, :, :, None], shape=(m_liq.shape[0], m_liq.shape[1], m_liq.shape[2], m_liq.shape[2]))
            dw_liq__dm_liq = - np.einsum("zsmn,zs->zsmn", m_liq_broadcast, 1 / m_liq_tot ** 2) + np.einsum("vw,zs->zsvw", np.eye(N=self.LiquidStream.num_of_species), 1 / m_liq_tot)
            df1__dw_liq = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStream)
            df1__dm_liq = np.einsum("zsfw,zswv->zsfv", df1__dw_liq, dw_liq__dm_liq)
            df2__dm_liq = np.broadcast_to(array=self.matrix_liq["A"][None, None, :, :], shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.temp_K.shape[1], self.matrix_liq["A"].shape[0], self.matrix_liq["A"].shape[1])).copy()
            df__dm_liq = np.concatenate((df1__dm_liq, df2__dm_liq), axis=2)
            H = np.linalg.inv(df__dm_liq)
            dm_liq__dB = H[:, :, :, self.LiquidStream.num_of_rxn_insta::]
            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStream)
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStream, Kw)
            Hr = H[:, :, :, :self.LiquidStream.num_of_rxn_insta:]
            dm_liq__dT = - np.einsum("zscr,zsr->zsc", Hr, (1 / Kw) * dKwdT)

            # Change in "B-Vector" for Liquid Phase
            dB__dz = np.einsum("bw,zsw->zsb", self.matrix_liq["A"], dm_liq__dz)

            # Direct Heat Transfer (kW/m3)
            q_dir = self.get_heat_transfer_kW_m3()

            # Note:
            # <dm_liq__dz> give the change in mass flow when equilibrium constraints are NOT taken into account
            # <dm_liq_at_eq__dz> take the equilibrium concentrations into account.
            # However, the latter also depend on the temperature and vice versa.
            dm_liq_at_eq__dz = dm_liq__dz
            for _ in range(10):
                dm_at_eq__dz = np.concatenate((dm_liq_at_eq__dz, dm_gas__dz), axis=2)
                q_lat = np.einsum("zsr,zsr->zs", h, np.einsum("rm,zsm->zsr", self.matrix_liq_heat["R+"], dm_at_eq__dz)) / self.cross_sectional_area_m2[:,None]
                dT_liq__dz = self.cross_sectional_area_m2[:,None] * (1 / (cp_liq * m_liq_tot)) * (q_dir + q_abs + q_lat)
                dm_liq_at_eq__dz = np.einsum("zswb,zsb->zsw", dm_liq__dB, dB__dz) + np.einsum("zsm,zs->zsm", dm_liq__dT, dT_liq__dz)
            dm_liq__dz = dm_liq_at_eq__dz


            """""""""
            GAS STREAM
            """""""""

            n_gas_tot = self.GasStream.get_gas_flow_kmol_h() / 3600
            y_gas = self.GasStream.__molar_fractions_dic2vec__()
            n_gas = y_gas * n_gas_tot[:,:,None]
            cp_gas = self.GasStream.get_gas_heat_capacity_kJ_kmolK()

            # Exothermic Heat [kJ/kmol rxn] of the following reactions
            # - Gas Insta Reactions
            # - Gas Reversibe Reactions
            h_gas = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1], self.GasStream.num_of_rxn_insta + self.GasStream.num_of_rxn_reversible),dtype=np.float64)
            for rxn_i, rxn_id in enumerate(self.GasStream.rxn_insta.keys()):
                h_gas[:,:, rxn_i] = self.GasStream.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.GasStream.rxn_reversible.keys()):
                h_gas[:,:,rxn_i + self.GasStream.num_of_rxn_insta] = self.GasStream.get_rxn_reversible_exothermic_heat_kJ_kmol(rxn_id)

            # Calculate Rate-Reactions in Gas Phase
            r_gas = self.GasStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

            # Differential Equations, Equilibrium Reactions not taken into Account
            dn_gas__dz = self.cross_sectional_area_m2[:,None,None] * (aJ_specie_gas + self.void_fraction_m3_m3[:,None,None] * (1 - liquid_holdup_m3_m3[:,:, None]) * r_gas)
            dn_gas_tot__dz = np.sum(dn_gas__dz, axis=2, keepdims=False)

            # Sensitivity in Gas Phase (dndB, dndT)
            n_gas_broadcast = np.broadcast_to(array=n_gas[:, :, :, None],shape=(n_gas.shape[0], n_gas.shape[1], n_gas.shape[2], n_gas.shape[2]))
            dy_gas__dn_gas = - np.einsum("zsmn,zs->zsmn", n_gas_broadcast, 1 / n_gas_tot ** 2) + np.einsum("vw,zs->zsvw", np.eye(N=self.GasStream.num_of_species), 1 / n_gas_tot)
            df1__dy_gas = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(self.GasStream)
            df1__dn_gas = np.einsum("zsfw,zswv->zsfv", df1__dy_gas, dy_gas__dn_gas)
            df2__dn_gas = np.broadcast_to(array=self.matrix_gas["A"][None, None, :, :], shape=(self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1], self.matrix_gas["A"].shape[0], self.matrix_gas["A"].shape[1])).copy()
            df__dn_gas = np.concatenate((df1__dn_gas, df2__dn_gas), axis=2)
            H = np.linalg.inv(df__dn_gas)
            dn_gas__dB = H[:, :, :, self.GasStream.num_of_rxn_insta::]
            Ky = self.__get_GasEquilibrium_Ky__(self.GasStream)
            dKydT = self.__get_GasEquilibrium_dKydT__(self.GasStream, Ky)
            Hr = H[:, :, :, :self.GasStream.num_of_rxn_insta:]
            dn_gas__dT = - np.einsum("zscr,zsr->zsc", Hr, (1 / Ky) * dKydT)

            # Change in "B-Vector" for Gas Phase
            dB__dz = np.einsum("bw,zsw->zsb", self.matrix_gas["A"], dn_gas__dz)

            # Temperature Gradient
            dn_gas_at_eq__dz = dn_gas__dz
            for _ in range(10):
                q_lat = np.einsum("zsr,zsr->zs", h_gas, np.einsum("rm,zsm->zsr", self.matrix_gas_heat["R+"], dn_gas_at_eq__dz)) / self.cross_sectional_area_m2[:, None]
                dT_gas__dz = self.cross_sectional_area_m2[:, None] * (1 / (cp_gas * n_gas_tot)) * (-q_dir + q_des + q_lat)
                dn_gas_at_eq__dz = np.einsum("zswb,zsb->zsw", dn_gas__dB, dB__dz) + np.einsum("zsm,zs->zsm", dn_gas__dT, dT_gas__dz)
            dn_gas__dz = dn_gas_at_eq__dz


            """""""""
            Now we have obtained following derivatives.
            1) dm_liq__dz
            2) dn_gas__dz
            3) dT_liq__dz
            4) dT_gas__dz
            """""""""

            dz = self.position_m[1] - self.position_m[0]

            dm_liq = dm_liq__dz * dz
            dn_gas = dn_gas__dz * dz

            dT_gas = dT_gas__dz * dz
            dT_liq = dT_liq__dz * dz

            dm_liq_tot = np.sum(dm_liq, axis=2, keepdims=False)
            dn_gas_tot = np.sum(dn_gas, axis=2, keepdims=False)

            y_gas_new = (n_gas + dn_gas) / (n_gas_tot + dn_gas_tot)[:,:,None]
            dy_gas = y_gas_new - y_gas

            w_liq_new = (m_liq + dm_liq) / (m_liq_tot + dm_liq_tot)[:,:,None]
            dw_liq = w_liq_new - w_liq

            dm_liq = - dm_liq
            dw_liq = - dw_liq
            dT_liq = - dT_liq

            # -----------------------

            self.lr = 0.9
            n = self.num_of_heights - 1

            # Gas
            n_gas_tot_target = self.GasStream.flow_kmol_h[:n:, :] + 3600 * 0.5 * (dn_gas_tot[:n:, :] + dn_gas_tot[1::, :])
            T_gas_target = self.GasStream.temp_K[:n:,:] + 0.5 * (dT_gas[:n:,:] + dT_gas[1::,:])
            y_gas_target = np.zeros(shape=(y_gas.shape[0] - 1, y_gas.shape[1], y_gas.shape[2]))
            for i, id in enumerate(self.GasStream.specie.keys()):
                y_gas_target[:,:,i] = self.GasStream.specie[id]["Molar Fraction"][:n:,:] + 0.5 * (dy_gas[:n:,:,i] + dy_gas[1::,:,i])

            self.GasStream.flow_kmol_h[1::,:] = self.GasStream.flow_kmol_h[1::,:] + self.lr * (n_gas_tot_target - self.GasStream.flow_kmol_h[1::,:])
            self.GasStream.temp_K[1::,:] = np.clip(a=self.GasStream.temp_K[1::,:] + self.lr * (T_gas_target - self.GasStream.temp_K[1::,:]), a_min=self.GasStream.temp_K[1::,:] - 5, a_max=self.GasStream.temp_K[1::,:] + 5)
            #self.GasStream.pressure_bara[:, 1::] = self.GasStream.pressure_bara[:, 1::] + self.lr * (self.GasStream.pressure_bara[:, :n:] + self.dz * self.dp_gas[:,:n:] - self.GasStream.pressure_bara[:, 1::])
            for i, id in enumerate(self.GasStream.specie.keys()):
                self.GasStream.specie[id]["Molar Fraction"][1::,:] = self.GasStream.specie[id]["Molar Fraction"][1::,:] + self.lr * (y_gas_target[:,:,i] - self.GasStream.specie[id]["Molar Fraction"][1::,:])
                self.GasStream.specie[id]["Molar Fraction"] = np.maximum(self.GasStream.specie[id]["Molar Fraction"], 0)
            self.GasStream.normalize_molar_fractions()

            # Liquid
            m_liq_tot_target = self.LiquidStream.flow_kg_h[1::,:] + 0.5 * 3600 * (dm_liq_tot[1::,:] + dm_liq_tot[:n:,:])
            T_liq_target = self.LiquidStream.temp_K[1::,:] - 0.5 * (dT_liq[1::,:] + dT_liq[:n:,:])
            w_liq_target = np.zeros(shape=(w_liq.shape[0] - 1, w_liq.shape[1], w_liq.shape[2]))
            for i, id in enumerate(self.LiquidStream.specie.keys()):
                w_liq_target[:, :, i] = self.LiquidStream.specie[id]["Mass Fraction"][1::,:] - 0.5 * (dw_liq[1::,:,i] + dw_liq[:n:,:,i])

            self.LiquidStream.flow_kg_h[:n:,:] = self.LiquidStream.flow_kg_h[:n:,:] + self.lr * (m_liq_tot_target - self.LiquidStream.flow_kg_h[:n:,:])
            self.LiquidStream.temp_K[:n:,:] = np.clip(a=self.LiquidStream.temp_K[:n:,:] + self.lr * (T_liq_target - self.LiquidStream.temp_K[:n:,:]), a_min=self.LiquidStream.temp_K[:n:,:] - 5, a_max=self.LiquidStream.temp_K[:n:,:] + 5)
            for i, id in enumerate(self.LiquidStream.specie.keys()):
                self.LiquidStream.specie[id]["Mass Fraction"][:n:,:] = self.LiquidStream.specie[id]["Mass Fraction"][:n:,:] + self.lr * (w_liq_target[:,:,i] - self.LiquidStream.specie[id]["Mass Fraction"][:n:,:])
                self.LiquidStream.specie[id]["Mass Fraction"] = np.maximum(self.LiquidStream.specie[id]["Mass Fraction"], 10**(-18))
            self.LiquidStream.normalize_mass_fractions()

            # Perform One Single (Isothermal) Equilibrium Step for the sake of Numerical Stability........
            w = self.LiquidStream.__mass_fractions_dic2vec__()
            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStream)
            f = self.__get_LiquidEquilibrium_f_rxn_insta_log__(self.LiquidStream, Kw)
            dfdw = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStream)
            dfdr = np.einsum("zsrc,cq->zsrq", dfdw, self.equilibrium_liquid_inlet.matrix["R"])
            #dr_newton = - np.linalg.solve(dfdr, f)
            dr_newton = - np.einsum("zscd,zsd->zsc", np.linalg.pinv(dfdr), f)
            dw_newton = np.einsum("ij,zsj->zsi", self.equilibrium_liquid_inlet.matrix["R"], dr_newton)
            dw = self.lr * dw_newton
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            w = w + dw * tau
            self.LiquidStream.__mass_fractions_vec2dic__(np.maximum(w, 10**(-18)))

            y = self.GasStream.__molar_fractions_dic2vec__()
            w = self.GasStream.__molefrac2massfrac__(y)
            Kp = self.__get_GasEquilibrium_Kp__(self.GasStream)
            f = self.__get_GasEquilibrium_f_rxn_insta_log__(self.GasStream, Kp)
            dfdy = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(self.GasStream)
            dydw = self.GasStream.__dydw__(w)
            dfdw = np.einsum("zsfy,zsyw->zsfw", dfdy, dydw)
            dfdr = np.einsum("zsrc,cq->zsrq", dfdw, self.equilibrium_gas_inlet.matrix["R"])
            #dr_newton = - np.linalg.solve(dfdr, f)
            dr_newton = - np.einsum("zscd,zsd->zsc", np.linalg.pinv(dfdr), f)
            dw_newton = np.einsum("ij,zsj->zsi", self.equilibrium_gas_inlet.matrix["R"], dr_newton)
            dw = self.lr * dw_newton
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            w = w + dw * tau
            self.GasStream.__molar_fractions_vec2dic__(y=self.GasStream.__massfrac2molefrac__(w))

        # Update Outlet Streams
        self.GasStreamOut.temp_K = self.GasStream.temp_K[-1,:]
        self.GasStreamOut.flow_kmol_h = self.GasStream.flow_kmol_h[-1,:]
        for id in self.GasStreamOut.specie.keys():
            self.GasStreamOut.specie[id]["Molar Fraction"] = self.GasStream.specie[id]["Molar Fraction"][-1,:]
        for id in self.GasStreamOut.info.keys():
            self.GasStreamOut.info[id] = self.GasStream.info[id][-1, :]
        self.GasStreamOut = self.equilibrium_gas_outlet.react(self.GasStreamOut, lr=lr)

        self.LiquidStreamOut.temp_K = self.LiquidStream.temp_K[0, :]
        self.LiquidStreamOut.flow_kg_h = self.LiquidStream.flow_kg_h[0, :]
        for id in self.LiquidStreamOut.specie.keys():
            self.LiquidStreamOut.specie[id]["Mass Fraction"] = self.LiquidStream.specie[id]["Mass Fraction"][0, :]
        for id in self.LiquidStreamOut.info.keys():
            self.LiquidStreamOut.info[id] = self.LiquidStream.info[id][0, :]
        self.LiquidStreamOut = self.equilibrium_liquid_outlet.react(self.LiquidStreamOut, lr=lr)

        return self.GasStreamOut, self.LiquidStreamOut


class GasLiquidContactor_PFR_CoCurrent(_Stochiometry, Serializer, _LiquidEquilibrium, _GasEquilibrium, _StreamFunctions, _GasLiquidContactor, _GasLiquidContactor_PFR):

    def __init__(self, position_m, cross_sectional_area_m2, void_fraction_m3_m3):

        self.info = {}

        self.mass_transfer_kmol_m3s = {}
        self.num_of_mass_transfer = 0

        self.firstscan = True
        self.position_m = position_m
        self.cross_sectional_area_m2 = cross_sectional_area_m2
        self.void_fraction_m3_m3 = void_fraction_m3_m3

        self.equilibrium_liquid_inlet = LiquidEquilibrium_Adiabatic()
        self.equilibrium_liquid_outlet = LiquidEquilibrium_Adiabatic()
        self.equilibrium_gas_inlet = GasEquilibrium_Adiabatic_Isobaric()
        self.equilibrium_gas_outlet = GasEquilibrium_Adiabatic_Isobaric()

    def get_superficial_gas_velocity_m_s(self):
        A = np.interp(x=self.z, xp=self.position_m, fp=self.cross_sectional_area_m2)
        T = self.GasStream.get_gas_temp_K()
        n = self.GasStream.get_gas_flow_kmol_h() / 3600
        p = self.GasStream.get_gas_pressure_bara()
        R = 0.08314
        V = n * R * T / p
        v = V / A
        return v

    def get_superficial_liquid_velocity_m_s(self):
        A = np.interp(x=self.z, xp=self.position_m, fp=self.cross_sectional_area_m2)
        rho = self.LiquidStream.get_solution_density_kg_m3()
        m = self.LiquidStream.get_solution_flow_kg_h() / 3600
        V = m / rho
        v = V / A
        return v

    def __react__(self, GasStreamIn, LiquidStreamIn, step_size_m, lr=0.75):

        # Step Size of Integrator
        self.step_size_m = step_size_m

        # Only Strictly Positive Mass Fractions Allowed
        for id in LiquidStreamIn.specie.keys():
            LiquidStreamIn.specie[id]["Mass Fraction"] = np.maximum(LiquidStreamIn.specie[id]["Mass Fraction"], 10 ** (-18))
        LiquidStreamIn.normalize_mass_fractions()
        for id in GasStreamIn.specie.keys():
            GasStreamIn.specie[id]["Molar Fraction"] = np.maximum(GasStreamIn.specie[id]["Molar Fraction"], 10 ** (-18))
        GasStreamIn.normalize_molar_fractions()

        # Process Liquid and Gas Inlet Streams
        self.LiquidStreamIn = self.equilibrium_liquid_inlet.react(LiquidStreamIn, lr=lr)
        self.GasStreamIn = self.equilibrium_gas_inlet.react(GasStreamIn, lr=lr)

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
            self.firstscan = False

        # Initiate as Inlet
        self.LiquidStream = self.__broadcast_to_LiquidProfile__(LiquidStreamIn=self.LiquidStreamIn, num_of_heights=1)
        self.GasStream = self.__broadcast_to_GasProfile__(GasStreamIn=self.GasStreamIn, num_of_heights=1)

        self.LP = deepcopy(self.LiquidStream)
        self.GP = deepcopy(self.GasStream)
        self.z = 0
        self.LP.position_m = [self.z]
        self.GP.position_m = [self.z]


        while self.z < self.position_m[-1]:

            A = np.interp(x=self.z, xp=self.position_m, fp=self.cross_sectional_area_m2)
            epsilon = np.interp(x=self.z, xp=self.position_m, fp=self.void_fraction_m3_m3)

            """""""""
            LIQUID STREAM
            """""""""

            m_liq_tot = self.LiquidStream.get_solution_flow_kg_h() / 3600
            w_liq = self.LiquidStream.__mass_fractions_dic2vec__()
            m_liq = w_liq * m_liq_tot[:, :, None]
            cp_liq = self.LiquidStream.get_solution_heat_capacity_kJ_kgK()
            T_liq = self.LiquidStream.get_solution_temp_K()

            # Exothermic Heat [kJ/kmol rxn] of the following reactions
            # - Liquid Insta Reactions
            # - Liquid Reversibe Reactions
            # - Mass Transfer
            h = np.zeros(shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.temp_K.shape[1], self.LiquidStream.num_of_rxn_insta + self.LiquidStream.num_of_rxn_reversible + self.num_of_mass_transfer), dtype=np.float64)
            for rxn_i, rxn_id in enumerate(self.LiquidStream.rxn_insta.keys()):
                h[:, :, rxn_i] = self.LiquidStream.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.LiquidStream.rxn_reversible.keys()):
                h[:, :,
                rxn_i + self.LiquidStream.num_of_rxn_insta] = self.LiquidStream.get_rxn_reversible_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.mass_transfer_kmol_m3s.keys()):
                h[:, :, rxn_i + self.LiquidStream.num_of_rxn_insta + self.LiquidStream.num_of_rxn_reversible] = self.mass_transfer_kmol_m3s[rxn_id]["Exothermic Heat [kJ/kmol]"](self)

            # Calculate Absorption Rates and some Related Heat Dissipation.
            # [aJ] = [kmol/m3.s]
            # [q] = [kW/m3]
            aJ_rxn, aJ_specie_gas, aJ_specie_liq, q_abs, q_des = self.__get_mass_transfer_vectors__()

            # Calculate Rate-Reactions in Liquid Phase
            # [r] = [kmol/m3.s]
            r_liq = self.LiquidStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

            # Packing Hydrodynamics
            liquid_holdup_m3_m3 = self.get_liquid_holdup_m3_m3()

            # Differential Equations (dm/dz)
            dm_liq__dz = A * (self.M_liq[None, None, :] * aJ_specie_liq + epsilon * liquid_holdup_m3_m3[:, :, None] * self.M_liq[None, None, :] * r_liq)
            dm_liq_tot__dz = A * np.sum(dm_liq__dz, axis=2, keepdims=False)
            dm_gas__dz = A * self.M_gas[None, None, :] * aJ_specie_gas
            dm_gas_tot__dz = - dm_liq_tot__dz

            # Sensitivity in Liquid Phase (dmdB, dmdT)
            m_liq_broadcast = np.broadcast_to(array=m_liq[:, :, :, None], shape=(m_liq.shape[0], m_liq.shape[1], m_liq.shape[2], m_liq.shape[2]))
            dw_liq__dm_liq = - np.einsum("zsmn,zs->zsmn", m_liq_broadcast, 1 / m_liq_tot ** 2) + np.einsum("vw,zs->zsvw", np.eye(N=self.LiquidStream.num_of_species), 1 / m_liq_tot)
            df1__dw_liq = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStream)
            df1__dm_liq = np.einsum("zsfw,zswv->zsfv", df1__dw_liq, dw_liq__dm_liq)
            df2__dm_liq = np.broadcast_to(array=self.matrix_liq["A"][None, None, :, :], shape=(self.LiquidStream.temp_K.shape[0], self.LiquidStream.temp_K.shape[1], self.matrix_liq["A"].shape[0], self.matrix_liq["A"].shape[1])).copy()
            df__dm_liq = np.concatenate((df1__dm_liq, df2__dm_liq), axis=2)
            H = np.linalg.inv(df__dm_liq)
            dm_liq__dB = H[:, :, :, self.LiquidStream.num_of_rxn_insta::]
            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStream)
            dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidStream, Kw)
            Hr = H[:, :, :, :self.LiquidStream.num_of_rxn_insta:]
            dm_liq__dT = - np.einsum("zscr,zsr->zsc", Hr, (1 / Kw) * dKwdT)

            # Change in "B-Vector" for Liquid Phase
            dB__dz = np.einsum("bw,zsw->zsb", self.matrix_liq["A"], dm_liq__dz)

            # Direct Heat Transfer (kW/m3)
            q_dir = self.get_heat_transfer_kW_m3()

            # Initial Guess
            dm_liq_at_eq__dz = dm_liq__dz

            # Iterate Until "Convergence"
            for _ in range(10):
                dm_at_eq__dz = np.concatenate((dm_liq_at_eq__dz, dm_gas__dz), axis=2)
                q_lat = np.einsum("zsr,zsr->zs", h, np.einsum("rm,zsm->zsr", self.matrix_liq_heat["R+"], dm_at_eq__dz)) / A
                dT_liq__dz = A * (1 / (cp_liq * m_liq_tot)) * ( q_dir + q_abs + q_lat)
                dm_liq_at_eq__dz = np.einsum("zswb,zsb->zsw", dm_liq__dB, dB__dz) + np.einsum("zsm,zs->zsm", dm_liq__dT, dT_liq__dz)
            dm_liq__dz = dm_liq_at_eq__dz

            """""""""
            GAS STREAM
            """""""""

            n_gas_tot = self.GasStream.get_gas_flow_kmol_h() / 3600
            y_gas = self.GasStream.__molar_fractions_dic2vec__()
            n_gas = y_gas * n_gas_tot[:, :, None]
            cp_gas = self.GasStream.get_gas_heat_capacity_kJ_kmolK()
            T_gas = self.GasStream.get_gas_temp_K()

            # Exothermic Heat [kJ/kmol rxn] of the following reactions
            # - Gas Insta Reactions
            # - Gas Reversibe Reactions
            h_gas = np.zeros(shape=(self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1], self.GasStream.num_of_rxn_insta + self.GasStream.num_of_rxn_reversible),dtype=np.float64)
            for rxn_i, rxn_id in enumerate(self.GasStream.rxn_insta.keys()):
                h_gas[:,:, rxn_i] = self.GasStream.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
            for rxn_i, rxn_id in enumerate(self.GasStream.rxn_reversible.keys()):
                h_gas[:,:,rxn_i + self.GasStream.num_of_rxn_insta] = self.GasStream.get_rxn_reversible_exothermic_heat_kJ_kmol(rxn_id)

            # Calculate Rate-Reactions in Gas Phase
            r_gas = self.GasStream.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

            # Differential Equations, Equilibrium Reactions not taken into Account
            dn_gas__dz = A * (aJ_specie_gas + epsilon * (1 - liquid_holdup_m3_m3[:, :, None]) * r_gas)
            dn_gas_tot__dz = np.sum(dn_gas__dz, axis=2, keepdims=False)

            # Differential Equations, Equilibrium Reactions not taken into Account
            #dn_gas__dz = self.cross_sectional_area_m2[:, None, None] * (aJ_specie_gas + self.void_fraction_m3_m3[:, None, None] * (1 - liquid_holdup_m3_m3[:, :, None]) * r_gas)
            #dn_gas_tot__dz = np.sum(dn_gas__dz, axis=2, keepdims=False)

            # Sensitivity in Gas Phase (dndB, dndT)
            n_gas_broadcast = np.broadcast_to(array=n_gas[:, :, :, None], shape=(n_gas.shape[0], n_gas.shape[1], n_gas.shape[2], n_gas.shape[2]))
            dy_gas__dn_gas = - np.einsum("zsmn,zs->zsmn", n_gas_broadcast, 1 / n_gas_tot ** 2) + np.einsum("vw,zs->zsvw", np.eye(N=self.GasStream.num_of_species), 1 / n_gas_tot)
            df1__dy_gas = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(self.GasStream)
            df1__dn_gas = np.einsum("zsfw,zswv->zsfv", df1__dy_gas, dy_gas__dn_gas)
            df2__dn_gas = np.broadcast_to(array=self.matrix_gas["A"][None, None, :, :], shape=(
            self.GasStream.temp_K.shape[0], self.GasStream.temp_K.shape[1], self.matrix_gas["A"].shape[0],
            self.matrix_gas["A"].shape[1])).copy()
            df__dn_gas = np.concatenate((df1__dn_gas, df2__dn_gas), axis=2)
            H = np.linalg.inv(df__dn_gas)
            dn_gas__dB = H[:, :, :, self.GasStream.num_of_rxn_insta::]
            Ky = self.__get_GasEquilibrium_Ky__(self.GasStream)
            dKydT = self.__get_GasEquilibrium_dKydT__(self.GasStream, Ky)
            Hr = H[:, :, :, :self.GasStream.num_of_rxn_insta:]
            dn_gas__dT = - np.einsum("zscr,zsr->zsc", Hr, (1 / Ky) * dKydT)

            # Change in "B-Vector" for Gas Phase
            dB__dz = np.einsum("bw,zsw->zsb", self.matrix_gas["A"], dn_gas__dz)

            # Temperature Gradient
            dn_gas_at_eq__dz = dn_gas__dz
            for _ in range(10):
                q_lat = np.einsum("zsr,zsr->zs", h_gas, np.einsum("rm,zsm->zsr", self.matrix_gas_heat["R+"], dn_gas_at_eq__dz)) / A
                dT_gas__dz = A * (1 / (cp_gas * n_gas_tot)) * (-q_dir + q_des + q_lat)
                dn_gas_at_eq__dz = np.einsum("zswb,zsb->zsw", dn_gas__dB, dB__dz) + np.einsum("zsm,zs->zsm", dn_gas__dT, dT_gas__dz)
            dn_gas__dz = dn_gas_at_eq__dz

            """""""""
            Now we have obtained following derivatives.
            1) dm_liq__dz
            2) dn_gas__dz
            3) dT_liq__dz
            4) dT_gas__dz
            """""""""

            # Calculate Step Size
            with np.errstate(divide='ignore', invalid='ignore'):
                dz_max_01 = 0.1 * np.min(np.nan_to_num(m_liq / np.abs(dm_liq__dz), nan=np.inf, posinf=np.inf, neginf=np.inf))
            with np.errstate(divide='ignore', invalid='ignore'):
                dz_max_02 = np.nan_to_num(5 / np.abs(dT_liq__dz), nan=np.inf, posinf=np.inf, neginf=np.inf)[0,0]
            with np.errstate(divide='ignore', invalid='ignore'):
                dz_max_03 = 0.1 * np.min(np.nan_to_num(n_gas / np.abs(dn_gas__dz), nan=np.inf, posinf=np.inf, neginf=np.inf))
            with np.errstate(divide='ignore', invalid='ignore'):
                dz_max_04 = np.nan_to_num(5 / np.abs(dT_gas__dz), nan=np.inf, posinf=np.inf, neginf=np.inf)[0,0]
            dz_max_05 = self.position_m[-1] - self.z
            dz = np.min(a=[dz_max_01, dz_max_02, dz_max_03, dz_max_04, dz_max_05, self.step_size_m])

            # Integrate
            dm_liq = dm_liq__dz * dz
            dn_gas = dn_gas__dz * dz

            dT_gas = dT_gas__dz * dz
            dT_liq = dT_liq__dz * dz

            dm_liq_tot = np.sum(dm_liq, axis=2, keepdims=False)
            dn_gas_tot = np.sum(dn_gas, axis=2, keepdims=False)

            y_gas_new = (n_gas + dn_gas) / (n_gas_tot + dn_gas_tot)[:, :, None]
            w_liq_new = (m_liq + dm_liq) / (m_liq_tot + dm_liq_tot)[:, :, None]

            self.LiquidStream.temp_K = self.LiquidStream.temp_K + dT_liq
            self.LiquidStream.flow_kg_h = self.LiquidStream.flow_kg_h + 3600 * dm_liq_tot
            self.LiquidStream.__mass_fractions_vec2dic__(w=w_liq_new)

            self.GasStream.temp_K = self.GasStream.temp_K + dT_gas
            self.GasStream.flow_kmol_h = self.GasStream.flow_kmol_h + 3600 * dn_gas_tot
            self.GasStream.__molar_fractions_vec2dic__(y=y_gas_new)

            self.z = self.z + dz

            # Perform One Single (Isothermal) Equilibrium Step for the sake of Numerical Stability........
            w = self.LiquidStream.__mass_fractions_dic2vec__()
            Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidStream)
            f = self.__get_LiquidEquilibrium_f_rxn_insta_log__(self.LiquidStream, Kw)
            dfdw = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidStream)
            dfdr = np.einsum("zsrc,cq->zsrq", dfdw, self.equilibrium_liquid_inlet.matrix["R"])
            #dr_newton = - np.linalg.solve(dfdr, f)
            dr_newton = - np.einsum("zscd,zsd->zsc", np.linalg.pinv(dfdr), f)
            dw_newton = np.einsum("ij,zsj->zsi", self.equilibrium_liquid_inlet.matrix["R"], dr_newton)
            dw = lr * dw_newton
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            w = w + dw * tau
            self.LiquidStream.__mass_fractions_vec2dic__(np.maximum(w, 10 ** (-18)))

            y = self.GasStream.__molar_fractions_dic2vec__()
            w = self.GasStream.__molefrac2massfrac__(y)
            Kp = self.__get_GasEquilibrium_Kp__(self.GasStream)
            f = self.__get_GasEquilibrium_f_rxn_insta_log__(self.GasStream, Kp)
            dfdy = self.__get_GasEquilibrium_dfdy_rxn_insta_log__(self.GasStream)
            dydw = self.GasStream.__dydw__(w)
            dfdw = np.einsum("zsfy,zsyw->zsfw", dfdy, dydw)
            dfdr = np.einsum("zsrc,cq->zsrq", dfdw, self.equilibrium_gas_inlet.matrix["R"])
            #dr_newton = - np.linalg.solve(dfdr, f)
            dr_newton = - np.einsum("zscd,zsd->zsc", np.linalg.pinv(dfdr), f)
            dw_newton = np.einsum("ij,zsj->zsi", self.equilibrium_gas_inlet.matrix["R"], dr_newton)
            dw = lr * dw_newton
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            w = w + dw * tau
            self.GasStream.__molar_fractions_vec2dic__(y=self.GasStream.__massfrac2molefrac__(w))

            # Append to Profile
            self.LP.temp_K = np.vstack((self.LP.temp_K, self.LiquidStream.temp_K))
            self.LP.flow_kg_h = np.vstack((self.LP.flow_kg_h, self.LiquidStream.flow_kg_h))
            for id in self.LiquidStream.specie.keys():
                self.LP.specie[id]["Mass Fraction"] = np.vstack((self.LP.specie[id]["Mass Fraction"], self.LiquidStream.specie[id]["Mass Fraction"]))
            for id in self.LiquidStream.info.keys():
                self.LP.info[id] = np.vstack((self.LP.info[id], self.LiquidStream.info[id]))
            self.LP.position_m.append(self.z)
            self.GP.temp_K = np.vstack((self.GP.temp_K, self.GasStream.temp_K))
            self.GP.flow_kmol_h = np.vstack((self.GP.flow_kmol_h, self.GasStream.flow_kmol_h))
            self.GP.pressure_bara = np.vstack((self.GP.pressure_bara, self.GasStream.pressure_bara))
            for id in self.GasStream.specie.keys():
                self.GP.specie[id]["Molar Fraction"] = np.vstack((self.GP.specie[id]["Molar Fraction"], self.GasStream.specie[id]["Molar Fraction"]))
            for id in self.GasStream.info.keys():
                self.GP.info[id] = np.vstack((self.GP.info[id], self.GasStream.info[id]))
            self.GP.position_m.append(self.z)

        # From List to Numpy Array
        self.LP.position_m = np.array(self.LP.position_m)
        self.GP.position_m = np.array(self.GP.position_m)

        # Return Result
        self.GasStreamOut = self.__get_slice_from_GasProfile__(GasStreamProfile=self.GasStream, height_index=0, GasStreamOut=self.GasStream)
        self.LiquidStreamOut = self.__get_slice_from_LiquidProfile__(LiquidStreamProfile=self.LiquidStream, height_index=0, LiquidStreamOut=self.LiquidStream)

        self.GasStream = deepcopy(self.GP)
        self.LiquidStream = deepcopy(self.LP)

        del self.GP
        del self.LP

        return self.GasStreamOut, self.LiquidStreamOut



# ---------------------------------------------------------------------------------------


class Column_StructuredPacking_CounterCurrent(GasLiquidContactor_PFR_CounterCurrent):

    def __init__(self, height_m, num_of_heights, cross_sectional_area_m2, void_fraction_m3_m3, packing_area_m2_m3, corrugation_angle_degree):
        super().__init__(position_m=np.linspace(0, height_m, num_of_heights),
                         cross_sectional_area_m2=cross_sectional_area_m2 * np.ones(shape=(num_of_heights,)),
                         void_fraction_m3_m3=void_fraction_m3_m3 * np.ones(shape=(num_of_heights,)))
        self.packing_area_m2_m3 = packing_area_m2_m3
        self.corrugation_angle_degree = corrugation_angle_degree
        self.height_m = self.position_m

    def react(self, GasStreamIn, LiquidStreamIn, epochs, lr):
        GasStreamOut, LiquidStreamOut = self.__react__(GasStreamIn=GasStreamIn, LiquidStreamIn=LiquidStreamIn, epochs=epochs, lr=lr)
        return GasStreamOut, LiquidStreamOut

    def get_packing_area_m2_m3(self):
        return self.packing_area_m2_m3

    def get_corrugation_angle_degree(self):
        return self.corrugation_angle_degree


class Column_StructuredPacking_CoCurrent(GasLiquidContactor_PFR_CoCurrent):

    def __init__(self, height_m, cross_sectional_area_m2, void_fraction_m3_m3, packing_area_m2_m3, corrugation_angle_degree):
        super().__init__(position_m=np.array([0, height_m]),
                         cross_sectional_area_m2=np.array([cross_sectional_area_m2, cross_sectional_area_m2]),
                         void_fraction_m3_m3= np.array([void_fraction_m3_m3, void_fraction_m3_m3]))

        self.packing_area_m2_m3 = packing_area_m2_m3
        self.corrugation_angle_degree = corrugation_angle_degree

    def react(self, GasStreamIn, LiquidStreamIn, step_size_m, lr):
        GasStreamOut, LiquidStreamOut = self.__react__(GasStreamIn=GasStreamIn, LiquidStreamIn=LiquidStreamIn, step_size_m=step_size_m, lr=lr)
        self.height_m = np.max(self.GasStream.position_m) - self.GasStream.position_m
        return GasStreamOut, LiquidStreamOut

    def get_packing_area_m2_m3(self):
        return self.packing_area_m2_m3

    def get_corrugation_angle_degree(self):
        return self.corrugation_angle_degree



# ---------------------------------------------------------------------------------------


class _Recycle_Bin:


    class LiquidHeatExchanger_CounterCurrent(_Stochiometry, Serializer, _LiquidEquilibrium, _StreamFunctions):

        def __init__(self, interface_area_m2, volume1_m3=None, volume2_m3=None):
            self.firstscan = True
            self.volume1_m3 = 1.8 * 10 ** (-3) * interface_area_m2 if volume1_m3 is None else volume1_m3
            self.volume2_m3 = 1.8 * 10 ** (-3) * interface_area_m2 if volume2_m3 is None else volume2_m3
            self.interface_area_m2 = interface_area_m2

        def add_heat_transfer_coefficient_kW_m2K(self, heat_transfer_coefficient_kW_m2K):
            self._heat_transfer_coefficient_kW_m2K_ = heat_transfer_coefficient_kW_m2K

        def react(self, LiquidStreamIn1, LiquidStreamIn2, lr=0.75, ntu_method=True):
            if ntu_method:
                LiquidStreamOut1, LiquidStreamOut2 = self.__react_ntu__(LiquidStreamIn1, LiquidStreamIn2)
            else:
                LiquidStreamOut1, LiquidStreamOut2 = self.__react__(LiquidStreamIn1, LiquidStreamIn2, lr=lr)
            return LiquidStreamOut1, LiquidStreamOut2

        def __react__(self, LiquidStreamIn1, LiquidStreamIn2, lr):

            _, _ = self.__react_ntu__(LiquidStreamIn1, LiquidStreamIn2)

            # Geometry

            self.length1_m = np.ones(shape=self.volume1_m3.shape)  # Assumption
            self.length2_m = np.ones(shape=self.volume2_m3.shape)  # Assumption
            self.cross_sectional_area1_m2 = self.volume1_m3 / self.length1_m
            self.cross_sectional_area2_m2 = self.volume2_m3 / self.length2_m
            self.dz = 1 / self.num_of_heights

            # Strictly Positive Mass Fractions
            for id in LiquidStreamIn1.specie.keys():
                LiquidStreamIn1.specie[id]["Mass Fraction"] = np.maximum(LiquidStreamIn1.specie[id]["Mass Fraction"],
                                                                         10 ** (-18))
            LiquidStreamIn1.normalize_mass_fractions()

            for id in LiquidStreamIn2.specie.keys():
                LiquidStreamIn2.specie[id]["Mass Fraction"] = np.maximum(LiquidStreamIn2.specie[id]["Mass Fraction"],
                                                                         10 ** (-18))
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
                    self.LiquidStream2 = self.__get_slice_from_LiquidProfile__(self.LP2, i - 1, self.LiquidStream2)

                    # Calculate Gradients
                    dwdz, dTdz = self.__get_gradients_liq1__()

                    # Integrate
                    # self.LiquidStream1.__mass_fractions_vec2dic__(self.LiquidStream1.__mass_fractions_dic2vec__() + dwdz * dz[:, None])
                    self.LiquidStream1.temp_K = self.LiquidStream1.temp_K + dTdz * self.dz
                    # self.LiquidStream1.normalize_mass_fractions()

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
                    # self.LiquidStream2.__mass_fractions_vec2dic__(self.LiquidStream2.__mass_fractions_dic2vec__() + dwdz * dz[:, None])
                    self.LiquidStream2.temp_K = self.LiquidStream2.temp_K + dTdz * self.dz
                    # self.LiquidStream2.normalize_mass_fractions()

                    # Update Liquid Profile
                    self.LP2 = self.__insert_slice_to_LiquidProfile__(self.LP2, i, self.LiquidStream2)

            return self.LiquidStream1, self.LiquidStream2

        def __react_ntu__(self, LiquidStreamIn1, LiquidStreamIn2):

            # self.LP1 = self.__broadcast_to_LiquidProfile__(LiquidStreamIn1, self.num_of_heights)
            # self.LP2 = self.__broadcast_to_LiquidProfile__(LiquidStreamIn2, self.num_of_heights)

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
                # self.LiquidStream1.temp_K = 1.0 * LiquidStreamIn1.temp_K
                # self.LiquidStream2.temp_K = 1.0 * LiquidStreamOut2.temp_K
                # kH_left = self._heat_transfer_coefficient_kW_m2K_(self)

                # self.LiquidStream1.temp_K = 1.0 * LiquidStreamOut1.temp_K
                # self.LiquidStream2.temp_K = 1.0 * LiquidStreamIn2.temp_K
                # kH_right = self._heat_transfer_coefficient_kW_m2K_(self)

                # Average
                # kH = (kH_left + kH_right) / 2

                kH = self._heat_transfer_coefficient_kW_m2K_(self)

                # Number of Transfer Units
                NTU = kH * self.interface_area_m2 / Cmin

                # Actual Heat Transferred
                effectiveness = (1 - np.exp(-NTU * (1 - Cr))) / (1 - Cr * np.exp(- NTU * (1 - Cr)))
                Q = effectiveness * Qmax

                # Calculate Outlet Temperatures
                LiquidStreamOut1.temp_K = LiquidStreamIn1.temp_K - Q / C1
                LiquidStreamOut2.temp_K = LiquidStreamIn2.temp_K + Q / C2
                # dT_left = np.abs(LiquidStreamIn1.temp_K - LiquidStreamOut2.temp_K)
                # dT_right = np.abs(LiquidStreamOut1.temp_K - LiquidStreamIn2.temp_K)
                # LMTD = (dT_left - dT_right) / (np.log(dT_left) - np.log(dT_right))

            # for i in range(LiquidStreamIn1.temp_K.shape[0]):
            #    self.LP1.temp_K[:,i] = np.linspace(LiquidStreamIn1.temp_K[i], LiquidStreamOut1.temp_K[i], self.num_of_heights)
            #    self.LP2.temp_K[:, i] = np.linspace(LiquidStreamOut2.temp_K[i], LiquidStreamIn2.temp_K[i], self.num_of_heights)

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


    class LiquidCSTR_Batch_DYNAMIC(_Stochiometry, Serializer, _LiquidEquilibrium):

        def __init__(self):
            self.info = {}
            self.firstscan = True
            self.equilibrium_liquid_inlet = LiquidEquilibrium_Adiabatic()
            self.equilibrium_liquid_outlet = LiquidEquilibrium_Adiabatic()

        def react(self, LiquidBulkInit, step_size_s, terminal_time_s, isothermal=False, lr=0.75):

            # Step Size of Integrator
            self.step_size_s = step_size_s
            self.terminal_time_s = terminal_time_s

            # Only Strictly Positive Mass Fractions Allowed
            for id in LiquidBulkInit.specie.keys():
                LiquidBulkInit.specie[id]["Mass Fraction"] = np.maximum(LiquidBulkInit.specie[id]["Mass Fraction"],
                                                                        10 ** (-18))
            LiquidBulkInit.normalize_mass_fractions()

            # Process Liquid and Gas Inlet Streams
            self.LiquidBulkInit = self.equilibrium_liquid_inlet.react(LiquidBulkInit, lr=lr)

            if self.firstscan:

                # Extract Molar Masses
                self.M_liq = np.zeros(shape=(LiquidBulkInit.num_of_species,))
                for i, id in enumerate(LiquidBulkInit.specie.keys()):
                    self.M_liq[i] = LiquidBulkInit.get_specie_molar_mass_kg_kmol(id)

                # Matrix or Reaction Stochiometry: Used to Calculate Sensitivity Matrices
                self.matrix_liq = self.__get_the_matrix__(LiquidStreamIn=LiquidBulkInit, liq_rxn_insta=True)

                # Matrix of Reaction Stochiometry: Used to Calculate Heat of Reactions
                self.matrix_liq_heat = self.__get_the_matrix__(LiquidStreamIn=LiquidBulkInit,
                                                               liq_rxn_insta=True,
                                                               liq_rxn_reversible=True)

                self.firstscan = False

            # Initiate as Inlet
            self.LiquidBulk = deepcopy(self.LiquidBulkInit)
            self.LiquidBulk.temp_K = self.LiquidBulkInit.temp_K[None, :]
            self.LiquidBulk.quantum_kg = self.LiquidBulkInit.quantum_kg[None, :]
            for id in self.LiquidBulk.specie.keys():
                self.LiquidBulk.specie[id]["Mass Fraction"] = self.LiquidBulkInit.specie[id]["Mass Fraction"][None, :]
            for id in self.LiquidBulk.info.keys():
                self.LiquidBulk.info[id] = self.LiquidBulkInit.info[id][None, :]

            self.LP = deepcopy(self.LiquidBulk)
            self.t = 0
            self.LP.time_s = [self.t]

            while self.t < self.terminal_time_s:

                m_liq_tot = self.LiquidBulk.get_solution_quantum_kg()
                w_liq = self.LiquidBulk.__mass_fractions_dic2vec__()
                m_liq = w_liq * m_liq_tot[:, :, None]
                cp_liq = self.LiquidBulk.get_solution_heat_capacity_kJ_kgK()
                T_liq = self.LiquidBulk.get_solution_temp_K()
                rho = self.LiquidBulk.get_solution_density_kg_m3()
                V = m_liq_tot / rho

                # Exothermic Heat [kJ/kmol rxn] of the following reactions
                # - Liquid Insta Reactions
                # - Liquid Reversibe Reactions
                h = np.zeros(shape=(self.LiquidBulk.temp_K.shape[0], self.LiquidBulk.temp_K.shape[1],
                                    self.LiquidBulk.num_of_rxn_insta + self.LiquidBulk.num_of_rxn_reversible),
                             dtype=np.float64)
                for rxn_i, rxn_id in enumerate(self.LiquidBulk.rxn_insta.keys()):
                    h[:, :, rxn_i] = self.LiquidBulk.get_rxn_insta_exothermic_heat_kJ_kmol(rxn_id)
                for rxn_i, rxn_id in enumerate(self.LiquidBulk.rxn_reversible.keys()):
                    h[:, :,
                    rxn_i + self.LiquidBulk.num_of_rxn_insta] = self.LiquidBulk.get_rxn_reversible_exothermic_heat_kJ_kmol(
                        rxn_id)

                # Calculate Rate-Reactions in Liquid Phase
                # [r] = [kmol/m3.s]
                r_liq = self.LiquidBulk.get_rxn_reversible_rate_kmol_m3s_as_specie_vector()

                # Differential Equations (dm/dz)
                dw_liq__dt = self.M_liq[None, None, :] * r_liq / rho[:, :, None]  # 1/s = kg/kmol * kmol/m3.s * m3/kg

                # Sensitivity in Liquid Phase (dwdB, dwdT)
                df1__dw_liq = self.__get_LiquidEquilibrium_dfdw_rxn_insta_loq__(self.LiquidBulk)
                df2__dw_liq = np.broadcast_to(array=self.matrix_liq["A"][None, None, :, :], shape=(
                self.LiquidBulk.temp_K.shape[0], self.LiquidBulk.temp_K.shape[1], self.matrix_liq["A"].shape[0],
                self.matrix_liq["A"].shape[1])).copy()
                df__dw_liq = np.concatenate((df1__dw_liq, df2__dw_liq), axis=2)
                H = np.linalg.inv(df__dw_liq)
                Kw = self.__get_LiquidEquilibrium_Kw__(self.LiquidBulk)
                dKwdT = self.__get_LiquidEquilibrium_dKwdT__(self.LiquidBulk, Kw)
                Hr = H[:, :, :, :self.LiquidBulk.num_of_rxn_insta:]
                dw_liq__db = H[:, :, :, self.LiquidBulk.num_of_rxn_insta::]
                dw_liq__dT = - np.einsum("zscr,zsr->zsc", Hr, (1 / Kw) * dKwdT)

                # Change in "b-Vector" for Liquid Phase
                db__dt = np.einsum("bw,zsw->zsb", self.matrix_liq["A"], dw_liq__dt)

                # Iterate Until "Convergence"
                if isothermal:
                    dw_liq_at_eq__dt = np.einsum("zswb,zsb->zsw", dw_liq__db, db__dt)
                    dT_liq__dt = np.zeros(shape=(cp_liq.shape))
                else:
                    dw_liq_at_eq__dt = dw_liq__dt
                    for _ in range(10):
                        dT_liq__dt = (1 / cp_liq) * np.einsum("zsr,zsr->zs", h,
                                                              np.einsum("rm,zsm->zsr", self.matrix_liq_heat["R+"],
                                                                        dw_liq_at_eq__dt))
                        dw_liq_at_eq__dt = np.einsum("zswb,zsb->zsw", dw_liq__db, db__dt) + np.einsum("zsm,zs->zsm",
                                                                                                      dw_liq__dT,
                                                                                                      dT_liq__dt)

                dw_liq__dt = dw_liq_at_eq__dt

                """""""""
                Now we have obtained following derivatives.
                1) dm_liq__dz
                2) dT_liq__dz
                """""""""

                # Calculate Step Size
                with np.errstate(divide='ignore', invalid='ignore'):
                    dt_max_01 = 0.1 * np.min(
                        np.nan_to_num(w_liq / np.abs(dw_liq__dt), nan=np.inf, posinf=np.inf, neginf=np.inf))
                with np.errstate(divide='ignore', invalid='ignore'):
                    dt_max_02 = np.nan_to_num(5 / np.abs(dT_liq__dt), nan=np.inf, posinf=np.inf, neginf=np.inf)[0, 0]
                dt_max_03 = self.terminal_time_s - self.t
                dt = np.min(a=[dt_max_01, dt_max_02, dt_max_03, self.step_size_s])

                # Integrate
                dw_liq = dw_liq__dt * dt
                dT_liq = dT_liq__dt * dt

                self.LiquidBulk.temp_K = self.LiquidBulk.temp_K + dT_liq
                self.LiquidBulk.quantum_kg = self.LiquidBulk.quantum_kg + 0
                self.LiquidBulk.__mass_fractions_vec2dic__(w=w_liq + dw_liq)
                self.t = self.t + dt

                # Append to Profile
                self.LP.temp_K = np.vstack((self.LP.temp_K, self.LiquidBulk.temp_K))
                self.LP.quantum_kg = np.vstack((self.LP.quantum_kg, self.LiquidBulk.quantum_kg))
                for id in self.LiquidBulk.specie.keys():
                    self.LP.specie[id]["Mass Fraction"] = np.vstack(
                        (self.LP.specie[id]["Mass Fraction"], self.LiquidBulk.specie[id]["Mass Fraction"]))
                for id in self.LiquidBulk.info.keys():
                    self.LP.info[id] = np.vstack((self.LP.info[id], self.LiquidBulk.info[id]))
                self.LP.time_s.append(self.t)

            # From List to Numpy Array
            self.LP.time_s = np.array(self.LP.time_s)
            self.LiquidBulk = deepcopy(self.LP)
            del self.LP

            # Return Result
            self.LiquidBulkOut = deepcopy(LiquidBulkInit)
            self.LiquidBulkOut.temp_K = self.LiquidBulk.temp_K[-1, :]
            self.LiquidBulkOut.quantum_kg = self.LiquidBulk.quantum_kg[-1, :]
            for id in self.LiquidBulkOut.specie.keys():
                self.LiquidBulkOut.specie[id]["Mass Fraction"] = self.LiquidBulk.specie[id]["Mass Fraction"][-1, :]
            for id in self.LiquidBulkOut.info.keys():
                self.LiquidBulkOut.info[id] = self.LiquidBulk.info[id][-1, :]

            return self.LiquidBulkOut


    def GasBulkSum(bulks):

        GasBulkOut = deepcopy(bulks[0])
        GasBulkOut.temp_K = 0
        GasBulkOut.quantum_kmol = 0

        for bulk in bulks:
            GasBulkOut.quantum_kmol = GasBulkOut.quantum_kmol + np.abs(bulk.quantum_kmol)

        for id in bulks[0].specie.keys():
            n = 0
            for bulk in bulks:
                n = n + bulk.get_specie_quantum_kmol(id=id)
            GasBulkOut.specie[id]["Molar Fraction"] = n / GasBulkOut.quantum_kmol

        for bulk in bulks:
            T = bulk.temp_K
            r = np.abs(bulk.quantum_kmol) / GasBulkOut.quantum_kmol
            GasBulkOut.temp_K = GasBulkOut.temp_K + r * T

        GasBulkOut.normalize_molar_fractions()
        return GasBulkOut


    class GasBulk(_Gas):

        def __init__(self):
            super().__init__()
            self.quantum_kmol = None

        def set_gas_quantum_kmol(self, value):
            self.amount_kmol = value

        def get_gas_quantum_kmol(self):
            return self.amount_kmol

        def get_gas_quantum_kg(self):
            m = 0
            for id in self.specie.keys():
                y = self.get_specie_molar_fraction(id=id)
                M = self.get_specie_molar_mass_kg_kmol(id=id)
                m = m + y * M * self.get_gas_quantum_kmol()
            return m

        def get_specie_quantum_kmol(self, id):
            return self.get_gas_quantum_kmol() * self.get_specie_molar_fraction(id=id)

        def get_specie_quantum_kg(self, id):
            return self.get_specie_quantum_kmol(id=id) * self.get_specie_molar_mass_kg_kmol(id=id)


    def LiquidStream_to_LiquidBulk(LiquidStreamIn, quantity_kg):

        bulk = LiquidBulk(stream_id=LiquidStreamIn.id, solvent_id=LiquidStreamIn.solvent_id)

        bulk.info = LiquidStreamIn.info
        bulk.function = LiquidStreamIn.function

        bulk.specie = LiquidStreamIn.specie
        bulk.rxn_insta = LiquidStreamIn.rxn_insta
        bulk.rxn_reversible = LiquidStreamIn.rxn_reversible
        bulk.rxn_irreversible = LiquidStreamIn.rxn_irreversible
        bulk.vapor_pressure_bara = LiquidStreamIn.vapor_pressure_bara

        bulk.num_of_species = LiquidStreamIn.num_of_species
        bulk.num_of_rxn_insta = LiquidStreamIn.num_of_rxn_insta
        bulk.num_of_rxn_reversible = LiquidStreamIn.num_of_rxn_reversible
        bulk.num_of_rxn_irreversible = LiquidStreamIn.num_of_rxn_irreversible
        bulk.num_of_vapor_pressure_bara = LiquidStreamIn.num_of_vapor_pressure_bara

        bulk.temp_K = LiquidStreamIn.temp_K
        bulk.set_solution_quantum_kg(value=quantity_kg)

        bulk.density_kg_m3 = LiquidStreamIn.density_kg_m3
        bulk.heat_capacity_kJ_kgK = LiquidStreamIn.heat_capacity_kJ_kgK
        bulk.viscosity_Pas = LiquidStreamIn.viscosity_Pas
        bulk.thermal_conductivity_kW_mK = LiquidStreamIn.thermal_conductivity_kW_mK
        bulk.activity_coefficient = LiquidStreamIn.activity_coefficient
        bulk.diffusivity_m2_s = LiquidStreamIn.diffusivity_m2_s
        bulk.surface_tension_N_m = LiquidStreamIn.surface_tension_N_m

        return bulk


    def LiquidBulk_to_LiquidStream(LiquidBulkIn, flow_kg_h):
        stream = LiquidStream(stream_id=LiquidBulkIn.id, solvent_id=LiquidBulkIn.solvent_id)

        stream.info = LiquidBulkIn.info
        stream.function = LiquidBulkIn.function

        stream.specie = LiquidBulkIn.specie
        stream.rxn_insta = LiquidBulkIn.rxn_insta
        stream.rxn_reversible = LiquidBulkIn.rxn_reversible
        stream.rxn_irreversible = LiquidBulkIn.rxn_irreversible
        stream.vapor_pressure_bara = LiquidBulkIn.vapor_pressure_bara

        stream.num_of_species = LiquidBulkIn.num_of_species
        stream.num_of_rxn_insta = LiquidBulkIn.num_of_rxn_insta
        stream.num_of_rxn_reversible = LiquidBulkIn.num_of_rxn_reversible
        stream.num_of_rxn_irreversible = LiquidBulkIn.num_of_rxn_irreversible
        stream.num_of_vapor_pressure_bara = LiquidBulkIn.num_of_vapor_pressure_bara

        stream.temp_K = LiquidBulkIn.temp_K
        stream.set_solution_flow_kg_h(value=flow_kg_h)

        stream.density_kg_m3 = LiquidBulkIn.density_kg_m3
        stream.heat_capacity_kJ_kgK = LiquidBulkIn.heat_capacity_kJ_kgK
        stream.viscosity_Pas = LiquidBulkIn.viscosity_Pas
        stream.thermal_conductivity_kW_mK = LiquidBulkIn.thermal_conductivity_kW_mK
        stream.activity_coefficient = LiquidBulkIn.activity_coefficient
        stream.diffusivity_m2_s = LiquidBulkIn.diffusivity_m2_s
        stream.surface_tension_N_m = LiquidBulkIn.surface_tension_N_m

        return stream


    class LiquidBulk(_Liquid):

        def __init__(self, stream_id, solvent_id):
            super().__init__(stream_id, solvent_id)
            self.quantum_kg = None

        def set_solution_quantum_kg(self, value):
            self.quantum_kg = value.astype(np.float64)

        def get_solution_quantum_kg(self):
            return self.quantum_kg

        def get_solution_quantum_m3(self):
            return self.get_solution_quantum_kg() / self.get_solution_density_kg_m3()

        def get_solution_quantum_kmol(self):
            return self.get_solution_molarity_kmol_m3() * self.get_solution_quantum_m3()

        def get_specie_quantum_kg(self, id):
            return self.get_specie_mass_fraction(id=id) * self.get_solution_quantum_kg()

        def get_specie_quantum_kmol(self, id):
            return self.get_specie_molarity_kmol_m3(id=id) * self.get_solution_quantum_m3()


    def LiquidBulkSum(bulks):

        LiquidBulkOut = deepcopy(streams[0])
        LiquidBulkOut.temp_K = 0
        LiquidBulkOut.quantum_kg = 0

        for bulk in bulks:
            LiquidStreamBulk.quantum_kg = LiquidBulkOut.quantum_kg + np.abs(bulk.quantum_kg)

        for id in bulks[0].specie.keys():
            m = 0
            for bulk in bulkss:
                if id in bulk.specie.keys():
                    m = m + bulk.get_specie_quantum_kg(id=id)
            LiquidBulkOut.specie[id]["Mass Fraction"] = m / LiquidBulkOut.quantum_kg

        for bulk in bulks:
            T = bulk.temp_K
            r = np.abs(bulk.quantum_kg) / LiquidBulkOut.quantum_kg
            LiquidBulkOut.temp_K = LiquidBulkOut.temp_K + r * T

        LiquidBulkOut.normalize_mass_fractions()
        return LiquidBulkOut






