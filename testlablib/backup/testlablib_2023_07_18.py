import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from sklearn import neighbors
from sklearn.neighbors import BallTree
from sklearn.neighbors import KDTree
from sklearn.neighbors import KNeighborsRegressor
import time
import scipy.special



class RecycleBin:

    def __Kw__(self, LiquidStreamIn):

        num_of_samples = LiquidStreamIn.temp_K.shape[0]

        Kw = np.zeros(shape=(num_of_samples, LiquidStreamIn.eq["Number of Reactions"]), dtype=np.float64)
        for i in range(LiquidStreamIn.eq["Number of Reactions"]):
            Kw[:,i] = LiquidStreamIn.equilibrium_K[i](LiquidStreamIn)

        den = 0
        for id in LiquidStreamIn.specie.keys():
            den = den + LiquidStreamIn.get_specie_mass_fraction(id) / LiquidStreamIn.get_specie_molar_mass_kg_kmol(id)

        c = LiquidStreamIn.get_solution_molarity_kmol_m3()
        rho = LiquidStreamIn.get_solution_density_kg_m3()
        for i in range(LiquidStreamIn.eq["Number of Reactions"]):

            for r in LiquidStreamIn.equilibrium_reactants_type[i].keys():
                Kw[:, i] = Kw[:, i] * (1 / LiquidStreamIn.get_specie_activity_coefficient(id=r)) ** LiquidStreamIn.equilibrium_reactants[i][r]
                if LiquidStreamIn.equilibrium_reactants_type[i][r] == "c":
                    Kw[:,i] = Kw[:,i] * (c / (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r))) ** LiquidStreamIn.equilibrium_reactants[i][r]
                elif LiquidStreamIn.equilibrium_reactants_type[i][r] == "x":
                    Kw[:, i] = Kw[:, i] * (1 / (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r))) ** LiquidStreamIn.equilibrium_reactants[i][r]
                elif LiquidStreamIn.equilibrium_reactants_type[i][r] == "m":
                    Kw[:, i] = Kw[:, i] * ((1000 * c)/(rho * LiquidStreamIn.get_specie_mass_fraction(id=LiquidStreamIn.solvent_id) * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r) * den)) ** LiquidStreamIn.equilibrium_reactants[i][r]
                elif LiquidStreamIn.equilibrium_reactants_type[i][r] == None:
                    Kw[:, i] = Kw[:, i] * (1 / LiquidStreamIn.get_specie_mass_fraction(id=r)) ** LiquidStreamIn.equilibrium_reactants[i][r]

            for p in LiquidStreamIn.equilibrium_products_type[i].keys():
                Kw[:, i] = Kw[:, i] * LiquidStreamIn.get_specie_activity_coefficient(id=r) ** LiquidStreamIn.equilibrium_products[i][p]
                if LiquidStreamIn.equilibrium_products_type[i][p] == "c":
                    Kw[:,i] = Kw[:,i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(p) / c) ** LiquidStreamIn.equilibrium_products[i][p]
                elif LiquidStreamIn.equilibrium_products_type[i][p] == "x":
                    Kw[:, i] = Kw[:, i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(p)) ** LiquidStreamIn.equilibrium_products[i][p]
                elif LiquidStreamIn.equilibrium_products_type[i][p] == "m":
                    Kw[:, i] = Kw[:, i] * ((rho * LiquidStreamIn.get_specie_mass_fraction(id=LiquidStreamIn.solvent_id) * LiquidStreamIn.get_specie_molar_mass_kg_kmol(p) * den) / (1000 * c)) ** equilibrium.equilibrium_products[i][p]
                elif LiquidStreamIn.equilibrium_products_type[i][p] == None:
                    Kw[:, i] = Kw[:, i] * (1 / LiquidStreamIn.get_specie_mass_fraction(id=p)) ** LiquidStreamIn.equilibrium_products[i][p]

        return Kw

    def __f__(self, LiquidStreamIn, Kw):
        fr = np.zeros(shape=(self.num_of_samples, LiquidStreamIn.eq["Number of Reactions"]), dtype=np.float64)
        fp = np.zeros(shape=(self.num_of_samples, LiquidStreamIn.eq["Number of Reactions"]), dtype=np.float64)
        f = np.zeros(shape=(self.num_of_samples, LiquidStreamIn.eq["Number of Reactions"]), dtype=np.float64)
        for i in range(LiquidStreamIn.eq["Number of Reactions"]):
            fr[:, i] = Kw[:, i]
            fp[:, i] = np.ones(shape=(self.num_of_samples,))
            for id in LiquidStreamIn.equilibrium_reactants[i].keys():
                nu = np.abs(LiquidStreamIn.equilibrium_reactants[i][id])
                fr[:, i] = fr[:, i] * LiquidStreamIn.get_specie_mass_fraction(id=id) ** nu
            for id in LiquidStreamIn.equilibrium_products[i].keys():
                nu = np.abs(LiquidStreamIn.equilibrium_products[i][id])
                fp[:, i] = fp[:, i] * LiquidStreamIn.get_specie_mass_fraction(id=id) ** nu
            f[:, i] = fr[:, i] - fp[:, i]
        return f

    def __dfdw__(self, LiquidStreamIn, Kw):
        dfdw = np.zeros(shape=(self.num_of_samples, LiquidStreamIn.eq["Number of Reactions"], LiquidStreamIn.eq["Number of Species"]), dtype=np.float64)
        for i in range(LiquidStreamIn.eq["Number of Reactions"]):
            for id in LiquidStreamIn.specie.keys():
                j = LiquidStreamIn.specie[id]["Index"]
                if id in LiquidStreamIn.equilibrium_reactants[i].keys():
                    nu = np.abs(LiquidStreamIn.equilibrium_reactants[i][id])
                    dfdw[:, i, j] = Kw[:, i] * nu * LiquidStreamIn.get_specie_mass_fraction(id=id) ** (nu - 1)
                    for jd in LiquidStreamIn.equilibrium_reactants[i].keys():
                        if id != jd:
                            nu = np.abs(LiquidStreamIn.equilibrium_reactants[i][jd])
                            dfdw[:, i, j] = dfdw[:, i, j] * LiquidStreamIn.get_specie_mass_fraction(id=jd) ** nu
                elif id in LiquidStreamIn.equilibrium_products[i].keys():
                    nu = np.abs(LiquidStreamIn.equilibrium_products[i][id])
                    dfdw[:, i, j] = - nu * LiquidStreamIn.get_specie_mass_fraction(id=id) ** (nu - 1)
                    for jd in LiquidStreamIn.equilibrium_products[i].keys():
                        if id != jd:
                            nu = np.abs(LiquidStreamIn.equilibrium_products[i][jd])
                            dfdw[:, i, j] = dfdw[:, i, j] * LiquidStreamIn.get_specie_mass_fraction(id=jd) ** nu
                else:
                    dfdw[:, i, j] = np.zeros(shape=(self.num_of_samples,))
        return dfdw

    def __recycle_bin__(self):

        AIb = np.einsum("cr,shr->shc", self._AI_, self._b_)
        self._w0_ = np.einsum("cd,shd->shc", self._P_, self._w_ - AIb) + AIb
        self._w0_ = np.maximum(self._w0_, 0)
        for id in LiquidStreamOut.specie.keys():
            i = LiquidStreamOut.specie[id]["Index"]
            LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.maximum(self._w0_[:, :, i], 0.0))
        self._p_ = np.zeros(shape=(self.num_of_samples, self.num_of_heights, self.num_of_reaction))

        # Conserved Quantities
        w = np.zeros(shape=(self.num_of_samples, self.num_of_heights, self.num_of_species), dtype=np.float64)
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            w[:, :, i] = LiquidStreamIn.get_specie_mass_fraction(id)
        self._b_ = np.einsum("ij,shj->shi", self._A_, w)

        # Generate Matrix containing all the Stochiometric Coefficents in the Conservation Laws
        # self._A_ = np.zeros(shape=(self.num_of_conservation, self.num_of_species), dtype=np.float64)
        # for i in range(self.num_of_conservation):
        #    for id in LiquidStreamIn.specie.keys():
        #        j = LiquidStreamIn.specie[id]["Index"]
        #        if id in LiquidStreamIn.conservation_law[i].keys():
        #            nu = LiquidStreamIn.conservation_law[i][id]
        #            self._A_[i, j] = nu / LiquidStreamIn.get_specie_molar_mass_kg_kmol(id=id)
        #        else:
        #            self._A_[i, j] = 0

        # Inverse and Null Space of A and Inverse of A
        # self._AI_ = np.linalg.pinv(self._A_)
        # self._N_ = scipy.linalg.null_space(self._A_)
        # self._NI_ = np.linalg.pinv(self._N_)

        # Projection Matrix
        # NTN = np.linalg.pinv(np.einsum("cr,cq->rq",self._N_,self._N_))
        # NNTN = np.einsum("cr,rq->cq", self._N_, NTN)
        # self._P_ = np.einsum("cq,dq->cd", NNTN, self._N_)

        # Matrix for map between Concentration and Null Space
        # self._Q_ = np.einsum("rc,cd->rd", self._NI_, self._P_)

        # -----------------------------

        # Calculate Distances Between Samples in Database and New Dataset
        # Next, apply argmin function, which should return the index of the sample in the database which is closest to the sample in new dataset.
        # The length of the array "i" is then equal to the number of samples in the new array.
        # distance_b = cdist(b_new, b, "sqeuclidean")
        # distance_T = cdist(T_new, T, "sqeuclidean")
        # distance = distance_b + 10 ** (-6) * distance_T
        # distance = cdist(x_new, x, "sqeuclidean")

        # print("cdist")
        # start_time = time.process_time()
        # x = np.concatenate((b, 10**(-3) * T), axis=1)
        # x_new = np.concatenate((b_new, 10**(-3) * T_new), axis=1)
        # distance = cdist(x_new, x, "euclidean") # sqeuclidean
        # i = np.argmin(distance, axis=1)[:, None]
        # print(time.process_time() - start_time, "seconds")
        # print("---")

        # ---------------------------------------------------

        # ii = self.tree.query(x_new, k=3, return_distance=False, dualtree=True, breadth_first=False)
        # ii.shape = (200,3)
        # w_broadcast.shape = (200, 12939, 37)
        # w.shape = (200,37)

        # b_new_broadcast = np.broadcast_to(b_new[:, None, :],shape=(num_of_samples_new, num_of_samples_database, num_of_conservations))
        # T_new_broadcast = np.broadcast_to(T_new[:, None, :],shape=(num_of_samples_new, num_of_samples_database, 1))
        # b_broadcast = np.broadcast_to(b[None, :, :], shape=(num_of_samples_new, num_of_samples_database, num_of_conservations))
        # T_broadcast = np.broadcast_to(T[None, :, :],shape=(num_of_samples_new, num_of_samples_database, 1))
        # w_broadcast = np.broadcast_to(w[None,:,:], shape=(num_of_samples_new, num_of_samples_database, num_of_species))
        # dwdb_broadcast = np.broadcast_to(dwdb[None, :, :, :], shape=(num_of_samples_new, num_of_samples_database, num_of_species, num_of_conservations))
        # dwdT_broadcast = np.broadcast_to(dwdT[None, :, :, :], shape=(num_of_samples_new, num_of_samples_database, num_of_species, 1))

        # print(w_broadcast.shape)
        # print("---")

        # print("Take Along Axis")
        # start_time = time.process_time()
        ## Values in New Array is Compared to the closest point in the Database
        # db = np.take_along_axis(b_new_broadcast - b_broadcast, i[:, :, None], axis=1)
        # dT = np.take_along_axis(T_new_broadcast - T_broadcast, i[:, :, None], axis=1)
        # dwdb = np.take_along_axis(dwdb_broadcast, i[:, :, None, None], axis=1)
        # dwdT = np.take_along_axis(dwdT_broadcast, i[:, :, None, None], axis=1)
        # w = np.take_along_axis(w_broadcast, i[:, :, None], axis=1)
        # print(time.process_time() - start_time, "seconds")
        # print("---")

    @staticmethod
    def __equilibrium_isothermal_f__(LiquidStreamIn, Kw):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        fr = np.zeros(shape=(num_of_samples, LiquidStreamIn.config_isothermal["Number of Reactions"]), dtype=np.float64)
        fp = np.zeros(shape=(num_of_samples, LiquidStreamIn.config_isothermal["Number of Reactions"]), dtype=np.float64)
        f = np.zeros(shape=(num_of_samples, LiquidStreamIn.config_isothermal["Number of Reactions"]), dtype=np.float64)
        for i in range(LiquidStreamIn.config_isothermal["Number of Reactions"]):
            fr[:, i] = Kw[:, i]
            fp[:, i] = np.ones(shape=(num_of_samples,))
            for id in LiquidStreamIn.equilibrium_reactants[i].keys():
                nu = np.abs(LiquidStreamIn.equilibrium_reactants[i][id])
                fr[:, i] = fr[:, i] * LiquidStreamIn.get_specie_mass_fraction(id=id) ** nu
            for id in LiquidStreamIn.equilibrium_products[i].keys():
                nu = np.abs(LiquidStreamIn.equilibrium_products[i][id])
                fp[:, i] = fp[:, i] * LiquidStreamIn.get_specie_mass_fraction(id=id) ** nu
            f[:, i] = fr[:, i] - fp[:, i]
        return f

    @staticmethod
    def __equilibrium_isothermal_dfdw__(LiquidStreamIn, Kw):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        dfdw = np.zeros(shape=(num_of_samples, LiquidStreamIn.config_isothermal["Number of Reactions"],
                               LiquidStreamIn.config_isothermal["Number of Species"]), dtype=np.float64)
        for i in range(LiquidStreamIn.config_isothermal["Number of Reactions"]):
            for id in LiquidStreamIn.specie.keys():
                j = LiquidStreamIn.specie[id]["Index"]
                if id in LiquidStreamIn.equilibrium_reactants[i].keys():
                    nu = np.abs(LiquidStreamIn.equilibrium_reactants[i][id])
                    dfdw[:, i, j] = Kw[:, i] * nu * LiquidStreamIn.get_specie_mass_fraction(id=id) ** (nu - 1)
                    for jd in LiquidStreamIn.equilibrium_reactants[i].keys():
                        if id != jd:
                            nu = np.abs(LiquidStreamIn.equilibrium_reactants[i][jd])
                            dfdw[:, i, j] = dfdw[:, i, j] * LiquidStreamIn.get_specie_mass_fraction(id=jd) ** nu
                elif id in LiquidStreamIn.equilibrium_products[i].keys():
                    nu = np.abs(LiquidStreamIn.equilibrium_products[i][id])
                    dfdw[:, i, j] = - nu * LiquidStreamIn.get_specie_mass_fraction(id=id) ** (nu - 1)
                    for jd in LiquidStreamIn.equilibrium_products[i].keys():
                        if id != jd:
                            nu = np.abs(LiquidStreamIn.equilibrium_products[i][jd])
                            dfdw[:, i, j] = dfdw[:, i, j] * LiquidStreamIn.get_specie_mass_fraction(id=jd) ** nu
                else:
                    dfdw[:, i, j] = np.zeros(shape=(num_of_samples,))
        return dfdw

    @staticmethod
    def __equilibrium_isothermal_dFdw__(LiquidStreamIn, Kw):
        dfdw = LiquidStreamIn.__equilibrium_isothermal_dfdw__(LiquidStreamIn, Kw)
        A = LiquidStreamIn.eq["A"]
        A = np.reshape(A, newshape=(1, A.shape[0], A.shape[1]))
        A = np.broadcast_to(A, shape=(dfdw.shape[0], A.shape[1], A.shape[2])).copy()
        dFdw = np.concatenate((dfdw, A), axis=1)
        return dFdw

    @staticmethod
    def __equilibrium_adiabatic_f__(LiquidStreamIn, LiquidStreamInit, Kw, dH):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        cp = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        dw = np.zeros(shape=(num_of_samples, LiquidStreamIn.config_isothermal["Number of Species"]))
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            dw[:, i] = LiquidStreamIn.get_specie_mass_fraction(id=id) - LiquidStreamInit.get_specie_mass_fraction(id=id)
        dT = LiquidStreamIn.temp_K - LiquidStreamInit.temp_K
        f1 = LiquidStreamIn.__equilibrium_isothermal_f__(LiquidStreamIn, Kw)
        f2 = (1 / cp) * np.einsum("sr,sr->s", dH,
                                  np.einsum("rw,sw->sr", LiquidStreamIn.config_isothermal["R+"], dw)) - dT
        f = np.concatenate((f1, f2[:, None]), axis=1)
        return f

    @staticmethod
    def __equilibrium_adiabatic_dfdr__(LiquidStreamIn, Kw, dH, dKdT):
        df1dw = LiquidStreamIn.__equilibrium_isothermal_dfdw__(LiquidStreamIn, Kw)
        df1dr = np.einsum("src,cq->srq", df1dw, LiquidStreamIn.config_isothermal["R"])
        cp = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        df2dr = dH / cp[:, None]
        dfdr = np.concatenate((df1dr, df2dr[:, None, :]), axis=1)
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        dfdT = - np.ones(shape=(num_of_samples, LiquidStreamIn.config_isothermal["Number of Reactions"] + 1),
                         dtype=np.float64)
        for i in range(LiquidStreamIn.config_isothermal["Number of Reactions"]):
            dfdT[:, i] = dKdT[:, i]
            for id in LiquidStreamIn.equilibrium_reactants[i].keys():
                nu = np.abs(LiquidStreamIn.equilibrium_reactants[i][id])
                dfdT[:, i] = dfdT[:, i] * LiquidStreamIn.get_specie_mass_fraction(id=id) ** nu
        dfdr = np.concatenate((dfdr, dfdT[:, :, None]), axis=2)
        return dfdr

    @staticmethod
    def __equilibrium_adiabatic_dfdr_old__(LiquidStreamIn, Kw, dH, dKdT):
        df1dw = LiquidStreamIn.__equilibrium_isothermal_dfdw__(LiquidStreamIn, Kw)
        df1dr = np.einsum("src,cq->srq", df1dw, LiquidStreamIn.RLL)
        cp = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        df2dr = dH / cp[:, None]
        dfdr = np.concatenate((df1dr, df2dr[:, None, :]), axis=1)
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        dfdT = - np.ones(shape=(num_of_samples, LiquidStreamIn.num_of_reactions + 1), dtype=np.float64)
        for i in range(LiquidStreamIn.num_of_reactions):
            dfdT[:, i] = (1 / Kw[:, i]) * dKdT[:, i]
        dfdr = np.concatenate((dfdr, dfdT[:, :, None]), axis=2)
        return dfdr


    # ------ FROM LIQUIDSTREAM ---------------

    @staticmethod
    def __Kw_vle__(LiquidStreamIn):

        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        Kw = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_volatile), dtype=np.float64)
        for i in range(LiquidStreamIn.num_of_volatile):
            Kw[:, i] = LiquidStreamIn.vle_K[i](LiquidStreamIn)
        den = 0

        for id in LiquidStreamIn.specie.keys():
            den = den + LiquidStreamIn.get_specie_mass_fraction(id) / LiquidStreamIn.get_specie_molar_mass_kg_kmol(id)
        c = LiquidStreamIn.get_solution_molarity_kmol_m3()
        rho = LiquidStreamIn.get_solution_density_kg_m3()

        for i in range(LiquidStreamIn.num_of_volatile):
            for p in LiquidStreamIn.vle_liq_products_type[i].keys():

                Kw[:, i] = Kw[:, i] * LiquidStreamIn.get_specie_activity_coefficient(id=p) ** \
                           LiquidStreamIn.vle_liq_products[i][p]

                if LiquidStreamIn.vle_liq_products_type[i][p] == "c":
                    Kw[:, i] = Kw[:, i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(p) / c) ** \
                               LiquidStreamIn.vle_liq_products[i][p]
                elif LiquidStreamIn.vle_liq_products_type[i][p] == "x":
                    Kw[:, i] = Kw[:, i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(p)) ** \
                               LiquidStreamIn.vle_liq_products[i][p]
                elif LiquidStreamIn.vle_liq_products_type[i][p] == "m":
                    Kw[:, i] = Kw[:, i] * ((rho * LiquidStreamIn.get_specie_mass_fraction(
                        id=LiquidStreamIn.solvent_id) * LiquidStreamIn.get_specie_molar_mass_kg_kmol(p) * den) / (
                                                       1000 * c)) ** LiquidStreamIn.vle_liq_products[i][p]
                elif LiquidStreamIn.vle_liq_products_type[i][p] == None:
                    Kw[:, i] = Kw[:, i] * (1 / LiquidStreamIn.get_specie_mass_fraction(id=p)) ** \
                               LiquidStreamIn.vle_liq_products[i][p]
        return Kw

    @staticmethod
    def __Kw_and_heat_of_reaction_kJ_kmol__(LiquidStreamIn):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        dT = 0.1

        Kw1 = LiquidStreamIn.__Kw__(LiquidStreamIn)
        ln_K1 = np.log(Kw1)
        T1 = np.broadcast_to(LiquidStreamIn.temp_K[:, None], shape=(num_of_samples, LiquidStreamIn.num_of_reactions))
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + dT

        Kw2 = LiquidStreamIn.__Kw__(LiquidStreamIn)
        ln_K2 = np.log(Kw2)
        T2 = np.broadcast_to(LiquidStreamIn.temp_K[:, None], shape=(num_of_samples, LiquidStreamIn.num_of_reactions))
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - dT

        dH = 8.314 * (ln_K2 - ln_K1) / ((1 / T2) - (1 / T1))
        dKdT = (Kw2 - Kw1) / dT

        return Kw1, dH, dKdT

    @staticmethod
    def __heat_of_vaporization_kJ_kmol__(LiquidStreamIn):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        T1 = LiquidStreamIn.temp_K
        p1 = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_volatile), dtype=np.float64)
        for i in range(LiquidStreamIn.num_of_volatile):
            p1[:, i] = LiquidStreamIn.vle_K[i](LiquidStreamIn)
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + 0.1

        T2 = LiquidStreamIn.temp_K
        p2 = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_volatile), dtype=np.float64)
        for i in range(LiquidStreamIn.num_of_volatile):
            p2[:, i] = LiquidStreamIn.vle_K[i](LiquidStreamIn)
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - 0.1

        dH = 8.314 * (np.log(p2) - np.log(p1)) / (1 / T2[:, None] - 1 / T1[:, None])
        return dH

    @staticmethod
    def __equilibrium_isothermal_heat_kW__(LiquidStreamIn, LiquidStreamOut):
        Kw1, dH, dKdT = LiquidStreamOut.__Kw_and_heat_of_reaction_kJ_kmol__(LiquidStreamOut)
        c = LiquidStreamOut.get_solution_molarity_kmol_m3()
        rho = LiquidStreamOut.get_solution_density_kg_m3()
        m = LiquidStreamOut.get_solution_flow_kg_h()
        dw_M = np.zeros(shape=(self.num_of_samples, LiquidStreamIn.eq["Number of Species"]))
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            dw_M[:, i] = (LiquidStreamOut.get_specie_mass_fraction(id=id) - LiquidStreamIn.get_specie_mass_fraction(
                id=id)) / LiquidStreamOut.get_specie_molar_mass_kg_kmol(id)
        den = 0
        for id in LiquidStreamOut.specie.keys():
            den = den + LiquidStreamOut.get_specie_mass_fraction(id=id) / LiquidStreamOut.get_specie_molar_mass_kg_kmol(
                id)
        Q = (m / 3600) * (1 / den) * (c / rho) * np.einsum("si,si->s",
                                                           np.einsum("sr,ri->si", dH, LiquidStreamIn.eq["R+"]), dw_M)
        return Q

    # ---------------------------------------------------------------------------

    def __eq_get_conserved_quantities__(self):
        num_of_samples = self.temp_K.shape[0]
        w = np.zeros(shape=(num_of_samples, self.num_of_species), dtype=np.float64)
        for id in self.specie.keys():
            i = self.specie[id]["Index"]
            w[:, i] = self.get_specie_mass_fraction(id=id)
        b = np.einsum("ic,sc->si", self.A, w).astype(np.float64)
        return b

    def __eq_get_solution_vector__(self):
        num_of_samples = self.temp_K.shape[0]
        w = np.zeros(shape=(num_of_samples, self.eq["Number of Species"]), dtype=np.float64)
        for id in self.specie.keys():
            i = self.specie[id]["Index"]
            w[:, i] = self.get_specie_mass_fraction(id=id)
        b = self.__eq_get_conserved_quantities__()
        Zb = np.einsum("ci,si->sc", self.eq["Z"], b)
        r = np.einsum("rc,sc->sr", self.eq["R+"], w - Zb).astype(np.float64)
        return r

    def __load_equilibrium_species__(self):

        # List all Species included in at least one equilibrium constraint
        reactive_species = []
        for j, jd in enumerate(d.keys()):
            for r in d[jd]["Stoch"].keys():
                if r not in reactive_species:
                    reactive_species.append(r)

        # List all species not included in any equilibrium constraints
        non_reactive_species = []
        for id in GasStreamIn.specie.keys():
            if id not in reactive_species:
                non_reactive_species.append(id)

        # One-Hot Array, Reactive Species
        reactive_species_flag = np.zeros(shape=(num_of_species,))
        for id in GasStreamIn.specie.keys():
            i = GasStreamIn.specie[id]["Index"]
            if id in reactive_species:
                reactive_species_flag[i] = 1.0


        # --------------------------------------------------------------------------------

        # List all Species included in at least one equilibrium constraint
        reactive_species = []
        for j, jd in enumerate(d.keys()):
            for r in d[jd]["Stoch"].keys():
                if r not in reactive_species:
                    reactive_species.append(r)

        # List all species not included in any equilibrium constraints
        non_reactive_species = []
        for id in LiquidStreamIn.specie.keys():
            if id not in reactive_species:
                non_reactive_species.append(id)

        # One-Hot Array, Reactive Species
        reactive_species_flag = np.zeros(shape=(num_of_species,))
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            if id in reactive_species:
                reactive_species_flag[i] = 1.0


class FeatureBooster():

    def __init__(self):
        pass

    def save(self, instance, filename):
        file = open(filename + ".obj","wb")
        pickle.dump(instance, file)
        file.close()

    def load(self, instance, filename):
        file = open(filename + ".obj", 'rb')
        instance = pickle.load(file)
        file.close()
        return instance


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

    def add_LiquidStream_rxn_reversible(self, id, stoch, rate_kmol_m3s):
        self.LiquidStream_rxn_reversible[id] = {}
        self.LiquidStream_rxn_reversible[id]["Stoch"] = stoch
        self.LiquidStream_rxn_reversible[id]["Rate [kmol/m3.s]"] = rate_kmol_m3s

    def add_LiquidStream_rxn_irreversible(self, id, stoch, rate_kmol_m3s, exothermic_heat_kJ_kmol):
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
        self.surface_tension[function_id] = function

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


class __Reactor_Common__(FeatureBooster):

    def __init__(self):
        super().__init__()
        self.w_min_init = 10**(-16)
        self.w_min = 10**(-16)
        self.dw_min_td = 10**(-16)

    # -----------------------

    def __load_stochiometry_liq_equilibrium__(self, LiquidStreamIn):
        matrix = {}
        e = LiquidStreamIn.equilibrium_liq
        num_of_species = len(LiquidStreamIn.specie.keys())
        num_of_equilibriums = len(e.keys())
        num_of_conserved_quantities = num_of_species - num_of_equilibriums
        matrix["R"] = np.zeros(shape=(num_of_species, num_of_equilibriums), dtype=np.float64)
        for j, jd in enumerate(e.keys()):
            for id in LiquidStreamIn.specie.keys():
                i = LiquidStreamIn.specie[id]["Index"]
                if id in e[jd]["Stoch"].keys():
                    nu = e[jd]["Stoch"][id]
                    matrix["R"][i, j] = nu * LiquidStreamIn.get_specie_molar_mass_kg_kmol(id=id)
        matrix["R+"] = np.linalg.pinv(matrix["R"]).astype(np.float64)
        matrix["P"] = np.einsum("cq,dq->cd", np.einsum("cr,rq->cq", matrix["R"], np.linalg.pinv(np.einsum("cr,cq->rq", matrix["R"], matrix["R"]))), matrix["R"]).astype(np.float64)
        U, D, _ = np.linalg.svd(matrix["R"], full_matrices=True)
        matrix["A"] = U[:, num_of_equilibriums::].T.astype(np.float64)
        return matrix

    def __load_stochiometry_liq_reaction__(self, LiquidStreamIn):
        matrix = {}
        r = LiquidStreamIn.reaction_liq
        num_of_species = len(LiquidStreamIn.specie.keys())
        num_of_reactions = len(r.keys())
        num_of_conserved_quantities = num_of_species - num_of_reactions
        matrix["R"] = np.zeros(shape=(num_of_species, num_of_reactions), dtype=np.float64)
        for j, jd in enumerate(r.keys()):
            for id in LiquidStreamIn.specie.keys():
                i = LiquidStreamIn.specie[id]["Index"]
                if id in r[jd]["Stoch"].keys():
                    nu = r[jd]["Stoch"][id]
                    matrix["R"][i, j] = nu * LiquidStreamIn.get_specie_molar_mass_kg_kmol(id=id)
        matrix["R+"] = np.linalg.pinv(matrix["R"]).astype(np.float64)
        matrix["P"] = np.einsum("cq,dq->cd", np.einsum("cr,rq->cq", matrix["R"], np.linalg.pinv( np.einsum("cr,cq->rq", matrix["R"], matrix["R"]))), matrix["R"]).astype(np.float64)
        U, D, _ = np.linalg.svd(matrix["R"], full_matrices=True)
        matrix["A"] = U[:, num_of_reactions::].T.astype(np.float64)
        return matrix

    def __load_stochiometry_liq_all__(self, LiquidStreamIn):
        R1 = self.__load_stochiometry_liq_equilibrium__(LiquidStreamIn)["R"]
        R2 = self.__load_stochiometry_liq_reaction__(LiquidStreamIn)["R"]

        num_of_equilibriums = R1.shape[1]
        num_of_reactions = R2.shape[1]

        matrix = {}
        matrix["R"] = np.concatenate((R1, R2), axis=1)
        matrix["R+"] = np.linalg.pinv(matrix["R"]).astype(np.float64)
        matrix["P"] = np.einsum("cq,dq->cd", np.einsum("cr,rq->cq", matrix["R"], np.linalg.pinv(np.einsum("cr,cq->rq", matrix["R"], matrix["R"]))), matrix["R"]).astype(np.float64)
        U, D, _ = np.linalg.svd(matrix["R"], full_matrices=True)
        matrix["A"] = U[:, num_of_equilibriums + num_of_reactions::].T.astype(np.float64)
        return matrix

    # ---

    def __load_stochiometry_gas_equilibrium__(self, GasStreamIn):
        matrix = {}
        e = GasStreamIn.equilibrium_gas
        num_of_species = len(GasStreamIn.specie.keys())
        num_of_equilibriums = len(e.keys())
        num_of_conserved_quantities = num_of_species - num_of_equilibriums
        matrix["R"] = np.zeros(shape=(num_of_species, num_of_equilibriums), dtype=np.float64)
        for j, jd in enumerate(e.keys()):
            for id in GasStreamIn.specie.keys():
                i = GasStreamIn.specie[id]["Index"]
                if id in e[jd]["Stoch"].keys():
                    nu = e[jd]["Stoch"][id]
                    matrix["R"][i, j] = nu * GasStreamIn.get_specie_molar_mass_kg_kmol(id=id)
        matrix["R+"] = np.linalg.pinv(matrix["R"]).astype(np.float64)
        matrix["P"] = np.einsum("cq,dq->cd", np.einsum("cr,rq->cq", matrix["R"], np.linalg.pinv(np.einsum("cr,cq->rq", matrix["R"], matrix["R"]))), matrix["R"]).astype(np.float64)
        U, D, _ = np.linalg.svd(matrix["R"], full_matrices=True)
        matrix["A"] = U[:, num_of_equilibriums::].T.astype(np.float64)
        return matrix

    def __load_stochiometry_gas_reaction__(self, GasStreamIn):
        matrix = {}
        r = GasStreamIn.reaction_liq
        num_of_species = len(GasStreamIn.specie.keys())
        num_of_reactions = len(r.keys())
        num_of_conserved_quantities = num_of_species - num_of_reactions
        matrix["R"] = np.zeros(shape=(num_of_species, num_of_reactions), dtype=np.float64)
        for j, jd in enumerate(r.keys()):
            for id in GasStreamIn.specie.keys():
                i = GasStreamIn.specie[id]["Index"]
                if id in r[jd]["Stoch"].keys():
                    nu = r[jd]["Stoch"][id]
                    matrix["R"][i, j] = nu * GasStreamIn.get_specie_molar_mass_kg_kmol(id=id)
        matrix["R+"] = np.linalg.pinv(matrix["R"]).astype(np.float64)
        matrix["P"] = np.einsum("cq,dq->cd", np.einsum("cr,rq->cq", matrix["R"], np.linalg.pinv(np.einsum("cr,cq->rq", matrix["R"], matrix["R"]))), matrix["R"]).astype(np.float64)
        U, D, _ = np.linalg.svd(matrix["R"], full_matrices=True)
        matrix["A"] = U[:, num_of_reactions::].T.astype(np.float64)
        return matrix

    def __load_stochiometry_gas_all__(self, GasStreamIn):
        R1 = self.__load_stochiometry_gas_equilibrium__(GasStreamIn)["R"]
        R2 = self.__load_stochiometry_gas_reaction__(GasStreamIn)["R"]

        num_of_equilibriums = R1.shape[1]
        num_of_reactions = R2.shape[1]

        matrix = {}
        matrix["R"] = np.concatenate((R1, R2), axis=1)
        matrix["R+"] = np.linalg.pinv(matrix["R"]).astype(np.float64)
        matrix["P"] = np.einsum("cq,dq->cd", np.einsum("cr,rq->cq", matrix["R"], np.linalg.pinv(np.einsum("cr,cq->rq", matrix["R"], matrix["R"]))), matrix["R"]).astype(np.float64)
        U, D, _ = np.linalg.svd(matrix["R"], full_matrices=True)
        matrix["A"] = U[:, num_of_equilibriums + num_of_reactions::].T.astype(np.float64)
        return matrix

    # ---

    def __load_stochiometry_vap_equilibrium__(self, GasStreamIn, LiquidStreamIn):
        matrix = {}
        d = LiquidStreamIn.equilibrium_vap
        num_of_species_liq = len(LiquidStreamIn.specie.keys())
        num_of_species_gas = len(GasStreamIn.specie.keys())
        num_of_reactions = len(d.keys())
        matrix["RL"] = np.zeros(shape=(num_of_species_liq, num_of_reactions), dtype=np.float64)
        matrix["RG"] = np.zeros(shape=(num_of_species_gas, num_of_reactions), dtype=np.float64)
        for j, jd in enumerate(d.keys()):
            for id in LiquidStreamIn.specie.keys():
                i = LiquidStreamIn.specie[id]["Index"]
                if id in d[jd]["Stoch Liq"].keys():
                    nu = d[jd]["Stoch Liq"][id]
                    matrix["RL"][i, j] = nu * LiquidStreamIn.get_specie_molar_mass_kg_kmol(id=id)
        for j, jd in enumerate(d.keys()):
            for id in GasStreamIn.specie.keys():
                i = GasStreamIn.specie[id]["Index"]
                if id in d[jd]["Stoch Gas"].keys():
                    nu = d[jd]["Stoch Gas"][id]
                    matrix["RG"][i, j] = nu * GasStreamIn.get_specie_molar_mass_kg_kmol(id=id)
        return matrix

    def __load_stochiometry_vap_reaction__(self, GasStreamIn, LiquidStreamIn):
        pass

    def __load_stochiometry_vap_all__(self, GasStreamIn, LiquidStreamIn):
        pass

    # ---

    def __load_stochiometry_all__(self, GasStreamIn, LiquidStreamIn):
        pass

    def __load_stochiometry_equilibrium__(self, GasStreamIn, LiquidStreamIn):

        matrix = {}
        matrix_liq = self.__load_stochiometry_liq_equilibrium__(LiquidStreamIn)
        matrix_vap = self.__load_stochiometry_vap_equilibrium__(GasStreamIn, LiquidStreamIn)
        matrix_gas = self.__load_stochiometry_gas_equilibrium__(GasStreamIn)

        num_of_equilibriums = len(LiquidStreamIn.equilibrium_liq.keys()) + len(LiquidStreamIn.equilibrium_vap.keys()) + len(GasStreamIn.equilibrium_gas.keys())

        O1 = np.zeros(shape=(matrix_liq["R"].shape[0], matrix_gas["R"].shape[1]))
        O2 = np.zeros(shape=(matrix_gas["R"].shape[0], matrix_liq["R"].shape[1]))
        Z1 = np.concatenate((matrix_liq["R"], matrix_vap["RL"], O1), axis=1)
        Z2 = np.concatenate((O2, matrix_vap["RG"], matrix_gas["R"]), axis=1)
        matrix["R"] = np.concatenate((Z1, Z2), axis=0)

        matrix["R+"] = np.linalg.pinv(matrix["R"]).astype(np.float64)
        self.Z_project = np.einsum("cq,dq->cd", np.einsum("cr,rq->cq", matrix["R"], np.linalg.pinv(np.einsum("cr,cq->rq", matrix["R"], matrix["R"]))), matrix["R"]).astype(np.float64)
        U, _, _ = np.linalg.svd(matrix["R"], full_matrices=True)
        matrix["A"] = U[:, num_of_equilibriums::].T.astype(np.float64)

        return matrix

    # -------------------------------------------------------------------------------------------

    def __get_equilibrium_liq_K__(self, LiquidStreamIn):

        d = LiquidStreamIn.equilibrium_liq
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = LiquidStreamIn.num_of_equilibrium_liq
        K = np.zeros(shape=(num_of_samples, num_of_reactions), dtype=np.float64)
        for i, id in enumerate(d.keys()):
            K[:, i] = d[id]["K"](LiquidStreamIn)

        den = 0
        for id in LiquidStreamIn.specie.keys():
            den = den + LiquidStreamIn.get_specie_mass_fraction(id) / LiquidStreamIn.get_specie_molar_mass_kg_kmol(id)
        c = LiquidStreamIn.get_solution_molarity_kmol_m3()
        rho = LiquidStreamIn.get_solution_density_kg_m3()

        for i, id in enumerate(d.keys()):
            for r in d[id]["Unit"].keys():
                power = d[id]["Stoch"][r]
                K[:, i] = K[:, i] * LiquidStreamIn.get_specie_activity_coefficient(id=r) ** (-power)
                if d[id]["Unit"][r] == "c":
                    K[:, i] = K[:, i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r) / c) ** power
                elif d[id]["Unit"][r] == "x":
                    K[:, i] = K[:, i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r)) ** power
                elif d[id]["Unit"][r] == "m":
                    K[:, i] = K[:, i] * ((rho * LiquidStreamIn.get_specie_mass_fraction(id=LiquidStreamIn.solvent_id) * LiquidStreamIn.get_specie_molar_mass_kg_kmol(p) * den) / (1000 * c)) ** power
                elif d[id]["Unit"][r] == None:
                    K[:, i] = K[:, i] * (1 / LiquidStreamIn.get_specie_mass_fraction(id=p)) ** power
                elif d[id]["Unit"][r] == "w":
                    pass
        return K

    def __get_equilibrium_liq_dKdT_q__(self, LiquidStreamIn, K):
        dT = 0.05
        K0 = K
        T0 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + 0.1
        K1 = self.__get_equilibrium_liq_K__(LiquidStreamIn)
        T1 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - 0.1
        dKdT = (K1 - K0) / dT
        exothermic_heat_kJ_kmol = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1[:,None]) - (1 / T0[:,None]))
        return dKdT, exothermic_heat_kJ_kmol

    def __get_equilibrium_liq_l__(self, LiquidStreamIn, K):
        d = LiquidStreamIn.equilibrium_liq
        l = np.log(K)
        for j, jd in enumerate(d.keys()):
            for id in d[jd]["Stoch"].keys():
                nu = - d[jd]["Stoch"][id]
                w = LiquidStreamIn.get_specie_mass_fraction(id=id)
                l[:, j] = l[:, j] + nu * np.log(w)
        return l

    def __get_equilibrium_liq_l2__(self, LiquidStreamIn, K):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_equilibriums = LiquidStreamIn.num_of_equilibrium_liq
        d = LiquidStreamIn.equilibrium_liq
        fr = 1.0 * K
        fp = np.ones(shape=(num_of_samples, num_of_equilibriums), dtype=np.float64)
        for j, jd in enumerate(d.keys()):
            for id in d[jd]["Stoch"].keys():
                nu = d[jd]["Stoch"][id]
                if nu > 0:
                    fp[:, j] = fp[:, j] * LiquidStreamIn.get_specie_mass_fraction(id=id) ** np.abs(nu)
                else:
                    fr[:, j] = fr[:, j] * LiquidStreamIn.get_specie_mass_fraction(id=id) ** np.abs(nu)
        f = fr - fp
        return f

    def __get_equilibrium_liq_dldw__(self, LiquidStreamIn):
        d = LiquidStreamIn.equilibrium_liq
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = LiquidStreamIn.num_of_equilibrium_liq
        dldw = np.zeros(shape=(num_of_samples, num_of_reactions, LiquidStreamIn.num_of_species), dtype=np.float64)
        for i, id in enumerate(d.keys()):
            for s in LiquidStreamIn.specie.keys():
                si = LiquidStreamIn.specie[s]["Index"]
                if s in d[id]["Stoch"].keys():
                    nu = - d[id]["Stoch"][s]
                    w = LiquidStreamIn.get_specie_mass_fraction(id=s)
                    dldw[:, i, si] = nu / w
        return dldw

    def __get_equilibrium_liq_dldw2__(self, LiquidStreamIn, K):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_equilibriums = LiquidStreamIn.num_of_equilibrium_liq
        num_of_species = LiquidStreamIn.num_of_species
        d = LiquidStreamIn.equilibrium_liq
        dldw = np.zeros(shape=(num_of_samples, num_of_equilibriums, num_of_species), dtype=np.float64)
        for equilibrium_i, equilibrium_id in enumerate(d.keys()):
            for specie_i, specie_id in enumerate(LiquidStreamIn.specie.keys()):
                if specie_id in d[equilibrium_id]["Stoch"].keys():
                    nu = d[equilibrium_id]["Stoch"][specie_id]
                    w = LiquidStreamIn.get_specie_mass_fraction(id=specie_id).astype(np.float64)
                    if nu < 0:
                        dldw[:, equilibrium_i, specie_i] = K[:, equilibrium_i] * np.abs(nu) * w ** (np.abs(nu) - 1)
                    else:
                        dldw[:, equilibrium_i, specie_i] = - np.abs(nu) * w ** (np.abs(nu) - 1)
                    for specie_jd in d[equilibrium_id]["Stoch"].keys():
                        if specie_id != specie_jd:
                            nu = d[equilibrium_id]["Stoch"][specie_jd]
                            w = LiquidStreamIn.get_specie_mass_fraction(id=specie_jd).astype(np.float64)
                            dldw[:, equilibrium_i, specie_i] = dldw[:, equilibrium_i, specie_i] * w ** np.abs(nu)
        return dldw

    def __get_equilibrium_liq_dldw2td__(self, LiquidStreamIn, K):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_equilibriums = LiquidStreamIn.num_of_equilibrium_liq
        num_of_species = LiquidStreamIn.num_of_species
        l = self.__get_equilibrium_liq_l2__(LiquidStreamIn, K)
        dldw = np.zeros(shape=(num_of_samples, num_of_equilibriums, num_of_species), dtype=np.float64)
        for specie_i, specie_id in enumerate(LiquidStreamIn.specie.keys()):
            w = LiquidStreamIn.get_specie_mass_fraction(id=specie_id)
            dw = np.maximum(0.01 * w, self.dw_min_td)
            LiquidStreamIn.set_specie_mass_fraction(id=specie_id, value=w + dw)
            l_new = self.__get_equilibrium_liq_l2__(LiquidStreamIn, K)
            LiquidStreamIn.set_specie_mass_fraction(id=specie_id, value=w)
            dldw[:, :, specie_i] = (l_new - l) / dw[:,None]
        return dldw

    def __get_equilibrium_liq_dldT__(self, LiquidStreamIn, K, dKdT):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = LiquidStreamIn.num_of_equilibrium_liq
        dldT = np.zeros(shape=(num_of_samples, num_of_reactions, 1), dtype=np.float64)
        for i in range(num_of_reactions):
            dldT[:, i, 0] = (1 / K[:, i]) * dKdT[:, i]
        return dldT

    def __get_equilibrium_liq_L__(self, LiquidStreamIn, K, b):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        w = np.zeros(shape=(num_of_samples, LiquidStreamIn.eq["Number of Species"]), dtype=np.float64)
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            w[:, i] = LiquidStreamIn.get_specie_mass_fraction(id=id)
        # b = np.einsum("ic,sc->si", self.R_conservation, w).astype(np.float64)
        l1 = self.__get_equilibrium_liq_l__(LiquidStreamIn, K)
        l2 = np.einsum("ic,sc->si", self.R_conservation, w) - b
        L = np.hstack((l1, l2))
        return L

    def __get_equilibrium_liq_dLdw__(self, LiquidStreamIn):
        dldw = self.__get_equilibrium_liq_dldw__(LiquidStreamIn)
        A = 1.0 * self.R_conservation
        A = np.reshape(A, newshape=(1, A.shape[0], A.shape[1]))
        A = np.broadcast_to(A, shape=(dldw.shape[0], A.shape[1], A.shape[2])).copy()
        dLdw = np.concatenate((dldw, A), axis=1)
        return dLdw

    # -------------------------------------------------------------------------------------------

    def __get_reaction_liq_q__(self, LiquidStreamIn):

        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = LiquidStreamIn.num_of_reaction_liq
        q = np.zeros(shape=(num_of_samples, num_of_reactions), dtype=np.float64)

        for i, id in enumerate(LiquidStreamIn.reaction_liq.keys()):
            r_forward, r_reverse = LiquidStreamIn.reaction_liq[id]["Rate [kmol/m3.s]"](LiquidStreamIn)
            K0 = r_forward / r_reverse
            T0 = LiquidStreamIn.get_solution_temp_K()
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + 0.05
            r_forward, r_reverse = LiquidStreamIn.reaction_liq[id]["Rate [kmol/m3.s]"](LiquidStreamIn)
            K1 = r_forward / r_reverse
            T1 = LiquidStreamIn.get_solution_temp_K()
            LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - 0.05
            q[:, i] = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1) - (1 / T0))

        return q

    def __get_reaction_gas_q__(self, LiquidStreamIn):
        pass

    # -------------------------------------------------------------------------------------------

    def __get_equilibrium_gas_K__(self, GasStreamIn):

        d = GasStreamIn.equilibrium_gas
        num_of_samples = GasStreamIn.temp_K.shape[0]
        num_of_reactions = GasStreamIn.num_of_equilibrium_gas
        K = np.zeros(shape=(num_of_samples, num_of_reactions), dtype=np.float64)
        for i, id in enumerate(d.keys()):
            K[:, i] = d[id]["K"](GasStreamIn)

        # FUGACITY MISSING....
        # for i, id in enumerate(d.keys()):
        #    for r in d[id]["Unit"].keys():
        #        K[:, i] = K[:, i] * GasStreamIn.get_specie_activity_coefficient(id=r) ** (-d[id]["Stoch"][r])

        return K

    def __get_equilibrium_gas_dKdT_q__(self, GasStreamIn, K):
        dT = 0.05
        K0 = K
        T0 = GasStreamIn.temp_K
        GasStreamIn.temp_K = GasStreamIn.temp_K + 0.1
        K1 = self.__get_equilibrium_gas_K__(GasStreamIn)
        T1 = GasStreamIn.temp_K
        GasStreamIn.temp_K = GasStreamIn.temp_K - 0.1
        dKdT = (K1 - K0) / dT
        exothermic_heat_kJ_kmol = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1[:, None]) - (1 / T0[:, None]))
        return dKdT, exothermic_heat_kJ_kmol

    def __get_equilibrium_gas_g__(self, GasStreamIn, K):
        d = GasStreamIn.equilibrium_gas
        l = np.log(K)
        for j, jd in enumerate(d.keys()):
            for id in d[jd]["Stoch"].keys():
                nu = - d[jd]["Stoch"][id]
                p = GasStreamIn.get_specie_pressure_bara(id=id)
                l[:, j] = l[:, j] + nu * np.log(p)
        return l

    def __get_equilibrium_gas_dgdy__(self, GasStreamIn):
        d = GasStreamIn.equilibrium_gas
        num_of_samples = GasStreamIn.temp_K.shape[0]
        num_of_reactions = GasStreamIn.num_of_equilibrium_gas
        dldw = np.zeros(shape=(num_of_samples, num_of_reactions, GasStreamIn.num_of_species), dtype=np.float64)
        for i, id in enumerate(d.keys()):
            for s in GasStreamIn.specie.keys():
                si = GasStreamIn.specie[s]["Index"]
                if s in d[id]["Stoch"].keys():
                    nu = - d[id]["Stoch"][s]
                    p = GasStreamIn.get_specie_pressure_bara(id=s)
                    dldw[:, i, si] = nu / p
        return dldw

    def __get_equilibrium_gas_dgdT__(self, GasStreamIn, K, dKdT):
        num_of_samples = GasStreamIn.temp_K.shape[0]
        num_of_reactions = GasStreamIn.num_of_equilibrium_gas
        dldT = np.zeros(shape=(num_of_samples, num_of_reactions, 1), dtype=np.float64)
        for i in range(num_of_reactions):
            dldT[:, i, 0] = (1 / K[:, i]) * dKdT[:, i]
        return dldT

    # -------------------------------------------------------------------------------------------

    def __get_equilibrium_vap_K__(self, GasStreamIn, LiquidStreamIn):

        d = LiquidStreamIn.equilibrium_vap
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = LiquidStreamIn.num_of_equilibrium_vap
        K = np.zeros(shape=(num_of_samples, num_of_reactions), dtype=np.float64)
        for i, id in enumerate(d.keys()):
            K[:, i] = d[id]["K"](LiquidStreamIn)

        den = 0
        for id in LiquidStreamIn.specie.keys():
            den = den + LiquidStreamIn.get_specie_mass_fraction(id) / LiquidStreamIn.get_specie_molar_mass_kg_kmol(id)
        c = LiquidStreamIn.get_solution_molarity_kmol_m3()
        rho = LiquidStreamIn.get_solution_density_kg_m3()

        for i, id in enumerate(d.keys()):
            for r in d[id]["Unit Liq"].keys():
                K[:, i] = K[:, i] * LiquidStreamIn.get_specie_activity_coefficient(id=r) ** (-d[id]["Stoch Liq"][r])
                if d[id]["Unit Liq"][r] == "c":
                    K[:, i] = K[:, i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r) / c) ** d[id]["Stoch Liq"][r]
                elif d[id]["Unit Liq"][r] == "x":
                    K[:, i] = K[:, i] * (den * LiquidStreamIn.get_specie_molar_mass_kg_kmol(r)) ** d[id]["Stoch Liq"][r]
                elif d[id]["Unit Liq"][r] == "m":
                    K[:, i] = K[:, i] * ((rho * LiquidStreamIn.get_specie_mass_fraction(id=LiquidStreamIn.solvent_id) * LiquidStreamIn.get_specie_molar_mass_kg_kmol(p) * den) / (1000 * c)) ** d[id]["Stoch Liq"][r]
                elif d[id]["Unit Liq"][r] == None:
                    K[:, i] = K[:, i] * (1 / LiquidStreamIn.get_specie_mass_fraction(id=p)) ** d[id]["Stoch Liq"][r]
                elif d[id]["Unit Liq"][r] == "w":
                    pass
        return K

    def __get_equilibrium_vap_dKdT_q__(self, GasStreamIn, LiquidStreamIn, K):
        dT = 0.05
        K0 = K
        T0 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + 0.1
        K1 = self.__get_equilibrium_vap_K__(GasStreamIn, LiquidStreamIn)
        T1 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - 0.1
        dKdT = (K1 - K0) / dT
        exothermic_heat_kJ_kmol = 8.314 * (np.log(K1) - np.log(K0)) / ((1 / T1[:, None]) - (1 / T0[:, None]))
        return dKdT, exothermic_heat_kJ_kmol

    def __get_equilibrium_vap_v__(self, GasStreamIn, LiquidStreamIn, K):
        d = LiquidStreamIn.equilibrium_vap
        l = np.log(K)
        for j, jd in enumerate(d.keys()):
            for id in d[jd]["Stoch Liq"].keys():
                nu = - d[jd]["Stoch Liq"][id]
                w = LiquidStreamIn.get_specie_mass_fraction(id=id)
                l[:, j] = l[:, j] + nu * np.log(w)
        for j, jd in enumerate(d.keys()):
            for id in d[jd]["Stoch Gas"].keys():
                nu = - d[jd]["Stoch Gas"][id]
                p = GasStreamIn.get_specie_pressure_bara(id=id)
                l[:, j] = l[:, j] + nu * np.log(p)
        return l

    def __get_equilibrium_vap_dvdw__(self, LiquidStreamIn):
        d = LiquidStreamIn.equilibrium_vap
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = LiquidStreamIn.num_of_equilibrium_vap
        dldw = np.zeros(shape=(num_of_samples, num_of_reactions, LiquidStreamIn.num_of_species), dtype=np.float64)
        for i, id in enumerate(d.keys()):
            for s in LiquidStreamIn.specie.keys():
                si = LiquidStreamIn.specie[s]["Index"]
                if s in d[id]["Stoch Liq"].keys():
                    nu = - d[id]["Stoch Liq"][s]
                    w = LiquidStreamIn.get_specie_mass_fraction(id=s)
                    dldw[:, i, si] = nu / w
        return dldw

    def __get_equilibrium_vap_dvdy__(self, GasStreamIn, LiquidStreamIn):
        d = LiquidStreamIn.equilibrium_vap
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = LiquidStreamIn.num_of_equilibrium_vap
        dldp = np.zeros(shape=(num_of_samples, num_of_reactions, GasStreamIn.num_of_species), dtype=np.float64)
        for i, id in enumerate(d.keys()):
            for s in GasStreamIn.specie.keys():
                si = GasStreamIn.specie[s]["Index"]
                if s in d[id]["Stoch Gas"].keys():
                    nu = - d[id]["Stoch Gas"][s]
                    p = GasStreamIn.get_specie_pressure_bara(id=s)
                    dldp[:, i, si] = nu / p
        return dldp

    def __get_equilibrium_vap_dvdT__(self, GasStreamIn, LiquidStreamIn, K, dKdT):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = LiquidStreamIn.num_of_equilibrium_vap
        dldT = np.zeros(shape=(num_of_samples, num_of_reactions, 1), dtype=np.float64)
        for i in range(num_of_reactions):
            dldT[:, i, 0] = (1 / K[:, i]) * dKdT[:, i]
        return dldT



# ---- STREAMS & GAS-LIQUID CONTACTOR -------------------------------------------------


class GasStream(FeatureBooster):

    def __init__(self):
        super().__init__()

        # Container for Various Info
        self.info = {}

        # Universal Gas Constant (bar.m3 / kmol.K)
        self.R = 8.314 * 10**(-2)

        # Concentration
        self.specie = {}
        self.equilibrium_gas = {}

        self.num_of_species = 0
        self.num_of_equilibrium_gas = 0

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
        for k in library.specie_gas[id].keys():
            self.specie[id][k] = library.specie_gas[id][k]
        self.specie[id]["Molar Fraction"] = None
        self.num_of_species = self.num_of_species + 1

    def add_equilibrium_gas(self, id, library):
        self.equilibrium_gas[id] = {}
        for k in library.equilibrium_gas[id].keys():
            self.equilibrium_gas[id][k] = library.equilibrium_gas[id][k]
        self.equilibrium_gas[id]["Index"] = self.num_of_equilibrium_gas
        self.num_of_equilibrium_gas = self.num_of_equilibrium_gas + 1

    # --------------------------------------------------------------------------------------------

    def load_viscosity(self, id, library):
        self.viscosity_Pas = library.viscosity[id]

    def load_thermal_conductivity_kW_mK(self, id, library):
        self.thermal_conductivity_kW_mK = library.thermal_conductivity[id]

    def load_heat_capacity_kJ_kmolK(self, id, library):
        self.heat_capacity_kJ_kmolK = library.heat_capacity[id]

    def load_diffusivity_m2_s(self, id, library):
        self.diffusivity_m2_s = library.diffusivity[id]

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
        return self.molarity_kmol_m3

    def get_gas_density_kg_m3(self):
        c = self.molarity_kmol_m3
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


class LiquidStream(FeatureBooster):

    def __init__(self, solvent_id):

        super().__init__()

        # Solvent (E.g. H2O) when calulating molality
        self.solvent_id = solvent_id

        # Container for Various Info
        self.info = {}

        # Concentration
        self.specie = {}
        self.equilibrium_liq = {}
        self.equilibrium_vap = {}
        self.reaction_liq = {}

        self.num_of_species = 0
        self.num_of_equilibrium_liq = 0
        self.num_of_equilibrium_vap = 0
        self.num_of_reaction_liq = 0

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

    def add_equilibrium_liq(self, id, library):
        self.equilibrium_liq[id] = {}
        for k in library.LiquidStream_rxn_insta[id].keys():
            self.equilibrium_liq[id][k] = library.LiquidStream_rxn_insta[id][k]
        self.equilibrium_liq[id]["Index"] = self.num_of_equilibrium_liq
        self.num_of_equilibrium_liq = self.num_of_equilibrium_liq + 1

    def add_equilibrium_vap(self, id, library):
        self.equilibrium_vap[id] = {}
        for k in library.equilibrium_vap[id].keys():
            self.equilibrium_vap[id][k] = library.LiquidStream_vapor_pressure_bara[id][k]
        self.equilibrium_vap[id]["Index"] = self.num_of_equilibrium_vap
        self.num_of_equilibrium_vap = self.num_of_equilibrium_vap + 1

    def add_reaction(self, id, library):
        self.reaction_liq[id] = {}
        for k in library.LiquidStream_rxn_reversible[id].keys():
            self.reaction_liq[id][k] = library.LiquidStream_rxn_reversible[id][k]
        self.reaction_liq[id]["Index"] = self.num_of_reaction_liq
        self.num_of_reaction_liq = self.num_of_reaction_liq + 1

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
        self.temp_K = value

    def set_solution_flow_kg_h(self, value):
        self.flow_kg_h = value

    def set_specie_mass_fraction(self, id, value):
        self.specie[id]["Mass Fraction"] = value

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

    # --------------------------------------------------------------------

    def get_specie_molar_mass_kg_kmol(self, id):
        return self.specie[id]["Molar Mass [kg/kmol]"]

    def get_specie_charge(self, id):
        return self.specie[id]["Charge"]

    def get_specie_molarity_kmol_m3(self, id):
        rho = self.get_solution_density_kg_m3()
        return rho * self.specie[id]["Mass Fraction"] / self.specie[id]["Molar Mass [kg/kmol]"]

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

        for id in self.equilibrium_vap.keys():
            if gas_id in self.equilibrium_vap[id]["Stoch Gas"].keys():
                reaction_id = id
        liq_id = list(self.equilibrium_vap[reaction_id]["Stoch Liq"].keys())[0]
        unit = self.equilibrium_vap[reaction_id]["Unit Liq"][liq_id]
        gamma = self.get_specie_activity_coefficient(id=liq_id)

        if self.equilibrium_vap[reaction_id]["H"] == None:
            K = 1 / self.equilibrium_vap[reaction_id]["p0"](self)
        else:
            K = self.equilibrium_vap[reaction_id]["H"](self)

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


def LiquidStream_Sum(streams):
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

    LiquidStreamOut.normalize()
    return LiquidStreamOut


def GasStream_Sum(streams):
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


class Equilibrium_Liquid_Isothermal(__Reactor_Common__):

    def __init__(self):
        super().__init__()
        self._first_training_iteration_ = True

    def approximate(self, LiquidStreamIn):

        LiquidStreamOut = deepcopy(LiquidStreamIn)

        # Approx Mass Fractions at Equilibrium
        T = LiquidStreamIn.get_solution_temp_K()
        #b = LiquidStreamIn.__eq_get_conserved_quantities__()

        w = LiquidStreamOut.__mass_fractions_dic2vec__()
        b = np.einsum("ic,sc->si", self.A, w).astype(np.float64)
        w_new, sample_accurate, sample_valid = self.__approximate__(b, T)

        # Return Result
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            if id in self.reactive_species:
                LiquidStreamOut.set_specie_mass_fraction(id=id,
                                                         value=np.array(np.maximum(w_new[:, i], 0.0), dtype=np.float32))
            else:
                LiquidStreamOut.set_specie_mass_fraction(id=id, value=LiquidStreamIn.get_specie_mass_fraction(id=id))
        LiquidStreamOut.normalize_mass_fractions()

        return LiquidStreamOut



    def approximate_react(self, LiquidStreamIn):

        if self._first_training_iteration_:
            LiquidStreamOut, _ = self.__react__(LiquidStreamIn)

        else:

            LiquidStreamOut = deepcopy(LiquidStreamIn)

            # Approx Mass Fractions at Equilibrium
            T = LiquidStreamIn.get_solution_temp_K()
            w = LiquidStreamOut.__mass_fractions_dic2vec__()
            b = np.einsum("ic,sc->si", self.A, w).astype(np.float64)

            #b = LiquidStreamIn.__eq_get_conserved_quantities__()
            w_new, sample_accurate, sample_valid = self.__approximate__(b, T)

            # Only Valid Samples...
            w_in = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_species))
            for id in LiquidStreamIn.specie.keys():
                i = LiquidStreamIn.specie[id]["Index"]
                w_in[:, i] = LiquidStreamIn.get_specie_mass_fraction(id=id)
            w_new = w_new * sample_valid[:, None] + w_in * np.invert(sample_valid[:, None])

            # Approximation...
            for id in LiquidStreamIn.specie.keys():
                i = LiquidStreamIn.specie[id]["Index"]
                if id in self.reactive_species:
                    LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.array(np.maximum(w_new[:, i], 0.0),
                                                                                   dtype=np.float64))
                else:
                    LiquidStreamOut.set_specie_mass_fraction(id=id,
                                                             value=LiquidStreamIn.get_specie_mass_fraction(id=id))
            LiquidStreamOut.normalize_mass_fractions()

            # Fine Tune...
            LiquidStreamOut, _ = self.__react__(LiquidStreamOut)

        return LiquidStreamOut

    def approximate_react_train(self, LiquidStreamIn):
        pass

    def react(self, LiquidStreamIn):
        LiquidStreamOut, _ = self.__react__(LiquidStreamIn)
        return LiquidStreamOut

    def train(self, LiquidStreamIn, max_error_pct=25, degrees_of_freedom=None):

        # Calculate Chemical Equilibrium aka Target Value
        LiquidStreamEq = self.react(LiquidStreamIn)

        # All samples in LiquidStreamOut Object is used for training
        # - If no training have been performed it is set equal to LiquidStreamEq
        # - If training have been performed the Target is compared to value from Approximator
        # --- If approximator is accurate sample is not used for training
        # --- If approximator is inaccurate sample is used for training

        if self._first_training_iteration_:
            LiquidStreamOut = deepcopy(LiquidStreamEq)
        else:
            LiquidStreamApprox = self.approximate(LiquidStreamIn)
            w_eq = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_species), dtype=np.float64)
            w_approx = np.zeros(shape=(LiquidStreamIn.temp_K.shape[0], LiquidStreamIn.num_of_species), dtype=np.float64)
            for id in LiquidStreamIn.specie.keys():
                i = LiquidStreamIn.specie[id]["Index"]
                w_eq[:, i] = LiquidStreamEq.get_specie_mass_fraction(id=id)
                w_approx[:, i] = LiquidStreamApprox.get_specie_mass_fraction(id=id)
            sample_need_training = ((w_approx > (1.0 + max_error_pct / 100) * w_eq) & (w_eq > 0)) | ((w_approx < (1 - max_error_pct / 100) * w_eq) & (w_eq > 0))
            sample_need_training = np.max(sample_need_training, axis=1, keepdims=False)
            LiquidStreamOut = self.__compress__(LiquidStreamIn=LiquidStreamEq, condition=sample_need_training)

        # Dimensions
        num_of_samples = LiquidStreamOut.temp_K.shape[0]
        num_of_reactions = len(self._liq_K_)

        if num_of_samples > 0:

            # Mass Fractions
            w = np.zeros(shape=(num_of_samples, LiquidStreamOut.num_of_species), dtype=np.float64)
            for id in LiquidStreamOut.specie.keys():
                i = LiquidStreamOut.specie[id]["Index"]
                w[:, i] = LiquidStreamOut.get_specie_mass_fraction(id=id)

            # Calculate Hessian
            Kw = self.__get_equilibrium_K_liq__(LiquidStreamOut)
            dLdw = self.__dLdw__(LiquidStreamOut, Kw)
            H = np.linalg.inv(dLdw)

            # Sensitivity Matrix with respect to "b"
            dwdb = H[:, :, num_of_reactions::]

            # Sensitivity Matrix with respect to "T"
            LiquidStreamOut2 = deepcopy(LiquidStreamOut)
            LiquidStreamOut2.temp_K = LiquidStreamOut.temp_K + 0.1
            Kw2 = self.__get_equilibrium_K_liq__(LiquidStreamOut2)

            dKwdT = (Kw2 - Kw) / 0.1
            Hr = H[:, :, :num_of_reactions:]
            dwdT = - np.einsum("scr,sr->sc", Hr, (1 / Kw) * dKwdT)[:, :, None]

            # Features (x)
            #b = LiquidStreamOut.__eq_get_conserved_quantities__()
            b = np.einsum("ic,sc->si", self.A, w).astype(np.float64)

            T = LiquidStreamOut.get_solution_temp_K()[:, None]

            LiquidStreamOut.info["dwdb"] = dwdb
            LiquidStreamOut.info["dwdT"] = dwdT
            LiquidStreamOut.info["b"] = b
            LiquidStreamOut.info["T"] = T
            LiquidStreamOut.info["w"] = w

            # Save for later...
            if self._first_training_iteration_:
                self.LiquidStreamTarget = deepcopy(LiquidStreamOut)
            else:
                self.LiquidStreamTarget.temp_K = np.concatenate(
                    (self.LiquidStreamTarget.temp_K, LiquidStreamOut.temp_K), axis=0)
                self.LiquidStreamTarget.flow_kg_h = np.concatenate(
                    (self.LiquidStreamTarget.flow_kg_h, LiquidStreamOut.flow_kg_h), axis=0)
                for id in self.LiquidStreamTarget.specie.keys():
                    self.LiquidStreamTarget.set_specie_mass_fraction(id=id, value=np.concatenate((
                                                                                                 self.LiquidStreamTarget.get_specie_mass_fraction(
                                                                                                     id=id),
                                                                                                 LiquidStreamOut.get_specie_mass_fraction(
                                                                                                     id=id)), axis=0))
                for id in self.LiquidStreamTarget.info.keys():
                    self.LiquidStreamTarget.info[id] = np.concatenate(
                        (self.LiquidStreamTarget.info[id], LiquidStreamOut.info[id]), axis=0)

            # Input (x) to Nearest Neighbor Algorithm
            mu = np.mean(self.LiquidStreamTarget.info["b"], axis=0, keepdims=True)
            U, D, _ = np.linalg.svd((self.LiquidStreamTarget.info["b"] - mu).T, full_matrices=False)
            dim = np.argmax(D / D[0] < 10 ** (-6)) if degrees_of_freedom is None else degrees_of_freedom
            W = U[0:U.shape[1], 0:dim]
            p = np.einsum("so,on->sn", self.LiquidStreamTarget.info["b"] - mu, W)
            dwdp = np.einsum("swb,bp->swp", self.LiquidStreamTarget.info["dwdb"], W)
            G = np.clip(a=np.mean(np.abs(self.LiquidStreamTarget.info["dwdT"])) / np.mean(np.abs(dwdp)),
                        a_min=10 ** (-12), a_max=1.0)
            x = np.concatenate((p, G * self.LiquidStreamTarget.info["T"]), axis=1)

            # Target
            s_size = self.LiquidStreamTarget.info["w"].shape[0]
            w_size = self.LiquidStreamTarget.info["w"].shape[1]
            b_size = self.LiquidStreamTarget.info["dwdb"].shape[2]
            T_size = 1

            y1 = self.LiquidStreamTarget.info["w"]
            y2 = self.LiquidStreamTarget.info["b"]
            y3 = np.reshape(self.LiquidStreamTarget.info["dwdb"], newshape=(s_size, w_size * b_size))
            y4 = np.reshape(self.LiquidStreamTarget.info["dwdT"], newshape=(s_size, w_size))
            y5 = self.LiquidStreamTarget.info["T"]
            y = np.concatenate((y1, y2, y3, y4, y5), axis=1)

            self.nn_mu = mu
            self.nn_W = W
            self.nn_G = G
            self.nn = KNeighborsRegressor(n_neighbors=1, metric="euclidean", algorithm="ball_tree").fit(x, y)

            # Wrap all up...
            self._first_training_iteration_ = False


    def __approximate__(self, b_new, T_new):

        # Database
        w = self.LiquidStreamTarget.info["w"]
        b = self.LiquidStreamTarget.info["b"]
        T = self.LiquidStreamTarget.info["T"]
        dwdb = self.LiquidStreamTarget.info["dwdb"]
        dwdT = self.LiquidStreamTarget.info["dwdT"]

        # New Samples
        b_new = 1.0 * b_new
        T_new = 1.0 * T_new[:, None]
        # p_new = np.einsum("so,on->sn", b_new - self.database_mu, self.database_W)
        # x_new = np.concatenate((p_new, self.database_G * T_new), axis=1)

        # Extract Dimensions
        num_of_samples_database = b.shape[0]
        num_of_samples_new = b_new.shape[0]
        num_of_species = w.shape[1]
        num_of_conservations = b.shape[1]
        num_of_reactions = num_of_species - num_of_conservations

        p_new = np.einsum("so,on->sn", b_new - self.nn_mu, self.nn_W)
        x_new = np.concatenate((p_new, self.nn_G * T_new), axis=1)
        y_new = self.nn.predict(x_new)

        s0 = 0
        s1 = s0 + num_of_species
        s2 = s1 + num_of_conservations
        s3 = s2 + num_of_species * num_of_conservations
        s4 = s3 + num_of_species
        s5 = s4 + 1

        # Unpack
        w = y_new[:, 0: s1: 1][:, None, :]
        db = (b_new - y_new[:, s1: s2: 1])[:, None, :]
        dwdb = np.reshape(a=y_new[:, s2: s3: 1], newshape=(num_of_samples_new, 1, num_of_species, num_of_conservations))
        dwdT = np.reshape(a=y_new[:, s3: s4: 1], newshape=(num_of_samples_new, 1, num_of_species, 1))
        dT = (T_new - y_new[:, s4: s5: 1])[:, None, :]

        # Extrapolate using Taylor First Order Expansion
        dw = np.einsum("smci,smi->smc", dwdb, db) + np.einsum("smci,smi->smc", dwdT, dT)

        # Update Mass Fractions...
        w_new = np.mean(w + dw, axis=1)
        w = np.mean(w, axis=1)
        dw = np.mean(dw, axis=1)

        # Check if Approximation are Valid (Not Extrapolated to Negative Mass Fractions)
        with np.errstate(divide='ignore', invalid='ignore'):
            tau = np.nan_to_num(- w / dw, nan=1.0, posinf=1.0, neginf=0.0)
        tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
        tau = tau * (w > 0) + 1.0 * (w <= 0)
        tau = np.min(tau, axis=1, keepdims=False)
        sample_valid = (tau >= 1)

        # Check if Approximation is "Accurate"
        epsilon = 4.0
        not_increased_to_much = (w_new < epsilon * w)
        not_decreased_to_much = (w_new > w / epsilon)
        dont_care = (w == 0)
        sample_accurate = (not_decreased_to_much & not_increased_to_much) | dont_care
        sample_accurate = np.min(sample_accurate, axis=1)

        # Return Result and some Info...
        return w_new, sample_accurate, sample_valid




    def __react__(self, LiquidStreamIn):

        matrix = self.__load_stochiometry_liq_equilibrium__(LiquidStreamIn)

        # Dimension
        self.num_of_samples = len(LiquidStreamIn.get_solution_temp_K())

        # Initiate Outlet Stream (Which will be in Equilibrium)
        LiquidStreamOut = deepcopy(LiquidStreamIn)
        for id in LiquidStreamOut.specie.keys():
            i = LiquidStreamOut.specie[id]["Index"]
            #if reactive_species_flag[i] == 1.0:
            LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.maximum(LiquidStreamOut.get_specie_mass_fraction(id=id), 10 ** (-18)))
        LiquidStreamOut.normalize_mass_fractions()

        # Starting Point: Initial Concentrations
        w0 = LiquidStreamOut.__mass_fractions_dic2vec__()
        w = 1.0 * w0

        r = np.zeros(shape=(self.num_of_samples, LiquidStreamOut.num_of_equilibrium_liq))

        Kw = self.__get_equilibrium_liq_K__(LiquidStreamOut)
        #l = self.__get_equilibrium_liq_l_option_2__(LiquidStreamOut, Kw)
        #dldw = self.__get_equilibrium_liq_dldw_option_2_td__(LiquidStreamOut, Kw)
        l = self.__get_equilibrium_liq_l__(LiquidStreamOut, Kw)
        dldw = self.__get_equilibrium_liq_dldw__(LiquidStreamOut)
        dldr = np.einsum("src,cq->srq", dldw, matrix["R"])

        # Preparations...
        converged = False
        epoch = 0
        error = []
        lr = 0.9
        iterations = np.zeros(shape=(self.num_of_samples,))

        # Iterate Until Convergence
        while converged == False:

            # Newton's Method
            dr_newton = - np.linalg.solve(dldr, l)
            dw_newton = np.einsum("ij,sj->si", matrix["R"], dr_newton)

            # Reducing Step Size (Slightly)
            dr = lr * dr_newton
            dw = np.einsum("ij,sj->si", matrix["R"], dr)

            # Backtrack to Ensure only Positive Concentrations
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)

            # Update Mass Fractions...
            r = r + tau * dr
            w = w0 + np.einsum("ij,sj->si", matrix["R"], r)

            # Project Onto Null Space
            #w = np.einsum("cd,sd->sc", self.R_project, w - w0) + w0

            # Update Outlet Stream
            LiquidStreamOut.__mass_fractions_vec2dic__(np.maximum(w, 10**(-18)))

            # Re-evaluate objective functions
            Kw = self.__get_equilibrium_liq_K__(LiquidStreamOut)
            #l = self.__get_equilibrium_liq_l_option_2__(LiquidStreamOut, Kw)
            #dldw = self.__get_equilibrium_liq_dldw_option_2_td__(LiquidStreamOut, Kw)
            l = self.__get_equilibrium_liq_l__(LiquidStreamOut, Kw)
            dldw = self.__get_equilibrium_liq_dldw__(LiquidStreamOut)
            dldr = np.einsum("src,cq->srq", dldw, matrix["R"])

            # Check if Algorithm have Converged
            specie_converged_1 = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)  # 0.005
            specie_converged_2 = np.array(np.abs(dw_newton) < 10 ** (-16), dtype=np.float32)  # 10 ** (-16)
            specie_converged = np.minimum(specie_converged_1 + specie_converged_2, 1.0)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        print("Epochs: \t" + str(epoch))
        return LiquidStreamOut, iterations

    def __compress__(self, LiquidStreamIn, condition):
        LiquidStreamOut = deepcopy(LiquidStreamIn)
        LiquidStreamOut.temp_K = np.compress(condition=condition, a=LiquidStreamIn.temp_K, axis=0)
        LiquidStreamOut.flow_kg_h = np.compress(condition=condition, a=LiquidStreamIn.flow_kg_h, axis=0)
        for id in LiquidStreamOut.specie.keys():
            LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.compress(condition=condition,
                                                                              a=LiquidStreamIn.get_specie_mass_fraction(
                                                                                  id=id), axis=0))
        for id in LiquidStreamOut.info.keys():
            LiquidStreamOut.info[id] = np.compress(condition=condition, a=LiquidStreamIn.info[id], axis=0)
        return LiquidStreamOut


class Equilibrium_Liquid_Adiabatic(__Reactor_Common__):

    def __init__(self):
        super().__init__()

    def react(self, LiquidStreamIn):
        LiquidStreamOut = self.__react__(LiquidStreamIn)
        return LiquidStreamOut

    def __react__(self, LiquidStreamIn):

        self.matrix = self.__load_stochiometry_liq_equilibrium__(LiquidStreamIn)

        LiquidStreamOut = deepcopy(LiquidStreamIn)
        for id in LiquidStreamOut.specie.keys():
            i = LiquidStreamOut.specie[id]["Index"]
            LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.maximum(LiquidStreamOut.get_specie_mass_fraction(id=id), 10 ** (-18)))
        LiquidStreamOut.normalize_mass_fractions()
        LiquidStreamInit = deepcopy(LiquidStreamOut)

        # Dimension
        num_of_samples = len(LiquidStreamIn.get_solution_temp_K())

        # Starting Point: Initial Concentrations
        w0 = LiquidStreamOut.__mass_fractions_dic2vec__()
        w = 1.0 * w0

        K = self.__get_equilibrium_liq_K__(LiquidStreamOut)
        dKdT, dH = self.__get_equilibrium_liq_dKdT_q__(LiquidStreamIn, K)

        #         | l |
        # f =     |   |
        #         | e |

        #         | dl/dr  dl/dT |
        # df/dX = |              |
        #         | de/dr  de/dT |

        l = self.__get_equilibrium_liq_l__(LiquidStreamOut, K)
        e = self.__get_equilibrium_liq_e__(LiquidStreamOut, LiquidStreamInit, dH)
        f = np.concatenate((l, e), axis=1)

        dldw = self.__get_equilibrium_liq_dldw__(LiquidStreamOut)
        dldr = np.einsum("src,cq->srq", dldw, self.matrix["R"])
        dedr = self.__get_equilibrium_liq_dedr__(LiquidStreamOut, dH)
        dfdr = np.concatenate((dldr, dedr), axis=1)

        dldT = self.__get_equilibrium_liq_dldT__(LiquidStreamOut, K, dKdT)
        dedT = - np.ones(shape=(num_of_samples, 1,1))
        dfdT = np.concatenate((dldT, dedT), axis=1)

        dfdX = np.concatenate((dfdr, dfdT), axis=2)

        # Preparations...
        converged = False
        epoch = 0
        error = []
        lr = 0.9
        iterations = np.zeros(shape=(num_of_samples,))

        # Iterate Until Convergence
        while converged == False:

            dX_newton = - np.linalg.solve(dfdX, f)
            dT_newton = dX_newton[:, dX_newton.shape[1] - 1]
            dr_newton = dX_newton[:, :dX_newton.shape[1] - 1:]
            dw_newton = np.einsum("ij,sj->si", self.matrix["R"], dr_newton)

            # Reducing Step Size (Slightly)
            dr = lr * dr_newton
            dT = lr * dT_newton
            dw = np.einsum("ij,sj->si", self.matrix["R"], dr)

            # Backtrack to Ensure only Positive Concentrations
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            tau = np.broadcast_to(tau, shape=(num_of_samples, LiquidStreamIn.num_of_species))

            # Update Mass Fractions...
            w = w + tau * dw

            # Project Onto Null Space
            #w = np.einsum("cd,sd->sc", LiquidStreamIn.R_project, w - w0) + w0
            #for id in LiquidStreamOut.specie.keys():
            #    i = LiquidStreamOut.specie[id]["Index"]
            #    if LiquidStreamOut.reactive_species_flag[i] == 1.0:
            #        w[:, i] = np.maximum(w[:, i], 10 ** (-18))

            # Update Outlet Stream
            LiquidStreamOut.__mass_fractions_vec2dic__(w)
            LiquidStreamOut.temp_K = LiquidStreamOut.temp_K + tau[:, 0] * dT

            # Re-evaluate objective functions
            l = self.__get_equilibrium_liq_l__(LiquidStreamOut, K)
            e = self.__get_equilibrium_liq_e__(LiquidStreamOut, LiquidStreamInit, dH)
            f = np.concatenate((l, e), axis=1)

            dldw = self.__get_equilibrium_liq_dldw__(LiquidStreamOut)
            dldr = np.einsum("src,cq->srq", dldw, self.matrix["R"])
            dedr = self.__get_equilibrium_liq_dedr__(LiquidStreamOut, dH)
            dfdr = np.concatenate((dldr, dedr), axis=1)

            dldT = self.__get_equilibrium_liq_dldT__(LiquidStreamOut, K, dKdT)
            dedT = - np.ones(shape=(num_of_samples, 1, 1))
            dfdT = np.concatenate((dldT, dedT), axis=1)

            dfdX = np.concatenate((dfdr, dfdT), axis=2)

            # Check if Algorithm have Converged
            specie_converged_1 = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)  # 0.005
            specie_converged_2 = np.array(np.abs(dw_newton) < 10 ** (-16), dtype=np.float32)  # 10 ** (-16)
            specie_converged = np.minimum(specie_converged_1 + specie_converged_2, 1.0)
            sample_converged = np.min(specie_converged, axis=1)
            sample_converged = sample_converged * (np.abs(dT_newton) < 10 ** (-2))
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) and (epoch > 0)) or (epoch > 198)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        print("Epochs: \t" + str(epoch))

        return LiquidStreamOut

    def __get_equilibrium_liq_e__(self, LiquidStreamIn, LiquidStreamInit, dH):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        cp = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        dw = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_species))
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            dw[:, i] = LiquidStreamIn.get_specie_mass_fraction(id=id) - LiquidStreamInit.get_specie_mass_fraction(id=id)
        dT = LiquidStreamIn.temp_K - LiquidStreamInit.temp_K
        e = (1 / cp) * np.einsum("sr,sr->s", dH, np.einsum("rw,sw->sr", self.matrix["R+"], dw)) - dT
        return e[:, None]

    def __get_equilibrium_liq_dedr__(self, LiquidStreamIn, dH):
        cp = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        dedr = dH / cp[:, None]
        return dedr[:, None, :]


class Equilibrium_Liquid_QPFlash(FeatureBooster):

    def __init__(self):
        pass

    def react(self, LiquidStreamIn, heat_kW, pressure_bara):

        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_volatile = len(LiquidStreamIn.vle_gas_reactant)

        # Make New Instance
        LiquidStreamOut = deepcopy(LiquidStreamIn)

        # Heat Up Solution
        cp = LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()
        m = LiquidStreamOut.get_solution_flow_kg_h() / 3600
        LiquidStreamOut.temp_K = LiquidStreamOut.temp_K + heat_kW / (cp * m)

        # Calculate Equilibrium (Assuming no Flashing)
        LiquidStreamOut = Equilibrium_Liquid_Adiabatic().react(LiquidStreamOut)

        # Calculate Vapor Pressure from Solution
        p_vap = 0
        p_vap_specie = {}
        y_vap_specie = {}
        for i in range(num_of_volatile):
            id = LiquidStreamOut.vle_gas_reactant[i]
            p_vap_specie[id] = LiquidStreamOut.get_specie_vapor_pressure_bara(gas_id=id)
            p_vap = p_vap + p_vap_specie[id]
        for i in range(num_of_volatile):
            id = LiquidStreamOut.vle_gas_reactant[i]
            y_vap_specie[id] = p_vap_specie[id] / p_vap

        # GasStream Object - Set Molar fractions based on Vapor Pressure
        for i in range(num_of_volatile):
            id = LiquidStreamOut.vle_gas_reactant[i]
            LiquidStreamOut.Vapor.set_specie_molar_fraction(id=id, value=y_vap_specie[id])

        # GasStream Object - Set Temperature, Molarity and Flow. (Flow is set to zero for now...)
        LiquidStreamOut.Vapor.set_gas_temp_K(value=LiquidStreamOut.temp_K)
        LiquidStreamOut.Vapor.set_gas_molarity_kmol_m3(value=pressure_bara / (0.08314 * LiquidStreamOut.Vapor.temp_K))
        LiquidStreamOut.Vapor.set_gas_flow_kmol_h(value=np.zeros(shape=(num_of_samples,)))

        # Extract samples where flash will occur
        flash = (p_vap > 1.00001 * pressure_bara)
        LiquidStreamFlash = self.__compress_liquid__(LiquidStreamOut, flash)
        LiquidStreamFlash.Vapor = self.__compress_gas__(LiquidStreamOut.Vapor, flash)

        # Calculate Flash on selected samples
        LiquidStreamFlash = self.__react__(LiquidStreamFlash, 0 * np.compress(condition=flash, a=heat_kW, axis=0), np.compress(condition=flash, a=pressure_bara, axis=0))

        # Combine...
        LiquidStreamOut.temp_K[flash] = LiquidStreamFlash.temp_K
        LiquidStreamOut.flow_kg_h[flash] = LiquidStreamFlash.flow_kg_h
        for id in LiquidStreamOut.specie.keys():
            LiquidStreamOut.specie[id]["Mass Fraction"][flash] = LiquidStreamFlash.specie[id]["Mass Fraction"]
        for id in LiquidStreamOut.info.keys():
            LiquidStreamOut.info[id][flash] = LiquidStreamFlash.info[id]

        # Combine...
        LiquidStreamOut.Vapor.temp_K[flash] = LiquidStreamFlash.Vapor.temp_K
        LiquidStreamOut.Vapor.flow_kmol_h[flash] = LiquidStreamFlash.Vapor.flow_kmol_h
        for id in LiquidStreamOut.Vapor.specie.keys():
            LiquidStreamOut.Vapor.specie[id]["Molar Fraction"][flash] = LiquidStreamFlash.Vapor.specie[id]["Molar Fraction"]
        for id in LiquidStreamOut.Vapor.info.keys():
            LiquidStreamOut.Vapor.info[id][flash] = LiquidStreamFlash.Vapor.info[id]

        # Return Result...
        return LiquidStreamOut

    def __react__(self, LiquidStreamIn, heat_kW, pressure_bara):

        num_of_volatile = len(LiquidStreamIn.vle_gas_reactant)
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = len(LiquidStreamIn.equilibrium_K)
        num_of_objectives = num_of_reactions + num_of_volatile + 1
        num_of_species = LiquidStreamIn.num_of_species

        # Init...
        LiquidStreamInit = deepcopy(LiquidStreamIn)

        # The Object "LiquidStreamIn" have zero Gas Flow.
        # The flow should be non-zero for sake of numerical stability.
        w = self.__mass_fraction_dictionary_to_array__(LiquidStreamIn)
        y_vap = np.zeros(shape=(num_of_samples, num_of_volatile))
        for i in range(num_of_volatile):
            id = LiquidStreamIn.vle_gas_reactant[i]
            y_vap[:, i] = LiquidStreamIn.Vapor.get_specie_molar_fraction(id=id)
        epsilon = np.einsum("wy,sy->sw", LiquidStreamIn.PL, y_vap)

        with np.errstate(divide='ignore', invalid='ignore'):
            dp_max = np.nan_to_num(w / epsilon, nan=np.inf, posinf=np.inf, neginf=np.inf)
        dp_max = np.min(dp_max, axis=1)
        dp = np.einsum("s,sy->sy", - 0.05 * dp_max, y_vap)
        dr = np.zeros(shape=(num_of_samples, num_of_reactions))
        dz = np.concatenate((dr, dp), axis=1)

        # Initial Flow (Gas Flow is Zero)
        m0 = np.zeros(shape=(num_of_samples, num_of_species + num_of_volatile), dtype=np.float64)
        for i, id in enumerate(LiquidStreamIn.specie.keys()):
            m0[:, i] = LiquidStreamIn.get_specie_flow_kg_h(id=id)
        for i, id in enumerate(LiquidStreamIn.Vapor.specie.keys()):
            m0[:, i + num_of_species] = LiquidStreamIn.Vapor.get_specie_flow_kg_h(id=id)

        # Initial Flow (Gas Flow non-Zero)
        dm = np.einsum("s,sm->sm", LiquidStreamIn.get_solution_flow_kg_h(), np.einsum("iz,mz->im", dz, LiquidStreamIn.Z))
        m = m0 + dm

        # Update Liquid
        m_liq = np.sum(m[:, 0:num_of_species:1], axis=1)
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            w = m[:, i] / m_liq
            LiquidStreamIn.set_specie_mass_fraction(id=id, value=w)

        # GasStream Object - Calculate average Molar Mass
        #M_vapor_avg = 0
        #for i in range(num_of_volatile):
        #    id = LiquidStreamIn.vle_gas_reactant[i]
        #    M_vapor_avg = M_vapor_avg + y_vap_specie[id] * LiquidStreamOut.Vapor.get_specie_molar_mass_kg_kmol(id=id)

        # Update Gas
        #m_gas = np.sum(m[:, num_of_species: num_of_species + num_of_volatile: 1] + dm[:,num_of_species: num_of_species + num_of_volatile: 1], axis=1)
        n_gas = 0
        n_specie = {}
        for id in LiquidStreamIn.Vapor.specie.keys():
            i = LiquidStreamIn.Vapor.specie[id]["Index"]
            n_specie[id] = m[:, i + num_of_species] / LiquidStreamIn.Vapor.get_specie_molar_mass_kg_kmol(id=id)
            n_gas = n_gas + n_specie[id]

        for id in LiquidStreamIn.Vapor.specie.keys():
            i = LiquidStreamIn.Vapor.specie[id]["Index"]
            LiquidStreamIn.Vapor.set_specie_molar_fraction(id=id, value=n_specie[id] / n_gas)

        LiquidStreamIn.Vapor.set_gas_flow_kmol_h(value=n_gas)
        LiquidStreamIn.Vapor.set_gas_temp_K(value=LiquidStreamIn.temp_K)
        LiquidStreamIn.Vapor.set_gas_molarity_kmol_m3(value=pressure_bara / (0.08314 * LiquidStreamIn.Vapor.temp_K))

        # Initiating Liquid Streams
        LiquidStreamOut = deepcopy(LiquidStreamIn)

        # Objective Functions
        lamda = LiquidStreamOut.__heat_of_vaporization_kJ_kmol__(LiquidStreamOut)
        Kw, dH, dKdT = LiquidStreamOut.__Kw_and_heat_of_reaction_kJ_kmol__(LiquidStreamOut)
        Kw_vle = LiquidStreamOut.__Kw_vle__(LiquidStreamOut)
        f = self.__f__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)

        # Partial Derivatives
        dldz = self.__dldz__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        dldT = self.__dldT__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        dvdz = self.__dvdz__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        dvdT = self.__dvdT__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        dedz = self.__dedz__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        dedT = self.__dedT__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)

        dldx = np.concatenate((dldz, dldT), axis=2)
        dvdx = np.concatenate((dvdz, dvdT), axis=2)
        dedx = np.concatenate((dedz, dedT), axis=2)
        dfdx = np.concatenate((dldx, dvdx, dedx), axis=1)

        # Preparations...
        converged = False
        epoch = 0
        error = []
        lr = 0.2
        iterations = np.zeros(shape=(num_of_samples,))

        # Iterate Until Convergence
        while converged == False:

            dx_newton = - np.linalg.solve(dfdx, f)

            #dr_newton = dx_newton[:, 0: num_of_reactions: 1]
            #dp_newton = dx_newton[:, num_of_reactions: num_of_reactions + num_of_volatile : 1]
            dz_newton = dx_newton[:, 0: num_of_reactions + num_of_volatile: 1]
            dT_newton = dx_newton[:, num_of_reactions + num_of_volatile: num_of_reactions + num_of_volatile + 1: 1]
            #dw_newton = np.einsum("ij,sj->si", LiquidStreamIn.RLL, dr_newton) + np.einsum("ij,sj->si", LiquidStreamIn.RVL, dp_newton)
            #dy_newton = np.einsum("ij,sj->si", LiquidStreamIn.RVG, dp_newton)
            dm_newton = np.einsum("i,im->im", LiquidStreamInit.get_solution_flow_kg_h(), np.einsum("iz,mz->im", dz_newton, LiquidStreamOut.Z))

            # Reducing Step Size (Slightly)
            #dr = lr * dr_newton
            #dp = lr * dp_newton
            dT = lr * dT_newton
            dm = lr * dm_newton

            # Backtrack to Ensure only Positive Concentrations
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * m / dm, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (m > 0) + 1.0 * (m <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            tau = np.broadcast_to(tau, shape=(num_of_samples, num_of_species + num_of_volatile))

            dm = tau * dm
            dT = tau * dT

            # Update Mass Fractions in Liquid Phase and Molar Fractions in Gas Phase...
            m = m + dm

            # Project Onto Null Space
            #m = np.einsum("cd,sd->sc", LiquidStreamIn.Z_project, m - m0) + m0

            m_liq = np.sum(m[:, 0:num_of_species:1], axis=1)
            m_gas = np.sum(m[:, num_of_species: num_of_species + num_of_volatile: 1], axis=1)
            m_tot = np.sum(m, axis=1)
            #m[:, 0:num_of_species:1] = np.maximum(m[:, 0:num_of_species:1], 10**(-18) * m_liq[:,None])
            #m[:, num_of_species: num_of_species + num_of_volatile: 1] = np.maximum(m[:, num_of_species: num_of_species + num_of_volatile: 1], 10**(-18) * m_gas[:,None])

            #for id in LiquidStreamOut.specie.keys():
            #    i = LiquidStreamOut.specie[id]["Index"]
            #    if LiquidStreamOut.reactive_species_flag[i] == 1.0:
            #        m[:, i] = np.maximum(m[:, i], 10**(-18) * m_liq)

            #for id in LiquidStreamOut.Vapor.specie.keys():
            #    i = LiquidStreamOut.Vapor.specie[id]["Index"]
            #    m[:, i + num_of_species] = np.maximum(m[:, i + num_of_species], 10 ** (-18) * m_gas)

            # Update Liquid Stream
            for id in LiquidStreamOut.specie.keys():
                i = LiquidStreamOut.specie[id]["Index"]
                w = m[:, i] / m_liq
                LiquidStreamOut.set_specie_mass_fraction(id=id, value=w)
            LiquidStreamOut.normalize_mass_fractions()
            LiquidStreamOut.temp_K = LiquidStreamOut.temp_K + dT[:,0]
            LiquidStreamOut.set_solution_flow_kg_h(value=m_liq)

            # Update Vapor Stream
            n_gas = 0
            for id in LiquidStreamOut.Vapor.specie.keys():
                i = LiquidStreamOut.Vapor.specie[id]["Index"]
                n_gas = n_gas + m[:, i + num_of_species] / LiquidStreamOut.Vapor.get_specie_molar_mass_kg_kmol(id=id)
            for id in LiquidStreamOut.Vapor.specie.keys():
                i = LiquidStreamOut.Vapor.specie[id]["Index"]
                y = m[:, i + num_of_species] / LiquidStreamOut.Vapor.get_specie_molar_mass_kg_kmol(id=id)
                LiquidStreamOut.Vapor.set_specie_molar_fraction(id=id, value=y)
            LiquidStreamOut.Vapor.set_gas_flow_kmol_h(value=n_gas)
            LiquidStreamOut.Vapor.set_gas_temp_K(value=LiquidStreamOut.temp_K)
            LiquidStreamOut.Vapor.set_gas_molarity_kmol_m3(value=pressure_bara / (0.08314 * LiquidStreamOut.Vapor.temp_K))

            # Re-evaluate objective functions
            lamda = LiquidStreamOut.__heat_of_vaporization_kJ_kmol__(LiquidStreamOut)
            Kw, dH, dKdT = LiquidStreamOut.__Kw_and_heat_of_reaction_kJ_kmol__(LiquidStreamOut)
            Kw_vle = LiquidStreamOut.__Kw_vle__(LiquidStreamOut)
            f = self.__f__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)

            # Partial Derivatives
            dldz = self.__dldz__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
            dldT = self.__dldT__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
            dvdz = self.__dvdz__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
            dvdT = self.__dvdT__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
            dedz = self.__dedz__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
            dedT = self.__dedT__(LiquidStreamOut.Vapor, LiquidStreamOut, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)

            dldx = np.concatenate((dldz, dldT), axis=2)
            dvdx = np.concatenate((dvdz, dvdT), axis=2)
            dedx = np.concatenate((dedz, dedT), axis=2)
            dfdx = np.concatenate((dldx, dvdx, dedx), axis=1)

            """""""""
            # Check if Algorithm have Converged
            specie_converged_1 = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)  # 0.005
            specie_converged_2 = np.array(np.abs(dw_newton) < 10 ** (-16), dtype=np.float32)  # 10 ** (-16)
            specie_converged = np.minimum(specie_converged_1 + specie_converged_2, 1.0)
            sample_converged = np.min(specie_converged, axis=1)
            sample_converged = sample_converged * (np.abs(dT_newton) < 10 ** (-2))
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) and (epoch > 0)) or (epoch > 198)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1
            """""""""

            epoch = epoch + 1

            converged = (epoch > 150)


        print("Epochs: \t" + str(epoch))
        return LiquidStreamOut

    def __compress_liquid__(self, LiquidStreamIn, condition):
        LiquidStreamOut = deepcopy(LiquidStreamIn)
        LiquidStreamOut.temp_K = np.compress(condition=condition, a=LiquidStreamIn.temp_K, axis=0)
        LiquidStreamOut.flow_kg_h = np.compress(condition=condition, a=LiquidStreamIn.flow_kg_h, axis=0)
        for id in LiquidStreamOut.specie.keys():
            LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.compress(condition=condition, a=LiquidStreamIn.get_specie_mass_fraction(id=id), axis=0))
        for id in LiquidStreamOut.info.keys():
            LiquidStreamOut.info[id] = np.compress(condition=condition, a=LiquidStreamIn.info[id], axis=0)
        return LiquidStreamOut

    def __compress_gas__(self, GasStreamIn, condition):
        GasStreamOut = deepcopy(GasStreamIn)
        GasStreamOut.temp_K = np.compress(condition=condition, a=GasStreamIn.temp_K, axis=0)
        GasStreamOut.flow_kmol_h = np.compress(condition=condition, a=GasStreamIn.flow_kmol_h, axis=0)
        for id in GasStreamOut.specie.keys():
            GasStreamOut.set_specie_molar_fraction(id=id, value=np.compress(condition=condition, a=GasStreamIn.get_specie_molar_fraction(id=id), axis=0))
        for id in GasStreamOut.info.keys():
            GasStreamOut.info[id] = np.compress(condition=condition, a=GasStreamIn.info[id], axis=0)
        return GasStreamOut

    def __l__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        l = np.log(Kw)
        for i in range(LiquidStreamIn.num_of_reactions):
            for id in LiquidStreamIn.equilibrium_reactants[i].keys():
                nu = LiquidStreamIn.equilibrium_reactants[i][id]
                w = LiquidStreamIn.get_specie_mass_fraction(id=id)
                l[:, i] = l[:, i] + nu * np.log(w)
            for id in LiquidStreamIn.equilibrium_products[i].keys():
                nu = - LiquidStreamIn.equilibrium_products[i][id]
                w = LiquidStreamIn.get_specie_mass_fraction(id=id)
                l[:, i] = l[:, i] + nu * np.log(w)
        return l

    def __v__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_volatile = len(LiquidStreamIn.vle_gas_reactant)
        v = np.zeros(shape=(num_of_samples, num_of_volatile), dtype=np.float64)
        n_tot = GasStreamIn.get_gas_flow_kmol_h()
        for i in range(num_of_volatile):
            gas_id = LiquidStreamIn.vle_gas_reactant[i]
            n = GasStreamIn.get_specie_flow_kmol_h(id=gas_id)
            #v[:, i] = np.log(GasStreamIn.get_specie_molar_fraction(id=gas_id)) - np.log(LiquidStreamIn.get_specie_vapor_pressure_bara(gas_id=gas_id) / pressure_bara)
            #v[:, i] = GasStreamIn.get_specie_molar_fraction(id=gas_id) - LiquidStreamIn.get_specie_vapor_pressure_bara(gas_id=gas_id) / pressure_bara
            #v[:, i] = GasStreamIn.get_specie_molar_fraction(id=gas_id) * pressure_bara - LiquidStreamIn.get_specie_vapor_pressure_bara(gas_id=gas_id)
            #v[:, i] = n * pressure_bara - LiquidStreamIn.get_specie_vapor_pressure_bara(gas_id=gas_id) * n_tot
            v[:, i] = GasStreamIn.get_specie_molar_fraction(id=gas_id) - LiquidStreamIn.get_specie_vapor_pressure_bara(gas_id=gas_id) / pressure_bara
            #v[:, i] = np.exp(-GasStreamIn.get_specie_molar_fraction(id=gas_id)) - np.exp(-LiquidStreamIn.get_specie_vapor_pressure_bara(gas_id=gas_id) / pressure_bara)
        return v

    def __e__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_volatile = len(LiquidStreamIn.vle_gas_reactant)
        e = np.zeros(shape=(num_of_samples, 1))
        cp = (LiquidStreamIn.get_solution_heat_capacity_kJ_kgK() + LiquidStreamInit.get_solution_heat_capacity_kJ_kgK()) / 2
        m_init = LiquidStreamInit.get_solution_flow_kg_h() / 3600
        T = LiquidStreamIn.get_solution_temp_K()
        T_init = LiquidStreamInit.get_solution_temp_K()
        Q_sens = m_init * cp * (T - T_init)
        dm = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_species + num_of_volatile))

        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            dm[:, i] = (LiquidStreamIn.get_specie_flow_kg_h(id=id) - LiquidStreamInit.get_specie_flow_kg_h(
                id=id)) / 3600
        for id in GasStreamIn.specie.keys():
            i = GasStreamIn.specie[id]["Index"]
            dm[:, i + LiquidStreamIn.num_of_species] = GasStreamIn.get_specie_flow_kg_h(id=id) / 3600

        dH = np.concatenate((dH, lamda), axis=1)
        Q_exo = np.einsum("sr,sr->s", dH, np.einsum("rw,sw->sr", LiquidStreamIn.Z_inverse, dm))
        e[:, 0] = (heat_kW + Q_exo - Q_sens) / m_init
        return e

    def __f__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        l = self.__l__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        v = self.__v__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        e = self.__e__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        f = np.concatenate((l, v, e), axis=1)
        return f

    def __dldw__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        dldw = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_reactions, LiquidStreamIn.num_of_species), dtype=np.float64)
        for i in range(LiquidStreamIn.num_of_reactions):
            for id in LiquidStreamIn.specie.keys():
                j = LiquidStreamIn.specie[id]["Index"]
                if id in LiquidStreamIn.equilibrium_reactants[i].keys():
                    nu = LiquidStreamIn.equilibrium_reactants[i][id]
                    w = LiquidStreamIn.get_specie_mass_fraction(id=id)
                    dldw[:, i, j] = nu / w
                elif id in LiquidStreamIn.equilibrium_products[i].keys():
                    nu = - LiquidStreamIn.equilibrium_products[i][id]
                    w = LiquidStreamIn.get_specie_mass_fraction(id=id)
                    dldw[:, i, j] = nu / w
        return dldw

    def __dldz__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):

        # dl/dr = dl/dw * dw/dr
        # dl/dp = dl/dw * dw/dp

        num_of_samples = len(LiquidStreamIn.temp_K)
        dldw = self.__dldw__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        m_in = LiquidStreamInit.get_solution_flow_kg_h()
        m_liq = LiquidStreamIn.get_solution_flow_kg_h()
        dwdr = np.einsum("s,ab->sab", m_in / m_liq, LiquidStreamIn.R)
        dldr = np.einsum("slw,swr->slr", dldw, dwdr)

        w = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_species), dtype=np.float64)
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            w[:, i] = LiquidStreamIn.get_specie_mass_fraction(id=id)
        p = np.sum(LiquidStreamIn.PL, axis=0)
        wp = np.einsum("sw,x->swx", w, p)
        dwdp = np.einsum("s,swx->swx", m_in / m_liq, (LiquidStreamIn.PL - wp))
        dldp = np.einsum("slw,swr->slr", dldw, dwdp)

        dldz = np.concatenate((dldr, dldp), axis=2)
        return dldz

    def __dldT__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        dldT = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_reactions), dtype=np.float64)
        for i in range(LiquidStreamIn.num_of_reactions):
            dldT[:, i] = (1 / Kw[:, i]) * dKdT[:, i]
        return dldT[:, :, None]

    def __dvdw_log__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        dvdw = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_volatile, LiquidStreamIn.num_of_species), dtype=np.float64)
        for i in range(LiquidStreamIn.num_of_volatile):
            for id in LiquidStreamIn.specie.keys():
                j = LiquidStreamIn.specie[id]["Index"]
                if id in LiquidStreamIn.vle_liq_products[i].keys():
                    nu = - LiquidStreamIn.vle_liq_products[i][id]
                    dvdw[:, i, j] = nu / LiquidStreamIn.get_specie_mass_fraction(id=id)
        return dvdw

    def __dvdy_log__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        dvdy = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_volatile, LiquidStreamIn.num_of_volatile), dtype=np.float64)
        for i, id in enumerate(GasStreamIn.specie.keys()):
            dvdy[:, i, i] = 1 / GasStreamIn.get_specie_molar_fraction(id=id)
        return dvdy

    def __dvdz_log__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):

        # Line 1:   dv/dr = dv/dy * dy/dr + dv/dw * dw/dr
        # Line 2:   dv/dp = dv/dy * dy/dp + dv/dw * dw/dp

        num_of_samples = len(LiquidStreamIn.temp_K)
        m_in = LiquidStreamInit.get_solution_flow_kg_h()
        m_liq = LiquidStreamIn.get_solution_flow_kg_h()
        n_gas = GasStreamIn.get_gas_flow_kmol_h()

        # Line 1
        dvdy = self.__dvdy__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        dydr = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_volatile, LiquidStreamIn.num_of_reactions))
        dvdw = self.__dvdw__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        dwdr = np.einsum("s,ab->sab", m_in / m_liq, LiquidStreamIn.R)
        dvdr = np.einsum("slw,swr->slr", dvdy, dydr) + np.einsum("slw,swr->slr", dvdw, dwdr)

        # Line 2
        dydp = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_volatile, LiquidStreamIn.num_of_volatile), dtype=np.float64)
        for i, id in enumerate(GasStreamIn.specie.keys()):
            for j, jd in enumerate(GasStreamIn.specie.keys()):
                dydp[:, i, j] = (m_in / n_gas) * (GasStreamIn.get_specie_molar_fraction(id) - 1) if i == j else m_in / n_gas

        w = np.zeros(shape=(num_of_samples, LiquidStreamIn.num_of_species), dtype=np.float64)
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            w[:, i] = LiquidStreamIn.get_specie_mass_fraction(id=id)
        p = np.sum(LiquidStreamIn.PL, axis=0)
        wp = np.einsum("sw,x->swx", w, p)
        dwdp = np.einsum("s,swx->swx", m_in / m_liq, (LiquidStreamIn.PL - wp))

        dvdp = np.einsum("slw,swr->slr", dvdy, dydp) + np.einsum("slw,swr->slr", dvdw, dwdp)
        dvdz = np.concatenate((dvdr, dvdp), axis=2)
        return dvdz * 0

    def __dvdw__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        num_of_volatile = len(LiquidStreamIn.vle_gas_reactant)
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = len(LiquidStreamIn.equilibrium_K)
        num_of_objectives = num_of_reactions + num_of_volatile + 1
        num_of_species = LiquidStreamIn.num_of_species
        w0 = self.__mass_fraction_dictionary_to_array__(LiquidStreamIn)
        dvdw = np.zeros(shape=(num_of_samples, num_of_volatile, num_of_species))
        v0 = self.__v__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        for i in range(num_of_species):
            dw = np.maximum(10**(-5) * w0[:,i], 10**(-15))
            w = 1.0 * w0
            w[:, i] = w[:, i] + dw
            LiquidStreamIn = self.__mass_fraction_array_to_dictionary__(LiquidStreamIn, w)
            lamda = LiquidStreamIn.__heat_of_vaporization_kJ_kmol__(LiquidStreamIn)
            Kw, dH, dKdT = LiquidStreamIn.__Kw_and_heat_of_reaction_kJ_kmol__(LiquidStreamIn)
            Kw_vle = LiquidStreamIn.__Kw_vle__(LiquidStreamIn)
            v = self.__v__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)

            # v.shape = (num_of_samples_, num_of_volatile)
            # dw.shape = (num_of_samples,)
            # dvdw = (num_of_samples, num_of_volatile, num_of_species)

            dvdw[:, :, i] = (v - v0) / dw[:, None]

        return dvdw

    def __dvdy__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        num_of_volatile = len(LiquidStreamIn.vle_gas_reactant)
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = len(LiquidStreamIn.equilibrium_K)
        num_of_objectives = num_of_reactions + num_of_volatile + 1
        num_of_species = LiquidStreamIn.num_of_species
        w0 = self.__mass_fraction_dictionary_to_array__(LiquidStreamIn)
        dvdw = np.zeros(shape=(num_of_samples, num_of_volatile, num_of_species))
        v0 = self.__v__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara,
                        Kw_vle)
        for i in range(num_of_species):
            dw = np.maximum(10 ** (-5) * w0[:, i], 10 ** (-15))
            w = 1.0 * w0
            w[:, i] = w[:, i] + dw
            LiquidStreamIn = self.__mass_fraction_array_to_dictionary__(LiquidStreamIn, w)
            lamda = LiquidStreamIn.__heat_of_vaporization_kJ_kmol__(LiquidStreamIn)
            Kw, dH, dKdT = LiquidStreamIn.__Kw_and_heat_of_reaction_kJ_kmol__(LiquidStreamIn)
            Kw_vle = LiquidStreamIn.__Kw_vle__(LiquidStreamIn)
            v = self.__v__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara,
                           Kw_vle)

            # v.shape = (num_of_samples_, num_of_volatile)
            # dw.shape = (num_of_samples,)
            # dvdw = (num_of_samples, num_of_volatile, num_of_species)

            dvdw[:, :, i] = (v - v0) / dw[:, None]

        return dvdw

    def __dvdz__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        num_of_volatile = len(LiquidStreamIn.vle_gas_reactant)
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_reactions = len(LiquidStreamIn.equilibrium_K)
        num_of_objectives = num_of_reactions + num_of_volatile + 1
        num_of_species = LiquidStreamIn.num_of_species

        w = self.__mass_fraction_dictionary_to_array__(LiquidStreamIn)
        y = self.__molar_fraction_dictionary_to_array__(GasStreamIn)

        mm = LiquidStreamInit.get_solution_flow_kg_h() / LiquidStreamIn.get_solution_flow_kg_h()
        mn = LiquidStreamInit.get_solution_flow_kg_h() / GasStreamIn.get_gas_flow_kmol_h()

        dvdw = self.__dvdw__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        dwdr = np.einsum("s,wr->swr", mm, LiquidStreamIn.R)
        dvdy = self.__dvdy__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        dydr = np.zeros(shape=(num_of_samples, num_of_volatile, num_of_reactions))
        dwdp = np.einsum("s,wp->swp", mm, LiquidStreamIn.PL) + np.einsum("s,swp->swp", mm, np.einsum("sw,p->swp", w, np.diag(LiquidStreamIn.PG)))
        dydp = np.zeros(shape=(num_of_samples, num_of_volatile, num_of_volatile))
        for i in range(num_of_volatile):
            for j in range(num_of_volatile):
                dydp[:, i, j] = mn * (y[:, i] - 1) if i==j else mn * y[:, i]
        #dydp = np.swapaxes(a=dydp, axis1=1, axis2=2)

        dvdr = np.einsum("svw,swr->svr", dvdw, dwdr) + np.einsum("svy,syr->svr", dvdy, dydr)
        dvdp = np.einsum("svw,swp->svp", dvdw, dwdp) + np.einsum("svy,syp->svp", dvdy, dydp)
        dvdz = np.concatenate((dvdr, dvdp), axis=2)
        return dvdz

    def __dvdT__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        # K0 = LiquidStreamIn.__Kw_vle__(LiquidStreamIn)
        # LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + 0.05
        # K1 = LiquidStreamIn.__Kw_vle__(LiquidStreamIn)
        # LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - 0.05
        # dvdT = (1 / K0) * (K1 - K0) / 0.05
        # return dvdT[:,:,None]

        v0 = self.__v__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + 0.02
        lamda = LiquidStreamIn.__heat_of_vaporization_kJ_kmol__(LiquidStreamIn)
        Kw, dH, dKdT = LiquidStreamIn.__Kw_and_heat_of_reaction_kJ_kmol__(LiquidStreamIn)
        Kw_vle = LiquidStreamIn.__Kw_vle__(LiquidStreamIn)
        v1 = self.__v__(GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle)
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - 0.02
        dvdT = (v1 - v0) / 0.02
        return dvdT[:, :, None]

    def __dedz__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        dedz = np.concatenate((dH, lamda), axis=1)
        return dedz[:, None, :]

    def __dedT__(self, GasStreamIn, LiquidStreamIn, LiquidStreamInit, Kw, dH, dKdT, lamda, heat_kW, pressure_bara, Kw_vle):
        cp = LiquidStreamIn.get_solution_heat_capacity_kJ_kgK()
        dedT = - cp
        return dedT[:, None, None]

    def __mass_fraction_dictionary_to_array__(self, LiquidStreamIn):
        num_of_samples = LiquidStreamIn.temp_K.shape[0]
        num_of_species = LiquidStreamIn.num_of_species
        w = np.zeros(shape=(num_of_samples, num_of_species), dtype=np.float64)
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            w[:, i] = LiquidStreamIn.get_specie_mass_fraction(id=id)
        return w

    def __mass_fraction_array_to_dictionary__(self, LiquidStreamIn, w):
        for id in LiquidStreamIn.specie.keys():
            i = LiquidStreamIn.specie[id]["Index"]
            LiquidStreamIn.set_specie_mass_fraction(id=id, value=w[:,i])
        return LiquidStreamIn

    def __molar_fraction_dictionary_to_array__(self, GasStreamIn):
        num_of_samples = GasStreamIn.temp_K.shape[0]
        num_of_species = GasStreamIn.num_of_species
        y = np.zeros(shape=(num_of_samples, num_of_species), dtype=np.float64)
        for id in GasStreamIn.specie.keys():
            i = GasStreamIn.specie[id]["Index"]
            y[:, i] = GasStreamIn.get_specie_molar_fraction(id=id)
        return y

    def __molar_fraction_array_to_dictionary__(self, GasStreamIn, y):
        for id in GasStreamIn.specie.keys():
            i = GasStreamIn.specie[id]["Index"]
            GasStreamIn.set_specie_molar_fraction(id=id, value=y[:,i])
        return GasStreamIn


class Equilibrium_Liquid_TPFlash(FeatureBooster):

    def __init__(self):
        pass


class Equilibrium_Gas_TP(__Reactor_Common__):

    def __init__(self):
        super().__init__()

    def react(self, GasStreamIn):
        GasStreamOut, _ = self.__react__(GasStreamIn)
        return GasStreamOut

    def __react__(self, GasStreamIn):

        self.matrix = self.__load_stochiometry_gas_equilibrium__(GasStreamIn)

        # Dimension
        self.num_of_samples = len(GasStreamIn.get_gas_temp_K())

        # Initiate Outlet Stream (Which will be in Equilibrium)
        GasStreamOut = deepcopy(GasStreamIn)
        for id in GasStreamOut.specie.keys():
            i = GasStreamOut.specie[id]["Index"]
            GasStreamOut.set_specie_molar_fraction(id=id, value=np.maximum(GasStreamOut.get_specie_molar_fraction(id=id), 10 ** (-18)))
        GasStreamOut.normalize_molar_fractions()

        # Starting Point: Initial Concentrations
        y0 = GasStreamOut.__molar_fractions_dic2vec__()
        w0 = GasStreamOut.__molefrac2massfrac__(y0)
        y = 1.0 * y0
        w = 1.0 * w0

        K = self.__get_equilibrium_gas_K__(GasStreamOut)
        l = self.__get_equilibrium_gas_g__(GasStreamOut, K)
        dldp = self.__get_equilibrium_gas_dgdy__(GasStreamOut)
        dldq = np.einsum("src,cq->srq", dldp, self.matrix["R"])

        # Preparations...
        converged = False
        epoch = 0
        error = []
        lr = 0.9
        iterations = np.zeros(shape=(self.num_of_samples, ))

        # Iterate Until Convergence
        while converged == False:

            # Newton's Method
            dq_newton = - np.linalg.solve(dldq, l)
            dw_newton = np.einsum("ij,sj->si", self.matrix["R"], dq_newton)

            # Reducing Step Size (Slightly)
            dq = lr * dq_newton
            dw = lr * dw_newton

            # Backtrack to Ensure only Positive Concentrations
            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(- 0.9 * w / dw, nan=1.0, posinf=1.0, neginf=0.0)
            tau = 1.0 * (tau <= 0.0) + np.minimum(tau, 1.0) * (tau > 0.0)
            tau = tau * (w > 0) + 1.0 * (w <= 0)
            tau = np.min(tau, axis=1, keepdims=True)
            tau = np.broadcast_to(tau, shape=(self.num_of_samples, GasStreamIn.num_of_species)).copy()

            # Update Mass Fractions...
            w = w + tau * dw

            # Update Outlet Stream
            y = GasStreamOut.__massfrac2molefrac__(w)
            GasStreamOut.__molar_fractions_vec2dic__(y)
            GasStreamOut.normalize_molar_fractions()

            # Re-evaluate objective functions
            K = self.__get_equilibrium_gas_K__(GasStreamOut)
            l = self.__get_equilibrium_gas_g__(GasStreamOut, K)
            dldp = self.__get_equilibrium_gas_dgdy__(GasStreamOut)
            dldq = np.einsum("src,cq->srq", dldp, self.matrix["R"])

            # Check if Algorithm have Converged
            specie_converged_1 = np.array(np.abs(dw_newton) < 0.005 * np.abs(w), dtype=np.float32)      # 0.005
            specie_converged_2 = np.array(np.abs(dw_newton) < 10 ** (-18), dtype=np.float32)            # 10 ** (-16)
            specie_converged = np.minimum(specie_converged_1 + specie_converged_2, 1.0)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 198)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        print("Epochs: \t" + str(epoch))
        return GasStreamOut, iterations


class Equilibrium_VaporLiquid(__Reactor_Common__):

    def __init__(self, area_m2, height_m, void_fraction_m3_m3):
        super().__init__(area_m2, height_m, void_fraction_m3_m3)

    def react(self, GasStreamIn, LiquidStreamIn):
        GasStreamOut, LiquidStreamOut = self.__react__(GasStreamIn, LiquidStreamIn)
        return GasStreamOut, LiquidStreamOut

    def __react__(self, GasStreamIn, LiquidStreamIn):

        num_of_samples = GasStreamIn.temp_K.shape[0]

        self.matrix = self.__load_stochiometry_equilibrium__(GasStreamIn, LiquidStreamIn)

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
        lr = 0.9
        iterations = np.zeros(shape=(num_of_samples,))

        # Iterate Until Convergence
        while converged == False:

            # Mass Flow
            m_liq_tot = LiquidStreamOut.get_solution_flow_kg_h()
            m_gas_tot = GasStreamOut.get_gas_flow_kg_h()
            m_liq = m_liq_tot[:, None] * w_liq
            m_gas = m_gas_tot[:, None] * w_gas

            # Equilibrium Constants
            K_liq = self.__get_equilibrium_liq_K__(LiquidStreamOut)
            K_vap = self.__get_equilibrium_vap_K__(GasStreamOut, LiquidStreamOut)
            K_gas = self.__get_equilibrium_gas_K__(GasStreamOut)

            # Exothermic Heat from Reactions
            dKdT_liq, dH_liq = self.__get_equilibrium_liq_dKdT_q__(LiquidStreamOut, K_liq)
            dKdT_vap, dH_vap = self.__get_equilibrium_vap_dKdT_q__(GasStreamOut, LiquidStreamOut, K_vap)
            dKdT_gas, dH_gas = self.__get_equilibrium_gas_dKdT_q__(GasStreamOut, K_gas)
            dH = np.concatenate((dH_liq, dH_vap, dH_gas), axis=1)

            # Objective Function
            l = self.__get_equilibrium_liq_l__(LiquidStreamOut, K_liq)
            v = self.__get_equilibrium_vap_v__(GasStreamOut, LiquidStreamOut, K_vap)
            g = self.__get_equilibrium_gas_g__(GasStreamOut, K_gas)

            dm_liq = m_liq - m_in_liq
            dm_gas = m_gas - m_in_gas
            dm = np.concatenate((dm_liq, dm_gas), axis=1)
            e = (m_in_liq_tot / 3600) * cp_in_liq * (T - T_in_liq) + (n_in_gas_tot / 3600) * cp_in_gas * (T - T_in_gas) - np.einsum("sm,sm->s", np.einsum("sr,rm->sm", dH, self.matrix["R+"]), dm / 3600)
            f = np.concatenate((l, v, g, e[:, None]), axis=1)

            # Partial Derivatives
            dldw = self.__get_equilibrium_liq_dldw__(LiquidStreamOut)
            dvdw = self.__get_equilibrium_vap_dvdw__(LiquidStreamOut)
            dvdy = self.__get_equilibrium_vap_dvdy__(GasStreamOut, LiquidStreamOut)
            dgdy = self.__get_equilibrium_gas_dgdy__(GasStreamOut)

            dmdz = np.einsum("s,wr->swr", m_in_tot, self.matrix["R"])
            dmdz_liq = dmdz[:, 0:LiquidStreamOut.num_of_species:1,:]
            dmdz_gas = dmdz[:, LiquidStreamOut.num_of_species::,:]
            dmdz_liq_tot = np.sum(dmdz_liq, axis=1, keepdims=False)
            dmdz_gas_tot = np.sum(dmdz_gas, axis=1, keepdims=False)

            dwdz_liq = np.einsum("s,smz->smz", m_liq_tot ** (-2), np.einsum("s,smz->smz", m_liq_tot, dmdz_liq) - np.einsum("sm,sz->smz", m_liq, dmdz_liq_tot))
            dwdz_gas = np.einsum("s,smz->smz", m_gas_tot ** (-2), np.einsum("s,smz->smz", m_gas_tot, dmdz_gas) - np.einsum("sm,sz->smz", m_gas, dmdz_gas_tot))

            dydw_gas = GasStreamOut.__dydw__(w_gas)

            dldz = np.einsum("slw,swz->slz", dldw, dwdz_liq)
            dvdz = np.einsum("slw,swz->slz", dvdw, dwdz_liq) + np.einsum("slw,swz->slz", np.einsum("sly,syw->slw", dvdy, dydw_gas), dwdz_gas)
            dgdz = np.einsum("slw,swz->slz", np.einsum("sly,syw->slw", dgdy, dydw_gas), dwdz_gas)
            dedz = - np.einsum("s,sr->sr", m_in_tot / 3600, dH)[:, None, :]

            dldT = self.__get_equilibrium_liq_dldT__(LiquidStreamOut, K_liq, dKdT_liq)
            dvdT = self.__get_equilibrium_vap_dvdT__(GasStreamOut, LiquidStreamOut, K_vap, dKdT_vap)
            dgdT = self.__get_equilibrium_gas_dgdT__(GasStreamOut, K_gas, dKdT_gas)
            dedT = (m_in_liq_tot / 3600) * cp_in_liq + (n_in_gas_tot / 3600) * cp_in_gas

            dfdz = np.concatenate((dldz, dvdz, dgdz, dedz), axis=1)
            dfdT = np.concatenate((dldT, dvdT, dgdT, dedT[:, None, None]), axis=1)
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

            # Mass Fractions and Temp
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

        print("Epochs: \t" + str(epoch))
        return GasStreamOut, LiquidStreamOut


# ---------------------------------------------------------------------------------------


class LiquidReactor_HeatExchanger_CounterFlow(FeatureBooster):

    def __init__(self, id, num_of_height, area_m2, heat_transfer_coefficient_kW_m2K):
        super().__init__()
        self.num_of_heights = num_of_height
        self.area_m2 = area_m2
        self.zeta = np.linspace(0, 1, self.num_of_heights)
        self.heat_transfer_coefficient_kW_m2K = heat_transfer_coefficient_kW_m2K

    def react(self, LiquidStreamIn1, LiquidStreamIn2):
        LiquidStreamIn2.flow_kg_h = - np.abs(LiquidStreamIn2.flow_kg_h)
        self.num_of_samples = len(LiquidStreamIn1.temp_K)
        shape = np.ones(shape=(self.num_of_heights,))
        self.LiquidStreamIn1 = deepcopy(LiquidStreamIn1)
        self.LiquidStreamIn2 = deepcopy(LiquidStreamIn2)
        self.LiquidStreamOut1 = deepcopy(LiquidStreamIn1)
        self.LiquidStreamOut2 = deepcopy(LiquidStreamIn2)
        # ------------------------------------------------------------
        for self.i in range(self.num_of_samples):
            # ------------------------------------------------------------
            self.LiquidStream1 = deepcopy(LiquidStreamIn1)
            for key in LiquidStreamIn1.specie.keys():
                self.LiquidStream1.specie[key]["Mass Fraction"] = LiquidStreamIn1.specie[key]["Mass Fraction"][self.i] * shape
            self.LiquidStream1.temp_K = LiquidStreamIn1.temp_K[self.i] * shape
            self.LiquidStream1.flow_kg_h = LiquidStreamIn1.flow_kg_h[self.i] * shape
            # ------------------------------------------------------------
            self.LiquidStream2 = deepcopy(LiquidStreamIn2)
            for key in LiquidStreamIn2.specie.keys():
                self.LiquidStream2.specie[key]["Mass Fraction"] = LiquidStreamIn2.specie[key]["Mass Fraction"][self.i] * shape
            self.LiquidStream2.temp_K = LiquidStreamIn2.temp_K[self.i] * shape
            self.LiquidStream2.flow_kg_h = LiquidStreamIn2.flow_kg_h[self.i] * shape
            # ------------------------------------------------------------
            y = np.vstack((self.LiquidStream1.temp_K, self.LiquidStream2.temp_K))
            problem = solve_bvp(fun=self.__system_dynamics__, bc=self.__boundary_conditions__, x=self.zeta, y=y)
            # ------------------------------------------------------------
            x = np.array([1])
            self.LiquidStreamOut1.temp_K[self.i] = problem.sol(x)[0,0]
            x = np.array([0])
            self.LiquidStreamOut2.temp_K[self.i] = problem.sol(x)[1,0]
            # ------------------------------------------------------------
        return self.LiquidStreamOut1, self.LiquidStreamOut2

    def __system_dynamics__(self, x, y):

        N = y.shape[1]
        shape = np.ones(shape=(N,))

        self.LiquidStream1 = deepcopy(self.LiquidStreamIn1)
        self.LiquidStream1.temp_K = y[0, :]
        self.LiquidStream1.flow_kg_h = self.LiquidStreamIn1.flow_kg_h[self.i] * shape
        for key in self.LiquidStreamIn1.specie.keys():
            self.LiquidStream1.specie[key]["Mass Fraction"] = self.LiquidStreamIn1.specie[key]["Mass Fraction"][self.i] * shape

        self.LiquidStream2 = deepcopy(self.LiquidStreamIn2)
        self.LiquidStream2.temp_K = y[1, :]
        self.LiquidStream2.flow_kg_h = self.LiquidStreamIn2.flow_kg_h[self.i] * shape
        for key in self.LiquidStreamIn2.specie.keys():
            self.LiquidStream2.specie[key]["Mass Fraction"] = self.LiquidStreamIn2.specie[key]["Mass Fraction"][self.i] * shape

        # Heat Transfer Coefficient, Heat Transfer and Heat Capacity
        self.kH_kW_m2K = self.heat_transfer_coefficient_kW_m2K(self)
        qH = self.kH_kW_m2K * (self.LiquidStream1.temp_K - self.LiquidStream2.temp_K)
        cp1 = self.LiquidStream1.get_solution_heat_capacity_kJ_kgK()
        cp2 = self.LiquidStream2.get_solution_heat_capacity_kJ_kgK()

        # Conservation of Energy
        dT1 = - 3600 * (1 / cp1) * self.area_m2 * (qH / self.LiquidStream1.get_solution_flow_kg_h())
        dT2 = 3600 * (1 / cp2) * self.area_m2 * (qH / self.LiquidStream2.get_solution_flow_kg_h())
        dydx = np.vstack((dT1, dT2))
        return dydx

    def __boundary_conditions__(self, ya, yb):
        con = np.zeros(shape=(2,))
        con[0] = ya[0] - self.LiquidStreamIn1.temp_K[self.i]
        con[1] = yb[1] - self.LiquidStreamIn2.temp_K[self.i]
        return con


class LiquidReactor_HeatExchanger_CounterFlow_NTU(FeatureBooster):

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


class LiquidReactor_CSTR_Isothermal(__Reactor_Common__):

    def __init__(self, train=True, max_error_pct=25, degrees_of_freedom=None):
        super().__init__()
        self._train_ = train
        self._max_error_pct_ = max_error_pct
        self._degrees_of_freedom_ = degrees_of_freedom

    def react(self, LiquidStreamIn, volume_m3):

        self.volume_m3 = volume_m3

        self._num_of_samples_ = LiquidStreamIn.temp_K.shape[0]
        self._num_of_reactions_ = LiquidStreamIn.num_of_reaction_liq
        self._num_of_equilibriums_ = LiquidStreamIn.num_of_equilibrium_liq
        self._num_of_species_ = LiquidStreamIn.num_of_species
        self._num_of_conservations_ = self._num_of_species_ - self._num_of_equilibriums_ - self._num_of_reactions_

        self._matrix_ = self.__load_stochiometry_liq_all__(LiquidStreamIn)
        self._matrix_["A broadcast"] = np.broadcast_to(array=self._matrix_["A"][None, :, :], shape=(self._num_of_samples_, self._matrix_["A"].shape[0], self._matrix_["A"].shape[1]))
        self._matrix_["R+ broadcast"] = np.broadcast_to(array=self._matrix_["R+"][None, :, :], shape=(self._num_of_samples_, self._matrix_["R+"].shape[0], self._matrix_["R+"].shape[1]))

        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        self._m0_ = self.LiquidStreamIn.get_solution_flow_kg_h() / 3600
        self._w0_ = self.LiquidStreamIn.__mass_fractions_dic2vec__()
        self._b0_ = np.einsum("cw,sw->sc", self._matrix_["A"], self._w0_)
        self._omega_ = self.volume_m3 / self._m0_

        self.LiquidStreamOut = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamOut.specie.keys():
            i = self.LiquidStreamOut.specie[id]["Index"]
            self.LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.maximum(self.LiquidStreamOut.get_specie_mass_fraction(id=id), self.w_min_init))
        self.LiquidStreamOut.normalize_mass_fractions()
        self._w_ = self.LiquidStreamOut.__mass_fractions_dic2vec__()

        # Preparations...
        converged = False
        epoch = 0
        lr = 0.9
        iterations = np.zeros(shape=(self._num_of_samples_,))

        # Iterate Until Convergence
        while converged == False:

            # Objective Function and Partial Derivatives
            f = self.__f__()
            dfdw = self.__dfdw__()

            # Damped Newton's Method
            dw_newton = - np.linalg.solve(dfdw, f)
            self._w_ = self._w_ + lr * dw_newton
            self._w_ = np.maximum(self._w_, self.w_min)

            # Update Outlet Stream
            self.LiquidStreamOut.__mass_fractions_vec2dic__(self._w_)

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dw_newton) < 0.005 * np.abs(self._w_), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 498)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        print("Epochs: \t" + str(epoch))
        return self.LiquidStreamOut

    def __f__(self):
        self.forward_rate_kmol_m3s = np.zeros(shape=(self._num_of_samples_, self._num_of_reactions_))
        self.backward_rate_kmol_m3s = np.zeros(shape=(self._num_of_samples_, self._num_of_reactions_))
        for i, id in enumerate(self.LiquidStreamOut.reaction_liq.keys()):
            self.forward_rate_kmol_m3s[:, i], self.backward_rate_kmol_m3s[:, i] = self.LiquidStreamOut.reaction_liq[id]["Rate [kmol/m3.s]"](self.LiquidStreamOut)
        self.rate_kmol_m3s = self.forward_rate_kmol_m3s - self.backward_rate_kmol_m3s

        self.K = self.__get_equilibrium_liq_K__(self.LiquidStreamOut)
        l1 = self.__get_equilibrium_liq_l2__(self.LiquidStreamOut, self.K)
        l2 = np.einsum("rw,sw->sr", self._matrix_["R+"], self._w_ - self._w0_)[:, self._num_of_equilibriums_::] - self._omega_[:,None] * self.rate_kmol_m3s
        l3 = np.einsum("cw,sw->sc", self._matrix_["A"], self._w_) - self._b0_
        l = np.concatenate((l1, l2, l3), axis=1)
        return l

    def __dfdw__(self):

        drdw = np.zeros(shape=(self._num_of_samples_, self._num_of_reactions_, self._num_of_species_), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(self.LiquidStreamOut.reaction_liq.keys()):
            for specie_id in self.LiquidStreamOut.specie.keys():
            #for _, specie_id in enumerate(self.LiquidStreamOut.reaction_liq[rxn_id]["Rate Dependencies"]):
                specie_i = self.LiquidStreamOut.specie[specie_id]["Index"]
                w_probe = self.LiquidStreamOut.get_specie_mass_fraction(id=specie_id)
                dw_probe = np.maximum(0.01 * w_probe, self.dw_min_td)
                self.LiquidStreamOut.set_specie_mass_fraction(id=specie_id, value=w_probe + dw_probe)
                forward_rate_pertubation, backward_rate_pertubation = self.LiquidStreamOut.reaction_liq[rxn_id]["Rate [kmol/m3.s]"](self.LiquidStreamOut)
                rate_pertubation = forward_rate_pertubation - backward_rate_pertubation
                self.LiquidStreamOut.set_specie_mass_fraction(id=specie_id, value=w_probe)
                drdw[:, rxn_i, specie_i] = (rate_pertubation - self.rate_kmol_m3s[:, rxn_i]) / dw_probe

        dldw1 = self.__get_equilibrium_liq_dldw2td__(self.LiquidStreamOut, self.K)
        dldw2 = self._matrix_["R+ broadcast"][:, self._num_of_equilibriums_::, :] - np.einsum("s,srw->srw", self._omega_, drdw)
        dldw3 = self._matrix_["A broadcast"]
        dldw = np.concatenate((dldw1, dldw2, dldw3), axis=1)
        return dldw


class LiquidReactor_CSTR_Adiabatic(__Reactor_Common__):

    def __init__(self, train=True, max_error_pct=25, degrees_of_freedom=None):
        super().__init__()
        self._train_ = train
        self._max_error_pct_ = max_error_pct
        self._degrees_of_freedom_ = degrees_of_freedom

    def react(self, LiquidStreamIn, volume_m3):

        self.volume_m3 = volume_m3

        self._num_of_samples_ = LiquidStreamIn.temp_K.shape[0]
        self._num_of_reactions_ = LiquidStreamIn.num_of_reaction_liq
        self._num_of_equilibriums_ = LiquidStreamIn.num_of_equilibrium_liq
        self._num_of_species_ = LiquidStreamIn.num_of_species
        self._num_of_conservations_ = self._num_of_species_ - self._num_of_equilibriums_ - self._num_of_reactions_

        self._matrix_ = self.__load_stochiometry_liq_all__(LiquidStreamIn)
        self._matrix_["A broadcast"] = np.broadcast_to(array=self._matrix_["A"][None, :, :], shape=(self._num_of_samples_, self._matrix_["A"].shape[0], self._matrix_["A"].shape[1]))
        self._matrix_["R+ broadcast"] = np.broadcast_to(array=self._matrix_["R+"][None, :, :], shape=(self._num_of_samples_, self._matrix_["R+"].shape[0], self._matrix_["R+"].shape[1]))

        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        self._m0_ = self.LiquidStreamIn.get_solution_flow_kg_h() / 3600
        self._w0_ = self.LiquidStreamIn.__mass_fractions_dic2vec__()
        self._b0_ = np.einsum("cw,sw->sc", self._matrix_["A"], self._w0_)
        self._T0_ = self.LiquidStreamIn.temp_K
        self._omega_ = self.volume_m3 / self._m0_

        self.LiquidStreamOut = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamOut.specie.keys():
            i = self.LiquidStreamOut.specie[id]["Index"]
            self.LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.maximum(self.LiquidStreamOut.get_specie_mass_fraction(id=id), self.w_min_init))
        self.LiquidStreamOut.normalize_mass_fractions()
        self._w_ = self.LiquidStreamOut.__mass_fractions_dic2vec__()
        self._T_ = self.LiquidStreamOut.get_solution_temp_K()

        # Preparations...
        converged = False
        epoch = 0
        lr = 0.75
        iterations = np.zeros(shape=(self._num_of_samples_,))

        # Iterate Until Convergence
        while converged == False:

            # Objective Function and Partial Derivatives
            f = self.__f__()
            dfdX = self.__dfdX__()

            # Damped Newton's Method
            dX_newton = - np.linalg.solve(dfdX, f)
            dw_newton = dX_newton[:, :dX_newton.shape[1] - 1:]
            dT_newton = dX_newton[:, dX_newton.shape[1] - 1]
            dw = lr * dw_newton
            dT = lr * dT_newton

            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(x=np.minimum(5.0 / np.abs(dT), 1.0), nan=1.0, posinf=1.0, neginf=1.0)

            dT = tau * dT
            dw = tau[:, None] * dw

            self._w_ = self._w_ + dw
            self._w_ = np.maximum(self._w_, self.w_min)
            self._T_ = self._T_ + dT

            # Update Outlet Stream
            self.LiquidStreamOut.__mass_fractions_vec2dic__(self._w_)
            self.LiquidStreamOut.set_solution_temp_K(value=self._T_)

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dw_newton) < 0.005 * np.abs(self._w_), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 498)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        print("Epochs: \t" + str(epoch))
        return self.LiquidStreamOut

    def __f__(self):

        self.forward_rate_kmol_m3s = np.zeros(shape=(self._num_of_samples_, self._num_of_reactions_))
        self.backward_rate_kmol_m3s = np.zeros(shape=(self._num_of_samples_, self._num_of_reactions_))

        self.K = self.__get_equilibrium_liq_K__(self.LiquidStreamOut)

        self._dKdT_, self._q_equ_ = self.__get_equilibrium_liq_dKdT_q__(self.LiquidStreamOut, self.K)
        self._q_rxn_ = self.__get_reaction_liq_q__(self.LiquidStreamOut)
        self._q_ = np.concatenate((self._q_equ_, self._q_rxn_), axis=1)
        self._cp_ = self.LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()

        for i, id in enumerate(self.LiquidStreamOut.reaction_liq.keys()):
            self.forward_rate_kmol_m3s[:, i], self.backward_rate_kmol_m3s[:, i] = self.LiquidStreamOut.reaction_liq[id]["Rate [kmol/m3.s]"](self.LiquidStreamOut)
        self.rate_kmol_m3s = self.forward_rate_kmol_m3s - self.backward_rate_kmol_m3s

        f1 = self.__get_equilibrium_liq_l2__(self.LiquidStreamOut, self.K)
        f2 = np.einsum("rw,sw->sr", self._matrix_["R+"], self._w_ - self._w0_)[:, self._num_of_equilibriums_::] - self._omega_[:,None] * self.rate_kmol_m3s
        f3 = np.einsum("cw,sw->sc", self._matrix_["A"], self._w_) - self._b0_
        f4 = np.einsum("sw,sw->s", np.einsum("sr,rw->sw", self._q_, self._matrix_["R+"]), self._w_ - self._w0_) / self._cp_ - (self._T_ - self._T0_)

        f = np.concatenate((f1, f2, f3, f4[:, None]), axis=1)
        return f

    def __dfdX__(self):
        drdw = np.zeros(shape=(self._num_of_samples_, self._num_of_reactions_, self._num_of_species_), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(self.LiquidStreamOut.reaction_liq.keys()):
            for specie_id in self.LiquidStreamOut.specie.keys():
            #for _, specie_id in enumerate(self.LiquidStreamOut.reaction_liq[rxn_id]["Rate Dependencies"]):
                specie_i = self.LiquidStreamOut.specie[specie_id]["Index"]
                w_probe = self.LiquidStreamOut.get_specie_mass_fraction(id=specie_id)
                dw_probe = np.maximum(0.01 * w_probe, self.dw_min_td)
                self.LiquidStreamOut.set_specie_mass_fraction(id=specie_id, value=w_probe + dw_probe)
                forward_rate_pertubation, backward_rate_pertubation = self.LiquidStreamOut.reaction_liq[rxn_id]["Rate [kmol/m3.s]"](self.LiquidStreamOut)
                rate_pertubation = forward_rate_pertubation - backward_rate_pertubation
                self.LiquidStreamOut.set_specie_mass_fraction(id=specie_id, value=w_probe)
                drdw[:, rxn_i, specie_i] = (rate_pertubation - self.rate_kmol_m3s[:, rxn_i]) / dw_probe

        drdT = np.zeros(shape=(self._num_of_samples_, self._num_of_reactions_, 1), dtype=np.float64)

        for rxn_i, rxn_id in enumerate(self.LiquidStreamOut.reaction_liq.keys()):
            self.LiquidStreamOut.set_solution_temp_K(value=self.LiquidStreamOut.get_solution_temp_K() + 0.01)
            forward_rate_pertubation, backward_rate_pertubation = self.LiquidStreamOut.reaction_liq[rxn_id]["Rate [kmol/m3.s]"](self.LiquidStreamOut)
            rate_pertubation = forward_rate_pertubation - backward_rate_pertubation
            self.LiquidStreamOut.set_solution_temp_K(value=self.LiquidStreamOut.get_solution_temp_K() - 0.01)
            drdT[:, rxn_i, 0] = (rate_pertubation - self.rate_kmol_m3s[:, rxn_i]) / 0.01

        dfdw1 = self.__get_equilibrium_liq_dldw2td__(self.LiquidStreamOut, self.K)
        dfdw2 = self._matrix_["R+ broadcast"][:, self._num_of_equilibriums_::, :] - np.einsum("s,srw->srw", self._omega_, drdw)
        dfdw3 = self._matrix_["A broadcast"]
        dfdw4 = (np.einsum("sr,rw->sw", self._q_, self._matrix_["R+"]) / self._cp_[:, None])[:, None, :]


        dfdT1 = self._dKdT_[:, :, None]
        for rxn_i, rxn_id in enumerate(self.LiquidStreamOut.equilibrium_liq.keys()):
            for specie_id in self.LiquidStreamOut.equilibrium_liq[rxn_id]["Stoch"].keys():
                nu = self.LiquidStreamOut.equilibrium_liq[rxn_id]["Stoch"][specie_id]
                if nu < 0:
                    dfdT1[:,rxn_i,0] = dfdT1[:,rxn_i,0] * self.LiquidStreamOut.get_specie_mass_fraction(id=specie_id) ** (-nu)


        dfdT2 = - np.einsum("s,srw->srw", self._omega_, drdT)
        dfdT3 = np.zeros(shape=(self._num_of_samples_, self._num_of_conservations_, 1), dtype=np.float64)
        dfdT4 = - np.ones(shape=(self._num_of_samples_, 1, 1), dtype=np.float64)

        dfdw = np.concatenate((dfdw1, dfdw2, dfdw3, dfdw4), axis=1)
        dfdT = np.concatenate((dfdT1, dfdT2, dfdT3, dfdT4), axis=1)

        dfdX = np.concatenate((dfdw, dfdT), axis=2)
        return dfdX


# ---------------------------------------------------------------------------------------


class __GasLiquidContactor__(__Reactor_Common__):

    def __init__(self):

        self.info = {}
        self.mass_transfer_kmol_m3s = {}
        self.num_of_mass_transfer = 0
        self.heat_transfer_kW_m3 = None
        self.liquid_holdup_m3_m3 = None
        self.GasStream = None
        self.LiquidStream = None

    def add_info(self, key, value):
        self.info[key] = value

    def add_mass_transfer_insta(self, id, stoch_gas, stoch_liq, equilibrium_constant):
        pass

    def add_mass_transfer_kmol_m3s(self, id, stoch_gas, stoch_liq, rate_kmol_m3s, exothermic_heat_kJ_kmol=None):
        self.mass_transfer_kmol_m3s[id] = {}
        self.mass_transfer_kmol_m3s[id]["Stoch Gas"] = stoch_gas
        self.mass_transfer_kmol_m3s[id]["Stoch Liq"] = stoch_liq
        self.mass_transfer_kmol_m3s[id]["Rate [kmol/m3.s]"] = rate_kmol_m3s
        self.mass_transfer_kmol_m3s[id]["Exothermic Heat [kJ/kmol]"] = exothermic_heat_kJ_kmol
        self.mass_transfer_kmol_m3s[id]["Index"] = self.num_of_mass_transfer
        self.num_of_mass_transfer = self.num_of_mass_transfer + 1

    def add_heat_transfer_kW_m3(self, id, heat_transfer_kW_m3):
        self.heat_transfer_kW_m3 = heat_transfer_kW_m3

    def add_liquid_holdup_m3_m3(self, id, liquid_holdup_m3_m3):
        self.liquid_holdup_m3_m3 = liquid_holdup_m3_m3

    def add_pressure_drop_Pa_m(self, id, presure_drop_Pa_m):
        pass


class __GasLiquidContactor_Column__(__GasLiquidContactor__):

    def __init__(self, area_m2, height_m, void_fraction_m3_m3):
        self.area_m2 = area_m2
        self.height_m = height_m
        self.void_fraction_m3_m3 = void_fraction_m3_m3
        self.volume_m3 = area_m2 * height_m * void_fraction_m3_m3

    def get_superficial_gas_velocity_m_s(self):
        pass

    def get_superficial_liquid_velocity_m_s(self):
        pass


class __GasLiquidContactor_Centrifuge__(__GasLiquidContactor__):

    def __init__(self, area_m2, height_m, void_fraction_m3_m3):
        self.area_m2 = area_m2
        self.height_m = height_m
        self.void_fraction_m3_m3 = void_fraction_m3_m3
        self.volume_m3 = area_m2 * height_m * void_fraction_m3_m3

    def get_superficial_gas_velocity_m_s(self):
        pass

    def get_superficial_liquid_velocity_m_s(self):
        pass


class GasLiquidContactor_Column_Stirred(__GasLiquidContactor_Column__):

    def __init__(self, area_m2, height_m, void_fraction_m3_m3):
        super().__init__(area_m2, height_m, void_fraction_m3_m3)

    def react(self, GasStreamIn, LiquidStreamIn):

        self.volume_m3 = area_m2 * height_m * void_fraction_m3_m3

        self._num_of_samples_ = LiquidStreamIn.temp_K.shape[0]
        # ----
        self._num_of_liq_reactions_ = LiquidStreamIn.num_of_reaction_liq
        self._num_of_liq_equilibriums_ = LiquidStreamIn.num_of_equilibrium_liq
        self._num_of_liq_species_ = LiquidStreamIn.num_of_species
        self._num_of_liq_conservations_ = self._num_of_liq_species_ - self._num_of_liq_equilibriums_ - self._num_of_liq_reactions_
        # ----
        self._num_of_gas_reactions_ = GasStreamIn.num_of_reaction_gas
        self._num_of_gas_equilibriums_ = GasStreamIn.num_of_equilibrium_gas
        self._num_of_gas_species_ = GasStreamIn.num_of_species
        self._num_of_gas_conservations_ = self._num_of_gas_species_ - self._num_of_gas_equilibriums_ - self._num_of_gas_reactions_
        # ---

        self._matrix_ = self.__load_stochiometry_liq_all__(LiquidStreamIn)
        self._matrix_["A broadcast"] = np.broadcast_to(array=self._matrix_["A"][None, :, :], shape=(self._num_of_samples_, self._matrix_["A"].shape[0], self._matrix_["A"].shape[1]))
        self._matrix_["R+ broadcast"] = np.broadcast_to(array=self._matrix_["R+"][None, :, :], shape=(self._num_of_samples_, self._matrix_["R+"].shape[0], self._matrix_["R+"].shape[1]))

        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        self._m0_ = self.LiquidStreamIn.get_solution_flow_kg_h() / 3600
        self._w0_ = self.LiquidStreamIn.__mass_fractions_dic2vec__()
        self._b0_ = np.einsum("cw,sw->sc", self._matrix_["A"], self._w0_)
        self._T0_ = self.LiquidStreamIn.temp_K
        self._omega_ = self.volume_m3 / self._m0_

        self.LiquidStreamOut = deepcopy(LiquidStreamIn)
        for id in self.LiquidStreamOut.specie.keys():
            i = self.LiquidStreamOut.specie[id]["Index"]
            self.LiquidStreamOut.set_specie_mass_fraction(id=id, value=np.maximum(self.LiquidStreamOut.get_specie_mass_fraction(id=id), self.w_min_init))
        self.LiquidStreamOut.normalize_mass_fractions()
        self._w_ = self.LiquidStreamOut.__mass_fractions_dic2vec__()
        self._T_ = self.LiquidStreamOut.get_solution_temp_K()

        # Preparations...
        converged = False
        epoch = 0
        lr = 0.9
        iterations = np.zeros(shape=(self._num_of_samples_,))

        # Iterate Until Convergence
        while converged == False:

            # Objective Function and Partial Derivatives
            f = self.__f__()
            dfdX = self.__dfdX__()

            # Damped Newton's Method
            dX_newton = - np.linalg.solve(dfdX, f)
            dw_newton = dX_newton[:, :dX_newton.shape[1] - 1:]
            dT_newton = dX_newton[:, dX_newton.shape[1] - 1]
            dw = lr * dw_newton
            dT = lr * dT_newton

            with np.errstate(divide='ignore', invalid='ignore'):
                tau = np.nan_to_num(x=np.minimum(5 / np.abs(dT), 1.0), nan=1.0, posinf=1.0, neginf=1.0)

            dT = tau * dT
            dw = tau[:, None] * dw

            self._w_ = self._w_ + dw
            self._w_ = np.maximum(self._w_, self.w_min)
            self._T_ = self._T_ + dT

            # Update Outlet Stream
            self.LiquidStreamOut.__mass_fractions_vec2dic__(self._w_)
            self.LiquidStreamOut.set_solution_temp_K(value=self._T_)

            # Check if Algorithm have Converged
            specie_converged = np.array(np.abs(dw_newton) < 0.005 * np.abs(self._w_), dtype=np.float32)
            sample_converged = np.min(specie_converged, axis=1)
            converged = np.min(sample_converged, axis=0)
            converged = (bool(converged) or (epoch > 498)) and (epoch > 0)
            iterations = iterations + (1 - sample_converged)
            epoch = epoch + 1

        print("Epochs: \t" + str(epoch))
        return self.LiquidStreamOut

    def __f__(self):

        self.rate_kmol_m3s = np.zeros(shape=(self._num_of_samples_, self._num_of_reactions_))
        self.K = self.__get_equilibrium_liq_K__(self.LiquidStreamOut)

        self._dKdT_, self._q_equ_ = self.__get_equilibrium_liq_dKdT_q__(self.LiquidStreamOut, self.K)
        self._q_rxn_ = self.__get_reaction_liq_q__(self.LiquidStreamOut)
        self._q_ = np.concatenate((self._q_equ_, self._q_rxn_), axis=1)
        self._cp_ = self.LiquidStreamOut.get_solution_heat_capacity_kJ_kgK()

        for i, id in enumerate(self.LiquidStreamOut.reaction_liq.keys()):
            self.rate_kmol_m3s[:, i] = self.LiquidStreamOut.reaction_liq[id]["Rate [kmol/m3.s]"](self.LiquidStreamOut)

        f1 = self.__get_equilibrium_liq_l_option_2__(self.LiquidStreamOut, self.K)
        f2 = np.einsum("rw,sw->sr", self._matrix_["R+"], self._w_ - self._w0_)[:, self._num_of_equilibriums_::] - self._omega_[:,None] * self.rate_kmol_m3s
        f3 = np.einsum("cw,sw->sc", self._matrix_["A"], self._w_) - self._b0_
        f4 = np.einsum("sw,sw->s", np.einsum("sr,rw->sw", self._q_, self._matrix_["R+"]), self._w_ - self._w0_) / self._cp_ - (self._T_ - self._T0_)

        f = np.concatenate((f1, f2, f3, f4[:, None]), axis=1)
        return f

    def __dfdX__(self):
        drdw = np.zeros(shape=(self._num_of_samples_, self._num_of_reactions_, self._num_of_species_), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(self.LiquidStreamOut.reaction_liq.keys()):
            for _, specie_id in enumerate(self.LiquidStreamOut.reaction_liq[rxn_id]["Rate Dependencies"]):
                specie_i = self.LiquidStreamOut.specie[specie_id]["Index"]
                w_probe = self.LiquidStreamOut.get_specie_mass_fraction(id=specie_id)
                dw_probe = np.maximum(0.01 * w_probe, self.dw_min_td)
                self.LiquidStreamOut.set_specie_mass_fraction(id=specie_id, value=w_probe + dw_probe)
                rate_pertubation = self.LiquidStreamOut.reaction_liq[rxn_id]["Rate [kmol/m3.s]"](self.LiquidStreamOut)
                self.LiquidStreamOut.set_specie_mass_fraction(id=specie_id, value=w_probe)
                drdw[:, rxn_i, specie_i] = (rate_pertubation - self.rate_kmol_m3s[:, rxn_i]) / dw_probe

        drdT = np.zeros(shape=(self._num_of_samples_, self._num_of_reactions_, 1), dtype=np.float64)
        for rxn_i, rxn_id in enumerate(self.LiquidStreamOut.reaction_liq.keys()):
            self.LiquidStreamOut.set_solution_temp_K(value=self.LiquidStreamOut.get_solution_temp_K() + 0.1)
            rate_pertubation = self.LiquidStreamOut.reaction_liq[rxn_id]["Rate [kmol/m3.s]"](self.LiquidStreamOut)
            self.LiquidStreamOut.set_solution_temp_K(value=self.LiquidStreamOut.get_solution_temp_K() - 0.1)
            drdT[:, rxn_i, 0] = (rate_pertubation - self.rate_kmol_m3s[:, rxn_i]) / 0.1

        dfdw1 = self.__get_equilibrium_liq_dldw_option_2_td__(self.LiquidStreamOut, self.K)
        dfdw2 = self._matrix_["R+ broadcast"][:, self._num_of_equilibriums_::, :] - np.einsum("s,srw->srw", self._omega_, drdw)
        dfdw3 = self._matrix_["A broadcast"]
        dfdw4 = (np.einsum("sr,rw->sw", self._q_, self._matrix_["R+"]) / self._cp_[:, None])[:, None, :]


        dfdT1 = self._dKdT_[:, :, None]
        for rxn_i, rxn_id in enumerate(self.LiquidStreamOut.equilibrium_liq.keys()):
            for specie_id in self.LiquidStreamOut.equilibrium_liq[rxn_id]["Stoch"].keys():
                nu = self.LiquidStreamOut.equilibrium_liq[rxn_id]["Stoch"][specie_id]
                if nu < 0:
                    dfdT1[:,rxn_i,0] = dfdT1[:,rxn_i,0] * self.LiquidStreamOut.get_specie_mass_fraction(id=specie_id) ** (-nu)

        dfdT2 = - np.einsum("s,srw->srw", self._omega_, drdT)
        dfdT3 = np.zeros(shape=(self._num_of_samples_, self._num_of_conservations_, 1), dtype=np.float64)
        dfdT4 = - np.ones(shape=(self._num_of_samples_, 1, 1), dtype=np.float64)

        dfdw = np.concatenate((dfdw1, dfdw2, dfdw3, dfdw4), axis=1)
        dfdT = np.concatenate((dfdT1, dfdT2, dfdT3, dfdT4), axis=1)

        dfdX = np.concatenate((dfdw, dfdT), axis=2)
        return dfdX


class GasLiquidContactor_Column_PlugFlow_CounterCurrent(__GasLiquidContactor_Column__):

    def __init__(self, area_m2, height_m, void_fraction_m3_m3):
        pass

    def add_info(self, key, value):
        pass

    def react(self, GasStreamIn, LiquidStreamIn):
        pass


class GasLiquidContactor_Column_PlugFlow_CoCurrent(__GasLiquidContactor_Column__):

    def __init__(self):
        pass


class GasLiquidContactor_Column_GasStirred_LiquidPlugFlow(__GasLiquidContactor_Column__):

    def __init__(self):
        pass


class GasLiquidContactor_Column_GasPlugFlow_LiquidStirred(__GasLiquidContactor_Column__):

    def __init__(self):
        pass


class GasLiquidContactor_Centrifuge_PlugFlow_CounterCurrent(__GasLiquidContactor_Column__):

    def __init__(self, radius_inner_m, radius_outer_m, axial_height_m):
        pass

    def react(self, GasStreamIn, LiquidStreamin):
        pass


class GasLiquidContactor_Centrifuge_PlugFlow_CoCurrent(__GasLiquidContactor_Column__):

    def __init__(self):
        pass




# ---- COMPRESSORS (IDEAL GAS, NO CHEMICAL REACTIONS) ---------------------------------

class Compressor_Isentropic(FeatureBooster):

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


class Compressor_Isothermal(FeatureBooster):

    def __init__(self):
        super().__init__()



# ---- GAS-LIQUID CONTACTORS -----------------------------------------------------------

class PackedColumn_CoFlow(FeatureBooster):

    def __init__(self):
        super().__init__()

    def react(self):
        pass


class PackedColumn_CounterFlow(FeatureBooster):

    def __init__(self, id, direct_heat_transfer_kW_m3, num_of_heights, height_m, area_m2, packing_area_m2_m3, corrugation_angle_degree, voidspace_m3_m3, liquid_holdup_m3_m3):
        super().__init__()
        self.num_of_samples = None
        self.direct_heat_transfer_kW_m3 = direct_heat_transfer_kW_m3
        self.num_of_heights = num_of_heights
        self.info = {}
        self.height_m = height_m
        self.area_m2 = area_m2
        self.packing_area_m2_m3 = packing_area_m2_m3
        self.corrugation_angle_degree = corrugation_angle_degree
        self.voidspace_m3_m3 = voidspace_m3_m3
        self.theta = self.corrugation_angle_degree
        self.liquid_holdup_m3_m3 = liquid_holdup_m3_m3
        self.print_info = False
        self.id = id
        self.R = 8.314 * 10 ** (-2)
        self.g = 9.81
        self.absorption = {}
        self.z = np.linspace(0, self.height_m, self.num_of_heights)

    def add_absorption(self, id, gas_reactant, liq_products, absorption_rate_kmol_m3s):
        self.absorption[id] = {}
        self.absorption[id]["Gas Reactant"] = gas_reactant
        self.absorption[id]["Liquid Products"] = liq_products
        self.absorption[id]["Absorption Rate [kmol/m3*s]"] = absorption_rate_kmol_m3s
        self.absorption[id]["Stoch Ratio Gas"] = {}
        self.absorption[id]["Stoch Ratio Liquid"] = {}
        self.absorption[id]["Stoch Ratio Gas"][gas_reactant] = - 1.0
        for key in liq_products.keys():
            self.absorption[id]["Stoch Ratio Liquid"][key] = liq_products[key]

    def react(self, GasStreamIn, LiquidStreamIn):

        self.num_of_samples = len(LiquidStreamIn.temp_K)
        shape = np.ones(shape=(self.num_of_heights,))
        # ------------------------------------------------------------
        self.LiquidStreamIn = deepcopy(LiquidStreamIn)
        self.GasStreamIn = deepcopy(GasStreamIn)
        # ------------------------------------------------------------
        self.LiquidStreamIn.v_sup = - np.abs((self.LiquidStreamIn.flow_kg_h / (3600 * self.area_m2 * self.LiquidStreamIn.get_solution_density_kg_m3())))
        self.GasStreamIn.v_sup = (self.GasStreamIn.flow_kmol_h * self.R * self.GasStreamIn.temp_K) / (self.GasStreamIn.pressure_bara * self.area_m2 * 3600)
        # ------------------------------------------------------------
        self.LiquidStreamOut = deepcopy(self.LiquidStreamIn)
        self.GasStreamOut = deepcopy(self.GasStreamIn)

        # ------------------------------------------------------------
        for self.i in range(self.num_of_samples):
            # ------------------------------------------------------------
            y = []
            self.LiquidStream = deepcopy(LiquidStreamIn)
            for key in LiquidStreamIn.specie.keys():
                self.LiquidStream.specie[key]["Mass Fraction"] = LiquidStreamIn.specie[key]["Mass Fraction"][self.i] * shape
                y.append(self.LiquidStream.specie[key]["Mass Fraction"])
            self.LiquidStream.temp_K = LiquidStreamIn.temp_K[self.i] * shape
            self.LiquidStream.flow_kg_h = LiquidStreamIn.flow_kg_h[self.i] * shape
            self.LiquidStream.v_sup = self.LiquidStreamIn.v_sup[self.i] * shape
            y.append(self.LiquidStream.temp_K)
            y.append(self.LiquidStream.v_sup)
            # ------------------------------------------------------------
            self.GasStream = deepcopy(GasStreamIn)
            for key in GasStreamIn.specie.keys():
                self.GasStream.specie[key]["Molar Fraction"] = GasStreamIn.specie[key]["Molar Fraction"][self.i] * shape
                y.append(self.GasStream.specie[key]["Molar Fraction"])
            self.GasStream.pressure_bara = GasStreamIn.pressure_bara[self.i] * shape
            self.GasStream.temp_K = GasStreamIn.temp_K[self.i] * shape
            self.GasStream.flow_kmol_h = GasStreamIn.flow_kmol_h[self.i] * shape
            self.GasStream.v_sup = self.GasStreamIn.v_sup[self.i] * shape
            self.GasStream.pressure_bara = GasStreamIn.pressure_bara[self.i] * shape
            y.append(self.GasStream.temp_K)
            y.append(self.GasStream.v_sup)
            # ------------------------------------------------------------
            y = np.vstack(y)
            problem = solve_bvp(fun=self.__system_dynamics__, bc=self.__boundary_conditions__, x=self.z, y=y)
            # ------------------------------------------------------------
            NL = len(self.LiquidStreamIn.specie)
            NG = len(self.GasStreamIn.specie)
            x = np.array([0])
            for i, key in enumerate(self.LiquidStreamOut.specie.keys()):
                self.LiquidStreamOut.specie[key]["Mass Fraction"][self.i] = problem.sol(x)[i, 0]
            self.LiquidStreamOut.temp_K[self.i] = problem.sol(x)[NL, 0]
            self.LiquidStreamOut.v_sup[self.i] = problem.sol(x)[NL + 1, 0]
            self.LiquidStreamOut.flow_kg_h[self.i] = 3600 * np.abs(self.LiquidStreamOut.v_sup[self.i]) * self.area_m2 * self.LiquidStreamOut.get_solution_density_kg_m3()
            # ------------------------------------------------------------
            x = np.array([self.height_m])
            for i, key in enumerate(self.GasStreamOut.specie.keys()):
                self.GasStreamOut.specie[key]["Molar Fraction"][self.i] = problem.sol(x)[NL + 2 + i, 0]
            self.GasStreamOut.temp_K[self.i] = problem.sol(x)[NL + 2 + NG, 0]
            self.GasStreamOut.v_sup[self.i] = problem.sol(x)[NL + 2 + NG + 1, 0]
            self.GasStreamOut.flow_kmol_h[self.i] = 3600 * self.GasStreamOut.pressure_bara[self.i] * self.GasStreamOut.v_sup[self.i] * self.area_m2 / (self.R * self.GasStreamOut.temp_K[self.i])
            # ------------------------------------------------------------
        return self.GasStreamOut, self.LiquidStreamOut

    def __get_absorption_rate_kmol_m3s__(self):
        absorption_rate_kmol_m3s = {}
        absorption_rate_tot_kmol_m3s = 0
        absorption_rate_tot_kg_m3s = 0
        for id in self.GasStream.specie.keys():
            absorption_rate_kmol_m3s[id] = 0
        for a in self.absorption.keys():
            id = self.absorption[a]["Gas Reactant"]
            absorption_rate_kmol_m3s[id] = self.absorption[a]["Absorption Rate [kmol/m3*s]"](self)
            absorption_rate_tot_kmol_m3s = absorption_rate_tot_kmol_m3s + absorption_rate_kmol_m3s[id]
            absorption_rate_tot_kg_m3s = absorption_rate_tot_kg_m3s + absorption_rate_kmol_m3s[id] * self.GasStream.specie[id]["Molar Mass [kg/kmol]"]
        return absorption_rate_kmol_m3s, absorption_rate_tot_kmol_m3s, absorption_rate_tot_kg_m3s

    def __get_q_rxn_gas_kW_m3__(self):
        q_rxn_gas = 0
        n_rnx_gas = len(self.GasStream.reaction_rate_kmol_m3s)
        c_gain_gas = 0
        for i in range(n_rnx_gas):
            rate = self.GasStream.reaction_rate_kmol_m3s[i](self.GasStream)
            enthalpy = self.GasStream.reaction_enthalpy_kJ_kmol[i](self.GasStream)
            q_rxn_gas = q_rxn_gas + rate * enthalpy  # kW/m3 = kmol/m3.s * kJ/kmol
        return q_rxn_gas

    def __get_q_rxn_liq_kW_m3__(self):
        q_rxn_liq = 0
        n_rnx_liq = len(self.LiquidStream.reaction_rate_kmol_m3s)
        for i in range(n_rnx_liq):
            rate = self.LiquidStream.reaction_rate_kmol_m3s[i](self.LiquidStream)
            enthalpy = self.LiquidStream.reaction_enthalpy_kJ_kmol[i](self.LiquidStream)
            q_rxn_liq = q_rxn_liq + rate * enthalpy
        return q_rxn_liq

    def __get_q_lat_kW_m3__(self, rate):
        q_lat = 0
        q_abs = 0
        q_des = 0
        dT = self.GasStream.get_gas_temp_K() - self.LiquidStream.get_solution_temp_K()
        for a in self.absorption.keys():
            id = self.absorption[a]["Gas Reactant"]
            cp = self.GasStream.get_specie_heat_capacity_kJ_kmolK(id=self.absorption[a]["Gas Reactant"])
            q_des = q_des + np.minimum(rate[id], 0) * dT * cp
            q_abs = q_abs + np.maximum(rate[id], 0) * dT * cp
            enthalpy = self.absorption[a]["Enthalpy [kJ/kmol]"](self.LiquidStream)
            q_lat = q_lat + rate[id] * enthalpy  # kW/m3 = kmol/m3.s * kJ/kmol
        return q_lat, q_abs, q_des

    def __get_c_gas_gain_kmol_m3s__(self):
        n_rnx_gas = len(self.GasStream.reaction_rate_kmol_m3s)
        c_gain_gas = 0
        for i in range(n_rnx_gas):
            # Keeping track to whether the total moles increase or decrease
            rate = self.GasStream.reaction_rate_kmol_m3s[i](self.GasStream)
            stoch = 0
            for id in self.GasStream.reaction_reactants[i].keys():
                stoch = stoch - self.GasStream.reaction_reactants[i][id]
            for id in self.GasStream.reaction_products[i].keys():
                stoch = stoch + self.GasStream.reaction_products[i][id]
            c_gain_gas = c_gain_gas + rate * stoch
        return c_gain_gas

    def __system_dynamics__(self, x, y):

        N = y.shape[1]
        NL = len(self.LiquidStreamIn.specie)
        NG = len(self.GasStreamIn.specie)

        # Unpack
        #shape = np.ones(shape=(N,))
        self.LiquidStream = deepcopy(self.LiquidStreamIn)
        for i, key in enumerate(self.LiquidStreamIn.specie.keys()):
            self.LiquidStream.specie[key]["Mass Fraction"] = np.maximum(y[i, :], 10**(-18))
        self.LiquidStream.temp_K = np.clip(a=y[NL, :], a_min=273.15, a_max=150 + 273.15)
        self.LiquidStream.v_sup = np.minimum(y[NL + 1, :], -10**(-9))
        self.LiquidStream.flow_kg_h = np.zeros(shape=(N,))

        self.GasStream = deepcopy(self.GasStreamIn)
        for i, key in enumerate(self.GasStreamIn.specie.keys()):
            self.GasStream.specie[key]["Molar Fraction"] = np.maximum(y[i + NL + 2, :], 10**(-18))
        self.GasStream.temp_K = np.clip(a=y[NL + 2 + NG, :], a_min=273.15, a_max=150 + 273.15)
        self.GasStream.v_sup = np.maximum(y[NL + 2 + NG + 1, :], 10**(-9))
        self.GasStream.pressure_bara = self.GasStreamIn.pressure_bara[self.i] * np.ones(shape=(N,))

        # --------------------------------------------------------------------------------------
        self.rho_gas = self.GasStream.get_gas_density_kg_m3()
        self.mu_gas = self.GasStream.get_gas_viscosity_Pas()
        self.T_gas = self.GasStream.get_gas_temp_K()
        self.p_gas = self.GasStream.get_gas_pressure_bara()
        # --------------------------------------------------------------------------------------
        self.rho_liq = self.LiquidStream.get_solution_density_kg_m3() * np.ones(shape=(N,))
        self.T_liq = self.LiquidStream.get_solution_temp_K()
        # --------------------------------------------------------------------------------------

        # Packing Hydraulics
        # Absorption / Desorption
        # Direct Heat Transfer and Heat due to Chemical Reactions and Absorption / Desorption
        # Overall increase in moles in gas phase due to (non-instantaneous) chemical reactions
        self.dp_gas = 0
        self.h_liq = self.liquid_holdup_m3_m3(self)
        absorption_rate_kmol_m3s, absorption_rate_tot_kmol_m3s, absorption_rate_tot_kg_m3s = self.__get_absorption_rate_kmol_m3s__()
        q_rxn_liq = self.__get_q_rxn_liq_kW_m3__()
        q_rxn_gas = self.__get_q_rxn_gas_kW_m3__()
        q_lat, q_abs, q_des = self.__get_q_lat_kW_m3__(rate=absorption_rate_kmol_m3s)
        q_dir = self.direct_heat_transfer_kW_m3(self)
        c_gain_gas = self.__get_c_gas_gain_kmol_m3s__()

        # Liquid Density
        drho_liq = np.zeros(shape=self.rho_liq.shape)
        #drho_liq[0] = (self.rho_liq[1] - self.rho_liq[0]) / (self.z[1] - self.z[0])
        #drho_liq[N - 1] = (self.rho_liq[N - 1] - self.rho_liq[N - 2]) / (self.z[N - 1] - self.z[N - 2])
        #drho_liq[1:N - 1:] = (self.rho_liq[2:N:] - self.rho_liq[0:N - 2:]) / (self.z[2:N:] - self.z[0:N - 2:])

        # Gas Temperature
        cp_gas = self.GasStream.get_gas_heat_capacity_kJ_kmolK()
        dT_gas = (self.R * self.T_gas / (self.p_gas * cp_gas * self.GasStream.v_sup)) * (-q_dir + q_rxn_gas + q_des)

        # Liquid Temperature
        cp_liq = self.LiquidStream.get_solution_heat_capacity_kJ_kgK()
        dT_liq = (1 / self.LiquidStream.v_sup) * (1 / (self.rho_liq * cp_liq)) * (q_dir + q_rxn_liq + q_lat + q_abs)

        # Superficial Velocities
        dv_gas_sup = self.GasStream.v_sup * (dT_gas / self.T_gas - self.dp_gas / self.p_gas) - self.R * self.T_gas * (absorption_rate_tot_kmol_m3s - self.voidspace_m3_m3 * (1 - self.h_liq) * c_gain_gas) / self.p_gas
        dv_liq_sup = - self.LiquidStream.v_sup * drho_liq / self.rho_liq + (1 / self.rho_liq) * absorption_rate_tot_kg_m3s

        # Molar Fraction in Gas Phase
        dy = {}
        for id in self.GasStream.specie.keys():
            _y_ = self.GasStream.specie[id]["Molar Fraction"]
            dy[id] = _y_ * (dT_gas / self.T_gas - self.dp_gas / self.p_gas - dv_gas_sup / self.GasStream.v_sup)
            dy[id] = dy[id] - (self.R * self.T_gas * absorption_rate_kmol_m3s[id]) / (self.p_gas * self.GasStream.v_sup)
        n_rnx_gas = len(self.GasStream.reaction_rate_kmol_m3s)

        for i in range(n_rnx_gas):
            rate = self.GasStream.reaction_rate_kmol_m3s[i](self.GasStream)
            for id in self.GasStream.reaction_reactants[i].keys():
                stoch = - self.GasStream.reaction_reactants[i][id]
                dy[id] = dy[id] + self.voidspace_m3_m3 * (1 - self.h_liq)(self.R * self.T_gas * rate * stoch) / (self.p_gas * self.GasStream.v_sup)
            for id in self.GasStream.reaction_products[i].keys():
                stoch = self.GasStream.reaction_products[i][id]
                dy[id] = dy[id] + self.voidspace_m3_m3 * (1 - self.h_liq) * (self.R * self.T_gas * rate * stoch) / (self.p_gas * self.GasStream.v_sup)

        # Mass Fraction in Liquid Phase
        dw = {}
        for id in self.LiquidStream.specie.keys():
            w = self.LiquidStream.specie[id]["Mass Fraction"]
            dw[id] = w * (- dv_liq_sup / self.LiquidStream.v_sup - drho_liq / self.rho_liq)
        for a in self.absorption.keys():
            id_absorbent = self.absorption[a]["Gas Reactant"]
            for id in self.absorption[a]["Liquid Products"].keys():
                M = self.LiquidStream.specie[id]["Molar Mass [kg/kmol]"]
                stoch = self.absorption[a]["Stoch Ratio Liquid"][id]
                dw[id] = dw[id] + M * absorption_rate_kmol_m3s[id_absorbent] * stoch / (self.rho_liq * self.LiquidStream.v_sup)

        n_rnx_liq = len(self.LiquidStream.reaction_rate_kmol_m3s)
        for i in range(n_rnx_liq):
            rate = self.LiquidStream.reaction_rate_kmol_m3s[i](self.LiquidStream)
            for id in self.LiquidStream.reaction_reactants[i].keys():
                M = self.LiquidStream.specie[id]["Molar Mass [kg/kmol]"]
                stoch = - self.LiquidStream.reaction_reactants[i][id]
                dw[id] = dw[id] + self.voidspace_m3_m3 * self.h_liq * M * rate * stoch / (self.rho_liq * self.LiquidStream.v_sup)
            for id in self.LiquidStream.reaction_products[i].keys():
                M = self.LiquidStream.specie[id]["Molar Mass [kg/kmol]"]
                stoch = self.LiquidStream.reaction_products[i][id]
                dw[id] = dw[id] + self.voidspace_m3_m3 * self.h_liq * M * rate * stoch / (self.rho_liq * self.LiquidStream.v_sup)

        # Return Result
        dydx = np.zeros(shape=y.shape)
        for i, id in enumerate(self.LiquidStream.specie.keys()):
            dydx[i, :] = dw[id]
        dydx[NL, :] = dT_liq
        dydx[NL + 1, :] = dv_liq_sup
        for i, id in enumerate(self.GasStream.specie.keys()):
            dydx[i + NL + 2, :] = dy[id]
        dydx[NL + 2 + NG, :] = dT_gas
        dydx[NL + 2 + NG + 1, :] = dv_gas_sup
        return dydx

    def __boundary_conditions__(self, ya, yb):
        NL = len(self.LiquidStreamIn.specie)
        NG = len(self.GasStreamIn.specie)
        con = np.zeros(shape=(NL + NG + 4,))
        # --------------------------------------------------------------------------------------
        for i, key in enumerate(self.LiquidStreamIn.specie.keys()):
            con[i] = yb[i] - self.LiquidStreamIn.specie[key]["Mass Fraction"][self.i]
        con[NL] = yb[NL] - self.LiquidStreamIn.temp_K[self.i]
        con[NL + 1] = yb[NL + 1] - self.LiquidStreamIn.v_sup[self.i]
        # --------------------------------------------------------------------------------------
        for i, key in enumerate(self.GasStreamIn.specie.keys()):
            con[i + NL + 2] = ya[i + NL + 2] - self.GasStreamIn.specie[key]["Molar Fraction"][self.i]
        con[NL + 2 + NG] = ya[NL + 2 + NG] - self.GasStreamIn.temp_K[self.i]
        con[NL + 2 + NG + 1] = ya[NL + 2 + NG + 1] - self.GasStreamIn.v_sup[self.i]
        return con


class OpenSpray_CoFlow(FeatureBooster):

    def __init__(self):
        super().__init__()


class OpenSpray_CounterFlow(FeatureBooster):

    def __init__(self):
        super().__init__()


class RPB_CoFlow(FeatureBooster):

    def __init__(self):
        super().__init__()


class RPB_CounterFlow(FeatureBooster):

    def __init__(self):
        super().__init__()





# ---------------------------------------------------------------------------------------


class Calculator():

    def __init__(self):
        pass

    def get_equilibrium_constant_from_enthalpy_gibbs(self, stoch_gas, stoch_liq, stoch_sol):
        reactH = 0
        reactG = 0

        prodH = 0
        prodG = 0

        for id in stoch_gas.keys():
            if stoch_gas[id] < 0:
                reactH = reactH + 1000 * self.specie_gas[id]["Enthalpy of Formation [kJ/mol]"] * np.abs(stoch_gas[id])
                reactG = reactG + 1000 * self.specie_gas[id]["Gibbs free Energy of Formation [kJ/mol]"] * np.abs(stoch_gas[id])
            if stoch_gas[id] > 0:
                prodH = prodH + 1000 * self.specie_gas[id]["Enthalpy of Formation [kJ/mol]"] * np.abs(stoch_gas[id])
                prodG = prodG + 1000 * self.specie_gas[id]["Gibbs free Energy of Formation [kJ/mol]"] * np.abs(stoch_gas[id])

        for id in stoch_liq.keys():
            if stoch_liq[id] < 0:
                reactH = reactH + 1000 * self.specie_liq[id]["Enthalpy of Formation [kJ/mol]"] * np.abs(stoch_liq[id])
                reactG = reactG + 1000 * self.specie_liq[id]["Gibbs free Energy of Formation [kJ/mol]"] * np.abs(stoch_liq[id])
            if stoch_liq[id] > 0:
                prodH = prodH + 1000 * self.specie_liq[id]["Enthalpy of Formation [kJ/mol]"] * np.abs(stoch_liq[id])
                prodG = prodG + 1000 * self.specie_liq[id]["Gibbs free Energy of Formation [kJ/mol]"] * np.abs(stoch_liq[id])

        for id in stoch_sol.keys():
            if stoch_sol[id] < 0:
                reactH = reactH + 1000 * self.specie_sol[id]["Enthalpy of Formation [kJ/mol]"] * np.abs(stoch_sol[id])
                reactG = reactG + 1000 * self.specie_sol[id]["Gibbs free Energy of Formation [kJ/mol]"] * np.abs(stoch_sol[id])
            if stoch_sol[id] > 0:
                prodH = prodH + 1000 * self.specie_sol[id]["Enthalpy of Formation [kJ/mol]"] * np.abs(stoch_sol[id])
                prodG = prodG + 1000 * self.specie_sol[id]["Gibbs free Energy of Formation [kJ/mol]"] * np.abs(stoch_sol[id])


        dH = prodH - reactH
        dG = prodG - reactG
        K298 = np.exp(-dG / (8.314 * 298.15))
        TempCoeff = - dH / 8.314
        rxnHeat = reactH - prodH

        result = "Gas Phase \t" + str(stoch_gas) + "\n"
        result = result + "Liq Phase \t" + str(stoch_liq) + "\n"
        result = result + "K(T) \t\t= " + str(K298) + " * exp(" + str(TempCoeff) + " *(1/T - 1/298) ) \n"
        result = result + "pKa(298K) \t= " + str(-np.log10(K298)) + "\n"
        result = result + "q \t\t\t=" + " " + str(rxnHeat) + " kJ/kmol exothermic heat" +  "\n"
        return result

    def get_heat_of_vaporization_kJ_kmol(self, LiquidStreamIn, gas_id):
        dT = 0.5
        LiquidStreamOut1 = self.eq.react(LiquidStreamIn)
        T1 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K + dT
        LiquidStreamOut2 = self.eq.react(LiquidStreamIn)
        T2 = LiquidStreamIn.temp_K
        LiquidStreamIn.temp_K = LiquidStreamIn.temp_K - dT
        p1 = LiquidStreamOut1.get_specie_vapor_pressure_bara(gas_id=gas_id)
        p2 = LiquidStreamOut2.get_specie_vapor_pressure_bara(gas_id=gas_id)
        q = - 8.314 * np.log(p2 / p1) * (1/T2 - 1/T1) ** (-1)
        return q





