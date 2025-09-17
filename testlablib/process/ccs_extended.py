import testlablib as lab
from testlablib.library.library import library




class ExhaustGas(lab.GasStream):

    def __init__(self, flow_Nm3_h_dry, pressure_bara, temp_K, SO2_ppm_dry, CO2_pct_dry, NO_ppm_dry, NO2_ppm_dry, O2_pct_dry, H2O_pct):

        super().__init__()

        self.load_viscosity(id="Exhaust Gas", library=gases)
        self.load_diffusivity_m2_s(id="Exhaust Gas", library=gases)
        self.load_thermal_conductivity_kW_mK(id="Exhaust Gas", library=gases)
        self.load_heat_capacity_kJ_kmolK(id="Exhaust Gas", library=gases)

        for id in species.specie_gas.keys():
            self.add_specie(id=id, library=species)
        #self.add_specie(id="CO2", library=species)
        #self.add_specie(id="SO2", library=species)
        #self.add_specie(id="NO", library=species)
        #self.add_specie(id="NO2", library=species)
        #self.add_specie(id="O2", library=species)
        #self.add_specie(id="H2O", library=species)
        #self.add_specie(id="N2", library=species)
        #self.add_specie(id="N2O2", library=species)
        #self.add_specie(id="N2O3", library=species)
        #self.add_specie(id="N2O4", library=species)
        #self.add_specie(id="HNO2", library=species)

        self.add_equilibrium_gas(id="O2 + N2 = 2NO", library=reactions)
        self.add_equilibrium_gas(id="2NO + O2 = 2NO2", library=reactions)
        self.add_equilibrium_gas(id="2NO2 = N2O4", library=reactions)
        self.add_equilibrium_gas(id="NO2 + NO = N2O3", library=reactions)
        self.add_equilibrium_gas(id="2NO = N2O2", library=reactions)
        self.add_equilibrium_gas(id="N2O3 + H2O = 2HNO2", library=reactions)

        self.set_gas_temp_K(value=temp_K)
        self.set_gas_pressure_bara(value=pressure_bara)

        shape = np.ones(shape=temp_K.shape)
        for id in self.specie.keys():
            self.set_specie_molar_fraction(id=id, value=0*shape)

        y_SO2 = SO2_ppm_dry * 10**(-6)
        y_CO2 = CO2_pct_dry / 100
        y_NO = NO_ppm_dry * 10**(-6)
        y_NO2 = NO2_ppm_dry * 10**(-6)
        y_O2 = O2_pct_dry / 100
        y_N2 = 1 - y_SO2 - y_CO2 - y_NO - y_NO2 - y_O2

        self.set_specie_molar_fraction(id="SO2", value=y_SO2)
        self.set_specie_molar_fraction(id="CO2", value=y_CO2)
        self.set_specie_molar_fraction(id="NO", value=y_NO)
        self.set_specie_molar_fraction(id="NO2", value=y_NO2)
        self.set_specie_molar_fraction(id="O2", value=y_O2)
        self.set_specie_molar_fraction(id="N2", value=y_N2)

        if H2O_pct is None:
            p_H2O = self.H2O_vapor_pressure_bara(temp_K=temp_K)
            y_H2O = p_H2O/ pressure_bara
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






