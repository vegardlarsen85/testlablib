import numpy as np
import testlablib as lab

library = lab.Library()

library.add_GasStream_specie_gas(id="CH4", molar_mass_kg_kmol=16, charge=0, enthalpy_of_formation_kJ_mol=-74.81, entropy_J_molK=186.264, gibbs_free_energy_of_formation_kJ_mol=-50.752)
library.add_GasStream_specie_gas(id="CO2", molar_mass_kg_kmol=44, charge=0, enthalpy_of_formation_kJ_mol=-393.509, entropy_J_molK=213.74, gibbs_free_energy_of_formation_kJ_mol=-394.359)
library.add_GasStream_specie_gas(id="H2", molar_mass_kg_kmol=2, charge=0, enthalpy_of_formation_kJ_mol=0, entropy_J_molK=130.684, gibbs_free_energy_of_formation_kJ_mol=0)
library.add_GasStream_specie_gas(id="H2O", molar_mass_kg_kmol=18, charge=0, enthalpy_of_formation_kJ_mol=-241.818, entropy_J_molK=188.25, gibbs_free_energy_of_formation_kJ_mol=-228.572)
library.add_GasStream_specie_gas(id="HNO2", molar_mass_kg_kmol=47, charge=0, enthalpy_of_formation_kJ_mol=-78.6592, entropy_J_molK=249.1572, gibbs_free_energy_of_formation_kJ_mol=-43.932)
library.add_GasStream_specie_gas(id="HNO3", molar_mass_kg_kmol=63, charge=0, enthalpy_of_formation_kJ_mol=-135.06, entropy_J_molK=266.38, gibbs_free_energy_of_formation_kJ_mol=-74.72)
library.add_GasStream_specie_gas(id="MDEA", molar_mass_kg_kmol=119, charge=0, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_GasStream_specie_gas(id="MEA", molar_mass_kg_kmol=61, charge=0, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_GasStream_specie_gas(id="N2", molar_mass_kg_kmol=28, charge=0, enthalpy_of_formation_kJ_mol=0, entropy_J_molK=191.61, gibbs_free_energy_of_formation_kJ_mol=0)
library.add_GasStream_specie_gas(id="NO", molar_mass_kg_kmol=30, charge=0, enthalpy_of_formation_kJ_mol=90.25, entropy_J_molK=210.76, gibbs_free_energy_of_formation_kJ_mol=86.55)
library.add_GasStream_specie_gas(id="NO2", molar_mass_kg_kmol=46, charge=0, enthalpy_of_formation_kJ_mol=33.18, entropy_J_molK=240.06, gibbs_free_energy_of_formation_kJ_mol=51.31)
library.add_GasStream_specie_gas(id="NO3", molar_mass_kg_kmol=62, charge=0, enthalpy_of_formation_kJ_mol=70.9188, entropy_J_molK=252.54624, gibbs_free_energy_of_formation_kJ_mol=114.47424)
library.add_GasStream_specie_gas(id="N2O4", molar_mass_kg_kmol=92, charge=0, enthalpy_of_formation_kJ_mol=9.16, entropy_J_molK=304.29, gibbs_free_energy_of_formation_kJ_mol=97.89)
library.add_GasStream_specie_gas(id="N2O3", molar_mass_kg_kmol=76, charge=0, enthalpy_of_formation_kJ_mol=83.72184, entropy_J_molK=312.16824, gibbs_free_energy_of_formation_kJ_mol=139.41088)
library.add_GasStream_specie_gas(id="N2O2", molar_mass_kg_kmol=60, charge=0, enthalpy_of_formation_kJ_mol=170.37248, entropy_J_molK=287.52448, gibbs_free_energy_of_formation_kJ_mol=202.88216)
library.add_GasStream_specie_gas(id="N2O", molar_mass_kg_kmol=44, charge=0, enthalpy_of_formation_kJ_mol=82.05, entropy_J_molK=219.85, gibbs_free_energy_of_formation_kJ_mol=104.2)
library.add_GasStream_specie_gas(id="NH3", molar_mass_kg_kmol=17, charge=0, enthalpy_of_formation_kJ_mol=-46.11, entropy_J_molK=192.45, gibbs_free_energy_of_formation_kJ_mol=-16.45)
library.add_GasStream_specie_gas(id="O2", molar_mass_kg_kmol=32, charge=0, enthalpy_of_formation_kJ_mol=0, entropy_J_molK=205.138, gibbs_free_energy_of_formation_kJ_mol=0)
library.add_GasStream_specie_gas(id="O3", molar_mass_kg_kmol=48, charge=0, enthalpy_of_formation_kJ_mol=142.7, entropy_J_molK=238.93, gibbs_free_energy_of_formation_kJ_mol=163.2)
library.add_GasStream_specie_gas(id="PZ", molar_mass_kg_kmol=86, charge=0, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_GasStream_specie_gas(id="SO2", molar_mass_kg_kmol=64, charge=0, enthalpy_of_formation_kJ_mol=-296.83, entropy_J_molK=248.21, gibbs_free_energy_of_formation_kJ_mol=-300.09)
library.add_GasStream_specie_gas(id="SO3", molar_mass_kg_kmol=80, charge=0, enthalpy_of_formation_kJ_mol=-395.72, entropy_J_molK=256.76, gibbs_free_energy_of_formation_kJ_mol=-371.06)


