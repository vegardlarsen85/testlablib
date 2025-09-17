import numpy as np
import testlablib as lab



library = lab.Library()




"""""""""

Species

The specie K'+ is Potassium originating (same as K+), but originating from K2CO3. Separation of the two is necessary to be able to define CO2 Load.

AMP
Ci      Citric Acid
Gly     Glycine
Lys     Lysine
MDEA    xxx
MEA     Monoethanolamine
PZ      Piperazine
Sar     Sarcosine

"""""""""



library.add_LiquidStream_specie_liq(id="AMP", molar_mass_kg_kmol=89, charge=0, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="AMPCOO-", molar_mass_kg_kmol=132, charge=-1, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="AMPH+", molar_mass_kg_kmol=90, charge=1, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Ca+2", molar_mass_kg_kmol=40, charge=2, enthalpy_of_formation_kJ_mol=-542.96, entropy_J_molK=-55.2, gibbs_free_energy_of_formation_kJ_mol=-553.04)
library.add_LiquidStream_specie_liq(id="Ci", molar_mass_kg_kmol=192, charge=0, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Ci-", molar_mass_kg_kmol=191, charge=-1, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Ci-2", molar_mass_kg_kmol=190, charge=-2, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Ci-3", molar_mass_kg_kmol=189, charge=-3, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Cl-", molar_mass_kg_kmol=35.5, charge=-1, enthalpy_of_formation_kJ_mol=-167.16, entropy_J_molK=56.5, gibbs_free_energy_of_formation_kJ_mol=-131.26)
library.add_LiquidStream_specie_liq(id="CO2", molar_mass_kg_kmol=44, charge=0, enthalpy_of_formation_kJ_mol=-413.8, entropy_J_molK=117.6, gibbs_free_energy_of_formation_kJ_mol=-386.0)
library.add_LiquidStream_specie_liq(id="CO3-2", molar_mass_kg_kmol=60, charge=-2, enthalpy_of_formation_kJ_mol=-677.14, entropy_J_molK=-56.9, gibbs_free_energy_of_formation_kJ_mol=-527.81)
library.add_LiquidStream_specie_liq(id="Gly", molar_mass_kg_kmol=75, charge=0, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Gly-", molar_mass_kg_kmol=74, charge=-1, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Gly+", molar_mass_kg_kmol=76, charge=1, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="GlyCOO-2", molar_mass_kg_kmol=117, charge=-2, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="H2O", molar_mass_kg_kmol=18, charge=0, enthalpy_of_formation_kJ_mol=-285.830, entropy_J_molK=69.91, gibbs_free_energy_of_formation_kJ_mol=-237.129)
library.add_LiquidStream_specie_liq(id="H2O2", molar_mass_kg_kmol=34, charge=0, enthalpy_of_formation_kJ_mol=-191.17, entropy_J_molK=143.9, gibbs_free_energy_of_formation_kJ_mol=-134.03)
library.add_LiquidStream_specie_liq(id="H+", molar_mass_kg_kmol=1, charge=1, enthalpy_of_formation_kJ_mol=0, entropy_J_molK=0, gibbs_free_energy_of_formation_kJ_mol=0)
library.add_LiquidStream_specie_liq(id="HCO3-", molar_mass_kg_kmol=61, charge=-1, enthalpy_of_formation_kJ_mol=-691.99, entropy_J_molK=91.2, gibbs_free_energy_of_formation_kJ_mol=-586.77)
library.add_LiquidStream_specie_liq(id="HSO3-", molar_mass_kg_kmol=81, charge=-1, enthalpy_of_formation_kJ_mol=-627.41, entropy_J_molK=134.17, gibbs_free_energy_of_formation_kJ_mol=-527.14)
library.add_LiquidStream_specie_liq(id="HSO4-", molar_mass_kg_kmol=97, charge=-1, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="HNO2", molar_mass_kg_kmol=47, charge=0, enthalpy_of_formation_kJ_mol=-119.2, entropy_J_molK=135.6, gibbs_free_energy_of_formation_kJ_mol=-50.6)
library.add_LiquidStream_specie_liq(id="HNO3", molar_mass_kg_kmol=63, charge=0, enthalpy_of_formation_kJ_mol=-207.36, entropy_J_molK=146.4, gibbs_free_energy_of_formation_kJ_mol=-111.25)
library.add_LiquidStream_specie_liq(id="HOOCPZCOO-", molar_mass_kg_kmol=173, charge=-1, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)     # Zwitterion
library.add_LiquidStream_specie_liq(id="HPZCOO", molar_mass_kg_kmol=130, charge=0,enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)           # Zwitterion
library.add_LiquidStream_specie_liq(id="HMEACOO", molar_mass_kg_kmol=105, charge=0,enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)          # Zwitterion
library.add_LiquidStream_specie_liq(id="HAMPCOO", molar_mass_kg_kmol=133, charge=0,enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)          # Zwitterion
library.add_LiquidStream_specie_liq(id="K'+", molar_mass_kg_kmol=39, charge=1, enthalpy_of_formation_kJ_mol=-252.38, entropy_J_molK=102.5, gibbs_free_energy_of_formation_kJ_mol=-283.27)
library.add_LiquidStream_specie_liq(id="K+", molar_mass_kg_kmol=39, charge=1, enthalpy_of_formation_kJ_mol=-252.38, entropy_J_molK=102.5, gibbs_free_energy_of_formation_kJ_mol=-283.27)
library.add_LiquidStream_specie_liq(id="Lys", molar_mass_kg_kmol=146, charge=0, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Lys-", molar_mass_kg_kmol=145, charge=-1, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Lys+", molar_mass_kg_kmol=147, charge=1, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Lys+2", molar_mass_kg_kmol=148, charge=2, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="LysCOO-2", molar_mass_kg_kmol=188, charge=-2, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="LysCOO-", molar_mass_kg_kmol=189, charge=-1, enthalpy_of_formation_kJ_mol=None, entropy_J_molK=None, gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="MDEA", molar_mass_kg_kmol=119, charge=0, enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="MDEAH+",molar_mass_kg_kmol=120,charge=1,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="MEA",molar_mass_kg_kmol=61,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="MEAH+",molar_mass_kg_kmol=62,charge=1,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="MEACOO-",molar_mass_kg_kmol=104,charge=-1,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Mg+2",molar_mass_kg_kmol=24.3,charge=2,enthalpy_of_formation_kJ_mol=-466.85,entropy_J_molK=-138.1,gibbs_free_energy_of_formation_kJ_mol=-454.8)
library.add_LiquidStream_specie_liq(id="Na+",molar_mass_kg_kmol=23,charge=1,enthalpy_of_formation_kJ_mol=-240.12,entropy_J_molK=59.0,gibbs_free_energy_of_formation_kJ_mol=-261.91)
library.add_LiquidStream_specie_liq(id="N2",molar_mass_kg_kmol=28,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="NH3",molar_mass_kg_kmol=17,charge=0,enthalpy_of_formation_kJ_mol=-80.29,entropy_J_molK=111,gibbs_free_energy_of_formation_kJ_mol=-26.6)
library.add_LiquidStream_specie_liq(id="NH4+",molar_mass_kg_kmol=18,charge=0,enthalpy_of_formation_kJ_mol=-132.51,entropy_J_molK=113.4,gibbs_free_energy_of_formation_kJ_mol=-79.31)
library.add_LiquidStream_specie_liq(id="NO2",molar_mass_kg_kmol=46,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="N2O4",molar_mass_kg_kmol=92,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="NO",molar_mass_kg_kmol=30,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="N2O2",molar_mass_kg_kmol=60,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="N2O3",molar_mass_kg_kmol=76,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="NO2-",molar_mass_kg_kmol=46,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="O2",molar_mass_kg_kmol=32,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="NO3-",molar_mass_kg_kmol=62,charge=0,enthalpy_of_formation_kJ_mol=-205.0,entropy_J_molK=146.4,gibbs_free_energy_of_formation_kJ_mol=-108.74)
library.add_LiquidStream_specie_liq(id="OH-",molar_mass_kg_kmol=17,charge=-1,enthalpy_of_formation_kJ_mol=-229.994,entropy_J_molK=-10.75,gibbs_free_energy_of_formation_kJ_mol=-157.244)
library.add_LiquidStream_specie_liq(id="PZ",molar_mass_kg_kmol=86,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="PZH+",molar_mass_kg_kmol=87,charge=1,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="PZ(2H)+",molar_mass_kg_kmol=88,charge=2,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="PZCOO-",molar_mass_kg_kmol=129,charge=-1,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="PZ(COO)2-2",molar_mass_kg_kmol=172,charge=-2,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Sar",molar_mass_kg_kmol=89,charge=0,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Sar-",molar_mass_kg_kmol=88,charge=-1,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="Sar+",molar_mass_kg_kmol=90,charge=1,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="SarCOO-2",molar_mass_kg_kmol=131,charge=-2,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_liq(id="SO2",molar_mass_kg_kmol=64,charge=0,enthalpy_of_formation_kJ_mol=-323.78,entropy_J_molK=159.48,gibbs_free_energy_of_formation_kJ_mol=-300.60)
library.add_LiquidStream_specie_liq(id="SO3-2",molar_mass_kg_kmol=80,charge=-2,enthalpy_of_formation_kJ_mol=-631.06,entropy_J_molK=-15.4,gibbs_free_energy_of_formation_kJ_mol=-486.20)
library.add_LiquidStream_specie_liq(id="SO4-2",molar_mass_kg_kmol=96,charge=-2,enthalpy_of_formation_kJ_mol=-909.27,entropy_J_molK=20.1,gibbs_free_energy_of_formation_kJ_mol=-744.53)

library.add_LiquidStream_specie_sol(id="CaCO3",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-1206.92,entropy_J_molK=92.9,gibbs_free_energy_of_formation_kJ_mol=-1128.79)
library.add_LiquidStream_specie_sol(id="Ca(NO3)2",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-938.39,entropy_J_molK=193.3,gibbs_free_energy_of_formation_kJ_mol=-743.07)
library.add_LiquidStream_specie_sol(id="CaSO3",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_sol(id="CaSO4",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-1434.11,entropy_J_molK=106.7,gibbs_free_energy_of_formation_kJ_mol=-1321.79)
library.add_LiquidStream_specie_sol(id="MgCO3",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-1095.8,entropy_J_molK=65.7,gibbs_free_energy_of_formation_kJ_mol=-1012.1)
library.add_LiquidStream_specie_sol(id="MgSO4",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-1284.9,entropy_J_molK=91.6,gibbs_free_energy_of_formation_kJ_mol=-1170.6)
library.add_LiquidStream_specie_sol(id="Mg(OH)2",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-924.54,entropy_J_molK=63.18,gibbs_free_energy_of_formation_kJ_mol=-833.58)
library.add_LiquidStream_specie_sol(id="NaCl",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-411.153,entropy_J_molK=72.13,gibbs_free_energy_of_formation_kJ_mol=-384.138)
library.add_LiquidStream_specie_sol(id="Na2SO3",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_sol(id="Na2SO4",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-1387.08,entropy_J_molK=149.58,gibbs_free_energy_of_formation_kJ_mol=-1270.16)
library.add_LiquidStream_specie_sol(id="NaHSO3",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_sol(id="NaHSO4",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=None,entropy_J_molK=None,gibbs_free_energy_of_formation_kJ_mol=None)
library.add_LiquidStream_specie_sol(id="Na2CO3",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-1130.68,entropy_J_molK=134.98,gibbs_free_energy_of_formation_kJ_mol=-1044.44)
library.add_LiquidStream_specie_sol(id="NaHCO3",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-950.81,entropy_J_molK=101.7,gibbs_free_energy_of_formation_kJ_mol=-851.0)
library.add_LiquidStream_specie_sol(id="KNO3",molar_mass_kg_kmol=None,enthalpy_of_formation_kJ_mol=-494.63,entropy_J_molK=133.05,gibbs_free_energy_of_formation_kJ_mol=-394.86)


