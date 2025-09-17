import numpy as np
import testlablib as lab


library = lab.Library()


# ---------------------------------------------------------------------------------------------------

from testlablib.library.library_gases.thermo.heat_capacity import *
from testlablib.library.library_gases.thermo.viscosity import *
from testlablib.library.library_gases.thermo.thermal_conductivity import *
from testlablib.library.library_gases.thermo.diffusivity import *

library.add_GasStream_heat_capacity_kJ_kmolK(id="Exhaust Gas", function=exhaust_gas_heat_capacity_kJ_kmolK)
library.add_GasStream_viscosity_Pas(id="Exhaust Gas", function=exhaust_gas_viscosity_Pas)
library.add_GasStream_thermal_conductivity_kW_mK(id="Exhaust Gas", function=exhaust_gas_thermal_conductivity_kW_mK)
library.add_GasStream_diffusivity(id="Exhaust Gas", function=exhaust_gas_diffusivity_m2_s)


# ---------------------------------------------------------------------------------------------------


from testlablib.library.library_liquids.thermo.density import *
from testlablib.library.library_liquids.thermo.activity import *
from testlablib.library.library_liquids.thermo.heat_capacity import *
from testlablib.library.library_liquids.thermo.viscosity import *
from testlablib.library.library_liquids.thermo.diffusivity import *
from testlablib.library.library_liquids.thermo.surface_tension import *

library.add_LiquidStream_density_kg_m3(id="Water", function=density_H2O_kg_m3)
library.add_LiquidStream_density_kg_m3(id="Aqueous Solution w/Amines and Amino Acids", function=density_amino_kg_m3)
library.add_LiquidStream_activity_coefficient(id="Truesdell Jones", function=activity_coefficient_truesdell_jones)
library.add_LiquidStream_activity_coefficient(id="Ideal Solution", function=activity_coefficient_ideal)
library.add_LiquidStream_heat_capacity_kJ_kgK(id="Water", function=heat_capacity_H2O_kJ_kgK)
library.add_LiquidStream_heat_capacity_kJ_kgK(id="Aqueous Solution w/Amines and Amino Acids", function=heat_capacity_amino_kJ_kgK)
library.add_LiquidStream_viscosity_Pas(id="Water", function=viscosity_H2O_Pas)
library.add_LiquidStream_viscosity_Pas(id="Aqueous Solution w/Amines and Amino Acids", function=viscosity_amino_Pas)
library.add_LiquidStream_diffusivity_m2_s(id="All", function=diffusivity_m2_s)
library.add_LiquidStream_surface_tension_N_m(id="Water", function=surface_tension_H2O_N_m)


# ---------------------------------------------------------------------------------------------------

from testlablib.library.library_gases.species.species import library as species
library.append(species)

from testlablib.library.library_liquids.species.species import library as species
library.append(species)

# ---------------------------------------------------------------------------------------------------

from testlablib.library.library_gases.reactions.nitrogen import library as nitrogen
library.append(nitrogen)

# ---------------------------------------------------------------------------------------------------

from testlablib.library.library_liquids.reactions.amines import library as amines
from testlablib.library.library_liquids.reactions.amino_acids import library as amino_acids
from testlablib.library.library_liquids.reactions.ammonia import library as ammonia
from testlablib.library.library_liquids.reactions.ascorbic_acid import library as ascorbic_acid
from testlablib.library.library_liquids.reactions.carbonic_acid import library as carbonic_acid
from testlablib.library.library_liquids.reactions.citric_acid import library as citric_acid
from testlablib.library.library_liquids.reactions.nitrogen import library as nitrogen
from testlablib.library.library_liquids.reactions.sulfur import library as sulfur
from testlablib.library.library_liquids.reactions.vapor_pressure import library as vapor_pressure
from testlablib.library.library_liquids.reactions.water import library as water


library.append(amines)
library.append(amino_acids)
library.append(ammonia)
library.append(ascorbic_acid)
library.append(carbonic_acid)
library.append(citric_acid)
library.append(nitrogen)
library.append(sulfur)
library.append(vapor_pressure)
library.append(water)



#library.print_info()


