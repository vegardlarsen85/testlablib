# testlablib

## 1. About


testlablib is a Python library for simulating Physical and Chemical Processes.  

Some keywords are  
- Contimuous Stirred Tank Reactors (CSTR)
- Plug Flow Reactors (PFR)
- Flash Calculations
- Heat Exchangers
- Compressors
- Gas Liquid Contactors
- Equilibrium Concentrations
- Vapor Liquid Equilibrium
- Vapor Pressure
- Heat of Absorption
- Distillation Columns
- Reaction Kinetics  


Note, The reactors (CSTR, PFR) only calculate Steady State. The library are therefore not well suited for simulating Batch Processes. In addition, the Gases are Modelled using Ideal Gas Law, and are thus only accurate at relatively low pressures. However, the plan is to update both these issues eventually.  

Hope you find the library useful  

Best Regards  
Vegard Larsen  
vegardlarsen85@gmail.com  



## 2. Installation

```
pip install git+https://github.com/vegardlarsen85/testlablib.git
```


## 3. Abbriviations

| Description                                       | Symbol     | Unit                               |
|:--------------------------------------------------|:-----------|:-----------------------------------|
| *Molar Flow*                                      | $\dot{n}$  | $kmol \cdot h^{-1}$                |
| *Mass Flow*                                       | $\dot{m}$  | $kg \cdot h^{-1}$                  |
| *Volume Flow*                                     | $\dot{V}$  | $m^3 \cdot h^{-1}$                 |
| *Gas Phase Molar Fraction*                        | $y_i$      | $kmol \cdot kmol^{-1}$             |
| *Liquid Phase Molar Fraction*                     | $x_i$      | $kmol \cdot kmol^{-1}$             |
| *Mass Fraction*                                   | $w_i$      | $kmol \cdot kmol^{-1}$             |
| *Molarity*                                        | $c_i$      | $kmol \cdot m^{-3}$                |
| *Molality of solute i*                            | $m_i$      | $mol \cdot kg^{-1}$                |
| *Gas Pressure*                                    | $p$        | $bar(a)$                           |
| *Gas Partial Pressure of Specie i*                | $p_i$      | $bar(a)$                           |
| *Vapor Pressure of Specie i*                      | $p_i^*$    | $bar(a)$                           |
| *Vapor Pressure over a pure solution of Specie i* | $p_i^0$    | $bar(a)$                           |
| *Temperature*                                     | $T$        | $K$                                |
| *Henry's Law Coefficient of solute i*             | $H_i$      | $mol \cdot kg^{-1} \cdot bar^{-1}$ |
| *Density*                                         | $\rho$     | $kg \cdot m^{-3}$                  |
| *Heat Capacity of Gas*                            | $c_p$      | $kJ\cdot kmol^{-1} \cdot K^{-1}$   |
| *Heat Capacity of Gas Specie i*                   | $c_{p,i}$  | $kJ\cdot kmol^{-1} \cdot K^{-1}$   |
| *Heat Capacity of Liquid*                         | $c_p$      | $kJ\cdot kg^{-1} \cdot K^{-1}$     |
| *Molar Mass of Specie i*                          | $M_i$      | $kg \cdot kmol^{-1}$               |
| *Activity Coefficient of Specie i in Liquid*      | $\gamma_i$ | $1$                                |

## 4. Formula Cheat Sheets

Concentrations in Liquid Phase.

| Dilute Aqueous Solution                      | Arbitrary Solution                                                                                                 | Note       |
|:---------------------------------------------|:-------------------------------------------------------------------------------------------------------------------|:-----------|
| $c=\rho/18$                                  | $c=\rho \sum_i \frac{w_i}{M_i}$                                                                                    |            |
| $c=\rho/18$                                  | $c=\rho \frac{1}{\sum_k M_k x_k}$                                                                                  |            |
| $c=c_{solvent}$                              | $c=\sum_i c_i$                                                                                                     |            |
| $c_i=\rho \frac{w_i}{M_i}$                   | $c_i=\rho \frac{w_i}{M_i}$                                                                                         |            |
| $c_i=\rho \frac{x_i}{18}$                    | $c_i=\rho \frac{x_i}{\sum_k M_k x_k}$                                                                              |            |
| $c_i=\frac{1}{1000} \cdot m_i \cdot \rho$    | $c_i = \frac{x_{solvent}}{1000} \cdot m_i \cdot M_{solvent} \cdot c$                                               | *i=solute* |
| $m_i=1000 \left( \frac{w_i}{M_i} \right)$    | $m_i= \left( \frac{1000}{w_{solvent}} \right) \left( \frac{w_i}{M_i} \right) $                                     | *i=solute* |
| $m_i=1000 \left( \frac{x_i}{18} \right)$     | $m_i = \left( \frac{1000}{x_{solvent}} \right) \left( \frac{x_i}{M_{solvent}} \right) $                            | *i=solute* |
| $m_i = \left( \frac{1000}{\rho} \right) c_i$ | $m_i = \left( \frac{1000}{x_{solvent}} \right) \left( \frac{1}{M_{solvent}} \right) \left( \frac{c_i}{c} \right) $ | *i=solute* |
| $x_i = \left( \frac{18}{M_i} \right) w_i$    | $x_i = \left( \frac{1/M_i}{ \sum_k w_k/M_k} \right) w_i$                                                           |            |
| $x_i = \frac{c_i}{c}$                        | $x_i = \frac{c_i}{c}$                                                                                              |            |
| $x_{solvent} = 1$                            | $x_{solvent}=\frac{1000}{1000+M_{solvent} \sum_i m_i}$                                                             | *i=solute* |
| $x_{solvent} = 1$                            | $x_{solvent}=\frac{1000}{1000+M_{solvent} \sum_i m_i}$                                                             | *i=solute* |


---


Two Film Theory, Flux of Solutes (Henry's Law)


| Flux                      | Mass Transfer Coefficient                                                  | Note                                                   |
|:--------------------------|:---------------------------------------------------------------------------|:-------------------------------------------------------|
| $J_i=K_G (p_i - p_i^*)$   | $K_G = \left( \frac{RT}{k_G} + \frac{1}{H_i} \frac{1}{E k_L} \right)^{-1}$ |                                                        |
| $J_i=K_L (H_i p_i - c_i)$ | $K_L = \left( \frac{RTH_i}{k_G} + \frac{1}{E k_L} \right)^{-1}$            | where $c_i$ refer to liquid phase molarity of specie i |


---

Two Film Theory, Flux of Solvents (Raoult's Law)


| Flux                                    | Mass Transfer Coefficient                                                      | Note                                               |
|:----------------------------------------|:-------------------------------------------------------------------------------|:---------------------------------------------------|
| $J_i=K_G (p_i - p_i^*)$                 | $K_G = \left( \frac{RT}{k_G} + \frac{p^0_i}{c_L} \frac{1}{E k_L} \right)^{-1}$ |                                                    |
| $J_i=K_L c_L (\frac{p_i}{p_i^0} - x_i)$ | $K_L = \left( \frac{c_L}{p_i^0} \frac{RT}{k_G} + \frac{1}{E k_L} \right)^{-1}$ | where $c_L$ refer to overall liquid phase molarity |





## 5. Getting Started


### 5.1 LiquidStream

To get started, consider an simple example with Carbonic Acid.

First we import the necessary libraries.  

```python
import numpy as np
import matplotlib.pyplot as plt
import testlablib as lab
```

Next we create a LiquidStream Object which will eventually contain all the equations for configuring the solvent.
The argument *solvent_id* is used by the library to calculate molality of solutes.

```python
sol = lab.LiquidStream(solvent_id="H2O")
```

Next we add all the species in the liquid.  
Let's say our liquid mixture is a result of mixing K2CO3 and CO2 into water.
```python
sol.add_specie(id="CO2", molar_mass_kg_kmol=44, charge=0)
sol.add_specie(id="CO3-2", molar_mass_kg_kmol=60, charge=-2)
sol.add_specie(id="HCO3-", molar_mass_kg_kmol=61, charge=-1)
sol.add_specie(id="H2O", molar_mass_kg_kmol=18, charge=0)
sol.add_specie(id="H+", molar_mass_kg_kmol=1, charge=1)
sol.add_specie(id="OH-", molar_mass_kg_kmol=17, charge=-1)
sol.add_specie(id="K+", molar_mass_kg_kmol=39, charge=1)
```

We add functions defining Heat Capacity, Density and Activity Coefficients. The Heat Capacity and Density 
take a *LiquidStream* Object as argument, while Activity Coefficient take both a *LiquidStream* Object and an *id*
(string) as arguments, since the various species in the solution have varying activity coefficients.  

We also define the function *CO2Load* for convenience indicating the ratio CO2/K2CO3 added into the mixture.

```python
def CO2Load(LiquidStream):
    x_CO2 = 0
    x_Amine = 0
    C = {"CO2": 1, "HCO3-": 1, "CO3-2": 1, "K+": -0.5}
    A = {"K+": 0.5}
    for c in C.keys():
        x_CO2 = x_CO2 + C[c] * LiquidStream.get_specie_molar_fraction(id=c)
    for a in A.keys():
        x_Amine = x_Amine + A[a] * LiquidStream.get_specie_molar_fraction(id=a)
    alpha = x_CO2 / x_Amine
    return alpha

def heat_capacity_kJ_kgK(LiquidStream):
    return 4.2 * np.ones(shape=LiquidStream.get_solution_temp_K().shape)

def density_kg_m3(LiquidStream):
    return 1050 * np.ones(shape=LiquidStream.get_solution_temp_K().shape)

def activity_coefficient(LiquidStream, id):
    I = LiquidStream.get_solution_ionic_strength_mol_kg()
    z = LiquidStream.get_specie_charge(id)
    log10_gamma = - 0.51 * z ** 2 * np.sqrt(I) / (1 + 1.5 * np.sqrt(I))
    gamma = 10 ** log10_gamma
    return gamma

sol.load_heat_capacity_kJ_kgK(function=heat_capacity_kJ_kgK)
sol.load_density_kg_m3(function=density_kg_m3)
sol.load_activity_coefficient(function=activity_coefficient)
```

---
**Remark Regarding Density**  
Note that it's weight fractions that are stored as numerical values in a LiquidStream Object.  
When for example molarity are needed it is calculated from weight fraction using formulas shown in section 4.  
Consequently, in order to calculate molarity density are required.  
The density cannot therefore be set as a function of molarities as a catch 22 case would emerge.  
Instead weight fraction, molar fraction or molality should be used.  

**Example of a well defined density function**  
```python
def density_kg_m3(LiquidStream):
    m = sol.get_specie_molality_mol_kg(id="OH-")
    rho = 1000 + 0.1 * m
    return rho

sol.load_density_kg_m3(function=density_kg_m3)
```

**Example of a NOT well defined density function**  
```python
def density_kg_m3(LiquidStream):
    c = sol.get_specie_molarity_kmol_m3(id="OH-")
    rho = 1000 + 0.1 * c
    return rho

sol.load_density_kg_m3(function=density_kg_m3)
```

---

Next the equilibrium constants of all the chemical reactions in liquid phase is defined. The functions returning 
the equilibrium constants take a LiquidStream Object as argument. Adding the reactions to the Solvent is 
achieved using the *add_rxn_insta method*. The *id* argument must be an unique string. The *stoch* argument 
define the stochiometric coefficients of the reactions. It must be a negative number for reactants and positive 
for products. The argument *unit* define wich concentration is used. The options are the following  
- Molality = “m”  
- Molarity = “c”  
- Molar Fraction = “x”  
- Weight Fraction = “w”  
- None = None  

As an example, the constraint given by the water self-dissociation reaction is given by following equation.  
$K \cdot \gamma_{H_2O} x_{H_2O} = \gamma_{H^+} m_{H^+} \cdot \gamma_{OH^-} m_{OH^-} $  

Were,  
- $x_{H_2O}$ is the molar fraction of water 
- $m_{H^+}$ is the molality of hydronium $[mol / kg \ H_2O]$
- $m_{H^+}$ is the molality of hydroxide $[mol / kg \ H_2O]$

If the constraint for water self-dissociation instead had been defined without reference to the water molar 
fraction as shown in equation below the unit of H2O would have been 
set to *"None"* instead of “x” in the code below  
$K \cdot \gamma_{H_2O} = \gamma_{H^+} m_{H^+} \cdot \gamma_{OH^-} m_{OH^-} $  

```python
def water_autoprotolysis_eq_constant(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    Kw = 10 ** (-14) * np.exp(-13445.9 * (1 / T - 1 / 298) - 22.48 * np.log(T / 298))
    return Kw

def carbonic_acid_rxn_1_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-6.32) * np.exp(5139 * (1 / T - 1 / 298) + 14.5258479 * np.log(T / 298))
    return K

def carbonic_acid_rxn_2_eq_const(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    K = 10 ** (-10.33) * np.exp(22062 * (1 / T - 1 / 298) + 67.264072 * np.log(T / 298))
    return K

sol.add_rxn_insta(id="H2O = H+ + OH-",
                  stoch={"H2O": -1, "H+": 1, "OH-": 1},
                  unit={"H2O": "x", "H+": "m", "OH-": "m"},
                  equilibrium_constant=water_autoprotolysis_eq_constant)

sol.add_rxn_insta(id="CO2 + H2O = HCO3- + H+",
                  stoch={"H2O": -1, "CO2": -1, "H+": 1, "HCO3-": 1},
                  unit={"H2O": "x", "CO2": "m", "H+": "m", "HCO3-": "m"},
                  equilibrium_constant=carbonic_acid_rxn_1_eq_const)

sol.add_rxn_insta(id="HCO3- = CO3-2 + H+",
                  stoch={"HCO3-": -1, "CO3-2": 1, "H+": 1},
                  unit={"HCO3-": "m", "CO3-2": "m", "H+": "m"},
                  equilibrium_constant=carbonic_acid_rxn_2_eq_const)
```


Vapor pressure of volatile species is defined using Henry’s or Raoult’s Law using one of the methods
*add_vapor_pressure_bara_henry* or *add_vapor_pressure_bara_raoult*.
In this example CO2 and H2O are assumed to be volatile.
In the below code the argument *liq_unit* for CO2 is set to molality and the CO2 vapor pressure is therefore calculated using following equation.  
$p_{CO_2}^* = \frac{1}{H_{CO2}} \cdot \gamma_{CO_2} \cdot m_{CO_2}$  

For Raoult's Law the unit is always molar fraction, and need not be specified.  
$p_{H_2O}^* = \gamma_{H_2O} \cdot x_{H_2O} \cdot p^0_{H_2O}$  


```python
def CO2_henrys_constant(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    H_CO2 = 1.153 * np.exp((-T * (1713 * (1 - 0.0015453 * T) ** (1 / 3) + 3680) + 1198506) / T ** 2)
    return H_CO2

def H2O_vapor_pressure_bara(LiquidStream):
    T = LiquidStream.get_solution_temp_K()
    pc = 220.64
    Tc = 647.096
    tau = 1 - T / Tc
    a1, a2, a3, a4, a5, a6 = -7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502
    p = pc * np.exp((Tc / T) * (a1 * tau + a2 * tau ** 1.5 + a3 * tau ** 3 + a4 * tau ** 3.5 + a5 * tau ** 4 + a6 * tau ** 7.5))
    return p

sol.add_vapor_pressure_bara_henry(id="CO2(g) = CO2(aq)",
                                  gas_id="CO2",
                                  liq_id="CO2",
                                  liq_unit="m",
                                  henrys_coefficient=CO2_henrys_constant)

sol.add_vapor_pressure_bara_raoult(id="H2O(g) = H2O(l)",
                                   gas_id="H2O",
                                   liq_id="H2O",
                                   pure_vapor_pressure_bara=H2O_vapor_pressure_bara)
```

The Solvent is with the code above completely configured, which enable the following calculations to be made.
- Liquid Equilibrium (Isothermal)
- Liquid Equilibrium (Adiabatic)
- Flash Calculation
- Vapor Pressure from Solution
- Heat of Absorption for Volatile Species

As an example, let say we want the equilibrium concentration at a
temperature of 25°C for various CO2 Loads for a 20% K2CO3 Solution. Temperature and Flow of the solution is set using
*set_solution_temp_K* and *set_solution_flow_kg_h* methods. The argument have to be 1D Numpy Arrays.
Initial concentrations can either be defined using mass fractions as shown below, or alternatively using
*set_species_molar_fractions* or *set_species_molality*.

```python
N = 100
shape = np.ones(shape=(N,))
sol.set_solution_temp_K(value=273.15 + 40 * shape)
sol.set_solution_flow_kg_h(value=shape)
sol.set_specie_mass_fraction(id="CO2", value=np.linspace(0.001, 0.15, N))
sol.set_specie_mass_fraction(id="CO3-2", value=0.2 * shape * 1 * 60 / (2*39 + 1*60))
sol.set_specie_mass_fraction(id="HCO3-", value=0.0 * shape)
sol.set_specie_mass_fraction(id="H2O", value=0.8 * shape)
sol.set_specie_mass_fraction(id="H+", value=0.0 * shape)
sol.set_specie_mass_fraction(id="OH-", value=0.0 * shape)
sol.set_specie_mass_fraction(id="K+", value=0.2 * shape * 2 * 39 / (2*39 + 1*60))
sol.normalize_mass_fractions()
```


In the case *set_species_molar_fractions* or *set_species_molality* are used the argument must be a dictionary
defining the concentration of all species at once, as shown in example below. Note that if molality is used the
solvent specie (typically H2O) need not be specified.
```python
"""""""""
solutes_molality_mol_kg = {"CO2": np.linspace(0.01, 10, N),
                           "CO3-2": 5.0 * shape,
                           "HCO3-": 0.0 * shape,
                           "H+": 0.0 * shape,
                           "OH-": 0.0 * shape,
                           "K+" : 5.0 * shape}
sol.set_species_molality(solutes_molality_mol_kg=solutes_molality_mol_kg)
"""""""""
```

With initial concentration set equilibrium is calculated using the Class *LiquidEquilibrium_Isothermal*. The
argument *lr* refer to the learning rate of the algorithm and should be an number between zero and one. The
Algorithm used to find the equilibrium concentrations is Newton’s Method, and a learning rate equal to one correspond to one "Newton Step"".
If the code fail to converge the
learning rate must be reduced.

```python
sol = lab.LiquidEquilibrium_Isothermal().react(sol, lr=0.75)

plt.figure(1)
plt.plot(CO2Load(sol), sol.get_specie_molality_mol_kg(id="CO3-2"), label="CO3-2")
plt.plot(CO2Load(sol), sol.get_specie_molality_mol_kg(id="HCO3-"), label="HCO3-")
plt.plot(CO2Load(sol), sol.get_specie_molality_mol_kg(id="CO2"), label="CO2")
plt.xlabel("CO2 Load")
plt.ylabel("Molality")
plt.legend()
plt.grid(True)

plt.figure(2)
pH = - np.log10(sol.get_specie_molality_mol_kg(id="H+"))
plt.plot(CO2Load(sol), pH)
plt.xlabel("CO2 Load")
plt.ylabel("pH")
plt.grid(True)
```

With equilibrium concentrations calculated the vapor pressure of the volatile species can be estimated using
the *get_specie_vapor_pressure_bara* method. Note that equilibrium must be calculated before trying to
obtain the vapor pressure.

```python
plt.figure(3)
plt.plot(CO2Load(sol), sol.get_specie_vapor_pressure_bara(gas_id="CO2"), label="CO2")
plt.plot(CO2Load(sol), sol.get_specie_vapor_pressure_bara(gas_id="H2O"), label="H2O")
plt.xlabel("CO2 Load")
plt.ylabel("Vapor Pressure [bara]")
plt.legend()
plt.grid(True)
plt.yscale("log")
```

Heat of absorption is calculated using the Clausius–Clapeyron relation via the
*get_heat_of_vaporization_kJ_kmol* method.

```python
plt.figure(4)
eq = lab.LiquidEquilibrium_Isothermal()
heat_of_vap_CO2 = eq.get_heat_of_vaporization_kJ_kmol(sol, gas_id="CO2", lr=0.75)
heat_of_vap_H2O = eq.get_heat_of_vaporization_kJ_kmol(sol, gas_id="H2O", lr=0.75)
plt.plot(CO2Load(sol), heat_of_vap_CO2/1000, label="CO2")
plt.plot(CO2Load(sol), heat_of_vap_H2O/1000, label="H2O")
plt.xlabel("CO2 Load")
plt.ylabel("Heat of Vaporization [kJ/mol]")
plt.legend()
plt.grid(True)
```


Reboiler can be modelled using a Flash Calculation via the *LiquidEquilibrium_QPFlash* method. Let’s say
our initial solvent is sent to a reboiler at 40°C and a flow rate of 5000 kg/h. The Heat Input to the reboiler is
set to 300 kW. Next we assume a pressure of 2.0 bar(a). We varying the initial CO2 Load of the solution.  


The result from such a system is shown in graph
below. Note that below a certaint CO2 Load there is zero boiloff which imply that the total
vapor pressure from the solution is less than 2.0 bar(a) and consequently no flashing occur. Note, the algorithm run faster if the argument *flash_always_occur* is set to True. However, if it is set to True and flashing don’t
occur for every input the algorithm may get stuck.

```python
plt.figure(5)
sol.set_solution_temp_K(value=273.15 + 40 * shape)
sol.set_solution_flow_kg_h(value=5000 * shape)

reboiler = lab.LiquidEquilibrium_QPFlash()

boiloff, lean = reboiler.react(LiquidStreamIn=sol,
                               heat_kW=300 * shape,
                               pressure_bara=2.0 * shape,
                               lr=0.75,
                               flash_always_occur=False)

plt.plot(CO2Load(sol), boiloff.get_specie_flow_kg_h(id="CO2"), label="CO2")
plt.plot(CO2Load(sol), boiloff.get_specie_flow_kg_h(id="H2O"), label="H2O")
plt.xlabel("CO2 Load")
plt.ylabel("Boiloff [kg/h]")
plt.grid(True)
plt.legend()

plt.show()
```


### 5.2 GasStream

GasStream Objects are configured much the same way as a LiquidStream Objects. First an Instance of
GasStream is made, and then all species in the gas is added.


```python
import numpy as np
import matplotlib.pyplot as plt
import testlablib as lab

gas = lab.GasStream()

gas.add_specie(id="CO2", molar_mass_kg_kmol=44, charge=0)
gas.add_specie(id="O2", molar_mass_kg_kmol=32, charge=0)
gas.add_specie(id="H2O", molar_mass_kg_kmol=18, charge=0)
gas.add_specie(id="N2", molar_mass_kg_kmol=28, charge=0)
```


Next Heat Capacity is defined. Note that the units are in $kJ \cdot kmol^{-1} \cdot K^{-1}$ in contrast to liquids that use *kg*
instead of *kmol*. The function for Heat Capacity must take both a GasStream Object and an id (string) as
argument.


```python
def heat_capacity_kJ_kmolK(GasStream, id):
    T = GasStream.get_gas_temp_K()
    if id == "O2":
        A, B, C, D, E = 29.103, 10.040, 2526.5, 9.356, 1153.8
        Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
    elif id == "N2":
        A, B, C, D, E = 29.105, 8.6149, 1701.6, 0.10347, 909.79
        Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
    elif id == "H2O":
        A, B, C, D, E = 33.363, 26.790, 2610.5, 8.896, 1169
        Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
    elif id == "CO2":
        A, B, C, D, E = 29.370, 34.540, 1428, 26.4, 588
        Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
    return Cp

gas.load_heat_capacity_kJ_kmolK(function=heat_capacity_kJ_kmolK)

```


Configuration of the GasStream are with the code above complete. What remains are then to
set pressure, flow, temperature and molar fractions as shown below.


```python
N = 100
shape = np.ones(shape=(N,))
t = np.linspace(0,200,N)
T = t + 273.15

gas.set_gas_pressure_bara(value=1.0 * shape)
gas.set_gas_flow_kmol_h(value=shape)
gas.set_gas_temp_K(value=T)
gas.set_specie_molar_fraction(id="CO2", value=0.065 * shape)
gas.set_specie_molar_fraction(id="H2O", value=0 * shape)
gas.set_specie_molar_fraction(id="O2", value=0.123 * shape)
gas.set_specie_molar_fraction(id="N2", value=0.812 * shape)

plt.figure(6)
plt.plot(t, gas.get_gas_density_kg_m3())
plt.xlabel("Temp [C]")
plt.ylabel("Density [kg/m3]")
plt.grid(True)

plt.figure(7)
plt.plot(t, gas.get_specie_heat_capacity_kJ_kmolK(id="CO2"), label="CO2")
plt.plot(t, gas.get_specie_heat_capacity_kJ_kmolK(id="O2"), label="O2")
plt.plot(t, gas.get_specie_heat_capacity_kJ_kmolK(id="N2"), label="N2")
plt.plot(t, gas.get_specie_heat_capacity_kJ_kmolK(id="H2O"), label="H2O")
plt.xlabel("Temp [C]")
plt.ylabel("Heat Capacity [kJ/kmol.K]")
plt.grid(True)
plt.legend()
plt.show()
```


### 5.3 Gas Liquid Contactors

To get started with Gas-Liquid-Contactors consider an example of Air Contaminated with NH3, and to remove the NH3 from the Air it is directed through a Packed Column containing Acidic Water.
The Air and Acidic water is flowing Counter-Current.  

We have the following reactions;  

Gas-Liquid Phase  
$H_2O(g) \rightleftarrows H_2O(l)$  
$NH_3(g) \rightleftarrows NH_3(aq)$  

Liquid Phase  
$NH_4^+ \rightleftarrows NH_3 + H^+$  
$H_2O \rightleftarrows OH^- + H^+$    


To avoid having to re-configure the LiquidStream and GasStream Objects every time one need another instance they can be
wrapped inside a Class as shown below. 



```python
import numpy as np
import matplotlib.pyplot as plt
import testlablib as lab


class Air(lab.GasStream):

    def __init__(self, flow_Nm3_h_dry, pressure_bara, temp_K, humidity_pct, NH3_parts_per_million):
        super().__init__()

        self.load_heat_capacity_kJ_kmolK(function=self.__heat_capacity_kJ_kmolK__)

        self.add_specie(id="NH3", molar_mass_kg_kmol=17, charge=0)
        self.add_specie(id="O2", molar_mass_kg_kmol=32, charge=0)
        self.add_specie(id="H2O", molar_mass_kg_kmol=18, charge=0)
        self.add_specie(id="N2", molar_mass_kg_kmol=28, charge=0)

        NH3_molar_fraction_dry = 10**(-6) * NH3_parts_per_million
        O2_molar_fraction_dry = 0.21 * shape
        N2_molar_fraction_dry = 1 - NH3_molar_fraction_dry - O2_molar_fraction_dry

        H2O_molar_fraction = (humidity_pct / 100) * self.__H2O_vapor_pressure_bara__(temp_K) / pressure_bara

        NH3_molar_fraction_wet = NH3_molar_fraction_dry * (1 - H2O_molar_fraction)
        O2_molar_fraction_wet = O2_molar_fraction_dry * (1 - H2O_molar_fraction)
        N2_molar_fraction_wet = N2_molar_fraction_dry * (1 - H2O_molar_fraction)

        self.set_specie_molar_fraction(id="NH3", value=NH3_molar_fraction_wet)
        self.set_specie_molar_fraction(id="O2", value=O2_molar_fraction_wet)
        self.set_specie_molar_fraction(id="N2", value=N2_molar_fraction_wet)
        self.set_specie_molar_fraction(id="H2O", value=H2O_molar_fraction)

        self.normalize_molar_fractions()

        self.set_gas_pressure_bara(value=pressure_bara)
        self.set_gas_temp_K(value=temp_K)

        flow_kmol_h_dry = flow_Nm3_h_dry / (0.08314 * 273.15)
        flow_kmol_h_H2O = flow_kmol_h_dry * H2O_molar_fraction
        flow_kmol_h_wet = flow_kmol_h_dry + flow_kmol_h_H2O
        self.set_gas_flow_kmol_h(value=flow_kmol_h_wet)

    def __heat_capacity_kJ_kmolK__(self, GasStream, id):
        T = GasStream.get_gas_temp_K()
        if id == "O2":
            A, B, C, D, E = 29.103, 10.040, 2526.5, 9.356, 1153.8
            Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
        elif id == "N2":
            A, B, C, D, E = 29.105, 8.6149, 1701.6, 0.10347, 909.79
            Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
        elif id == "H2O":
            A, B, C, D, E = 33.363, 26.790, 2610.5, 8.896, 1169
            Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
        elif id == "CO2":
            A, B, C, D, E = 29.370, 34.540, 1428, 26.4, 588
            Cp = A + B * ((C / T) / (np.sinh(C / T))) ** 2 + D * ((E / T) / (np.cosh(E / T))) ** 2
        else:
            Cp = 30 * np.ones(shape=T.shape)
        return Cp

    def __H2O_vapor_pressure_bara__(self, T):
        pc = 220.64
        Tc = 647.096
        tau = 1 - T / Tc
        a1, a2, a3, a4, a5, a6 = -7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502
        p = pc * np.exp((Tc / T) * (a1 * tau + a2 * tau ** 1.5 + a3 * tau ** 3 + a4 * tau ** 3.5 + a5 * tau ** 4 + a6 * tau ** 7.5))
        return p


class MySolvent(lab.LiquidStream):

    def __init__(self, temp_K, flow_kg_h, HCl_molality, NH3_molality):

        super().__init__(solvent_id="H2O")

        shape = np.ones(shape=temp_K.shape)

        self.load_activity_coefficient(function=self.__activity_coefficient__)
        self.load_density_kg_m3(function=self.__density_kg_m3__)
        self.load_heat_capacity_kJ_kgK(function=self.__heat_capacity_kJ_kgK__)

        self.add_specie(id="H2O", molar_mass_kg_kmol=18, charge=0)
        self.add_specie(id="H+", molar_mass_kg_kmol=1, charge=1)
        self.add_specie(id="OH-", molar_mass_kg_kmol=17, charge=-1)
        self.add_specie(id="Cl-", molar_mass_kg_kmol=35.5, charge=-1)
        self.add_specie(id="NH3", molar_mass_kg_kmol=17, charge=0)
        self.add_specie(id="NH4+", molar_mass_kg_kmol=18, charge=1)

        solutes_molality_mol_kg = {"H+" : HCl_molality,
                                   "OH-" : 0 * shape,
                                   "Cl-" : HCl_molality,
                                   "NH3" : NH3_molality,
                                   "NH4+" : 0 * shape}

        self.set_species_molality(solutes_molality_mol_kg=solutes_molality_mol_kg)

        self.add_rxn_insta(id="H2O = H+ + OH-",
                          stoch={"H2O": -1, "H+": 1, "OH-": 1},
                          unit={"H2O": "x", "H+": "m", "OH-": "m"},
                          equilibrium_constant=self.__water_autoprotolysis_eq_constant__)

        self.add_rxn_insta(id="NH4+ = NH3 + H+",
                          stoch={"NH4+": -1, "H+": 1, "NH3": 1},
                          unit={"NH4+": "m", "H+": "m", "NH3": "m"},
                          equilibrium_constant=self.__nh3_dissociation_constant__)

        self.add_vapor_pressure_bara_henry(id="NH3(g) = NH3(aq)",
                                           gas_id="NH3",
                                           liq_id="NH3",
                                           liq_unit="m",
                                           henrys_coefficient=self.__NH3_henrys_constant__)

        self.add_vapor_pressure_bara_raoult(id="H2O(g) = H2O(l)",
                                           gas_id="H2O",
                                           liq_id="H2O",
                                           pure_vapor_pressure_bara=self.__H2O_vapor_pressure_bara__)

        self.set_solution_temp_K(value=temp_K)
        self.set_solution_flow_kg_h(value=flow_kg_h)

    def __activity_coefficient__(self, LiquidStream, id):
        I = LiquidStream.get_solution_ionic_strength_mol_kg()
        z = LiquidStream.get_specie_charge(id)
        log10_gamma = - 0.51 * z ** 2 * np.sqrt(I) / (1 + 1.5 * np.sqrt(I))
        gamma = 10 ** log10_gamma
        return gamma

    def __heat_capacity_kJ_kgK__(self, LiquidStream):
        return 4.2 * np.ones(shape=LiquidStream.get_solution_temp_K().shape)

    def __density_kg_m3__(self, LiquidStream):
        return 1050 * np.ones(shape=LiquidStream.get_solution_temp_K().shape)

    def __water_autoprotolysis_eq_constant__(self, LiquidStream):
        T = LiquidStream.get_solution_temp_K()
        Kw = 10 ** (-14) * np.exp(-13445.9 * (1 / T - 1 / 298) - 22.48 * np.log(T / 298))
        return Kw

    def __nh3_dissociation_constant__(self, LiquidStream):
        T = LiquidStream.get_solution_temp_K()
        K = 10 ** (-10.34) * np.exp(-7988 * (1/T - 1/298))
        return K

    def __NH3_henrys_constant__(self, LiquidStream):
        T = LiquidStream.get_solution_temp_K()
        H_NH3 = 56 * np.exp(4100 * (1/T - 1/298))
        return H_NH3

    def __H2O_vapor_pressure_bara__(self, LiquidStream):
        T = LiquidStream.get_solution_temp_K()
        pc = 220.64
        Tc = 647.096
        tau = 1 - T / Tc
        a1, a2, a3, a4, a5, a6 = -7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502
        p = pc * np.exp((Tc / T) * (a1 * tau + a2 * tau ** 1.5 + a3 * tau ** 3 + a4 * tau ** 3.5 + a5 * tau ** 4 + a6 * tau ** 7.5))
        return p
```

There are two ways to simulate a Packed Column
1) Equilibrium Stages
2) Rate-Based Method


***Equilibrium Stages***  
Example of Equilibrium Stage Implementation is shown in below code.  
Gas flow is varied while keeping all other parameters stable.  

For a single Vapor-Liquid Equilibrium Calculation (One Single Stage) the Class
*VaporLiquidEquilibrium_Isothermal* or *VaporLiquidEquilibrium_Adiabatic* can be used as shown in below
code.  

Note that the initial solvent must be processed through a *LiquidEquilibrium_Isothermal* Object first if we
want to make sure that the solvent is 25 °C and at equilibrium when entering the absorber.
Recall that we have defined NH3 and H2O to be volatile species in the LiquidStream Object. The algorithm
will automatically assure that these two species are in equilibrium with the gas, while all other species are
considered non-volatile.




```python
N = 25
shape = np.ones(shape=(N,))
flow_Nm3_h_dry = np.linspace(200, 5000, N)

gas_in = Air(flow_Nm3_h_dry=flow_Nm3_h_dry,
             pressure_bara=1.0 * shape,
             temp_K=298.15*shape,
             humidity_pct=100*shape,
             NH3_parts_per_million=25000*shape)

liq_in = MySolvent(temp_K=298.15*shape,
                   flow_kg_h=5000*shape,
                   HCl_molality=0.5*shape,
                   NH3_molality=0.0*shape)
liq_in = lab.LiquidEquilibrium_Isothermal().react(liq_in, lr=0.75)


plt.figure(1)
absorber = lab.VaporLiquidEquilibrium_Adiabatic()
gas_out, liq_out = absorber.react(gas_in, liq_in, lr=0.75)
CO2_cap = 100 * (1 - gas_out.get_specie_flow_kg_h(id="NH3") / gas_in.get_specie_flow_kg_h(id="NH3"))
plt.plot(flow_Nm3_h_dry, CO2_cap, label="1 Stage")

absorber = lab.VaporLiquidEquilibrium_EquilibriumStages(num_of_stages=2)
gas_out, liq_out = absorber.react(gas_in, liq_in)
CO2_cap = 100 * (1 - gas_out.get_specie_flow_kg_h(id="NH3") / gas_in.get_specie_flow_kg_h(id="NH3"))
plt.plot(flow_Nm3_h_dry, CO2_cap, label="2 Stages")

absorber = lab.VaporLiquidEquilibrium_EquilibriumStages(num_of_stages=3)
gas_out, liq_out = absorber.react(gas_in, liq_in)
CO2_cap = 100 * (1 - gas_out.get_specie_flow_kg_h(id="NH3") / gas_in.get_specie_flow_kg_h(id="NH3"))
plt.plot(flow_Nm3_h_dry, CO2_cap, label="3 Stages")

plt.legend()
plt.grid(True)
plt.xlabel("Exhaust Gas Flow [Nm3/h]")
plt.ylabel("NH3 Captured [%]")
plt.show()
```

***Rate-Based Method***    
Simulating Absorption Column using Rate-Based Method is shown in below code. We then need to specify:
- The Liquid Holdup
- The Pressure Drop
- The Heat Transfer between the Gas and Liquid
- The Mass Transfer of the Volatile Species
- The Heat of Absorption of the Volatile Species


The calculated Mass and Heat Transfer Coefficients are for illustration only, and is not in any way accurate.


```python
absorber = lab.Column_StructuredPacking_CounterCurrent(height_m=5.6,
                                                       num_of_heights=100,
                                                       cross_sectional_area_m2=0.5,
                                                       void_fraction_m3_m3=0.98,
                                                       packing_area_m2_m3=350,
                                                       corrugation_angle_degree=60)


def pressure_drop_Pa_m(Column):
    return 0.0 * np.ones(shape=Column.LiquidStream.temp_K.shape)


def liquid_holdup_m3_m3(Column):
    return 0.1 * np.ones(shape=Column.LiquidStream.temp_K.shape)


def Heat_Transfer_kW_m3(Column):

    # Features
    ap = Column.get_packing_area_m2_m3()
    theta = Column.get_corrugation_angle_degree() * (np.pi / 180)
    M = (3 * ap ** 3 * np.sin(theta) * np.cos(theta)) / (16 * (np.sin(theta) ** 2 + 1) ** (3 / 2))
    Mi = M / ap ** 3
    m = Column.LiquidStream.get_solution_flow_kg_h()
    A = Column.get_cross_sectional_area_m2()
    T_gas = Column.GasStream.get_gas_temp_K()
    T_liq = Column.LiquidStream.get_solution_temp_K()
    v_gas = Column.get_superficial_gas_velocity_m_s()

    # Heat Transfer Coefficient
    kHa = 11.14 * ((m / 6000) * (0.5 / A)) ** 0.15 * (v_gas / 2.4) ** 0.54 * (Mi / 0.035) **0.29 * (ap / 350) ** 1.22

    # Heat Transfer [kW/m3]
    q = kHa * (T_gas - T_liq)
    return q


def Mass_Transfer_NH3_kJ_kmol(Column):
    T0 = Column.LiquidStream.temp_K
    T1 = Column.LiquidStream.temp_K + 0.05
    H0 = 56 * np.exp(4100 * (1 / T0 - 1 / 298))
    H1 = 56 * np.exp(4100 * (1 / T1 - 1 / 298))
    q = 8.314 * (np.log(H1) - np.log(H0)) / ((1 / T1) - (1 / T0))
    return q


def Mass_Transfer_NH3_kmol_m3s(Column):

    # Features
    ap = Column.get_packing_area_m2_m3()
    theta = Column.get_corrugation_angle_degree() * (np.pi / 180)
    v_gas = Column.get_superficial_gas_velocity_m_s()
    v_liq = Column.get_superficial_liquid_velocity_m_s()
    m = Column.LiquidStream.get_solution_flow_kg_h()
    A = Column.get_cross_sectional_area_m2()
    T_liq = Column.LiquidStream.get_solution_temp_K()

    # Driving Force
    p_NH3 = Column.GasStream.get_specie_pressure_bara(id="NH3")
    p_NH3_vap = Column.LiquidStream.get_specie_vapor_pressure_bara(gas_id="NH3")

    # Henry's Law Coefficient
    H_NH3 = 56 * np.exp(4100 * (1 / T_liq - 1 / 298))

    # Mass Transfer Coefficient
    KGa = 0.1

    #Absorption Rate
    r = KGa * (p_NH3 - p_NH3_vap)
    return r


def Mass_Transfer_H2O_kJ_kmol(Column):
 return 44000


def Mass_Transfer_H2O_kmol_m3s(Column):

    # Features
    ap = Column.get_packing_area_m2_m3()
    theta = Column.get_corrugation_angle_degree() * (np.pi / 180)
    M = (3 * ap ** 3 * np.sin(theta) * np.cos(theta)) / (16 * (np.sin(theta) ** 2 + 1) ** (3 / 2))
    Mi = M / ap ** 3
    v_gas = Column.get_superficial_gas_velocity_m_s()
    v_liq = Column.get_superficial_liquid_velocity_m_s()
    m = Column.LiquidStream.get_solution_flow_kg_h()
    A = Column.get_cross_sectional_area_m2()

    # Driving Force
    p_H2O = Column.GasStream.get_specie_pressure_bara(id="H2O")
    p_H2O_vap = Column.LiquidStream.get_specie_vapor_pressure_bara(gas_id="H2O")

    # Overall Volumetric Mass Transfer Coefficient
    KGa = 0.325 * ((m / 6000) * (0.5 / A)) ** 0.15 * (v_gas / 2.4) ** 0.54 * (Mi / 0.035) ** 0.29 * (ap / 350) ** 1.22
    r = KGa * (p_H2O - p_H2O_vap)
    return r


absorber.add_pressure_drop_Pa_m(pressure_drop_Pa_m=pressure_drop_Pa_m)

absorber.add_liquid_holdup_m3_m3(liquid_holdup_m3_m3=liquid_holdup_m3_m3)

absorber.add_heat_transfer_kW_m3(heat_transfer_kW_m3=Heat_Transfer_kW_m3)

absorber.add_mass_transfer_kmol_m3s(id="NH3(g) -> NH3(aq)",
                                    stoch_gas={"NH3": -1},
                                    stoch_liq={"NH3": 1},
                                    rate_kmol_m3s=Mass_Transfer_NH3_kmol_m3s,
                                    exothermic_heat_kJ_kmol=Mass_Transfer_NH3_kJ_kmol)

absorber.add_mass_transfer_kmol_m3s(id="H2O(g) -> H2O(aq)",
                                    stoch_gas={"H2O": -1},
                                    stoch_liq={"H2O": 1},
                                    rate_kmol_m3s=Mass_Transfer_H2O_kmol_m3s,
                                    exothermic_heat_kJ_kmol=Mass_Transfer_H2O_kJ_kmol)
```

With the above code the Absorber is Configured, and ready to be deployed.  
In below example Air Flow at Absorber inlet is varied, while all other variables are kept constant.  
A plot of the NH3 Capture rate as a function of Air Flow is then generated.  
Profile plots along the height of the Absorber is available. In the below example the temperature profile is plotted for the highest Air Flow Rate.

```python
N = 20
shape = np.ones(shape=(N,))
flow_Nm3_h_dry = np.linspace(1000, 5000, N)

gas_in = Air(flow_Nm3_h_dry=flow_Nm3_h_dry,
             pressure_bara=1.0 * shape,
             temp_K=298.15*shape,
             humidity_pct=100*shape,
             NH3_parts_per_million=25000 * shape)

liq_in = MySolvent(temp_K=298.15*shape,
                   flow_kg_h=5000*shape,
                   HCl_molality=0.5*shape,
                   NH3_molality=0.0*shape)

liq_in = lab.LiquidEquilibrium_Isothermal().react(liq_in, lr=0.75)

gas_out, liq_out = absorber.react(GasStreamIn=gas_in, LiquidStreamIn=liq_in, epochs=200, lr=0.25)


plt.figure(2)
CO2_cap = 100 * (1 - gas_out.get_specie_flow_kg_h(id="NH3") / gas_in.get_specie_flow_kg_h(id="NH3"))
plt.plot(flow_Nm3_h_dry, CO2_cap)
plt.grid(True)
plt.xlabel("Exhaust Gas Flow [Nm3/h]")
plt.ylabel("NH3 Captured [%]")
plt.show()

plt.figure(3)
plt.plot(absorber.GasStream.get_gas_temp_K()[:,19] - 273.15, absorber.height_m, label="Gas")
plt.plot(absorber.LiquidStream.get_solution_temp_K()[:,19] - 273.15, absorber.height_m, label="Solvent")
plt.legend()
plt.grid(True)
plt.xlabel("Temperature [C]")
plt.ylabel("Height [m]")
plt.show()
```

## 6. Custom Reactors

Custom reactors for both Gases and Liquids can be made.  
Below is an simple example of a reactor for drying an GasStream.  


```python
import numpy as np
import matplotlib.pyplot as plt
import testlablib as lab
from copy import deepcopy

class Dryer:
    
    def __init__(self):
        pass
    
    def react(self, GasStreamIn):
        shape = np.ones(shape=lab.GasStream().get_gas_temp_K())
        
        GasStreamOut = deepcopy(GasStreamIn)
        
        flow_H2O_kmol_h = GasStreamOut.get_gas_flow_kmol_h() * GasStreamOut.get_specie_molar_fraction(id="H2O")
    
        GasStreamOut.set_gas_flow_kmol_h(value=GasStreamOut.get_gas_flow_kmol_h - flow_H2O_kmol_h)    
        GasStreamOut.set_specie_molar_fraction(id="H2O", value=0 * shape)
        GasStreamOut.normalize_molar_fractions()
        
        return GasStreamOut
```


Note that for GasStreams, the following are stored as numerical values:
- Temperature
- Pressure
- Molar Fractions
- Molar Flow

Consequently, a Custom Reactor altering a GasStream must change the above values.  
Other values such as Molarity and Density are not stored as numerical values, but calculated using Ideal Gas Law.  

Similarly, for LiquidStreams, the following are stored as numerical values:
- Temperature
- Mass Fractions
- Mass Flow

