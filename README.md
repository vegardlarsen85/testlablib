<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

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

---

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
<img width="569" height="432" alt="01" src="https://github.com/user-attachments/assets/f3f2ed7a-37ef-4fb6-923b-6ef17f24cf5d" />

<img width="577" height="426" alt="02" src="https://github.com/user-attachments/assets/1ea35e03-d064-43f4-ab35-3ee569aa68fb" />

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
<img width="589" height="435" alt="03" src="https://github.com/user-attachments/assets/0812b219-21ab-4282-9aec-98e1c5aa36e6" />

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
<img width="570" height="431" alt="04" src="https://github.com/user-attachments/assets/f456f121-fe5c-4469-997d-6bce6a49c1ed" />


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
<img width="591" height="430" alt="05" src="https://github.com/user-attachments/assets/35c57308-ae2e-478e-a66a-ce68f82106c3" />


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
<img width="578" height="429" alt="06" src="https://github.com/user-attachments/assets/eda9d97f-e4d9-4e1f-affb-0595c024e22f" />

<img width="572" height="432" alt="07" src="https://github.com/user-attachments/assets/9c9700b2-971d-4a37-b9c9-d2523e0c3e9e" />


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
<img width="582" height="429" alt="08" src="https://github.com/user-attachments/assets/a6be7faf-4360-46d6-91a9-56d99e91a5a6" />

***Rate-Based Method***    
Example showing simulation of Absorption Column using Rate-Based Method is shown in below code. For such a simulation we need to specify:
- Liquid Holdup
- Pressure Drop
- Heat Transfer between the Gas and Liquid
- Mass Transfer of the Volatile Species
- Heat of Absorption of the Volatile Species


In addition, any Chemical Reactions occuring in one of the phases is automatically taken care of.  
In this example we have two reactions occuring in liquid phase.  

First we will creat an instance of the Absorber for which we will use the *Column_StructuredPacking_CounterCurrent* method.
Note, a method for simulating Structured Packing Co-Current also exists.
There also exist more generic classes such as *GasLiquidContactor_PFR_CounterCurrent* were it s
possible to simulate Columns with varying Cross Sectional Area. Examples could be Venturies or Rotary Packed Beds.  

The argument *num_of_heights* determine how many discrete height steps the absorber is split into.
The library generate the differential equations from mass and energy balance and solve them
numerically. When the flow is Co-Current solving the differential equations is relatively straight
forward, as they can be integrated forward from a fixed initial condition. In the case of
Counter-Current flow the problem is a TPBVP (Two Point Boundary Value Problem) which is harder to solve.
Converging of the algorithm therefore take longer time for Counter than Co-Current Columns.  

In addition, having equilibrium constraints (aka. instantanious chemical reactions) in gas and liquid phase
complicate things. To make the calculation converge fast the matemathics used in the article below was used.

*Ultra-Fast Reactive Transport Simulations When Chemical Reactions Meet Machine Learning: Chemical Equilibrium*  
*doi.org/10.48550/arXiv.1708.04825*

The result is a library that can solve the Differential Equations from Columns suprisingly fast.  
I would dare to say that this part of the code most likely superseed Aspen's Rate Based Models by a long streth
wich model such a system using a serie of CSTR's, while this library treat the Column as a true Plug Flow Reactor
by solving the differential equations directly enabling much more accuracy. 


```python
absorber = lab.Column_StructuredPacking_CounterCurrent(height_m=5.6,
                                                       num_of_heights=100,
                                                       cross_sectional_area_m2=0.5,
                                                       void_fraction_m3_m3=0.98,
                                                       packing_area_m2_m3=350,
                                                       corrugation_angle_degree=60)
```

Next we add functions for the Column Hydrodynamics.
For simplicity we neglect any pressure drop in the column, and set the liquid holdup fixed to 10%.

```python
def pressure_drop_Pa_m(Column):
    return 0.0 * np.ones(shape=Column.LiquidStream.temp_K.shape)


def liquid_holdup_m3_m3(Column):
    return 0.1 * np.ones(shape=Column.LiquidStream.temp_K.shape)

absorber.add_pressure_drop_Pa_m(pressure_drop_Pa_m=pressure_drop_Pa_m)
absorber.add_liquid_holdup_m3_m3(liquid_holdup_m3_m3=liquid_holdup_m3_m3)
```

Next we create a function defining the Heat transfer between gas and liquid.
The driving force is not suprisingly the temperature difference between the gas and liquid,
while the Heat Transfer Coefficient is set as a function of the flow rates and packing geometry.

```python
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

absorber.add_heat_transfer_kW_m3(heat_transfer_kW_m3=Heat_Transfer_kW_m3)
```

Next we add the mass transfer rates of the volatile species. In this example NH3 and H2O.
We also need to specify the exothermic heat of the reaction to proper simulate the temperature
profile of the column.

```python

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
<img width="585" height="422" alt="09" src="https://github.com/user-attachments/assets/2a99f3c3-2036-4ab0-a9da-b692076afd5a" />

<img width="565" height="431" alt="10" src="https://github.com/user-attachments/assets/52059a0a-57e5-48e8-9a7d-0cc25b27f3cd" />


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





## 7. Numerical Methods

In this section the algorithms used in the testlablib-code is derived and justified.  


### 7.1 Equilibrium

Consider SO2 Absorption into seawater, where the two reactions below occur in the liquid phase.  
$SO_2 + H_2O \rightleftarrows HSO_3^- + H^+$  
$CO_2 + H_2O \rightleftarrows HCO_3^- + H^+$  

The equilibrium constants are typically w.r.t, molality for solutes and molar fraction for the solvent initially.  
The activity coefficients is also not typically included in the constants.  
$K_1 \left( \frac{\gamma_{SO_2} \gamma_{H_2O}}{\gamma_{HSO_3^-}\gamma_{H^+}} \right) \cdot m_{SO_2} \cdot x_{H_2O} = m_{HSO_3^-} \cdot m_{H^+}$  
$K_2 \left( \frac{\gamma_{CO_2} \gamma_{H_2O}}{\gamma_{HCO_3^-}\gamma_{H^+}} \right) \cdot m_{CO_2} \cdot x_{H_2O} = m_{HCO_3^-} \cdot m_{H^+}$  


In the remaining of the section the equilibrium constants are assumed to be with respect to the mass fractions, as shown below. Recalculating the original equilibrium constants above into the new ones below can be achieved using formulas from section 4.  
$K_{1}  \cdot w_{SO_2} \cdot w_{H_2O} = w_{HSO_3^-} \cdot w_{H^+}$  
$K_{2}  \cdot w_{CO_2} \cdot w_{H_2O} = w_{HCO_3^-} \cdot w_{H^+}$  

In total there are six species and two reactions.
Defining the following vectors and matrices.  
```math
\mathbf{w}=
\left[
\begin{matrix}
w_0\\\  w_1\\\  w_2\\\ w_3\\\ w_4\\\ w_5
\end{matrix}
\right]
=
\left[
\begin{matrix}
w_{SO_2}\\\  w_{HSO_3^-}\\\  w_{CO_2}\\\ w_{HCO_3^-}\\\ w_{H_2O}\\\ w_{H^+}
\end{matrix}
\right]
```  

```math
\mathbf{v}=
\left[
\begin{matrix}
-1 & 0 \\
1 & 0 \\
0 & -1 \\
0 & 1 \\
-1 & -1 \\
1 & 1 \\
\end{matrix}
\right]
```  

```math
\mathbf{R}=
\left[
\begin{matrix}
-64 & 0 \\
81 & 0 \\
0 & -44 \\
0 & 61 \\
-18 & -18 \\
1 & 1 \\
\end{matrix}
\right]
```  

The vector $\mathbf{w}$ is the mass fractions.  
The matrix $\mathbf{v}$ is generated from the stochiometric coefficients of the reactions.  
The matrix $\mathbf{R}$ is equal to $\mathbf{v}$ but scaled with molar masses of the species via below formula.  
Note that $\mathbf{v}$ is unitless, and consequently the unit of $\mathbf{R}$ is kg⁄kmol.  
$R_{\alpha u} = M_{\alpha} \nu_{\alpha u}$  

Let $dr_u$ be a vector representing the increase or decrease in the reactions and represents *kmol reaction per kg solution*.  
A change in the mass fractions $(d \mathbf{w})$ caused by the reactions is then given by the following matrix product.  
$dw_{\alpha} = R_{\alpha u} \cdot dr_{u}$  

With 6 species in total and 2 reactions the next objective is to find 4 conserved quantities $(\mathbf{b})$ on the form.  
$b_i = A_{i \alpha} w_{\alpha}$  

Finding the matrix $\mathbf{A}$ could be done based on conservation of sulfur, carbon etc.  
However, it is possible to calculate $\mathbf{A}$ automatically.  

Note, we are searching for four quantities $(b_i)$ that is constant when adjusting $dw_{\alpha}$ according to the formula $d\mathbf{w} = \mathbf{R} d\mathbf{r}$.  
From below derivation we observe that if the four row vectors of the matrix $A_{i \alpha}$ is orthogonal to the two column vectors in the matrix $R_{\alpha u}$ the quantities $b_i$ become constants.  
$b_i = A_{i \alpha} w_{\alpha}$  
$db_i = A_{i \alpha} dw_{\alpha}$  
$db_i = A_{i \alpha} R_{\alpha u} dr_u$  
$db_i = 0$  

We are thus searching for a matrix A with the property $A_{i \alpha} R_{\alpha u} = 0$.  
To construct the matrix $\mathbf{A}$ several options are possible: most notably SVD and QR-Decomposition.  
In the below case Singular Value Decomposition is taken of the matrix $\mathbf{R}$.  
The first two column vectors of $\mathbf{U}$ span the same space as $\mathbf{R}$ while the remaining vectors are orthogonal to $\mathbf{R}$.  
The matrix $\mathbf{A}$ can thus be generated from the last four column vectors of $\mathbf{U}$.    
$\mathbf{U},\mathbf{D},\mathbf{V}^T = SVD(\mathbf{R})$  
$A_{i \alpha} = U_{\alpha, i+2}$

Were,  
$\alpha \in \{ 0,1,2,3,4,5 \}$  
$i \in \{ 0,1,2,3 \}$  


---

---

---


**Newton's Method (I)**  
Combining the two equilibrium constraints with the four conservation laws gives six equations which must be satisfied by adjusting the concentration of the six species.  
```math
\mathbf{f} =
\left[
\begin{matrix}
\ln K_1 + \ln w_{SO_2} + \ln w_{H_2O} - \ln w_{HSO_3^-} - \ln w_{H^+} \\
\ln K_2 + \ln w_{CO_2} + \ln w_{H_2O} - \ln w_{HCO_3^-} - \ln w_{H^+} \\
A_{0 \alpha} w_{\alpha} - b_{0,init} \\
A_{1 \alpha} w_{\alpha} - b_{1,init} \\
A_{2 \alpha} w_{\alpha} - b_{2,init} \\
A_{3 \alpha} w_{\alpha} - b_{3,init} \\
\end{matrix}
\right]
```    

Partial derivatives with respect to mass fraction is straight forward to calculate if the equilibrium constants are assumed to be true constants.  
```math
\frac{\partial \mathbf{f}}{\partial \mathbf{w}} =
\left[
\begin{matrix}
\frac{\partial f_0}{\partial w_0} &
\frac{\partial f_0}{\partial w_1} &
\frac{\partial f_0}{\partial w_2} &
\frac{\partial f_0}{\partial w_3} &
\frac{\partial f_0}{\partial w_4} &
\frac{\partial f_0}{\partial w_5} \\
\frac{\partial f_1}{\partial w_0} &
\frac{\partial f_1}{\partial w_1} &
\frac{\partial f_1}{\partial w_2} &
\frac{\partial f_1}{\partial w_3} &
\frac{\partial f_1}{\partial w_4} &
\frac{\partial f_1}{\partial w_5} \\
A_{00} & A_{01} & A_{02} & A_{03} & A_{04} & A_{05} \\
A_{10} & A_{11} & A_{12} & A_{13} & A_{14} & A_{15} \\
A_{20} & A_{21} & A_{22} & A_{23} & A_{24} & A_{25} \\
A_{30} & A_{31} & A_{32} & A_{33} & A_{34} & A_{35} \\
\end{matrix}
\right]
=
\left[
\begin{matrix}
\frac{1}{w_0} &
-\frac{1}{w_1} &
0 &
0 &
\frac{1}{w_4} &
-\frac{1}{w_5} \\
0 &
0 &
\frac{1}{w_2} &
-\frac{1}{w_3} &
\frac{1}{w_4} &
-\frac{1}{w_5} \\
A_{00} & A_{01} & A_{02} & A_{03} & A_{04} & A_{05} \\
A_{10} & A_{11} & A_{12} & A_{13} & A_{14} & A_{15} \\
A_{20} & A_{21} & A_{22} & A_{23} & A_{24} & A_{25} \\
A_{30} & A_{31} & A_{32} & A_{33} & A_{34} & A_{35} \\
\end{matrix}
\right]
```  

***Pseudocode***  
$b_{i,init} = A_{i \alpha} w_{\alpha}$  
$\epsilon \approx 10^{-16}$  
$lr \approx 0.75$  
*Repeat Until Convergence*  
--- *Calculate $\mathbf{K}$*  
--- *Calculate* $\mathbf{f}$  
--- *Calculate* $\frac{\partial \mathbf{f}}{\partial \mathbf{w}}$  
--- $\Delta w_{\alpha} = - lr \cdot \left( \frac{\partial f_{\alpha}}{\partial w_{\beta}} \right)^{-1} f_{\beta}$  
--- $w_{\alpha} = \max \{ \epsilon, w_{\alpha} + \Delta w_{\alpha}  \}$  


---

---

---

**Newton's Method (II)**  
In the above case the most resource demanding step is to invert the 6x6 matrix at every iteration.
The Algorithm can be speeded up by making some adjustments such that only inversion of a 2x2 matrix is necessary.  

To achieve this the function $\mathbf{f}$ is redefined as below.  
```math
\mathbf{f} =
\left[
\begin{matrix}
\ln K_1 + \ln w_{SO_2} + \ln w_{H_2O} - \ln w_{HSO_3^-} - \ln w_{H^+} \\
\ln K_2 + \ln w_{CO_2} + \ln w_{H_2O} - \ln w_{HCO_3^-} - \ln w_{H^+} \\
\end{matrix}
\right]
```  

We are looking to speed up the algorithm by reducing the dimensionality from 6 to 2 dimensions (6 species and 2 reactions in this problem).  
We will therefore only update along the 2-dimensional subspace given by the vector $\Delta \mathbf{r}$.  
$\Delta w_{\alpha} = R_{\alpha u} \Delta r_u$  

The 2-dimensional vector $\Delta \mathbf{r}$ is updated using Newton's Method.  
$\Delta r_u = -lr \cdot \left( \frac{\partial f_u}{\partial r_v} \right)^+ f_v$  

Obtaining the change in mass fraction can be made by matrix multiplication with the matrix $\mathbf{R}$.  
$\Delta w_{\alpha} = - lr \cdot R_{\alpha u} \left( \frac{\partial f_u}{\partial r_v} \right)^+ f_v$

Applying Chain Rule  
$\Delta w_{\alpha} = - lr \cdot R_{\alpha u} \left( \frac{\partial f_u}{\partial w_{\beta}} \frac{\partial w_{\beta}}{\partial r_v} \right)^+ f_v$  

Final Update Rule  
$\Delta w_{\alpha} = - lr \cdot R_{\alpha u} \left( \frac{\partial f_u}{\partial w_{\beta}} R_{\beta v} \right)^+ f_v$  

The downside by making this update is that the step size must be reduced such that it doesn’t try to impose negative values.
This is achieved by introducing the parameter $\tau$.
Let $w_{\alpha}$ be the current mass fractions and  $\Delta w_{\alpha}$ the requested update from Newton’s Method.
τ is then calculated via the following steps.  
$\tau_{\alpha} = - w_{\alpha} / \Delta w_{\alpha}$  
$\tau_{\alpha} = 1 \cdot (\tau_{\alpha} < 0) + \tau_{\alpha} (\tau_{\alpha} > 0)$  
$\tau = \min_{\alpha} \tau_{\alpha}$  
$\tau = \min \{ \tau ,1 \}$  


In addition, the mass fractions must initially obey the conservation laws as the update rule only can update along the subspace;   
$\Delta w_{\alpha} = R_{\alpha u} \cdot \Delta r_u$  

***Pseudocode***  
*Repeat Until Convergence*  
--- *Calculate* $\mathbf{K}$  
--- *Calculate* $\mathbf{f}$  
--- *Calculate* $\frac{\partial \mathbf{f}}{\partial \mathbf{w}}$  
--- $\Delta w_{\alpha} = - lr \cdot R_{\alpha u} \left( \frac{\partial f_u}{\partial w_{\beta}} R_{\beta v} \right)^+ f_v$  
--- *Calculate* $\tau$  
--- $w_{\alpha} = w_{\alpha} + \tau \cdot \Delta w_{\alpha} $  



---

---

---

***First Order Taylor Expansion***  
The isothermal equilibrium concentration is a state variable depending on the conserved quantities $\mathbf{b}$ and the temperature $T$,
were $\mathbf{b}$ is defined using the equation showed previous in this chapter.  
$w_{eq,\alpha} = f \left( \mathbf{b}, T \right)$  

In the contrary, the adiabatic equilibrium is not a state variable of $\mathbf{b}$ and $T$ as the initial concentration of all species must be known to properly calculate the exothermic heat from the reactions.

Anyway, back to the isothermal case; If the equilibrium concentrations are already found for a point $(\mathbf{b},T)$ the equilibrium concentration at a nearby point $(\mathbf{b'},T')$ can be approximated using First Order Taylor Expansion.  
$w'_{eq,\alpha} = w_{eq,\alpha} + \frac{\partial w_{eq,\alpha}}{\partial b_i} \left( b'_i - b_i\right) + \frac{\partial w_{eq,\alpha}}{\partial T} \left(T' - T \right)$  

To calculate Taylor Expansion the sensitivity matrices with respect to $\mathbf{b}$ and $T$ must be calculated.
The Objective Function $(\mathbf{f})$ shown below is the easiest to use for this case.
The equilibrium concentrations at point $(\mathbf{b},T)$ is assumed to be known and denoted by $w_{eq,\alpha}$.
Since the mass fractions are at equilibrium the objective functions $(\mathbf{f})$ are initially zero.
The temperature affects the equilibrium concentrations via its influence on the equilibrium constants.
The influence of a small perturbation $(d\mathbf{b}, dT)$ on the objective functions are shown in $\mathbf{f'}$.  

```math
\mathbf{f} =
\left[
\begin{matrix}
\ln K_1 + \ln w_{eq,SO_2} + \ln w_{eq,H_2O} - \ln w_{eq,HSO_3^-} - \ln w_{eq,H^+} \\
\ln K_2 + \ln w_{eq,CO_2} + \ln w_{eq,H_2O} - \ln w_{eq,HCO_3^-} - \ln w_{eq,H^+} \\
A_{0 \alpha} w_{eq,\alpha} - b_{0,init} \\
A_{1 \alpha} w_{eq,\alpha} - b_{1,init} \\
A_{2 \alpha} w_{eq,\alpha} - b_{2,init} \\
A_{3 \alpha} w_{eq,\alpha} - b_{3,init} \\
\end{matrix}
\right]
=
\left[
\begin{matrix}
0 \\
0 \\
0 \\
0 \\
0 \\
0 \\
\end{matrix}
\right]
```  

```math
\mathbf{f'} =
\left[
\begin{matrix}
\frac{1}{K_1} \frac{\partial K_1}{\partial T} dT \\
\frac{1}{K_2} \frac{\partial K_2}{\partial T} dT \\
-db_0 \\
-db_1 \\
-db_2 \\
-db_3 \\
\end{matrix}
\right]
```  

The partial derivatives $\partial \mathbf{f} / \partial \mathbf{w}$ (and therefore also the Hessian) only depend on $\mathbf{w}$ and have no direct dependence on $\mathbf{b}$ or $T$.
The Hessian at the point $(\mathbf{b'},T')$ will therefore be equal to the Hessian at point $\left( \mathbf{b}, T \right)$ as the mass fractions $(\mathbf{w})$ are left unchanged (for now).  

```math
\frac{\partial \mathbf{f}}{\partial \mathbf{w}} =
\frac{\partial \mathbf{f'}}{\partial \mathbf{w}} =
\left[
\begin{matrix}
\frac{\partial f_0}{\partial w_0} &
\frac{\partial f_0}{\partial w_1} &
\frac{\partial f_0}{\partial w_2} &
\frac{\partial f_0}{\partial w_3} &
\frac{\partial f_0}{\partial w_4} &
\frac{\partial f_0}{\partial w_5} \\
\frac{\partial f_1}{\partial w_0} &
\frac{\partial f_1}{\partial w_1} &
\frac{\partial f_1}{\partial w_2} &
\frac{\partial f_1}{\partial w_3} &
\frac{\partial f_1}{\partial w_4} &
\frac{\partial f_1}{\partial w_5} \\
A_{00} & A_{01} & A_{02} & A_{03} & A_{04} & A_{05} \\
A_{10} & A_{11} & A_{12} & A_{13} & A_{14} & A_{15} \\
A_{20} & A_{21} & A_{22} & A_{23} & A_{24} & A_{25} \\
A_{30} & A_{31} & A_{32} & A_{33} & A_{34} & A_{35} \\
\end{matrix}
\right]
```

Assuming we are at point (b,T) and make a perturbation (db,dT). Updating the mass fractions to be at equilibrium at the new point could be done using Newton’s Method, in which case the below expression is used.  
$dw_{\alpha} = - H_{\alpha \beta}' \cdot f_{\beta}'$  
$dw_{\alpha} = - H_{\alpha \beta} \cdot f_{\beta}'$  

Expanding into matrix form  
```math
  \left[
\begin{matrix}
dw_0\\ dw_1\\ dw_2\\ dw_3\\ dw_4\\ dw_5\\
\end{matrix}
\right] =
- \left[
\begin{matrix}
H_{00} & H_{01} & H_{02} & H_{03} & H_{04} & H_{05}\\\
H_{10} & H_{11} & H_{12} & H_{13} & H_{14} & H_{15}\\\
H_{20} & H_{21} & H_{22} & H_{23} & H_{24} & H_{25}\\\
H_{30} & H_{31} & H_{32} & H_{33} & H_{34} & H_{35}\\\
H_{40} & H_{41} & H_{42} & H_{43} & H_{44} & H_{45}\\\
H_{50} & H_{51} & H_{52} & H_{53} & H_{54} & H_{55}\\\
\end{matrix}
\right]
\left[
\begin{matrix}
\frac{1}{K_1} \frac{\partial K_1}{\partial T} dT\\
\frac{1}{K_2} \frac{\partial K_2}{\partial T} dT\\
-db_0\\
-db_1\\
-db_2\\
-db_3\\
\end{matrix}
\right]
```

From the above equation one observes that the sensitivity matrices can be extracted from the Hessian.  
$\frac{\partial w_{\alpha}}{\partial T} = - H_{\alpha u} \left( \frac{1}{K_u} \frac{\partial K_u}{\partial T} \right)$  
$\frac{\partial w_{\alpha}}{\partial b_{i}} = H_{\alpha, i+2}$  

With the above equations we are in a position to sketch a Nearest Neighbour inspired algorithm as used in the article *doi.org/10.48550/arXiv.1708.04825*;  

For every point where the equilibrium concentration $\left( \mathbf{w} \right)$ is found the feature vector $(\mathbf{b},T)$ and the partial derivatives $\left( \partial \mathbf{w}/\partial \mathbf{b} , \partial \mathbf{w} / \partial T \right)$ is stored into a database.
To calculate the equilibrium concentration at a new point $(\mathbf{b'}, T')$ one first finds the point closest in the database, and then extrapolate using First Order Taylor Expansion.  



