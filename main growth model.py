
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

time_step = 1 # time step, hours

# Molarities (mole/L):
kN = 3e-5  # Nitrogen
kP = 1e-5  # Phophorous
kCO2 = 2.5e-5  # Carbon dioxide

# Algae concentration:
C0 = 0.05  # Initial concentration, mg/L
max_conc = 2 #  mg/L

# Light parameters:
max_light = 800  # W/m^2
k = 0.002  # Light mu multiplier, W/m^2
mu_night = -0.1  # The rate at which the algae biomass will decrease at night due to cell division. 
number_of_days = 20  # days
hpd = 24  # hours per day

# Nutrients:
N0 = 0.1  # Nitrogen initial, moles/L
P0 = 0.01  # Phophorous initial, moles/L
CO20 = 0.1  # Carbon dioxide initial, moles/L


N_moles = 16
P_moles = 1
CO2_moles = 124

h = time_step * (1/hpd)  # time step --> days
time_steps = int(number_of_days / h)  # get time steps in the number of days
# g/mol of algae per mole of P or 16 moles of N
algae_conc = 12*106 + 1*263 + 16*110 + 14*16 + 31  # total g per mole 

time = np.linspace(0, (time_steps*h), time_steps+1)  # time vector

C = np.ones(time_steps+1) * C0 # mg/L
N = np.ones(time_steps+1) * N0
P = np.ones(time_steps+1) * P0
CO2 = np.ones(time_steps+1) * CO20
I0 = np.zeros(time_steps+1)
mu = np.ones(time_steps+1) * mu_night

for i in range(1, time_steps+1):
    # Available light intensity over the course of the entire day, per unit time (W/m^2)
    I0[i] = max_light * np.sin(2 * np.pi * (i-6)/hpd) 
    
    # Average light intensity for growth:
    I0_effective = I0[i] * (max_conc - C[i-1]) / max_conc 
    
    # max growth rate given available light
    mu[i] = I0_effective * k * N[i-1]/(kN + N[i-1]) * P[i-1]/(kP + P[i-1]) * CO2[i-1]/(kCO2 + CO2[i-1])
    
    if mu[i] < 0:
        mu[i] = mu_night
    else: mu[i]
   
    # Apply the change in concentration per time step:
    deltaC = C[i-1] * mu[i] * h
    C[i] = C[i-1] + deltaC
    
    # If the algae is growing, subtract the amount of nutrients used during growth process
    if mu[i] > 0:
        N[i] = N[i-1] - deltaC / algae_conc * N_moles
        P[i] = P[i-1] - deltaC / algae_conc * P_moles
        CO2[i] = CO2[i-1] - deltaC / algae_conc * CO2_moles
        
    # Otherwise, nutrients do not change:
    else:
        N[i] = N[i-1]
        P[i] = P[i-1]
        CO2[i] = CO2[i-1]

growth_script = pd.DataFrame(columns = ['Time', 'I0', 'mu', 'Conc', 'N', 'P', 'CO2'])
growth_script['Time'] = time
growth_script['I0'] = I0
growth_script['mu'] = mu
growth_script['Conc'] = C
growth_script['N'] = N
growth_script['P'] = P
growth_script['CO2'] = CO2

plt.figure(0)
plt.plot(growth_script['Time'], growth_script['Conc'], 'k-', label='C')
plt.plot(growth_script['Time'], growth_script['N'], 'b--', label='N')
plt.plot(growth_script['Time'], growth_script['P'], 'r-.', label='P')
plt.plot(growth_script['Time'], growth_script['CO2'], 'm:', label='$CO_2$')
plt.plot(growth_script['Time'], growth_script['mu'], 'y-', label='$\mu$')
plt.legend(loc='best', fancybox=True)
plt.xlabel('Time [days]')
plt.ylabel('Concentration $K,P,CO_2$ [mole/L], C [mg/L]')
plt.savefig('p1_algae_%ddays.png' % number_of_days, dpi=300, bbox_inches='tight')
plt.show()
