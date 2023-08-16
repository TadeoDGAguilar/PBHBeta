# Jupyter Notebook

The following examples are based on the calculations that were performed to obtain our results.

---

## Primordial Black Holes in Standard Cosmology

First step is import the PBHBeta library and necessary modules, 
our library containts the following basic modules: `constanst.py`, `constraints.py`, `functions.py`,
`constants.py`, `BfM.py`, `BfS.py`, `PfR.py`, `PfM.py` and `PfS.py`. You can find more
information about these modules in the [CLASS REFERENCE](../Classes.md), [MODULES AND COMPONENTS](../Module_1.md) sections.

On the other hand, it's necessary to import the dependencies (to more information see [Installing dependencies in python](https://pbhbeta.readthedocs.io/en/latest/Requeriments.html#installing-dependencies-in-python))

```python
from PBHBeta import *
import matplotlib.pyplot as plt
import numpy as np
```

Generate an array of masses using the put_M_array function: The initial function to execute is {py:func}`PBHBeta.functions.put_M_array`,
which generates an array of masses from a specified minimum and maximum mass value in grams.

**Note:** If you wanna run with our array data masses, you could import `constraints.data_M_tot`

```python
functions.put_M_array(2.17645e-5, 1e20)
#constraints.M_tot = constraints.data_M_tot
M_tot = np.array(constraints.M_tot)
omega=1/3
```

Once our array of masses in grams has been defined, we can now calculate the constraints for the masses and abundances.
Here, we will use all the main functions: {py:func}`PBHBeta.functions.Betas_BBN`, {py:func}`PBHBeta.functions.Betas_SD`,
{py:func}`PBHBeta.functions.CMB_AN`, {py:func}`PBHBeta.functions.Betas_GRB`, {py:func}`PBHBeta.functions.Betas_Reio`, 
and {py:func}`PBHBeta.functions.Betas_LSP`.

As most of the main functions return tuples, we can use the following syntax to assign each output to a variable 
(more information about the returns of each main function in the [Main functions](https://pbhbeta.readthedocs.io/en/latest/Module_1.html#main-functions) section.)


```hint
To ignore or not save the output of a main function, simply use an underscore (_).
```

```python
M_dm, betas_dm, M_relic, betas_relic ,_  = functions.Betas_DM(M_tot,omega)
```

```attention
To use the aforementioned hint, the number of variables must be equal to the number of outputs of each main function. 
Otherwise, you can use the notation described in the [Using Main Functions](https://pbhbeta.readthedocs.io/en/latest/Quickstart.html#using-main-functions) section.
```

Given the above, we execute all the main functions with the parameters `M_tot` and `omega`:

We save the results of masses and PBH abundances with constraints for: Dark Matter, Big Bang Nucleosynthesis,
CMB Anisotropies, Extragalactic Gamma-ray Background, Reionization, Lightest Supersymmetric Particles, and Planck Mass Relics, respectively.

```python
M_bbn = functions.Betas_BBN(M_tot,omega)[0]
M_an = functions.Betas_CMB_AN(M_tot,omega)[0]
M_grb1 = functions.Betas_GRB(M_tot,omega)[0]
M_grb2 = functions.Betas_GRB(M_tot,omega)[1]
M_reio = functions.Betas_Reio(M_tot,omega)[0]
M_lsp = functions.Betas_LSP(M_tot,omega)[0]
```

```python
betas_bbn = functions.Betas_BBN(M_tot,omega)[1]
betas_sd = functions.Betas_SD(M_tot,omega)[1]
betas_an = functions.Betas_CMB_AN(M_tot,omega)[1]
betas_grb1 = functions.Betas_GRB(M_tot,omega)[2]
betas_grb2 = functions.Betas_GRB(M_tot,omega)[3]
betas_reio = functions.Betas_Reio(M_tot,omega)[1]
betas_lsp = functions.Betas_LSP(M_tot,omega)[1]
betas_full = functions.get_Betas_full(M_tot)
```

The final instruction involves another main function called {py:func}`PBHBeta.functions.get_Betas_full`. The function calculates
a set of composite constraint values based on different scenarios involving PBHs, using the provided `M_tot`
parameter and pre-defined constraint arrays. The composite constraints represent the strongest constraints across different scenarios
for each specific mass value.


Within each main function, there are instructions to obtain the PBH abundance during BBN. 
At this epoch, PBHs with masses $M_{PBH} \leq 10^{9}$ g have already evaporated, and their constraints are imposed 
by following the evolution of the Planck mass remnants from their evaporation. 
This total abundance across various scenarios is stored within the `constraints.py` module to be used if needed. 
You can use them as follows:

```python
omegas_dm = constraints.Omega_DM_tot
omega_bbn = constraints.Omega_BBN_tot
omega_sd = constraints.Omega_SD_tot
omega_an = constraints.Omega_CMB_AN_tot
omega_grb = constraints.Omega_GRB_tot
omega_reio = constraints.Omega_Reio_tot
omega_lsp = constraints.Omega_LSP_tot
omegas_full = functions.get_Omegas_full(M_tot)
```

Al igual que la función `functions.get_Betas_full(M_tot)` de la celda anterior, 
la main-function {py:func}`PBHBeta.functions.get_Omegas_full` obtiene the most strongest constraints of PBHs.


Once the basic calculations have been performed, we can proceed to create the plot.

```python
plt.loglog(M_dm, betas_dm, "k:",label=r"$\Omega_{DM}$")
plt.fill_between(M_dm, betas_dm, betas_dm*0+10, color='grey',alpha=0.2)

plt.loglog(M_relic, betas_relic, "k--",label=r"$\Omega_{DM}^{relic}$")
plt.fill_between(M_relic, betas_relic, betas_relic*0+10, color='grey',alpha=0.2)

plt.loglog(M_bbn, betas_bbn, color = "darkkhaki",label = "BBN")
plt.fill_between(M_bbn, betas_bbn, betas_bbn*0+10, color='darkkhaki',alpha=0.2)

plt.loglog(M_an, betas_an, color='m', label = "CMB anisotropies")
plt.fill_between(M_an, betas_an, betas_an*0+10, color ="m",alpha=0.2)

plt.loglog(M_grb1, betas_grb1, color="g", label="EG GRB")
plt.fill_between(M_grb1, betas_grb1, betas_grb1*0+10, color='g',alpha=0.2)
plt.loglog(M_grb2, betas_grb2, color="g")
plt.fill_between(M_grb2, betas_grb2, betas_grb2*0+10, color='g',alpha=0.2)

plt.loglog(M_reio, betas_reio, color = "brown", label = "Reionization")
plt.fill_between(M_reio, betas_reio, betas_reio*0+10, color='brown',alpha=0.2)

plt.loglog(M_lsp, betas_lsp, color = "cadetblue",label = "LSP")
plt.fill_between(M_lsp, betas_lsp, betas_lsp*0+10, color='cadetblue',alpha=0.2)

plt.xlabel(r"$M_{\rm PBH}~[\rm{g}]$")
plt.ylabel(r"$\beta$")
plt.xlim([1,1e20])
plt.ylim([1e-30,1])

plt.legend(ncol=2,bbox_to_anchor=(1, 1.5))
plt.show()
```

Constraints on the abundance of monochromatic PBHs as a function of mass in the standard Big Bang scenario. The shaded
colored regions are excluded by observations.

![png](img/output_9_0.png)


To obtain the values of $k$ and the constraints in the power spectrum, $\mathcal{P}_{\zeta}(k)$, the function (class) {py:func}`PBHBeta.PfR.get_P_k_RD` will be used,
which performs the following operations:

\begin{equation*}
k = \left(\frac{7.1 \times 10^{2} * \gamma^{RD} * 1.5\times^{15} * H_{end}}{M_{PBH}}\right)**(1/2)\left(\frac{\rho_{0}^{RD}}{\rho_{end}^{inf}}\right)**(1/4)
\end{equation*}

Now, to calculate the constraints in the PPS by assuming that

\begin{equation*}
\beta = Erfc\left(\frac{\delta_c}{\sqrt{2}\sigma}\right), \ \ \ \ \Rightarrow \ \ \ \ \sigma = \frac{\delta_c}{\sqrt{2}Erfc^{-1}(\beta)},
\end{equation*}

with $\sigma^2\sim P(k)$ and $\delta_c = 0.41$ in the radiation case. 

Then, we will utilize the special function (class) {py:func}`PBHBeta.PfR.get_P_k_RD`, which takes the following parameters: an array of masses, 
an array of PBH abundances, and the value of $\delta_{c}$. This function returns two arrays: one containing the wave number $k$ 
and the other containing the constraints in the Power Spectrum for each $k$ value.

```{Important}
It is necessary for the arrays of mass and abundances to have the same dimensional size or number of elements.
```

```{note}
For a future version of `PBHBeta`, we might remove the dependency on the abundance array and solely use an array of masses.
```

We then proceed to utilize all the previously calculated and assigned values for each of the variables, 
including constrained masses and abundances of PBHs, assuming the standard evolution of the universe.

```python
delta_c = 0.41

k_DM, P_k_DM = PfR.get_P_k_RD(M_dm, betas_dm, delta_c)
k_BBN, P_k_BBN = PfR.get_P_k_RD(M_bbn, betas_bbn, delta_c)
k_AN, P_k_AN = PfR.get_P_k_RD(M_an, betas_an, delta_c)
k_GRB1, P_k_GRB1 = PfR.get_P_k_RD(M_grb1, betas_grb1, delta_c)
k_GRB2, P_k_GRB2 = PfR.get_P_k_RD(M_grb2, betas_grb2, delta_c)
k_Reio, P_k_Reio = PfR.get_P_k_RD(M_reio, betas_reio, delta_c)
k_LSP, P_k_LSP = PfR.get_P_k_RD(M_lsp, betas_lsp, delta_c)
k_DM_relic, P_k_DM_relic = PfR.get_P_k_RD(M_relic, betas_relic, delta_c)

```

```python
plt.loglog(k_DM, P_k_DM, "k:",label=r"$\Omega_{DM}$")
plt.fill_between(k_DM, P_k_DM, P_k_DM*0+10, color='grey',alpha=0.2)

plt.loglog(k_BBN, P_k_BBN, label="BBN",color = "darkkhaki")
plt.fill_between(k_BBN, P_k_BBN, P_k_BBN*0+10, color = "darkkhaki",alpha=0.2)

plt.loglog(k_AN, P_k_AN, "m", label="CMB anisotropies")
plt.fill_between(k_AN, P_k_AN, P_k_AN*0+10, color='m',alpha=0.2)

plt.loglog(k_GRB1, P_k_GRB1, color = "g", label="EG GRB")
plt.fill_between(k_GRB1, P_k_GRB1, P_k_GRB1*0+10, color='g',alpha=0.2)

plt.loglog(k_GRB2, P_k_GRB2, color = "g")
plt.fill_between(k_GRB2, P_k_GRB2, P_k_GRB2*0+10, color='g',alpha=0.2)

plt.loglog(k_Reio, P_k_Reio, color = "brown",label = "Reionization")
plt.fill_between(k_Reio, P_k_Reio, P_k_Reio*0+10, color='brown',alpha=0.2)

plt.loglog(k_LSP, P_k_LSP, color = "cadetblue",label = "LSP")
plt.fill_between(k_LSP, P_k_LSP, P_k_LSP*0+10, color='cadetblue',alpha=0.2)

plt.loglog(k_DM_relic, P_k_DM_relic, "k--",label=r"$\Omega_{DM}^{relic}$")
plt.fill_between(k_DM_relic, P_k_DM_relic, P_k_DM_relic*0+10, color='grey',alpha=0.2)

plt.xlabel(r"$\mathcal{k}$")
plt.ylabel(r"$\mathcal{P}_\zeta(k)$")
plt.xlim([1e13,1e22])
plt.ylim([6e-3,2e-2])
plt.show()

```

Plotting these values yields constraints imposed on the power spectrum, $\mathcal{P}_{\zeta}(k)$, 
by PBHs in the standard Big Bang scenario.
    
![png](img/output_11_0.png)


---


## Primordial Black Holes in non-standard cosmology

### Early Matter Dominated (MD) scenario

```python
k_rd, P_k_rd = PfR.get_P_k_RD(M_tot, betas_full, delta_c)

k_10s, P_k_10, beta_10 = PfM.get_P_k_MD(M_tot,10,0,1)
k_20s, P_k_20, beta_20 = PfM.get_P_k_MD(M_tot,20,0,1)
k_30s, P_k_30, beta_30 = PfM.get_P_k_MD(M_tot,30,0,1)
```

#### Abundance of PBHs ($\beta_{PBH}$) as function of their mass, where $N_{MD}$ is the total number of $e$-folds that the the MD epoch lasted.

```python
plt.loglog(M_tot,betas_full, label = "SBB", color='black')
plt.loglog(M_tot,beta_10, label = r"$N_{\rm MD}=10$",color='red')
plt.loglog(M_tot,beta_20, label = r"$N_{\rm MD}=20$",color='green')
plt.loglog(M_tot,beta_30, label = r"$N_{\rm MD}=30$",color='blue')
plt.loglog(M_tot,betas_full, color='black')
plt.ylim([1e-30,1])
plt.xlim([1e1,1e20])
plt.xlabel(r"$M_{\rm PBH}~[\rm{g}]$")
plt.ylabel(r"$\beta$")

plt.legend(ncol=2,bbox_to_anchor=(0.85, 1.23))
plt.show()
```

    
![png](img/output_13_0.png)
    
#### Constraints on $\mathcal{P}_{\zeta}(k)$ as a function of $k$ for different values of $N_{MD}$.

```python
plt.loglog(k_rd, P_k_rd,label = "SBB",color='black')
plt.fill_between(k_rd, P_k_rd, P_k_rd*0+10, color='black',alpha=0.2)

plt.loglog(k_10s[:567], P_k_10[:567],color='red' ,label = r"$N_{\rm MD}=10$")
plt.fill_between(k_10s[:567], P_k_10[:567], P_k_10[:567]*0+10, color='red',alpha=0.2)

plt.loglog(k_20s[:816],P_k_20[:816], label = r"$N_{\rm MD}=20$", color='green')
plt.fill_between(k_20s[:816], P_k_20[:816], P_k_10[:816]*0+10, color='green',alpha=0.2)

plt.loglog(k_30s[:1351], P_k_30[:1351], label = r"$N_{\rm MD}=30$", color='blue')
plt.fill_between(k_30s[:1351], P_k_30[:1351], P_k_30[:1351]*0+10, color='blue',alpha=0.2)

plt.loglog(k_rd,P_k_rd,color='black')

plt.xlabel(r"$k~[\rm{Mpc}^{-1}]$")
plt.ylabel(r"$\mathcal{P}_\zeta(k)$")
plt.ylim([1e-7,1])
plt.xlim(0.45e13,5.3e22)
plt.legend(ncol=2,bbox_to_anchor=(0.85, 1.23))
plt.show()
```


    
![png](img/output_14_0.png)


---    


### Stiff fluid Dominated (SD) Scenario

```python
k_phys_rad = np.array(PfS.k_rad(M_tot))
_, p_k_rad = PfR.get_P_k_RD(M_tot,constraints.betas_full,delta_c)
```

```python
k_5st, P_k_5s, beta_s5 = PfS.get_P_k_SD(M_tot,5,1,1)
k_10st, P_k_10s, beta_s10 = PfS.get_P_k_SD(M_tot,10,1,1)
k_15st, P_k_15s, beta_s15 = PfS.get_P_k_SD(M_tot,15,1,1)
```

#### Abundance of PBHs ($\beta_{PBH}$) as function of their mass, where $N_{SD}$ is the total number of $e$-folds that the stiff era lasted.

```python
plt.loglog(M_tot,betas_full, color='black',label = "SBB")
plt.loglog(M_tot,beta_s5,label = r"$N_{\rm SD} = 5$", color='red')
plt.loglog(M_tot,beta_s10,label = r"$N_{\rm SD} = 10$", color='green')
plt.loglog(M_tot,beta_s15,label = r"$N_{\rm SD} = 15$", color='blue')
plt.loglog(M_tot,betas_full, color ='black')

plt.xlabel(r"$M_{\rm PBH}~[\rm{g}]$")
plt.ylabel(r"$\beta$")
plt.ylim([1e-40,1e-5])
plt.xlim([1,1e20])

plt.legend(ncol=2,bbox_to_anchor=(0.85, 1.25))


plt.show()
```


    
![png](img/output_16_0.png)
    
#### Constraints on $\mathcal{P}_{\zeta}(k)$ as a function of $k$ for the extended SD scenario.

```python
plt.loglog(k_phys_rad, P_k_rd, color='black', label = "SBB")
plt.fill_between(k_phys_rad, P_k_rd, [10]*len(P_k_rd), color='black',alpha=0.2)

plt.loglog(k_5st[:567],P_k_5s[:567], color='red', label = r"$N_{\rm SD}=5$")
plt.fill_between(k_5st[:567], P_k_5s[:567], P_k_5s[:567]*0+10, color='red',alpha=0.2)

plt.loglog(k_10st[:818],P_k_10s[:818], color='green',label = r"$N_{\rm SD}=10$")
plt.fill_between(k_10st[:820], P_k_10s[:820], P_k_10s[:820]*0+10, color='green',alpha=0.2)

plt.loglog(k_15st[:1320],P_k_15s[:1320], color='blue',label = r"$N_{\rm SD}=15$")
plt.fill_between(k_15st[:1320], P_k_15s[:1320], P_k_15s[:1320]*0+10, color='blue',alpha=0.2)


plt.xlabel(r"$k~[\rm{Mpc}^{-1}]$")
plt.ylabel(r"$\mathcal{P}_\zeta(k)$")
plt.ylim([6.5e-3,2e-2])
plt.xlim([4.3e12,1e26])
plt.ylabel(r"$\mathcal{P}_\zeta(k)$")
plt.legend(ncol=2,bbox_to_anchor=(0.85, 1.25))

plt.show()
```

![png](img/output_17_0.png)
---





