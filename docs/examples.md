# Jupyter notebooks

The following examples are based on the calculations that were performed to obtain our results.

## Primordial Black Holes in non-standard cosmology

First step is import the PBHBeta library and necessary modules, our library containts the following modules: `constanst.py`, `classes.py`, `functions.py` and `constants.py`

```{code-block} python
from PBHBeta import *
import matplotlib.pyplot as plt
import numpy as np
```

Generate an array of masses using the put_M_array function: The initial function to execute is PBHBeta.functions.put_M_array(), which generates an array of masses from a specified minimum and maximum mass value in grams.

```{code-block} python
#functions.put_M_array(constants.M_pl_g, 1e20)
constraints.M_tot = constraints.data_M_tot
M_tot = np.array(constraints.M_tot)
omega = 1/3
```

Once our array of masses in grams has been defined, we can now calculate the constraints for the masses and abundances. Here, we will use all the main functions: `PBHBeta.functions.Betas_BBN`, `PBHBeta.functions.Betas_SD`, `PBHBeta.functions.CMB_AN`, `PBHBeta.functions.Betas_GRB`, `PBHBeta.functions.Betas_Reio`, and `PBHBeta.functions.Betas_LSP`.



```{code-block} python
M_n, betas, M_relic, betas_relic = functions.Betas_DM(M_tot)
betas_bbn = functions.Betas_BBN(M_tot,omega)[1]
betas_sd = functions.Betas_SD(M_tot,omega)[1]
betas_an = functions.Betas_CMB_AN(M_tot,omega)[1]
betas_grb1 = functions.Betas_GRB(M_tot,omega)[2]
betas_grb2 = functions.Betas_GRB(M_tot,omega)[3]
betas_reio = functions.Betas_Reio(M_tot,omega)[1]
betas_lsp = functions.Betas_LSP(M_tot,omega)[1]
betas_full = functions.get_Betas_full(M_tot)
```



```{code-block} python
M_bbn = functions.Betas_BBN(M_tot,omega)[0]
M_an = functions.Betas_CMB_AN(M_tot,omega)[0]
M_grb1 = functions.Betas_GRB(M_tot,omega)[0]
M_grb2 = functions.Betas_GRB(M_tot,omega)[1]
M_reio = functions.Betas_Reio(M_tot,omega)[0]
M_lsp = functions.Betas_LSP(M_tot,omega)[0]
```

```{code-block} python
omega_dm = functions.Omegas_DM(M_tot)
omega_bbn = constraints.Omega_BBN_tot
omega_sd = constraints.Omega_SD_tot
omega_an = constraints.Omega_CMB_AN_tot
omega_grb = constraints.Omega_GRB_tot
omega_reio = constraints.Omega_Reio_tot
omega_lsp = constraints.Omega_LSP_tot
functions.get_Omegas_full(M_tot)
```

Now, we calculate the constraints in the PPS by assuming that
\begin{equation}
\beta = Erfc\left(\frac{\delta_c}{\sqrt{2}\sigma}\right), \ \ \ \ \Rightarrow \ \ \ \ \sigma = \frac{\delta_c}{\sqrt{2}Erfc^{-1}(\beta)},
\end{equation} with $\sigma^2\sim P(k)$ and $\delta_c = 0.41$ in the radiation case.

```{code-block}
delta_c = 0.41
sigma = np.array(functions.inverse_error(betas,delta_c))
sigma_relic = np.array(functions.inverse_error(betas_relic,delta_c))
sigma_bbn = np.array(functions.inverse_error(betas_bbn,delta_c))
sigma_sd = np.array(functions.inverse_error(betas_sd,delta_c))
sigma_an = np.array(functions.inverse_error(betas_an,delta_c))
sigma_grb1 = np.array(functions.inverse_error(betas_grb1,delta_c))
sigma_grb2 = np.array(functions.inverse_error(betas_grb2,delta_c))
sigma_reio = np.array(functions.inverse_error(betas_reio,delta_c))
sigma_lsp = np.array(functions.inverse_error(betas_lsp,delta_c))
sigma_tot = np.array(functions.inverse_error(betas_full,delta_c))
```

### Early matter dominated scenario



```{code-block} python
k_md = KfN.k_MD(M_tot)
sigma_10k, k_10s, beta_10 = KfN.get_P_k_MD(M_tot,10,0,1)
sigma_20k, k_20s, beta_20 = KfN.get_P_k_MD(M_tot,20,0,1)
sigma_30k, k_30s, beta_30 = KfN.get_P_k_MD(M_tot,30,0,1)
```

```{code-block} python
plt.loglog(M_tot,betas_full,"k", label = "SBB", color='black')
plt.loglog(M_tot,beta_10,"k", label = r"$N_{\rm MD}=10$",color='red')
plt.loglog(M_tot,beta_20,"k", label = r"$N_{\rm MD}=20$",color='green')
plt.loglog(M_tot,beta_30,"k", label = r"$N_{\rm MD}=30$",color='blue')
plt.loglog(M_tot,betas_full,"k", color='black')
plt.ylim([1e-30,1])
plt.xlim([1,1e20])
plt.xlabel(r"$M_{\rm PBH}~[\rm{g}]$")
plt.ylabel(r"$\beta$")

plt.legend(ncol=2,bbox_to_anchor=(0.85, 1.23))
plt.show()
```

`out`
```{figure} /img/Beta_EMD_note1.png
:alt: Beta_emd
:class: bg-primary
:width: 500px
:align: center
:caption: Constraints on the abundance of PBHs as a function of their mass. $N_{MD}$ is the total number of e-folds that the MD epoch lasted.
```

```{code-block} python
B = 25/4

plt.loglog(k_md,B*sigma_tot**2,"k",label = "SBB",color='#0000FF')
plt.fill_between(k_md, B*sigma_tot**2, (B*sigma_tot**2)*0+10, color='#0000FF',alpha=0.2)

plt.loglog(k_10s[:567],(B*sigma_10k**2)[:567],"k",color='lime' ,label = r"$N_{\rm MD}=10$")
plt.fill_between(k_10s[:567], (B*sigma_10k**2)[:567], (B*sigma_10k**2)[:567]*0+10, color='lime',alpha=0.2)

plt.loglog(k_20s[:816],(B*sigma_20k**2)[:816],"k", label = r"$N_{\rm MD}=20$", color='red')
plt.fill_between(k_20s[:816], (B*sigma_20k**2)[:816], (B*sigma_20k**2)[:816]*0+10, color='red',alpha=0.2)

plt.loglog(k_30s[:1351],(B*sigma_30k**2)[:1351],"k", label = r"$N_{\rm MD}=30$", color='brown')
plt.fill_between(k_30s[:1351], (B*sigma_30k**2)[:1351], (B*sigma_30k**2)[:1351]*0+10, color='brown',alpha=0.2)

plt.loglog(k_md,B*sigma_tot**2,"k",color='#0000FF')

plt.xlabel(r"$k~[\rm{Mpc}^{-1}]$")
plt.ylabel(r"$\mathcal{P}_\zeta(k)$")
plt.legend(ncol=2,bbox_to_anchor=(0.85, 1.23))
plt.show()
```

#### Constraints on $\mathcal{P}_{\zeta}(k)$ as a function of $k$ for different values of $N_{MD}$.

`out`
```{figure} img/PofK_EMD_note1.png
:alt: Pofk_emd
:class: bg-primary
:width: 500px
:align: center
```

**Downloadable Notebook with this example:** 

---
### Stiff fluid dominated scenario

```{code-block} python
k_phys_rad = np.array(KfS.k_rad(M_tot))
sigma_5st, k_5st, beta_s5 = KfS.get_P_k_SD(M_tot,5,1,1)
sigma_10st, k_10st, beta_s10  = KfS.get_P_k_SD(M_tot,10,1,1)
sigma_15st, k_15st, beta_s15 = KfS.get_P_k_SD(M_tot,15,1,1)
```

```{code-block}
plt.loglog(M_tot,betas_full,"k", color='black',label = "SBB")
plt.loglog(M_tot,betas_s5,"k",label = r"$N_{\rm SFD} = 5$", color='red')
plt.loglog(M_tot,betas_s10,"k",label = r"$N_{\rm SFD} = 10$", color='green')
plt.loglog(M_tot,betas_s15,"k",label = r"$N_{\rm SFD} = 15$", color='blue')
plt.loglog(M_tot,betas_full,"k", color ='black')
plt.xlabel(r"$M_{\rm PBH}~[\rm{g}]$")
plt.ylabel(r"$\beta$")
plt.ylim([1e-40,1e-5])
plt.xlim([1,1e20])
plt.legend(ncol=2,bbox_to_anchor=(0.85, 1.25))
plt.show()
```

#### Abundance of PBHs as function of their mass, where NSD is the total number of $e$-folds that the stiff era lasted

`out`
```{figure} img/Beta_SFD_note2.png
:alt: img3
:class: bg-primary
:width: 500px
:align: center
```



```{code-block}
C = (8/2)**2/4

plt.loglog(k_phys_rad, A*sigma_tot**2, "k", color='black', label = "SBB")
plt.fill_between(k_phys_rad, A*sigma_tot**2, [10]*len(A*sigma_tot**2), color='black',alpha=0.2)

plt.loglog(k_5st[:567],(C*sigma_5st**2)[:567],"r", color='red', label = r"$N_{\rm SD}=5$")
plt.fill_between(k_5st[:567], (C*sigma_5st**2)[:567], (C*sigma_5st**2)[:567]*0+10, color='red',alpha=0.2)

plt.loglog(k_10st[:818],(C*sigma_10st**2)[:818],"g", color='green',label = r"$N_{\rm SD}=10$")
plt.fill_between(k_10st[:820], (C*sigma_10st**2)[:820], (C*sigma_10st**2)[:820]*0+10, color='green',alpha=0.2)

plt.loglog(k_15st[:1320],(C*sigma_15st**2)[:1320],"b", color='blue',label = r"$N_{\rm SD}=15$")
plt.fill_between(k_15st[:1320], (C*sigma_15st**2)[:1320], (C*sigma_15st**2)[:1320]*0+10, color='blue',alpha=0.2)

plt.loglog(k_5st[:567],(C*sigma_5st**2)[:567],"r", color='red')
plt.loglog(k_10st[:819],(C*sigma_10st**2)[:819],"g", color='green')
plt.loglog(k_phys_rad, A*sigma_tot**2, "k", color='black')

plt.xlabel(r"$k~[\rm{Mpc}^{-1}]$")
plt.ylabel(r"$\mathcal{P}_\zeta(k)$")
plt.ylim([6.5e-3,2e-2])
plt.xlim([4.3e12,1e26])
plt.ylabel(r"$\mathcal{P}_\zeta(k)$")
plt.legend(ncol=2,bbox_to_anchor=(0.85, 1.25))

plt.show()
```


#### Constraints on $\mathcal{P}_{\zeta}(k)$ as a function of $k$ for the extended SD scenario.

`out`
```{figure} img/PofK_SFD_note2.png
:alt: img4
:class: bg-primary
:width: 500px
:align: center
```

**Downloadable Notebook with this example:** 





