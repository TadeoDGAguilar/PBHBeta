# Jupyter notebooks

## Primordial Black Holes in non-standard cosmology


### Early matter dominated scenario

```{code-block} python
:lineno-start: 1
from PBHBeta import *
import matplotlib.pyplot as plt
import numpy as np
```

```{code-block} python
k_md = KfN.k_MD(M_tot)
sigma_10k, k_10s, beta_10 = KfN.get_k_Nreh(M_tot,10,0,1)
sigma_20k, k_20s, beta_20 = KfN.get_k_Nreh(M_tot,20,0,1)
sigma_30k, k_30s, beta_30 = KfN.get_k_Nreh(M_tot,30,0,1)
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

Constraints on the abundance of PBHs as a function of their mass.
$N_{MD}$ is the total number of e-folds that the MD epoch lasted.

`out`
```{figure} /img/Beta_EMD_note1.png
:alt: Beta_emd
:class: bg-primary
:width: 500px
:align: center
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

Constraints on $\mathcal{P}_{\zeta}(k)$ as a function of $k$ for different values of $N_{MD}$.

`out`
```{figure} img/PofK_EMD_note1.png
:alt: Pofk_emd
:class: bg-primary
:width: 500px
:align: center
```

**Downloadable Notebook with this example:** 


### Stiff fluid dominated scenario

```{code-block} python
:lineno-start: 1
from PBHBeta import *
import matplotlib.pyplot as plt
import numpy as np
```

Abundance of PBHs as function of their mass, where NSD is the total number of $e$-folds that the stiff era lasted

`out`
```{figure} img/Beta_SFD_note2.png
:alt: Pofk_emd
:class: bg-primary
:width: 500px
:align: center
```

Constraints on $\mathcal{P}_{\zeta}(k)$ as a function of $k$ for the extended SD scenario.

`out`
```{figure} img/PofK_SFD_note2.png
:alt: Pofk_emd
:class: bg-primary
:width: 500px
:align: center
```

**Downloadable Notebook with this example:** 





