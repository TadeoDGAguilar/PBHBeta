# Jupyter notebooks
## Obtain Betas in radiation dominated epoch


```{code-block} python
:lineno-start: 1
from betaPBH import functions, constants, constraints, BfN, BfS
import matplotlib.pyplot as plt
import numpy as np
```

```{code-block} python
:lineno-start: 4
functions.put_M_array(0.0123)
```

```{code-block} python
:lineno-start: 5
M_tot = np.array(constraints.M_tot)
omega=1/3

funtions.Betas_DM(M_tot)
funtions.Betas_BBN(M_tot,omega)
funtions.Betas_CMB_AN(M_tot,omega)
funtions.Betas_GRB(M_tot,omega)
funtions.Betas_Reio(M_tot,omega)
funtions.Betas_SD(M_tot,omega)
funtions.Betas_LSP(M_tot)
```

```{code-block} python
:lineno-start: 15
funtions.get_Betas_full(M_tot)
```

```{code-block} python
:lineno-start: 16
plt.loglog(M_tot,constraints.betas_full,"k:")
plt.ylim([1e-30,1e-13])
plt.xlim([1e-2,1e20])
plt.show()
```
`out:`
```{figure} img/betas_full.png
:alt: fishy
:class: bg-primary
:width: 500px
:align: center
```

## Computing abundances with efolds

$\beta_{PBH} = \beta(N_{\rm reh}, \omega, \gamma_{\rm reh})$

```{code-block} python
:lineno-start: 1
from betaPBH import functions, constants, constraints, BfN, BfS
import matplotlib.pyplot as plt
import numpy as np
```

```{code-block} python
:lineno-start: 4
functions.put_M_array(0.0123)
M_tot = np.array(constraints.M_tot)
```

```{code-block}
:lineno-start: 6
plt.loglog(M_tot,BfN.get_betas_reh_tot(10,0,1),label = r"$N_{\rm reh}=10$")
plt.loglog(M_tot,BfN.get_betas_reh_tot(20,0,1),label = r"$N_{\rm reh}=20$")
plt.loglog(M_tot,BfN.get_betas_reh_tot(30,0,1),label = r"$N_{\rm reh}=30$")
plt.ylim([1e-30,1])
plt.xlim([1,1e20])
plt.xlabel(r"$M_{\rm PBH}~[\rm{g}]$")
plt.ylabel(r"$\beta$")
plt.legend(ncol=2,bbox_to_anchor=(0.85, 1.5))
plt.show()
```

`out:`
```{figure} img/BfN.png
:alt: fishy
:class: bg-primary
:width: 500px
:align: center
```


## Computing abundances with Stiff

$\beta_{PBH} = \beta(N_{\rm stiff}, \omega, \gamma_{\rm stiff})$

```{code-block} python
:lineno-start: 1
from betaPBH import functions, constants, constraints, BfN, BfS
import matplotlib.pyplot as plt
import numpy as np
```

```{code-block} python
:lineno-start: 4
functions.put_M_array(0.0123)
M_tot = np.array(constraints.M_tot)
```

```{code-block}
:lineno-start: 6
plt.loglog(M_tot,BfS.get_betas_stiff_tot(5,1.,1.),"k",label = r"$N_{\rm stiff} = 5$",alpha = 0.8)
plt.loglog(M_tot,BfS.get_betas_stiff_tot(10,1.,1.),"k",label = r"$N_{\rm stiff} = 10$",alpha = 0.8)
plt.loglog(M_tot,BfS.get_betas_stiff_tot(15,1.,1.),"k",label = r"$N_{\rm stiff} = 15$",alpha = 0.8)
plt.xlabel(r"$M_{\rm PBH}~[\rm{g}]$")
plt.ylabel(r"$\beta$")
plt.ylim([1e-40,1e-5])
plt.xlim([1,1e20])
plt.legend(ncol=2,bbox_to_anchor=(0.85, 1.5))
plt.show()
```
`out:`
```{figure} img/BfS.png
:alt: fishy
:class: bg-primary
:width: 500px
:align: center
```

