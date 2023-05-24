# Quickstar

One way to verify that the betaPBH library has been installed correctly is to immediately import the library and all its modules.

```{code-block}
from PBHBeta import *
```

Additionally, it is recommended for the user to import numpy and matplotlib.pyplot, so the main code header would be written as follows:

```{code-block}
from PBHBeta import *
import matplotlib.pyplot as plt
import numpy as np
```

The first function we will execute is {py:func}`PBHBeta.functions.put_M_array`, which contains the instruction to generate an array of masses from a maximum value of mass in grams.

```{code-block}
PBHBeta.functions.put_M_array(1e20)
```

```{note}
The library contains a module called `constraints.py`, where the `return` of {py:func}`PBHBeta.functions.put_M_array` will be stored.
```

You can retrieve the values stored in `constraint.M_tot` and store them in a variable to perform operations within the PBHBeta functions. 

Here's an example:

```{code-block} python
:lineno-start: 1
M_tot = constraints.M_tot
```

To obtain the abundances (betas) considering the DM constraints, we will use `M_tot` into {py:func}`PBHBeta.functions.Betas_DM()`

```{code-block} python
:lineno-start: 2
functions.Betas_DM(M_tot)
```

```{note}
You can find in the documentation of the {py:func}`PBHBeta.functions.Betas_DM(M_tot)` that this function stores the values of total abundance in an array called `betas_DM_tot` within the `constraints.py` module. On the other hand, the function returns 4 arrays named:`M_n`, `betas`, `M_relic`, `betas_relic`.

If the user wants to manipulate those arrays, they can simply invoke `functions.Betas_DM(M_tot)[number]`, where `number` can be 0, 1, 2, or 3, corresponding to each output array respectively.
```

It is possible to plot these outputs using matplotlib.


```{code-block} python
:lineno-start: 3
plt.loglog(functions.Betas_DM(M_tot)[0],functions.Betas_DM(M_tot)[1], "k:",label=r"$\Omega_{DM}$")
plt.fill_between(functions.Betas_DM(M_tot)[0],functions.Betas_DM(M_tot)[1], functions.Betas_DM(M_tot)[1]*0+10, color='grey',alpha=0.2)
plt.loglog(functions.Betas_DM(M_tot)[2],functions.Betas_DM(M_tot)[3], "k:",label=r"$\Omega_{DM}$")
plt.fill_between(functions.Betas_DM(M_tot)[2],functions.Betas_DM(M_tot)[3], functions.Betas_DM(M_tot)[3]*0+10, color='grey',alpha=0.2)
plt.xlabel(r"$M_{\rm PBH}~[\rm{g}]$")
plt.ylabel(r"$\beta$")
plt.xlim([1,1e20])
plt.ylim([1e-30,1])

plt.legend()
plt.show()
```

`out:`
```{figure} img/IMG_DM.png
:alt: fishy
:class: bg-primary
:width: 500px
:align: center
```

If you successfully obtained the example image, congratulations! Now you are ready to use `PBHBeta`.
