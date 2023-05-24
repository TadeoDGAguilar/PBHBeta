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

The first function we will execute is {py:func}`PBHBeta.functions.put_M_array`, which contains the instruction to generate an array of masses

```{code-block}
PBHBeta.functions.put_M_array()
```


```{note}
The library contains a module called `constraints.py`, where the `return` of {py:func}`PBHBeta.functions.put_M_array` will be stored.
```

