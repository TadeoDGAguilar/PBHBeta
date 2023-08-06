# PBH-Beta
## Warning this project is a beta version

### Authors

[tadeodaguilar](https://www.linkedin.com/in/tadeodaguilar/) & [Luis E. Padilla](https://www.linkedin.com/in/luis-enrique-padilla-albores-052087199/)

[![Build Status](https://app.travis-ci.com/TadeoDGAguilar/PBHBeta.svg?branch=main)](https://app.travis-ci.com/TadeoDGAguilar/PBHBeta)

[![Documentation Status](https://readthedocs.org/projects/pbhbeta/badge/?version=latest)](https://pbhbeta.readthedocs.io/en/latest/?badge=latest)

[![arXiv](https://img.shields.io/badge/arXiv-1234.56789-f9f107.svg)](https://arxiv.org/abs/1234.56789)

## Prerequisites

The betaPBH library requires Python 3.10 or later to be installed on your system.

[Python 3](https://www.python.org/downloads/)

[pip package manager](https://pypi.org/project/pip/): The pip package manager is used to install betaPBH and its dependencies. It should be included with your Python installation by default.

In general, when you install `betaPBH`, the `setup.py` will install all dependences: `matplotlib` (v-3.7.1), `numpy` (v-1.22.4), `scipy` (v-1.10.1).
If this not happend, you need install manually to use `betaPBH`

  1. [Matplotlib](https://matplotlib.org/stable/users/installing/index.html)

  2. [NumPy](https://numpy.org/install/)

  3. [SciPy](https://scipy.org/install/)


**Note:** `betaPBH` runs both in `Python` 2.x and 3.x. However, we highly recommend Python 3.x

# Example 
## Abundances of PBHs with Number of e-folds
```{code}
  from betaPBH import functions, constants, constraints, BfN, BfS
  import matplotlib.pyplot as plt
  import numpy as np
```
```{code}
  functions.put_M_array()
  M_tot = np.array(constraints.M_tot)
```
```{code}
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

![Descripción de la imagen](https://pbhbeta.readthedocs.io/en/latest/_images/BfN.png)

# How to cite us

If you use $\beta$-PBH, please cite its pre-print, arXiv:.

## Regards

We would also like to extend our appreciation to [Encieh Erfani](https://globalyoungacademy.net/eerfani/) and [Juan Carlos Hidalgo](https://www.fis.unam.mx/directorio/540/juan-carlos-strong-hidalgo-strong-cuellar) for their collaboration and valuable contributions to the project. Furthermore, we would like to thanks our colleagues for their assistance, suggestions, and encouragement during the development of the PBHbeta library, with a special mention to [Isidro Gomez](https://igomezv.github.io/).
Lastly, we are grateful to the institutions that have supported this research, enabling us to pursue our scientific investigations and contribute to the field of cosmology. We would like to extend our appreciation to the Universidad Nacional Autónoma de México (UNAM), Instituto de Ciencias Físicas (ICF), International Centre For Theoretical Physics (ICTP), and the Consejo Nacional de Humanidades, Ciencias y Tecnologías (CONAHCYT) for their invaluable support and resources, which have been crucial in carrying out this research.

[![UNAM](docs/img/UNAM.png){.bg-warning w=100px h=130px}](https://www.unam.mx/)
[![ICF](docs/img/ICF.png){.bg-warning w=250px h=130px}](https://www.fis.unam.mx/)
[![ICTP](docs/img/ICTP.svg){.bg-warning w=130px h=130px}](https://www.ictp.it/)
[![CONAHCYT](docs/img/CONAHCYT.svg){.bg-warning w=190px h=130px}](https://conahcyt.mx/)

