# Requeriments

## Installing dependencies in python

Imperative libraries are:

1. **Python 3.6 or later:** The betaPBH library requires Python 3.6 or later to be installed on your system.

```{note}
This code runs in any Python 2x and 3x. However, we highly recommend Python 3x. In the GitHub repository and PyPI package for PBHBeta contains the basic libraries to run `betaPBH`, but this dependences need import into notebook or in the console to functions be available.
```

2. **pip package manager:** The pip package manager is used to install betaPBH and its dependencies. It should be included with your Python installation by default.

3. **NumPy and SciPy libraries:** The NumPy and SciPy libraries are required by betaPBH and can be installed using pip. These libraries provide support for numerical computations and scientific computing in Python.

4. **Matplotlib library:** The Matplotlib library is also required by betaPBH and can be installed using pip. This library provides support for creating various types of charts and graphs in Python.

To install 

```{block-code}
sudo pip install numpy matplotlib scipy
```

## Using Google Colab

For easy and immediate use, we recommend using Google Colab. Within it, it is necessary to install some libraries and `PBHBeta` in the first cells.

```{code}
%pip install numpy matplotlib scipy
!pip install -e git+https://github.com/TadeoDGAguilar/PBHBeta#egg=PBHBeta
%cd /content/src/pbhbeta
```

