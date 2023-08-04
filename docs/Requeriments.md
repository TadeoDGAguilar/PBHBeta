# Requeriments

## Installing dependencies in python

Imperative libraries are:

1. **NumPy:** NumPy is a Python library that provides support for working with multidimensional arrays and high-performance mathematical functions. It is widely used in data science and scientific computing due to its efficiency, linear algebra capabilities, and integration with other Python libraries. It is an essential tool for complex mathematical operations and numerical data manipulation in Python.

```{block-code}
pip install numpy
```

2. **SciPy:** SciPy is a Python library that extends the capabilities of NumPy and provides additional functionalities for scientific computing. It includes tools for integration, optimization, interpolation, signal and image processing, statistical functions, and more.

```{block-code}
pip install scipy
```

3. **Matplotlib library:** The Matplotlib library is also required by betaPBH and can be installed using pip. This library provides support for creating various types of charts and graphs in Python.

```{block-code}
pip install matplotlib
```

```{tip}
All in one copy-paste line:
---
emphasize-lines: 1
---
pip install numpy scipy matplotlib
```

## Using Google Colab

For easy and immediate use, we recommend using Google Colab. Within it, it is necessary to install some libraries and `PBHBeta` in the first cells.

```{code}
%pip install numpy matplotlib scipy
!pip install -e git+https://github.com/TadeoDGAguilar/PBHBeta#egg=PBHBeta
%cd /content/src/pbhbeta
```

1. **Python 3.6 or later:** The betaPBH library requires Python 3.6 or later to be installed on your system.

```{note}
This code runs in any Python 2x and 3x. However, we highly recommend Python 3x. In the GitHub repository and PyPI package for PBHBeta contains the basic libraries to run `betaPBH`, but this dependences need import into notebook or in the console to functions be available.
```

1. **pip package manager:** The pip package manager is used to install betaPBH and its dependencies. It should be included with your Python installation by default.
