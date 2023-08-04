# Requeriments

## Installing dependencies in python

Imperative libraries are:

1. **NumPy:** NumPy is a Python library that provides support for working with multidimensional arrays and high-performance mathematical functions. It is widely used in data science and scientific computing due to its efficiency, linear algebra capabilities, and integration with other Python libraries. It is an essential tool for complex mathematical operations and numerical data manipulation in Python.

To install copy-paste this:

```{block-code}
pip install numpy
```

2. **SciPy:** SciPy is a Python library that extends the capabilities of NumPy and provides additional functionalities for scientific computing. It includes tools for integration, optimization, interpolation, signal and image processing, statistical functions, and more.

To install copy-paste this:

```{block-code}
pip install scipy
```

3. **Matplotlib library:** This library provides support for creating various types of charts and graphs in Python.

To install copy-paste this:

```{block-code}
pip install matplotlib
```

If you prefer, you can use:

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
