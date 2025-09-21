# Structure Growth Numerical Integration

This repository contains Python code to numerically integrate the differential equation for the growth of density perturbations in an expanding universe, comparing different equations of state.

## Overview
The code is simulating the equations included into the document.

In the exmaple here, the code solves the linear growth equation for cosmological perturbations:

$\frac{d^2\delta}{da^2} + \bigl[\frac{d(ln H)}{da} + \frac{3}{a}\bigr]\frac{d\delta}{da} - \bigl[\frac{3}{2}\omega_m(a)\frac{(1-w)(1+3w)}{a^5H(a)^2/H_0^2}\bigr]\delta=0$
d²δ/da² + [d(lnH)/da + 3/a] dδ/da - [3/2 Ω_m(a) (1-w)(1+3w)/(a⁵ H(a)²/H₀²)] δ = 0

It specifically compares:
1. Standard Cold Dark Matter (CDM) with equation of state parameter w = 0
2. The proposed T-field gradient component with w = -1/3 from the paper "Quantum Time as a Scalar Field"

## Requirements

- Python 3.6+
- NumPy
- SciPy
- Matplotlib

Install requirements with:
```bash
pip install numpy scipy matplotlib
python QTSF_Figure_Growth_w.py
```
License
MIT License - feel free to use this code for research or educational purposes.

