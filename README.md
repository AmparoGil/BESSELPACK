# BESSELPACK

**BESSELPACK** is a MATLAB package for the accurate computation of Bessel functions, their zeros, and zeros of cross products of Bessel functions.

The main routines are:

* `bessel` — computation of Bessel functions $J_\nu(x)$, $Y_\nu(x)$, and their first derivatives.
* `besselz` — computation of the zeros of general cylinder functions.
* `besselzd` — computation of the zeros of a generalized derivative combination.
* `crosspz` — computation of the zeros of three cross products of Bessel functions.
* `BESSELPACKdemo` — demonstration examples.

---

## 1. Bessel functions

The routine

```matlab
bessel
```

computes

$$
J_\nu(x), \qquad Y_\nu(x),
$$

and their first derivatives

$$
J_\nu'(x), \qquad Y_\nu'(x).
$$

The implementation combines different numerical methods depending on the parameter region in order to provide accurate and efficient results over a wide range of orders and arguments.

---

## 2. Zeros of general cylinder functions

The routine

```matlab
besselz
```

computes the zeros of the general cylinder function

$$
C_{\nu,\alpha}(x)=\cos(\alpha)J_\nu(x)-
\sin(\alpha)Y_\nu(x).
$$

The zeros can be computed in a prescribed interval $(a,b)$, for given order $\nu$, parameter $\alpha$, and requested accuracy.

The routine returns a vector containing the computed zeros.

---

## 3. Zeros of generalized derivative combinations

The routine

```matlab
besselzd
```

computes the zeros of

$$
xC'_{\nu,\alpha}(x)
+
\gamma C_{\nu,\alpha}(x)=0.
$$

The routine is designed to provide reliable results also in parameter regions where zeros may become very close or nearly degenerate.

In addition to the computed zeros, information is provided for zeros for which the requested accuracy could not be achieved. This allows potentially problematic cases, such as nearly multiple zeros, to be identified.

---

## 4. Zeros of cross products

The routine

```matlab
crosspz
```

computes the zeros of three cross products of Bessel functions. For $\lambda>0, \lambda \neq 1$, the three functions considered are:

### A. Dirichlet–Dirichlet cross product

$$
J_\nu(x)Y_\nu(\lambda x)-
J_\nu(\lambda x)Y_\nu(x).
$$

### B. Neumann–Neumann cross product

$$
J_\nu'(x)Y_\nu'(\lambda x)-
J_\nu'(\lambda x)Y_\nu'(x).
$$

### C. Mixed cross product

$$
J_\nu'(x)Y_\nu(\lambda x)-
J_\nu(\lambda x)Y_\nu'(x).
$$

The cross product to be considered is selected by the input parameter `icho`.

A typical call has the form

```matlab
[xzer] = crosspz(icho,nu,lambda,a,b,eps);
```

where:

| Input    | Description                              |
| -------- | ---------------------------------------- |
| `icho`   | Selects the cross product to be computed |
| `nu`     | Bessel function order $\nu$              |
| `lambda` | Parameter $\lambda>0, \lambda \neq 1$    |
| `a`      | Lower endpoint of the interval           |
| `b`      | Upper endpoint of the interval           |
| `eps`    | Requested relative accuracy              |

The output is:

| Output | Description                          |
| ------ | ------------------------------------ |
| `xzer` | Vector containing the computed zeros |

---

## 5. Numerical validation

The routines have been tested over broad ranges of orders, parameters, and arguments.

For the zero-finding algorithms, the computed zeros are checked using independent numerical procedures. In particular, the number of computed zeros is compared with the number of sign changes of the corresponding function on a sufficiently fine independent grid, providing a check that no zeros have been missed within the prescribed interval.

The computed zeros are also evaluated in the corresponding defining functions to verify that the residuals are consistent with the requested accuracy.

---

## 6. Demonstration

The script

```matlab
BESSELPACKdemo
```

contains examples illustrating the use of the main routines and their application to the computation of Bessel functions and their zeros.

---

## 7. Requirements

* MATLAB
* No additional MATLAB toolboxes are required beyond the functionality used by the package.

---

## 8. Package structure

The main components of BESSELPACK are:

```text
BESSELPACK/
├── bessel.m
├── besselz.m
├── besselzd.m
├── crosspz.m
├── BESSELPACKdemo.m
└── ...
```

---

## 9. Reference

If you use BESSELPACK in academic work, please cite:

> *Numerical Software for Bessel Functions and Associated Values*
> Amparo Gil, Javier Segura, and Nico M. Temme
> Submitted

A BibTeX entry will be provided with the final release.

---

## 10. License

Please see the `LICENSE` file for the terms under which BESSELPACK is distributed.
