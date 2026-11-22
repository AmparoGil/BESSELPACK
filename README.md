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

```matlab
[J,Y,Jp,Yp,ier] = bessel(nu,x);
```

| Input | Description |
| ----- | ----------- |
| `nu`  | Order $\nu$ (real) |
| `x`   | Argument $x>0$ |

The error flag is `ier=0` for a successful computation, `ier=1` in case of underflow/overflow and `ier=2` if $x\le 0$. When `ier=1` and $\nu>0$, `J=Jp=0`, `Y=-realmax` and `Yp=realmax`.

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

```matlab
xzer = besselz(a,b,nu,alpha,tol);
```

The routine returns a vector containing the computed zeros in $(a,b)$, $a>0$.

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

```matlab
[xzer,xzbad] = besselzd(a,b,nu,gamma,alpha,tol);
```

In addition to the computed zeros `xzer`, the zeros for which the requested accuracy could not be achieved are returned in `xzbad` (they are also included in `xzer`). This allows potentially problematic cases, such as nearly multiple zeros, to be identified.

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
[xzer,xzbad] = crosspz(icho,nu,lambda,a,b,tol);
```

where:

| Input    | Description                              |
| -------- | ---------------------------------------- |
| `icho`   | Selects the cross product: 1 (A), 2 (B) or 3 (C) |
| `nu`     | Bessel function order $\nu$              |
| `lambda` | Parameter $\lambda>0, \lambda \neq 1$    |
| `a`      | Lower endpoint of the interval           |
| `b`      | Upper endpoint of the interval           |
| `tol`    | Requested relative accuracy (the routine uses `max(tol,10*eps/(1-lambda))`, with `lambda` replaced by `1/lambda` if `lambda>1`) |

The output is:

| Output | Description                          |
| ------ | ------------------------------------ |
| `xzer`  | Vector containing the computed zeros |
| `xzbad` | Zeros (also in `xzer`) for which the requested accuracy was not reached |

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

* MATLAB (tested with R2024b). No additional toolboxes are required.

---

## 8. Package structure

The package contains the following files:

```text
BESSELPACK/
├── bessel.m
├── besselz.m
├── besselzd.m
├── crosspz.m
├── BESSELPACKdemo.m
├── README.md
├── CITATION.cff
└── LICENSE
```

---

## 9. Reference

If you use BESSELPACK in academic work, please cite:

> *Numerical Software for Bessel Functions and Associated Values*
> Amparo Gil, Javier Segura, and Nico M. Temme
> Submitted to Computer Physics Communications

A BibTeX entry will be provided with the final release.

---

## 10. License

Please see the `LICENSE` file for the terms under which BESSELPACK is distributed.
