============
Introduction
============

The Problem of Spacetime Equivalence
------------------------------------

In the theory of general relativity, the gravitational field is represented by the metric of the spacetime, a rank two symmetric tensor of Lorentzian signature. One of the main features of this theory is that it is totally covariant, in the sense the equations of motion are written in a way that are invariant under the change of coordinates.

While this covariance is a desirable feature of the theory, it comes at a price: the difficulty of knowing whether two metrics are equivalent or not. Given two line elements expressed in different coordinate systems, it's generally very hard to check if they represent the same physical spacetime.

For instance, consider the following two line elements:

.. math::

   ds^2 = -\left(1-\frac{2M}{r} \right)dt^2 + \left(1-\frac{2M}{r} \right)^{-1}dr^2 + r^2 d\theta^2 + r^2 \sin^2\theta d\phi^2

.. math::

   d\tilde{s}^2 = - \frac{b^2}{ax^2}(a x^2 -c)e^{2\tau}d\tau^2 +  \frac{4a^3 x^4}{a x^2 -c} dx^2 +  \frac{a^2 c^2 x^4}{4 -c^2 y^2} dy^2
   +  \frac{a^2 c^2 x^4 y^2}{4(y^2 + z^2)^2} \left( z^2 dy^2 + y^2 dz^2 - 2 zy \,dydz  \right)

The first line element is the well-known Schwarzschild solution in spherical coordinates, which represents the gravitational field outside a static and spherically symmetric distribution of mass. The second line element appears completely different, but both metrics represent the same gravitational field whenever :math:`c = 2M`.

The Cartan-Karlhede Algorithm
-----------------------------

The Cartan-Karlhede algorithm provides a systematic procedure to determine whether two metrics represent the same physical spacetime, by analyzing the functional relations between curvature tensors and their covariant derivatives in appropriately chosen frames.

The algorithm consists of these main steps:

1. Define a constant metric :math:`\eta_{ab}` and find a frame for each metric.
2. Compute curvature components and their derivatives in these frames.
3. For each order of derivative, compute:
   - :math:`\tau_i`: the minimum number of functionally independent components
   - :math:`H_i`: the isotropy group that preserves the form of the curvature and its derivatives
4. Continue until :math:`\tau_q = \tau_{q-1}` and :math:`H_q = H_{q-1}`.
5. Compare the functional relations between the curvature components of both metrics.

The algorithm guarantees that two metrics are equivalent if and only if there exists a frame for each metric such that all the functional relations between the components of the curvature tensor and its derivatives are identical.

About This Implementation
------------------------

This Python implementation provides the tools to represent metrics, compute their curvature tensors, and apply the Cartan-Karlhede algorithm to determine equivalence. The implementation uses SymPy for symbolic tensor calculations, allowing for exact results rather than numerical approximations.

Key features of this implementation include:

- A flexible `Metric` class for representing spacetime metrics
- Support for both orthonormal (Lorentz) frames and null frames
- Calculation of curvature tensors (Riemann, Ricci, Weyl) and their derivatives
- Framework for comparing metrics according to the Cartan-Karlhede algorithm
- Example implementations from the paper, including the Schwarzschild metric in different coordinate systems 