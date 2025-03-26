===================
Mathematical Details
===================

This section provides the mathematical details of the Cartan-Karlhede algorithm.

The Equivalence Problem
----------------------

In general relativity, two spacetimes are considered equivalent if there exists a diffeomorphism between them that preserves the metric structure. In other words, if two metrics represent the same physical spacetime in different coordinate systems, there should exist a coordinate transformation between them.

However, finding such a transformation directly is often extremely difficult. The Cartan-Karlhede algorithm provides a systematic approach to determine whether two metrics are equivalent by analyzing their curvature properties.

Tensors and Frames
-----------------

Given a metric :math:`g_{\mu\\nu}`, we can compute:

- The Christoffel symbols:

.. math::

   \Gamma^{\lambda}_{\mu\\nu} = \\frac{1}{2} g^{\lambda\\sigma} (\\partial_{\mu} g_{\\nu\\sigma} + \\partial_{\\nu} g_{\\mu\\sigma} - \\partial_{\\sigma} g_{\\mu\\nu})

- The Riemann tensor:

.. math::

   R^{\\alpha}_{\\beta\\mu\\nu} = \\partial_{\\nu} \\Gamma^{\\alpha}_{\\beta\\mu} - \\partial_{\\mu} \\Gamma^{\\alpha}_{\\beta\\nu} + \\Gamma^{\\alpha}_{\\sigma\\nu} \\Gamma^{\\sigma}_{\\beta\\mu} - \\Gamma^{\\alpha}_{\\sigma\\mu} \\Gamma^{\\sigma}_{\\beta\\nu}

- The Ricci tensor:

.. math::

   R_{\\mu\\nu} = R^{\\alpha}_{\\mu\\alpha\\nu}

- The Ricci scalar:

.. math::

   R = g^{\\mu\\nu} R_{\\mu\\nu}

- The Weyl tensor (in dimension ≥ 3):

.. math::

   C_{\\mu\\nu\\rho\\sigma} = R_{\\mu\\nu\\rho\\sigma} - \\frac{1}{n-2}(g_{\\mu\\rho}R_{\\sigma\\nu} - g_{\\mu\\sigma}R_{\\rho\\nu} - g_{\\nu\\rho}R_{\\sigma\\mu} + g_{\\nu\\sigma}R_{\\rho\\mu}) + \\frac{R}{(n-1)(n-2)}(g_{\\mu\\rho}g_{\\sigma\\nu} - g_{\\mu\\sigma}g_{\\rho\\nu})

We can also introduce a frame :math:`e_a^\\mu` that transforms the metric to a constant form :math:`\\eta_{ab}`:

.. math::

   g_{\\mu\\nu} e_a^\\mu e_b^\\nu = \\eta_{ab}

In this frame, we can express the components of the Riemann tensor as:

.. math::

   R_{abcd} = R_{\\mu\\nu\\rho\\sigma} e_a^\\mu e_b^\\nu e_c^\\rho e_d^\\sigma

The Cartan-Karlhede Algorithm
----------------------------

The algorithm consists of the following steps:

1. Choose a constant frame metric :math:`\\eta_{ab}` (e.g., Minkowski).

2. For each metric, find a frame :math:`e_a^\\mu` such that:

   .. math::

      g_{\\mu\\nu} e_a^\\mu e_b^\\nu = \\eta_{ab}

3. Compute the components of the Riemann tensor and its covariant derivatives in the frame.

4. Define:
   - :math:`t_i`: The number of functionally independent components in the set :math:`\\{R_{abcd}, \\nabla_e R_{abcd}, \\ldots, \\nabla_{e_1} \\cdots \\nabla_{e_i} R_{abcd}\\}`
   - :math:`\\tau_i`: The minimum value of :math:`t_i` over all possible frame choices
   - :math:`H_i`: The isotropy group preserving the form of the curvature and its derivatives up to order :math:`i`
   
5. Continue computing derivatives until:
   - :math:`\\tau_q = \\tau_{q-1}` and :math:`H_q = H_{q-1}` for some :math:`q`

6. Compare the functional relations between the components of the curvature tensor and its derivatives for both metrics.

7. The metrics are equivalent if and only if there exists a frame for each metric such that all the functional relations are identical.

Special Case: Null Frames
------------------------

For 4D Lorentzian spacetimes, it's often useful to use null frames, particularly for classification via the Petrov type.

A null frame consists of a pair of real null vectors :math:`\\ell` and :math:`n` and a pair of complex conjugate null vectors :math:`m` and :math:`\\bar{m}` satisfying:

.. math::

   \\ell \\cdot \\ell = n \\cdot n = m \\cdot m = \\bar{m} \\cdot \\bar{m} = 0

.. math::

   \\ell \\cdot m = \\ell \\cdot \\bar{m} = n \\cdot m = n \\cdot \\bar{m} = 0

.. math::

   \\ell \\cdot n = -1, \\quad m \\cdot \\bar{m} = 1

In such a frame, the Weyl tensor can be characterized by five complex scalars:

.. math::

   \\Psi_0 = C_{\\mu\\nu\\rho\\sigma} \\ell^{\\mu} m^{\\nu} \\ell^{\\rho} m^{\\sigma}

.. math::

   \\Psi_1 = C_{\\mu\\nu\\rho\\sigma} \\ell^{\\mu} n^{\\nu} \\ell^{\\rho} m^{\\sigma}

.. math::

   \\Psi_2 = C_{\\mu\\nu\\rho\\sigma} \\ell^{\\mu} m^{\\nu} \\bar{m}^{\\rho} n^{\\sigma}

.. math::

   \\Psi_3 = C_{\\mu\\nu\\rho\\sigma} \\ell^{\\mu} n^{\\nu} \\bar{m}^{\\rho} n^{\\sigma}

.. math::

   \\Psi_4 = C_{\\mu\\nu\\rho\\sigma} n^{\\mu} \\bar{m}^{\\nu} n^{\\rho} \\bar{m}^{\\sigma}

The pattern of vanishing of these scalars determines the Petrov type of the spacetime.

Implementation Details
--------------------

In our implementation, we:

1. Represent the metric and compute its associated tensors using SymPy for symbolic calculations.

2. Implement both Lorentzian frames and null frames.

3. Compute the curvature components and their derivatives in these frames.

4. Analyze the functional relations between components (with some simplifications).

5. Compare metrics based on their curvature properties.

For practical reasons, we make some simplifications:

- We approximate :math:`\\tau_i` as :math:`t_i` (computing the minimum over all frames is computationally intensive).
- We use a simplified version of the isotropy group calculation.
- For covariant derivatives, we use a simplified approach that works well for coordinate-aligned frames.

For more details on the mathematical foundations of the algorithm, see the references section. 