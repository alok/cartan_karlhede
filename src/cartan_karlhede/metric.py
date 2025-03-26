"""Metric and frame definitions for the Cartan-Karlhede algorithm."""

from typing import Dict, List, Tuple, Union, Optional, Callable
import numpy as np
import sympy as sp
from sympy import Matrix, Symbol, Function, diff, simplify, solve


class Metric:
    """Class representing a spacetime metric.

    The metric is a symmetric rank 2 tensor with Lorentzian signature that defines
    the geometry of the spacetime in general relativity. In local coordinates, the
    metric is represented by a matrix of components :math:`g_{\mu\nu}`.

    According to the Cartan-Karlhede algorithm, the equivalence of two metrics can be
    determined by comparing the functional relations between the components of the
    curvature tensor and its covariant derivatives in appropriate frames.

    Attributes:
        components (Matrix): Symbolic matrix of metric components :math:`g_{\mu\nu}`
        coordinates (List[Symbol]): List of coordinate symbols
        name (str): Optional name for the metric
        dim (int): Dimension of the spacetime
        inverse (Matrix): Inverse metric components :math:`g^{\mu\nu}`

    References:
        Mergulhão, T. M., & Batista, C. (2023). "On the Equivalence of Spacetimes,
        the Cartan-Karlhede Algorithm". Section 2.
    """

    def __init__(
        self,
        components: Matrix,
        coordinates: List[Symbol],
        name: str = "Unnamed Metric",
    ):
        """Initialize a metric.

        Args:
            components: Symbolic matrix of metric components :math:`g_{\mu\nu}`
            coordinates: List of coordinate symbols
            name: Optional name for the metric
        """
        self.components = components
        self.coordinates = coordinates
        self.name = name
        self.dim = len(coordinates)

        # Validate dimensions
        if components.shape != (self.dim, self.dim):
            raise ValueError(
                f"Metric components must be a {self.dim}x{self.dim} matrix"
            )

        # Compute inverse metric
        self.inverse = components.inv()

        # Compute Christoffel symbols, Riemann tensor, etc. on demand
        self._christoffel_symbols = None
        self._riemann_tensor = None
        self._ricci_tensor = None
        self._ricci_scalar = None
        self._weyl_tensor = None

    def christoffel_symbols(self) -> Dict[Tuple[int, int, int], sp.Expr]:
        """Compute the Christoffel symbols of the metric.

        The Christoffel symbols are defined as:

        .. math::
            \Gamma^{\lambda}_{\mu\\nu} = \\frac{1}{2} g^{\lambda\\sigma}
            (\\partial_{\mu} g_{\\nu\\sigma} + \\partial_{\\nu} g_{\\mu\\sigma} -
            \\partial_{\\sigma} g_{\\mu\\nu})

        Returns:
            Dictionary mapping (μ, ν, λ) → Γ^λ_μν

        References:
            Mergulhão, T. M., & Batista, C. (2023). "On the Equivalence of Spacetimes,
            the Cartan-Karlhede Algorithm". Equation (2).
        """
        if self._christoffel_symbols is not None:
            return self._christoffel_symbols

        g = self.components
        g_inv = self.inverse
        coords = self.coordinates
        n = self.dim

        christoffel = {}

        for l in range(n):
            for i in range(n):
                for j in range(n):
                    christoffel[(i, j, l)] = 0
                    for k in range(n):
                        term = simplify(
                            g_inv[l, k]
                            * (
                                diff(g[i, k], coords[j])
                                + diff(g[j, k], coords[i])
                                - diff(g[i, j], coords[k])
                            )
                            / 2
                        )
                        christoffel[(i, j, l)] += term
                    christoffel[(i, j, l)] = simplify(christoffel[(i, j, l)])

        self._christoffel_symbols = christoffel
        return christoffel

    def riemann_tensor(self) -> Dict[Tuple[int, int, int, int], sp.Expr]:
        """Compute the Riemann curvature tensor of the metric.

        The Riemann curvature tensor is defined as:

        .. math::
            R^{\\alpha}_{\\beta\\mu\\nu} = \\partial_{\\nu} \\Gamma^{\\alpha}_{\\beta\\mu} -
            \\partial_{\\mu} \\Gamma^{\\alpha}_{\\beta\\nu} +
            \\Gamma^{\\alpha}_{\\sigma\\nu} \\Gamma^{\\sigma}_{\\beta\\mu} -
            \\Gamma^{\\alpha}_{\\sigma\\mu} \\Gamma^{\\sigma}_{\\beta\\nu}

        Returns:
            Dictionary mapping (α, β, μ, ν) → R^α_βμν

        References:
            Mergulhão, T. M., & Batista, C. (2023). "On the Equivalence of Spacetimes,
            the Cartan-Karlhede Algorithm". Equation (3).
        """
        if self._riemann_tensor is not None:
            return self._riemann_tensor

        coords = self.coordinates
        n = self.dim

        # Get Christoffel symbols
        gamma = self.christoffel_symbols()

        riemann = {}

        for a in range(n):
            for b in range(n):
                for m in range(n):
                    for n_idx in range(n):
                        # R^a_bμν = ∂_ν Γ^a_bμ - ∂_μ Γ^a_bν + Γ^a_σν Γ^σ_bμ - Γ^a_σμ Γ^σ_bν
                        term1 = diff(gamma[(b, m, a)], coords[n_idx])
                        term2 = diff(gamma[(b, n_idx, a)], coords[m])

                        term3 = 0
                        term4 = 0

                        for s in range(n):
                            term3 += gamma[(s, n_idx, a)] * gamma[(b, m, s)]
                            term4 += gamma[(s, m, a)] * gamma[(b, n_idx, s)]

                        riemann[(a, b, m, n_idx)] = simplify(
                            term1 - term2 + term3 - term4
                        )

        self._riemann_tensor = riemann
        return riemann

    def ricci_tensor(self) -> Dict[Tuple[int, int], sp.Expr]:
        """Compute the Ricci tensor of the metric.

        Returns:
            Dictionary mapping (μ, ν) -> R_μν
        """
        if self._ricci_tensor is not None:
            return self._ricci_tensor

        n = self.dim

        # Get Riemann tensor
        riemann = self.riemann_tensor()

        ricci = {}

        for mu in range(n):
            for nu in range(n):
                ricci[(mu, nu)] = 0
                for alpha in range(n):
                    ricci[(mu, nu)] += riemann[(alpha, mu, alpha, nu)]
                ricci[(mu, nu)] = simplify(ricci[(mu, nu)])

        self._ricci_tensor = ricci
        return ricci

    def ricci_scalar(self) -> sp.Expr:
        """Compute the Ricci scalar of the metric.

        Returns:
            The Ricci scalar R
        """
        if self._ricci_scalar is not None:
            return self._ricci_scalar

        g_inv = self.inverse
        ricci = self.ricci_tensor()
        n = self.dim

        r_scalar = 0
        for i in range(n):
            for j in range(n):
                r_scalar += g_inv[i, j] * ricci[(i, j)]

        self._ricci_scalar = simplify(r_scalar)
        return self._ricci_scalar

    def weyl_tensor(self) -> Dict[Tuple[int, int, int, int], sp.Expr]:
        """Compute the Weyl tensor of the metric (in dimension ≥ 3).

        Returns:
            Dictionary mapping (μ, ν, ρ, σ) -> C_μνρσ
        """
        if self._weyl_tensor is not None:
            return self._weyl_tensor

        if self.dim < 3:
            raise ValueError("Weyl tensor is only defined for dimension ≥ 3")

        n = self.dim
        g = self.components

        # Get Riemann tensor components (lowered indices)
        riemann_low = {}
        riemann = self.riemann_tensor()

        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        riemann_low[(i, j, k, l)] = 0
                        for a in range(n):
                            riemann_low[(i, j, k, l)] += g[i, a] * riemann[(a, j, k, l)]
                        riemann_low[(i, j, k, l)] = simplify(riemann_low[(i, j, k, l)])

        # Get Ricci tensor and scalar
        ricci = self.ricci_tensor()
        r_scalar = self.ricci_scalar()

        # Compute Weyl tensor
        weyl = {}
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        # C_μνρσ = R_μνρσ - (g_μ[ρ R_σ]ν - g_ν[ρ R_σ]μ) + (R/(n-1)(n-2)) g_μ[ρ g_σ]ν
                        term1 = riemann_low[(i, j, k, l)]

                        term2 = (
                            g[i, k] * ricci[(l, j)]
                            - g[i, l] * ricci[(k, j)]
                            - g[j, k] * ricci[(l, i)]
                            + g[j, l] * ricci[(k, i)]
                        ) / (n - 2)

                        term3 = (
                            r_scalar
                            * (g[i, k] * g[j, l] - g[i, l] * g[j, k])
                            / ((n - 1) * (n - 2))
                        )

                        weyl[(i, j, k, l)] = simplify(term1 - term2 + term3)

        self._weyl_tensor = weyl
        return weyl


class LorentzFrame:
    """Class representing a Lorentz frame for a metric."""

    def __init__(
        self, metric: Metric, frame_vectors: List[Matrix], eta: Optional[Matrix] = None
    ):
        """Initialize a Lorentz frame.

        Args:
            metric: The metric for which this is a frame
            frame_vectors: List of frame vector fields as column matrices
            eta: The constant frame metric (default: Minkowski)
        """
        self.metric = metric
        self.vectors = frame_vectors
        self.dim = metric.dim

        # Validate dimensions
        if len(frame_vectors) != self.dim:
            raise ValueError(f"Frame must have {self.dim} vectors")

        for v in frame_vectors:
            if v.shape != (self.dim, 1):
                raise ValueError(f"Frame vectors must be {self.dim}x1 column matrices")

        # Default to Minkowski metric if not provided
        if eta is None:
            eta = Matrix.diag(-1, *[1] * (self.dim - 1))

        self.eta = eta

        # Verify that the frame is orthonormal with respect to eta
        self._validate_frame()

        # Store frame components
        self.frame_components = Matrix.hstack(*frame_vectors)

        # Computed properties
        self._curvature_components = None
        self._curvature_derivatives = {}

    def _validate_frame(self):
        """Verify that the frame vectors satisfy the required inner products."""
        g = self.metric.components
        eta = self.eta

        for i, ei in enumerate(self.vectors):
            for j, ej in enumerate(self.vectors):
                # Compute inner product e_i · e_j using the metric
                inner_product = 0
                for mu in range(self.dim):
                    for nu in range(self.dim):
                        inner_product += ei[mu, 0] * ej[nu, 0] * g[mu, nu]

                # This should equal η_ij
                if not simplify(inner_product - eta[i, j]).is_zero:
                    raise ValueError(
                        f"Frame vectors do not satisfy orthonormality: "
                        f"e_{i} · e_{j} = {inner_product} ≠ {eta[i, j]}"
                    )

    def dual_frame(self) -> List[Matrix]:
        """Compute the dual frame 1-forms.

        Returns:
            List of dual frame 1-forms as row matrices
        """
        # For an orthonormal frame, the dual is just e^a_μ = η^ab g_μν e^ν_b
        dual_forms = []
        g = self.metric.components
        eta_inv = self.eta.inv()

        for a in range(self.dim):
            form = Matrix.zeros(1, self.dim)
            for b in range(self.dim):
                for mu in range(self.dim):
                    for nu in range(self.dim):
                        form[0, mu] += (
                            eta_inv[a, b] * g[mu, nu] * self.vectors[b][nu, 0]
                        )
            dual_forms.append(simplify(form))

        return dual_forms

    def curvature_components(self) -> Dict[Tuple[int, int, int, int], sp.Expr]:
        """Compute the components of the Riemann tensor in the frame.

        Returns:
            Dictionary mapping (a, b, c, d) -> R_abcd
        """
        if self._curvature_components is not None:
            return self._curvature_components

        # Get Riemann tensor components (lowered indices)
        riemann = self.metric.riemann_tensor()
        riemann_low = {}
        g = self.metric.components

        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        riemann_low[(i, j, k, l)] = 0
                        for a in range(self.dim):
                            riemann_low[(i, j, k, l)] += g[i, a] * riemann[(a, j, k, l)]
                        riemann_low[(i, j, k, l)] = simplify(riemann_low[(i, j, k, l)])

        # Project onto the frame
        frame_components = {}
        dual_forms = self.dual_frame()

        for a in range(self.dim):
            for b in range(self.dim):
                for c in range(self.dim):
                    for d in range(self.dim):
                        component = 0
                        for mu in range(self.dim):
                            for nu in range(self.dim):
                                for rho in range(self.dim):
                                    for sigma in range(self.dim):
                                        # R_abcd = R_μνρσ e^μ_a e^ν_b e^ρ_c e^σ_d
                                        component += (
                                            riemann_low[(mu, nu, rho, sigma)]
                                            * dual_forms[a][0, mu]
                                            * dual_forms[b][0, nu]
                                            * dual_forms[c][0, rho]
                                            * dual_forms[d][0, sigma]
                                        )

                        frame_components[(a, b, c, d)] = simplify(component)

        self._curvature_components = frame_components
        return frame_components

    def curvature_derivative(self, order: int = 1) -> Dict[Tuple, sp.Expr]:
        """Compute the components of the derivative of the Riemann tensor in the frame.

        Args:
            order: The order of the derivative to compute

        Returns:
            Dictionary mapping multi-indices to derivative components
        """
        if order in self._curvature_derivatives:
            return self._curvature_derivatives[order]

        if order == 0:
            return self.curvature_components()

        if order > 1 and (order - 1) not in self._curvature_derivatives:
            # Recursively compute previous order
            self.curvature_derivative(order - 1)

        # For simplicity in this implementation, we'll compute the covariant derivative
        # with respect to the coordinate basis and then project onto the frame
        # (This is a simplified version - a full implementation would involve connection forms)

        prev_order = self._curvature_derivatives.get(
            order - 1, self.curvature_components()
        )
        coords = self.metric.coordinates

        if order == 1:
            # First derivative
            derivatives = {}
            for a in range(self.dim):
                for b in range(self.dim):
                    for c in range(self.dim):
                        for d in range(self.dim):
                            for e in range(self.dim):
                                # ∇_e R_abcd
                                key = (e, a, b, c, d)
                                r_abcd = self._curvature_components[(a, b, c, d)]

                                # For simplicity, we assume the frame vectors are coordinate-aligned
                                # A full implementation would use the connection coefficients
                                derivatives[key] = simplify(diff(r_abcd, coords[e]))

            self._curvature_derivatives[order] = derivatives
            return derivatives
        else:
            # Higher derivatives (we just take additional coordinate derivatives for simplicity)
            derivatives = {}
            prev_derivatives = self._curvature_derivatives[order - 1]

            for indices, value in prev_derivatives.items():
                for e in range(self.dim):
                    new_indices = (e,) + indices
                    derivatives[new_indices] = simplify(diff(value, coords[e]))

            self._curvature_derivatives[order] = derivatives
            return derivatives


class NullFrame:
    """Class representing a null frame for a 4D Lorentzian metric."""

    def __init__(
        self,
        metric: Metric,
        l_vector: Matrix,
        n_vector: Matrix,
        m_vector: Matrix,
        mbar_vector: Matrix,
    ):
        """Initialize a null frame.

        Args:
            metric: The metric for which this is a frame
            l_vector: The first null vector field (ℓ)
            n_vector: The second null vector field (n)
            m_vector: The complex null vector field (m)
            mbar_vector: The complex conjugate of m (m̄)
        """
        if metric.dim != 4:
            raise ValueError("Null frame is currently only implemented for 4D metrics")

        self.metric = metric
        self.l = l_vector
        self.n = n_vector
        self.m = m_vector
        self.mbar = mbar_vector
        self.vectors = [l_vector, n_vector, m_vector, mbar_vector]

        # Validate dimensions
        for v in self.vectors:
            if v.shape != (4, 1):
                raise ValueError("Frame vectors must be 4x1 column matrices")

        # Verify that the frame satisfies the null frame conditions
        self._validate_frame()

        # Store frame components
        self.frame_components = Matrix.hstack(*self.vectors)

        # Computed properties
        self._curvature_components = None
        self._weyl_scalars = None
        self._curvature_derivatives = {}

    def _validate_frame(self):
        """Verify that the frame vectors satisfy the null frame conditions."""
        g = self.metric.components

        # l·l = n·n = m·m = m̄·m̄ = 0 (all vectors are null)
        # l·m = l·m̄ = n·m = n·m̄ = 0 (orthogonality)
        # l·n = -1, m·m̄ = 1 (normalization)

        def inner_product(v1, v2):
            result = 0
            for i in range(4):
                for j in range(4):
                    result += v1[i, 0] * v2[j, 0] * g[i, j]
            return simplify(result)

        conditions = [
            (inner_product(self.l, self.l), 0, "l·l = 0"),
            (inner_product(self.n, self.n), 0, "n·n = 0"),
            (inner_product(self.m, self.m), 0, "m·m = 0"),
            (inner_product(self.mbar, self.mbar), 0, "m̄·m̄ = 0"),
            (inner_product(self.l, self.m), 0, "l·m = 0"),
            (inner_product(self.l, self.mbar), 0, "l·m̄ = 0"),
            (inner_product(self.n, self.m), 0, "n·m = 0"),
            (inner_product(self.n, self.mbar), 0, "n·m̄ = 0"),
            (inner_product(self.l, self.n), -1, "l·n = -1"),
            (inner_product(self.m, self.mbar), 1, "m·m̄ = 1"),
        ]

        for actual, expected, message in conditions:
            if not simplify(actual - expected).is_zero:
                raise ValueError(
                    f"Null frame condition violated: {message}, got {actual}"
                )

    def weyl_scalars(self) -> Dict[str, sp.Expr]:
        """Compute the Weyl scalars for the null frame.

        Returns:
            Dictionary mapping scalar names to their values
        """
        if self._weyl_scalars is not None:
            return self._weyl_scalars

        # Get Weyl tensor
        weyl = self.metric.weyl_tensor()

        # Project onto null frame
        # Ψ₀ = C_μνρσ ℓ^μ m^ν ℓ^ρ m^σ
        # Ψ₁ = C_μνρσ ℓ^μ n^ν ℓ^ρ m^σ
        # Ψ₂ = C_μνρσ ℓ^μ m^ν m̄^ρ n^σ
        # Ψ₃ = C_μνρσ ℓ^μ n^ν m̄^ρ n^σ
        # Ψ₄ = C_μνρσ n^μ m̄^ν n^ρ m̄^σ

        l, n, m, mbar = self.l, self.n, self.m, self.mbar

        psi0 = 0
        psi1 = 0
        psi2 = 0
        psi3 = 0
        psi4 = 0

        for mu in range(4):
            for nu in range(4):
                for rho in range(4):
                    for sigma in range(4):
                        c = weyl.get((mu, nu, rho, sigma), 0)
                        psi0 += c * l[mu, 0] * m[nu, 0] * l[rho, 0] * m[sigma, 0]
                        psi1 += c * l[mu, 0] * n[nu, 0] * l[rho, 0] * m[sigma, 0]
                        psi2 += c * l[mu, 0] * m[nu, 0] * mbar[rho, 0] * n[sigma, 0]
                        psi3 += c * l[mu, 0] * n[nu, 0] * mbar[rho, 0] * n[sigma, 0]
                        psi4 += c * n[mu, 0] * mbar[nu, 0] * n[rho, 0] * mbar[sigma, 0]

        self._weyl_scalars = {
            "Psi0": simplify(psi0),
            "Psi1": simplify(psi1),
            "Psi2": simplify(psi2),
            "Psi3": simplify(psi3),
            "Psi4": simplify(psi4),
        }

        return self._weyl_scalars

    def petrov_type(self) -> str:
        """Determine the Petrov type of the spacetime.

        Returns:
            String indicating the Petrov type (I, II, D, III, N, or O)
        """
        scalars = self.weyl_scalars()

        # Check in order of increasing specialization
        if all(simplify(scalars[f"Psi{i}"]).is_zero for i in range(5)):
            return "O"  # Conformally flat

        if (
            simplify(scalars["Psi0"]).is_zero
            and simplify(scalars["Psi1"]).is_zero
            and simplify(scalars["Psi2"]).is_zero
            and simplify(scalars["Psi3"]).is_zero
        ):
            return "N"  # Psi4 is the only non-zero component

        if (
            simplify(scalars["Psi0"]).is_zero
            and simplify(scalars["Psi1"]).is_zero
            and simplify(scalars["Psi2"]).is_zero
        ):
            return "III"  # Psi3 and Psi4 are the only non-zero components

        if (
            simplify(scalars["Psi0"]).is_zero
            and simplify(scalars["Psi1"]).is_zero
            and simplify(scalars["Psi3"]).is_zero
            and simplify(scalars["Psi4"]).is_zero
        ):
            return "D"  # Psi2 is the only non-zero component

        if simplify(scalars["Psi0"]).is_zero and simplify(scalars["Psi1"]).is_zero:
            return "II"  # Psi2, Psi3, Psi4 are the only non-zero components

        return "I"  # Generic case

    def curvature_components(self) -> Dict[Tuple[int, int, int, int], sp.Expr]:
        """Compute the components of the Riemann tensor in the null frame.

        Returns:
            Dictionary mapping (a, b, c, d) -> R_abcd where a,b,c,d index the null frame
        """
        if self._curvature_components is not None:
            return self._curvature_components

        # Create an ordered list of frame vectors and compute components
        frame_vectors = [self.l, self.n, self.m, self.mbar]
        g = self.metric.components

        # Get Riemann tensor with lowered indices
        riemann = self.metric.riemann_tensor()
        riemann_low = {}

        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        riemann_low[(i, j, k, l)] = 0
                        for a in range(4):
                            riemann_low[(i, j, k, l)] += g[i, a] * riemann[(a, j, k, l)]
                        riemann_low[(i, j, k, l)] = simplify(riemann_low[(i, j, k, l)])

        # Project onto the null frame
        frame_components = {}

        for a in range(4):
            for b in range(4):
                for c in range(4):
                    for d in range(4):
                        component = 0
                        for mu in range(4):
                            for nu in range(4):
                                for rho in range(4):
                                    for sigma in range(4):
                                        component += (
                                            riemann_low[(mu, nu, rho, sigma)]
                                            * frame_vectors[a][mu, 0]
                                            * frame_vectors[b][nu, 0]
                                            * frame_vectors[c][rho, 0]
                                            * frame_vectors[d][sigma, 0]
                                        )

                        frame_components[(a, b, c, d)] = simplify(component)

        self._curvature_components = frame_components
        return frame_components

    def curvature_derivative(self, order: int = 1) -> Dict[Tuple, sp.Expr]:
        """Compute the components of the derivative of the Riemann tensor in the null frame.

        Args:
            order: The order of the derivative to compute

        Returns:
            Dictionary mapping multi-indices to derivative components
        """
        if order in self._curvature_derivatives:
            return self._curvature_derivatives[order]

        if order == 0:
            return self.curvature_components()

        if order > 1 and (order - 1) not in self._curvature_derivatives:
            # Recursively compute previous order
            self.curvature_derivative(order - 1)

        # For simplicity, we'll take a similar approach as in the LorentzFrame class
        # This is a simplified implementation that would need refinement in practice

        coords = self.metric.coordinates

        if order == 1:
            # First derivative
            derivatives = {}
            components = self.curvature_components()

            for a in range(4):
                for b in range(4):
                    for c in range(4):
                        for d in range(4):
                            for e in range(4):
                                # ∇_e R_abcd (simplified approach)
                                key = (e, a, b, c, d)
                                r_abcd = components[(a, b, c, d)]

                                # For simplicity, take coordinate derivative
                                # A complete implementation would use connection coefficients
                                derivatives[key] = simplify(diff(r_abcd, coords[e]))

            self._curvature_derivatives[order] = derivatives
            return derivatives
        else:
            # Higher derivatives
            derivatives = {}
            prev_derivatives = self._curvature_derivatives[order - 1]

            for indices, value in prev_derivatives.items():
                for e in range(4):
                    new_indices = (e,) + indices
                    derivatives[new_indices] = simplify(diff(value, coords[e]))

            self._curvature_derivatives[order] = derivatives
            return derivatives
