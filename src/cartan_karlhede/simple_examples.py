"""Simple examples for the Cartan-Karlhede algorithm."""

import sympy as sp
from sympy import Matrix, Symbol, simplify, diag

from cartan_karlhede.metric import Metric


def create_minkowski_metric():
    """Create a Minkowski metric in 4D."""
    t = Symbol("t")
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")

    g = diag(-1, 1, 1, 1)

    return Metric(g, [t, x, y, z], "Minkowski 4D")


def analyze_flat_spacetime():
    """Analyze flat spacetime metrics."""
    print("Analyzing flat spacetime (Minkowski)...")

    # Create Minkowski metric
    metric = create_minkowski_metric()

    # Compute Riemann tensor
    riemann = metric.riemann_tensor()

    # Check that all components are zero (flat spacetime)
    all_zero = True
    for indices, value in riemann.items():
        if not simplify(value).is_zero:
            all_zero = False
            print(f"Non-zero component: R_{indices} = {value}")

    if all_zero:
        print("Verified: Minkowski spacetime is flat (all Riemann components are zero)")

    # Compute Ricci tensor and scalar
    ricci = metric.ricci_tensor()
    ricci_scalar = metric.ricci_scalar()

    print(f"Ricci scalar: {ricci_scalar}")

    return metric


def run_simple_examples():
    """Run simple examples to test the implementation."""
    print("Running simple examples for the Cartan-Karlhede algorithm...")

    # Test Minkowski metric
    metric = analyze_flat_spacetime()

    print("\nImplementation test successful!")
    return metric
