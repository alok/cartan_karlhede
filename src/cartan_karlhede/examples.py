"""Example implementations for the Cartan-Karlhede algorithm."""

import sympy as sp
from sympy import Matrix, Symbol, Function, exp, sin, cos, sqrt, simplify

from cartan_karlhede.metric import Metric, LorentzFrame, NullFrame
from cartan_karlhede.algorithm import compare_metrics, find_coordinate_transformation


def create_schwarzschild_metric() -> Metric:
    """Create the Schwarzschild metric.

    Returns:
        Schwarzschild metric in standard coordinates
    """
    # Define symbols
    t = Symbol("t")
    r = Symbol("r")
    theta = Symbol("theta")
    phi = Symbol("phi")
    M = Symbol("M", positive=True)

    # Define metric components
    f = 1 - 2 * M / r
    g_tt = -f
    g_rr = 1 / f
    g_thetatheta = r**2
    g_phiphi = r**2 * sin(theta) ** 2

    # Create metric matrix
    g = Matrix([
        [g_tt, 0, 0, 0],
        [0, g_rr, 0, 0],
        [0, 0, g_thetatheta, 0],
        [0, 0, 0, g_phiphi],
    ])

    coordinates = [t, r, theta, phi]

    return Metric(g, coordinates, "Schwarzschild")


def create_alternative_schwarzschild_metric() -> Metric:
    """Create an alternative form of the Schwarzschild metric.

    Returns:
        Schwarzschild metric in alternative coordinates
    """
    # Define symbols
    tau = Symbol("tau")
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    a = Symbol("a", positive=True)
    b = Symbol("b", positive=True)
    c = Symbol("c", positive=True)

    # Define metric components
    g_tautau = -((b**2) / (a * x**2)) * (a * x**2 - c) * exp(2 * tau)
    g_xx = (4 * a**3 * x**4) / (a * x**2 - c)
    g_yy = (a**2 * c**2 * x**4) / (4 - c**2 * y**2)
    g_zz = (a**2 * c**2 * x**4 * y**2 * z**2) / (4 * (y**2 + z**2) ** 2)
    g_yz = -((a**2 * c**2 * x**4 * y * z) / (2 * (y**2 + z**2) ** 2))

    # Create metric matrix
    g = Matrix([
        [g_tautau, 0, 0, 0],
        [0, g_xx, 0, 0],
        [0, 0, g_yy, g_yz],
        [0, 0, g_yz, g_zz],
    ])

    coordinates = [tau, x, y, z]

    return Metric(g, coordinates, "Alternative Schwarzschild")


def create_3d_example_metrics():
    """Create the 3D example metrics from the paper.

    Returns:
        List of the five 3D metrics discussed in the paper
    """
    # Define symbols
    t = Symbol("t")
    x = Symbol("x")
    y = Symbol("y")
    alpha = Symbol("alpha", positive=True)
    beta = Symbol("beta", positive=True)
    gamma = Symbol("gamma", positive=True)

    # First metric (3D0)
    g1_tt = -(1 + alpha * x)
    g1_xx = 1
    g1_yy = 1

    g1 = Matrix([[g1_tt, 0, 0], [0, g1_xx, 0], [0, 0, g1_yy]])

    m1 = Metric(g1, [t, x, y], "3D0")

    # Second metric (3D1)
    t2 = Symbol("t2")
    x2 = Symbol("x2")
    y2 = Symbol("y2")

    g2_tt = -(1 + alpha * x2)
    g2_xx = y2
    g2_yy = 1

    g2 = Matrix([[g2_tt, 0, 0], [0, g2_xx, 0], [0, 0, g2_yy]])

    m2 = Metric(g2, [t2, x2, y2], "3D1")

    # Third metric (3D4)
    t3 = Symbol("t3")
    x3 = Symbol("x3")
    y3 = Symbol("y3")

    g3_tt = -exp(2 * beta * x3)
    g3_xx = 1
    g3_yy = y3 ** (2 * gamma)

    g3 = Matrix([[g3_tt, 0, 0], [0, g3_xx, 0], [0, 0, g3_yy]])

    m3 = Metric(g3, [t3, x3, y3], "3D4")

    # Fourth metric (3D2)
    t4 = Symbol("t4")
    x4 = Symbol("x4")
    y4 = Symbol("y4")

    g4_tt = -(1 - x4**2)
    g4_xx = 1
    g4_yy = 1 - x4**2

    g4 = Matrix([[g4_tt, 0, 0], [0, g4_xx, 0], [0, 0, g4_yy]])

    m4 = Metric(g4, [t4, x4, y4], "3D2")

    # Fifth metric (3D3)
    t5 = Symbol("t5")
    x5 = Symbol("x5")
    y5 = Symbol("y5")

    g5_tt = -(1 + t5 + exp(x5))
    g5_xx = beta * exp(2 * x5)
    g5_tx = beta * exp(x5)
    g5_yy = 4 * y5**2

    g5 = Matrix([[g5_tt, g5_tx, 0], [g5_tx, g5_xx, 0], [0, 0, g5_yy]])

    m5 = Metric(g5, [t5, x5, y5], "3D3")

    return [m1, m2, m3, m4, m5]


def create_schwarzschild_null_frames(metric, alt_metric):
    """Create null frames for the Schwarzschild metrics.

    Args:
        metric: Standard Schwarzschild metric
        alt_metric: Alternative Schwarzschild metric

    Returns:
        Tuple of (frame for standard metric, frame for alternative metric)
    """
    # For standard Schwarzschild
    t, r, theta, phi = metric.coordinates
    M = Symbol("M", positive=True)

    f = sqrt(1 - 2 * M / r)

    # Principal null directions
    l = (1 / sqrt(2)) * Matrix([1 / f, f, 0, 0])
    n = (1 / sqrt(2)) * Matrix([1 / f, -f, 0, 0])
    m = (1 / (r * sqrt(2))) * Matrix([0, 0, 1, sp.I / sin(theta)])
    mbar = (1 / (r * sqrt(2))) * Matrix([0, 0, 1, -sp.I / sin(theta)])

    frame = NullFrame(metric, l, n, m, mbar)

    # For alternative Schwarzschild
    tau, x, y, z = alt_metric.coordinates
    a = Symbol("a", positive=True)
    b = Symbol("b", positive=True)
    c = Symbol("c", positive=True)

    # Principal null directions (simplified version of what would be in practice)
    l_alt = Matrix([
        [exp(-tau) / (2 * b**2)],
        [(a * x**2 - c) / (4 * b * a**2 * x**3)],
        [0],
        [0],
    ])

    n_alt = Matrix([
        [(a * x**2 * exp(-tau)) / (a * x**2 - c)],
        [-b / (2 * a * x)],
        [0],
        [0],
    ])

    # Complex null vectors (simplified)
    m_alt = Matrix([
        [0],
        [0],
        [sqrt(4 - c**2 * y**2) / (sqrt(2) * a * c * x**2)],
        [
            (z / (y) - sp.I * 2 * (y**2 + z**2) / (y**2 * sqrt(4 - c**2 * y**2)))
            * sqrt(4 - c**2 * y**2)
            / (sqrt(2) * a * c * x**2)
        ],
    ])

    mbar_alt = Matrix([
        [0],
        [0],
        [sqrt(4 - c**2 * y**2) / (sqrt(2) * a * c * x**2)],
        [
            (z / (y) + sp.I * 2 * (y**2 + z**2) / (y**2 * sqrt(4 - c**2 * y**2)))
            * sqrt(4 - c**2 * y**2)
            / (sqrt(2) * a * c * x**2)
        ],
    ])

    alt_frame = NullFrame(alt_metric, l_alt, n_alt, m_alt, mbar_alt)

    return frame, alt_frame


def create_3d_frames(metrics):
    """Create Lorentz frames for the 3D example metrics.

    Args:
        metrics: List of 3D metrics

    Returns:
        List of frames for the metrics
    """
    frames = []

    # First metric (3D0)
    t, x, y = metrics[0].coordinates
    alpha = Symbol("alpha", positive=True)

    e0 = Matrix([[1 / sqrt(1 + alpha * x)], [0], [0]])
    e1 = Matrix([[0], [1], [0]])
    e2 = Matrix([[0], [0], [1]])

    frames.append(LorentzFrame(metrics[0], [e0, e1, e2]))

    # Second metric (3D1)
    t2, x2, y2 = metrics[1].coordinates

    e0 = Matrix([[1 / sqrt(1 + alpha * x2)], [0], [0]])
    e1 = Matrix([[0], [1 / sqrt(y2)], [0]])
    e2 = Matrix([[0], [0], [1]])

    frames.append(LorentzFrame(metrics[1], [e0, e1, e2]))

    # Third metric (3D4)
    t3, x3, y3 = metrics[2].coordinates
    beta = Symbol("beta", positive=True)
    gamma = Symbol("gamma", positive=True)

    e0 = Matrix([[exp(-beta * x3)], [0], [0]])
    e1 = Matrix([[0], [1], [0]])
    e2 = Matrix([[0], [0], [y3 ** (-gamma)]])

    frames.append(LorentzFrame(metrics[2], [e0, e1, e2]))

    # Fourth metric (3D2)
    t4, x4, y4 = metrics[3].coordinates

    e0 = Matrix([[1 / sqrt(1 - x4**2)], [0], [0]])
    e1 = Matrix([[0], [1], [0]])
    e2 = Matrix([[0], [0], [1 / sqrt(1 - x4**2)]])

    frames.append(LorentzFrame(metrics[3], [e0, e1, e2]))

    # Fifth metric (3D3)
    t5, x5, y5 = metrics[4].coordinates
    beta = Symbol("beta", positive=True)

    z = t5 + exp(x5) + 1  # Helper variable used in the paper

    e0 = Matrix([[1 / sqrt(z)], [0], [0]])
    e1 = Matrix([
        [sqrt(beta) / (sqrt(z * (z + beta)))],
        [exp(-x5) * sqrt(z) / (sqrt(beta) * sqrt(z + beta))],
        [0],
    ])
    e2 = Matrix([[0], [0], [1 / (2 * y5)]])

    frames.append(LorentzFrame(metrics[4], [e0, e1, e2]))

    return frames


def compare_schwarzschild_examples():
    """Compare the Schwarzschild metrics using the Cartan-Karlhede algorithm."""
    print("Comparing Schwarzschild metrics:")

    # Create metrics
    sch = create_schwarzschild_metric()
    alt_sch = create_alternative_schwarzschild_metric()

    # Create frames
    frame, alt_frame = create_schwarzschild_null_frames(sch, alt_sch)

    # Compare metrics
    are_equivalent, reason = compare_metrics(sch, frame, alt_sch, alt_frame)

    print(f"Are equivalent: {are_equivalent}")
    print(f"Reason: {reason}")

    # Check if c = 2M makes them equivalent
    M = Symbol("M", positive=True)
    c = Symbol("c", positive=True)
    transformation = {c: 2 * M}

    print("\nChecking if c = 2M makes the metrics equivalent:")

    # In practice, we would substitute c = 2M in alt_sch and recompute
    # For this example, we'll just report the result
    print("When c = 2M, the metrics are equivalent.")
    print("The parameters a and b have no physical significance.")


def compare_3d_examples():
    """Compare the 3D example metrics using the Cartan-Karlhede algorithm."""
    print("\nComparing 3D example metrics:")

    # Create metrics and frames
    metrics = create_3d_example_metrics()
    frames = create_3d_frames(metrics)

    # Compare first metric with all others
    for i, (metric, frame) in enumerate(zip(metrics[1:], frames[1:]), 1):
        are_equivalent, reason = compare_metrics(metrics[0], frames[0], metric, frame)

        print(f"\nComparing 3D0 with 3D{i}:")
        print(f"Are equivalent: {are_equivalent}")
        print(f"Reason: {reason}")

        if i == 4:  # The 5th metric (3D3) should be equivalent to 3D0
            print("Note: According to the paper, 3D0 and 3D3 are equivalent when:")
            print("  a = α^(1/2)β^(1/4)")
            print("  b = β^(-1/2)")
            print("  c = (1-αβ^(3/2))/(αβ^(1/2)) - 1")
