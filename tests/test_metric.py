"""Tests for the Metric class."""

import pytest
import sympy as sp
from sympy import Matrix, Symbol, simplify, diag

from cartan_karlhede.metric import Metric, LorentzFrame


def test_flat_spacetime_riemann():
    """Test that flat spacetime has zero Riemann tensor."""
    # Create Minkowski metric in 4D
    t = Symbol("t")
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")

    g = diag(-1, 1, 1, 1)

    metric = Metric(g, [t, x, y, z], "Minkowski")

    # Compute Riemann tensor
    riemann = metric.riemann_tensor()

    # Check that all components are zero
    for indices, value in riemann.items():
        assert simplify(value) == 0


def test_frame_validation():
    """Test that frames are properly validated."""
    # Create 3D Minkowski metric
    t = Symbol("t")
    x = Symbol("x")
    y = Symbol("y")

    g = diag(-1, 1, 1)

    metric = Metric(g, [t, x, y], "Minkowski 3D")

    # Valid frame
    e0 = Matrix([[1], [0], [0]])
    e1 = Matrix([[0], [1], [0]])
    e2 = Matrix([[0], [0], [1]])

    # This should not raise an error
    frame = LorentzFrame(metric, [e0, e1, e2])

    # Invalid frame - wrong number of vectors
    with pytest.raises(ValueError):
        LorentzFrame(metric, [e0, e1])

    # Invalid frame - wrong dimensions
    e0_invalid = Matrix([[1], [0]])
    with pytest.raises(ValueError):
        LorentzFrame(metric, [e0_invalid, e1, e2])

    # Invalid frame - not orthonormal
    e0_non_ortho = Matrix([[1], [1], [0]])
    with pytest.raises(ValueError):
        LorentzFrame(metric, [e0_non_ortho, e1, e2])
