========
Examples
========

This section provides examples of using the Cartan-Karlhede algorithm implementation.

Basic Examples
-------------

Flat Spacetime
^^^^^^^^^^^^^

The simplest example is analyzing flat (Minkowski) spacetime:

.. code-block:: python

    from cartan_karlhede.simple_examples import analyze_flat_spacetime
    
    # Analyze Minkowski metric
    metric = analyze_flat_spacetime()

This will compute the Riemann tensor for the Minkowski metric and verify that all components are zero, confirming that the spacetime is flat.

3D Examples
----------

The paper "On the Equivalence of Spacetimes, the Cartan-Karlhede Algorithm" by Mergulhão and Batista discusses several 3D metrics. This implementation includes these examples:

.. code-block:: python

    from cartan_karlhede.examples import compare_3d_examples
    
    # Compare the 3D example metrics
    compare_3d_examples()

This will:

1. Create five different 3D metrics from the paper
2. Construct appropriate frames for each metric
3. Compare the first metric with each of the others using the Cartan-Karlhede algorithm

The implementation shows that some of these seemingly different metrics actually represent the same geometric spacetime.

Schwarzschild Examples
--------------------

The implementation also includes the Schwarzschild metric in different coordinate systems:

.. code-block:: python

    from cartan_karlhede.examples import compare_schwarzschild_examples
    
    # Compare Schwarzschild metrics in different coordinate systems
    compare_schwarzschild_examples()

This example demonstrates that the standard Schwarzschild metric and an alternative representation are equivalent when :math:`c = 2M`.

Custom Examples
-------------

You can create your own examples by defining metrics and comparing them:

.. code-block:: python

    from cartan_karlhede.metric import Metric, LorentzFrame
    from cartan_karlhede.algorithm import compare_metrics
    from sympy import Matrix, Symbol, diag, exp
    
    # Create a simple 2D curved metric
    t = Symbol('t')
    x = Symbol('x')
    
    # First metric: ds^2 = -dt^2 + exp(2t) dx^2
    g1 = Matrix([[-1, 0], [0, exp(2*t)]])
    metric1 = Metric(g1, [t, x], "2D curved spacetime")
    
    # Create a frame for this metric
    e0 = Matrix([[1], [0]])
    e1 = Matrix([[0], [exp(-t)]])
    frame1 = LorentzFrame(metric1, [e0, e1])
    
    # Another representation of the same metric with a coordinate transformation
    T = Symbol('T')
    X = Symbol('X')
    
    # Second metric: ds^2 = -dT^2 + T^2 dX^2
    g2 = Matrix([[-1, 0], [0, T**2]])
    metric2 = Metric(g2, [T, X], "2D curved spacetime (alternative)")
    
    # Create a frame for this metric
    E0 = Matrix([[1], [0]])
    E1 = Matrix([[0], [1/T]])
    frame2 = LorentzFrame(metric2, [E0, E1])
    
    # Compare the metrics
    are_equivalent, reason = compare_metrics(metric1, frame1, metric2, frame2)
    
    print(f"Are the metrics equivalent? {are_equivalent}")
    print(f"Reason: {reason}")

Command Line Examples
-------------------

The package can also be run from the command line:

.. code-block:: bash

    # Run all examples
    python -m cartan_karlhede --verbose
    
    # Run only the 3D examples
    python -m cartan_karlhede --example 3d --verbose
    
    # Run only the Schwarzschild examples
    python -m cartan_karlhede --example schwarzschild --verbose

With the Justfile
^^^^^^^^^^^^^^^

If using the Justfile, you can run examples with:

.. code-block:: bash

    # Run all examples
    just run-examples
    
    # Run only the 3D examples
    just run-3d
    
    # Run only the Schwarzschild examples
    just run-schwarzschild 