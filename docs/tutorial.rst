========
Tutorial
========

This tutorial provides a step-by-step guide to using the Cartan-Karlhede algorithm implementation.

Getting Started
--------------

First, import the necessary modules:

.. code-block:: python

    from cartan_karlhede.metric import Metric, LorentzFrame
    from cartan_karlhede.algorithm import compare_metrics
    from sympy import Matrix, Symbol, diag

Creating a Metric
----------------

To create a metric, you need to define:

1. The coordinate symbols
2. The metric components as a Matrix
3. An optional name for the metric

Here's an example of creating a 4D Minkowski metric:

.. code-block:: python

    # Define coordinate symbols
    t = Symbol('t')
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    
    # Define metric components (using diag for a diagonal matrix)
    g = diag(-1, 1, 1, 1)
    
    # Create the metric object
    minkowski_metric = Metric(g, [t, x, y, z], "Minkowski 4D")

Computing Curvature Tensors
--------------------------

Once you have a metric, you can compute various curvature tensors:

.. code-block:: python

    # Compute the Christoffel symbols
    christoffel = minkowski_metric.christoffel_symbols()
    
    # Compute the Riemann tensor
    riemann = minkowski_metric.riemann_tensor()
    
    # Compute the Ricci tensor
    ricci = minkowski_metric.ricci_tensor()
    
    # Compute the Ricci scalar
    ricci_scalar = minkowski_metric.ricci_scalar()
    
    # Compute the Weyl tensor (in dimension ≥ 3)
    weyl = minkowski_metric.weyl_tensor()

Creating Lorentz Frames
----------------------

To apply the Cartan-Karlhede algorithm, you need to define frames for your metrics:

.. code-block:: python

    # Create a 3D Minkowski metric
    t = Symbol('t')
    x = Symbol('x')
    y = Symbol('y')
    g = diag(-1, 1, 1)
    metric = Metric(g, [t, x, y], "Minkowski 3D")
    
    # Define a frame (column vectors)
    e0 = Matrix([[1], [0], [0]])  # Timelike vector
    e1 = Matrix([[0], [1], [0]])  # Spacelike vector
    e2 = Matrix([[0], [0], [1]])  # Spacelike vector
    
    # Create a Lorentz frame
    frame = LorentzFrame(metric, [e0, e1, e2])
    
    # Compute curvature in this frame
    frame_curvature = frame.curvature_components()
    
    # Compute derivatives of the curvature
    first_derivatives = frame.curvature_derivative(order=1)
    second_derivatives = frame.curvature_derivative(order=2)

Comparing Metrics
----------------

To compare two metrics using the Cartan-Karlhede algorithm:

.. code-block:: python

    # Create two metrics and their frames
    # (assuming you have metric1, frame1, metric2, frame2 defined)
    
    # Compare the metrics
    are_equivalent, reason = compare_metrics(metric1, frame1, metric2, frame2)
    
    print(f"Are equivalent: {are_equivalent}")
    print(f"Reason: {reason}")

Complete Example
---------------

Here's a complete example that creates and analyzes a flat spacetime:

.. code-block:: python

    from cartan_karlhede.metric import Metric
    from sympy import Matrix, Symbol, simplify, diag
    
    # Create Minkowski metric in 4D
    t = Symbol('t')
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    
    g = diag(-1, 1, 1, 1)
    
    metric = Metric(g, [t, x, y, z], "Minkowski")
    
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

This tutorial covered the basics of using the Cartan-Karlhede algorithm implementation. For more advanced examples, including null frames and comparing metrics in different coordinate systems, see the :doc:`examples` section. 