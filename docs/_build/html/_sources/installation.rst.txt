============
Installation
============

Requirements
-----------

The Cartan-Karlhede algorithm implementation requires:

* Python 3.12 or newer
* sympy: For symbolic mathematics
* numpy: For numerical operations
* jax (optional): For acceleration of numerical computations
* matplotlib/plotly (optional): For visualization

This package uses ``uv`` for package management, which is a fast, reliable Python package installer and resolver.

Basic Installation
----------------

To install the latest release of the package:

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/yourusername/cartan-karlhede.git
    cd cartan-karlhede

    # Create a virtual environment and install dependencies
    uv venv .venv
    source .venv/bin/activate
    uv pip install -e .

Development Installation
----------------------

If you want to contribute to the project or run the tests, install with the development dependencies:

.. code-block:: bash

    # Install with development dependencies
    uv pip install -e ".[dev]"

Using the Justfile
-----------------

The project includes a ``Justfile`` with common commands. If you have ``just`` installed, you can use:

.. code-block:: bash

    # List available commands
    just

    # Install the package with development dependencies
    just install

    # Run the examples
    just run-examples

    # Run tests
    just test

    # Run type checking
    just typecheck

    # Format the code
    just format

    # Run linting checks
    just lint 