Advection of a square-shaped scalar filed
============================================


This example is designed to evaluate the performance of different numerical methods in handling discontinuous problems. A square scalar field is advected by a uniform velocity field, and the initial conditions are defined as follows:

Initial Conditions
-------------------

The scalar field :math:`C(x, y)` is given by:

.. math::

    C(x, y) =
    \begin{cases} 
    10, & 0 \leq |x - x_0| \leq \frac{a}{2}, \, 0 \leq |y - y_0| \leq \frac{a}{2}, \\
    0, & \text{otherwise}.
    \end{cases}

- **Center point of the square**: :math:`(x_0, y_0) = (-1.5, -1.5)`
- **Side length of the square**: :math:`a = 1.5`
- **Domain**: :math:`-3 \leq x, y \leq 3`
- **Velocity field**: Constant at :math:`1` in both :math:`x`- and :math:`y`-directions.

Simulation Details
-------------------

- **Time duration**: :math:`2.8 \, \text{s}`
- **Time step**: :math:`\Delta t = 0.005 \, \text{s}`
- **Final position of the square's center**: :math:`(1.3, 1.3)`

Results
-------

The simulation results include:

- **Maximum and minimum values** for a grid of :math:`64 \times 64` nodes (see Table 1).
- **Final scalar field distributions**, provided in Appendix 1.

This setup offers a benchmark for analyzing the accuracy and stability of numerical methods in capturing discontinuities within scalar fields.
