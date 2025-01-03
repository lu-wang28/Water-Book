Convection-Diffusion of a Gaussian Distribution
===============================================

This example is designed to evaluate the performance of different numerical methods in convection-diffusion problems. A two-dimensional Gaussian scalar field is diffused and advected by a uniform velocity field. The initial conditions and parameters are defined as follows:

Initial Conditions
-------------------

The scalar field :math:`C(x, y, 0)` is given by:

.. math::

    C(x, y, 0) = \exp\left(-\frac{(x - x_0)^2}{2\delta_x^2} - \frac{(y - y_0)^2}{2\delta_y^2}\right)

- **Center point of the Gaussian distribution**: :math:`(x_0, y_0) = (1300, 1300)`
- **Standard deviations**: :math:`\delta_x = \delta_y = 200`
- **Domain**: :math:`0 \leq x, y \leq 6400`
- **Diffusion coefficients**: :math:`\alpha = \beta = 0.1`
- **Velocity field**: Constant at :math:`1` in both :math:`x`- and :math:`y`-directions.

Simulation Details
-------------------

- **Time duration**: :math:`3500 \, \text{s}`
- **Time step**: :math:`\Delta t = 1 \, \text{s}`
- **Final position of the Gaussian's center**: :math:`(4800, 4800)`

Results
-------

The simulation results include:

- **Maximum and minimum values** for a grid of :math:`128 \times 128` nodes (see Table 2).
- **Final scalar field distributions** and **concentration profiles at :math:`y = 3200`**, provided in Appendix 2.

This setup provides a benchmark for analyzing the accuracy and stability of numerical methods in solving convection-diffusion problems.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   ../CSSL/CSSL
