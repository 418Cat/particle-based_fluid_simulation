# Particle-based fluid simulation


This project is based on the paper [<ins>*SPH Fluids in Computer Graphics*</ins>](https://cgl.ethz.ch/Downloads/Publications/Papers/2014/Sol14a/Sol14a.pdf) available on ETH Zurich's website.

In the current state, the fluid simulation part uses mainly the following equations (Pages 2 & 3):

___
### The Kernel function
$$W_{ij} = \frac{1}{h^d}f(q)$$
<br>
$f(q)$ being close to a gaussian with compact support: $f(q) \gt 0$ for $0 \le q \lt 1$ or $0 \le q \lt 2$ (__Eq 5__)
<br>
with $q = \frac{||\mathbf{x}_i - \mathbf{x}_j||}{h}$
<br>
<br>
This equation describes how much influence a neighbouring particle $j$ has on the current particle $i$'s properties, based on their distance


___
### Density $\rho_i$ computation
$$\rho_i = \sum_j m_j W_{ij}$$
with $j$ the neighbouring particles and $m_j$ their respective masses
<br>
<br>
Using the kernel function defined above, this equation computes the density $\rho$ at a particle $i$'s position

___
### Pressure $p_i$ computation
There are multiple ways to get the pressure $p_i$. Section 1.3 (p.2) shows __Eq 9__: $$p_i = k\left(\left(\frac{\rho_i}{\rho_0}\right)^7-1\right)$$
Section 3.1 (p.5) proposes a few alternatives:
$$p_i = k\left(\rho_i - \rho_0\right)$$
or
$$p_i = k\left(\frac{\rho_i}{\rho_0} - 1 \right)$$
with $\rho_0$ the rest density (or target density) and $k$ a stiffness constant, scaling the pressure

___
### Pressure force computation $F_i^{pressure}$
The pressure force acting on a particle $i$ can be obtained using the pressure gradient $\nabla p_i$ from __Eq 6__ (p.2)
$$F_i^{pressure} = -\frac{m_i}{\rho_i} \nabla p_i$$

<br>
<br>
The pressure gradient allows particles to be pushed away from high pressure zones towards low pressure ones

___
### Viscosity force $F_i^{viscosity}$
$$F_i^{viscosity} = m_i v \nabla^2 \mathbf{v}_i$$
$\nabla^2\mathbf{v}_i$ can be obtained using __Eq 8__ (p.2)

___
With $\mathbf{a}_i = \frac{\mathbf{F}_i}{m_i}$, time integration can be used to deduce velocity $\mathbf{v}_i$ and position $\mathbf{x}_i$


<br>
<br>
# Current state
The simulation implements:

- Collision physics
- Gravity
- Fluid physics

## Collision physics
Collision detection uses bounding volumes to optimize spatial lookup time
<br>
Collision physics are handled either using velocity (More stable, less accurate) or acceleration (More accurate, tends to be less stable)


## Gravity
Domain gravity can be axis based (individual $xyz$ axis have separate gravity) or radial (based on the norm of the gravity vector)
<br><br>
Particles gravity is still work in progress. Barnes-hut algorithm was used to speedup gravity computation, but the implementation isn't tested yet.

## Liquid sim
The liquid simulation part has been implemented, individual constants and settings need to be tweaked and tested in order to obtain a correct and visually pleasing result.


___
This project is under the MIT license
