---
editor_options: 
  markdown: 
    wrap: 72
---

# Finite-Sized Particle Interaction Model

This repository is a work in progress and contains the MATLAB scripts of
a novel model of particle-particle interactions, demonstrating how
**electrostatic**, **magnetic**, and **gravitational forces** emerge
from a Galilean framework of finite-sized particles emitting and
absorbing **vector corpuscles**. The model challenges conventional
physics by proposing that these forces arise from corpuscle exchanges,
potentially eliminating the need for dark matter in explaining galactic
rotation curves.

## Model Overview

The model assumes a universe filled with **vector
corpuscles**—hypothetical entities carrying position, velocity, and
orientation—that mediate interactions between finite-sized particles.
Each corpuscle is characterized by:

-   **Position**: $\mathbf{x_c}$,

-   **Velocity**: $\mathbf{v_c}$,

-   **Orientation**: $\mathbf{o_c}$ a unit vector.

The corpuscle velocity follows **Galilean relativity**:

$$ \mathbf{v_c} = \mathbf{v_p} + \mathbf{c},$$

where $\mathbf{v_p}$ is the velocity of the emitting particle, and
$\mathbf{c}$ is a vector with magnitude equal to the speed of light $c$.
The orientation, $\mathbf{o_c}$, is defined at the moment the corpuscle
is emitted by the particle.

## Finite-Sized Particles

The model defines two types of particles, each emitting $q$ corpuscles
per second:

-   **Positive particles** $\mathbf{p^+}$ : Emit corpuscles with
    orientation **parallel** to their velocity:

$$ \mathbf{o_c} \cdot \mathbf{c} = c.$$

-   **Negative particles** $\mathbf{p^-}$ : Emit corpuscles with
    orientation **anti-parallel** to their velocity:

$$ \mathbf{o_c} \cdot \mathbf{c} = -c.$$

The emission of the corpuscles is isotropic, and at the speed of light
$c$.

Particles **absorb** the corpuscles from the surrounding medium. The
absorption process generates an **action** that alters the particle's
velocity based on the corpuscles' density, velocity, and orientation:

-   For $\mathbf{p^+}$ , the action aligns with the absorbed corpuscle’s
    orientation.

-   For $\mathbf{p^-}$ , the action opposes the absorbed corpuscle’s
    orientation.

## Corpuscle Absorption Model

Here I develop the action model of an receiver particle acting on
corpuscles emitted by an independent particle. The basic model of
particle-particle interaction is based on an universe filled with
corpuscles and any particle that emits and absorbs corpuscles at a $q$
rate of corpuscles per second. An emitting particle sends the corpuscles
at the speed of light and any receiver particle absorb corpuscles at a
$\mu$ rate of particles per second as long as the emitting corpuscle
flow is inside the particle.

The coordinate system is centered at the observer. This system will be
used to describe how the corpuscle-particle interaction changes the
receiver velocity from the observer point of view. The model development
starts by assuming that there is a single moving emitting particle
$\mathbf{p_e}$ and a single moving receiver particle $\mathbf{P_r}$. At
emission time $t_o$ the emitting particle is at $\mathbf{x_e}(t_o)$
moving with a relative velocity $\mathbf{v_e}(t_o)$ and emitting
isotropic $q_e$ light speed corpuscles. At absorption time,
$\mathbf{P_r}$ is at $\mathbf{x_r}(t)$ and moving at $\mathbf{v_r}(t)$.

The density of corpuscles at the receiver is:

$$
\rho( \mathbf{x_r,t,t_o})= \frac{q_e}{4\pi (c \Delta t)^2 \| \mathbf{c}(t_o)+\mathbf{v_e}(t_o) \|};
$$

where $\Delta t$ is the time interval it took the corpuscle to reach the
receiver. The corpuscle will travel the distance
$r(t,t_o)=\| \mathbf{X_r}(t)-\mathbf{X_e}(t_o)\|$ at the Galilean speed
$v_c(t_o)=\| \mathbf{c}(t_o)+\mathbf{v_e}(t_o)\|$. Hence:

$$
\Delta t=\frac{r(t,t_o)}{\| \mathbf{c}(t_o)+\mathbf{v_e}(t_o)\|}
$$

Therefore the corpuscle density at $\mathbf{x_r(t)}$ is:

$$
\rho( \mathbf{x_r},t,t_o)= \frac{q_e \| \mathbf{c}(t_o)+\mathbf{v_e}(t_o) \|}{4\pi c^2 r^2(t,t_o)} ;
$$

The corpuscle vectorial flow is $\Phi = \rho v$ .Hence:

$$
\Phi( \mathbf{x_r},t,t_o)= \frac{q_e}{4\pi (c \Delta t)^2 \| \mathbf{c}(t_o)+\mathbf{v_e}(t_o) \|}
( \mathbf{c}(t_o)+\mathbf{v_e}(t_o)-\mathbf{v_r}(t) )
$$ or $$
\Phi( \mathbf{x_r},t,t_o)= \frac{q_e \| \mathbf{c}(t_o)+\mathbf{v_e}(t_o) \|}{4\pi c^2 r^2(t,t_o)} 
( \mathbf{c}(t_o)+\mathbf{v_e}(t_o)-\mathbf{v_r}(t) )
$$

The next step is to compute the total number of corpuscles absorbed by
$\mathbf{P_r}$. I will assume that the number of absorbed corpuscles at
a given time is proportional to the flow and the particle cross area
$\mathbf{A}=A \hat{o_e}$; where $\hat{o_e}$ is the orientation of the
interacting corpuscles at emission from $\mathbf{P_e}$ . The total of
absorbed corpuscles is:

$$
Abs(\mathbf{x_r},t,t_o)=\Phi( \mathbf{x_r},t,t_o) \cdot \mathbf{A} (1-e^{-\mu \delta t});  
$$

where $\delta t$ is the time required by the flow to cross
$\mathbf{P_r}$.

I will assume that the crossing time is:

$$
\delta t=\frac{r_p}{\| \mathbf{c}(t_o)+\mathbf{v_e}(t_o)-\mathbf{v_r}(t) \|};
$$

where $r_p$ is the particle radius. For very small $r_p$ the total
absorbed is:

$$
Abs(\mathbf{x_r},t,t_o)=\frac{\Phi( \mathbf{x_r},t,t_o) \cdot \mathbf{A} \mu r_p}{\| \mathbf{c}(t_o)+\mathbf{v_e}(t_o)-\mathbf{v_r}(t) \|}.
$$

The absorbed particles causes a change in the particle momentum in
direction of the corpuscle orientation:

$$
\frac{d\mathbf{mv_r}(t)}{dt} \propto  Abs(\mathbf{x_r},t,t_o) \left( \frac{(\mathbf{c}(t_o)+\mathbf{v_e}(t_o)-\mathbf{v_r}(t))}{\| \mathbf{c}(t_o)+\mathbf{v_e}(t_o)-\mathbf{v_r}(t) \|} \cdot \hat{o_e}(t_o) \right) \hat{o_e}(t_o);
$$

From the above equations we arrive to the final action formula:

$$
 \frac{d\mathbf{mv_r}(t)}{dt} \propto A r_p \frac{q_e \mu_r}{4 \pi r^2} \frac{\|\mathbf{c} + \mathbf{v_e}\|}{c^2}  \left( \frac{\hat{o_e} \cdot ( \mathbf{c} + \mathbf{v_e} - \mathbf{v_r} )}{ \|\mathbf{c} + \mathbf{v_e} - \mathbf{v_r}\|} \right)^2 \hat{o_e},
$$

I'll assume that the area is proportional to the square of the particle
radius ($A \propto r_p^2$); therefore:

$$
 \frac{d\mathbf{mv_r}(t)}{dt} \propto {r_p}^3 \frac{q_e \mu_r}{4 \pi r^2} \frac{\|\mathbf{c} + \mathbf{v_e}\|}{c^2}  \left( \frac{\hat{o_e} \cdot ( \mathbf{c} + \mathbf{v_e} - \mathbf{v_r} )}{\|\mathbf{c} + \mathbf{v_e} - \mathbf{v_r}\|} \right)^2 \hat{o_e}
$$

Dimensional Analysis

$r^2: m^2$

${r_p}^3: m^3$

$q_e: p/s$

$\mu_r: 1/s$

$v:(m/s)$

$c: (m/s)$ $$
 (p/s) = {m}^3 \frac{(p/s)(1/s)}{m^2} \frac{m/s}{{(m/s)}^2}
$$

The action is proportional to the number of absorbed corpuscles per
second.

## Electrostatic Force

The model predicts that the net force $\mathbf{f_2}$ of the particle
absorbing corpuscles $q_2$ (Absorbing) of an emitting charge $q_1$
(Emitting) and (e.g., electron charges $e$ separated by distance $r$ (at
emitting-absorbing) is:

$$ \mathbf{f_2} = \frac{k q_1 q_2}{4 \pi r^2} \frac{\|\mathbf{c} + \mathbf{v_1}\|}{\|c\|}  \left( \frac{\mathbf{c} \cdot ( \mathbf{c} + \mathbf{v_1} - \mathbf{v_2} )}{\|\mathbf{c}\| \|\mathbf{c} + \mathbf{v_1} - \mathbf{v_2}\|} \right)^2 \hat{o_1},$$

where:

-   $k$: Coulomb constant,

-   $q$ is positive ($>0$) for $\mathbf{p^+}$ and negative ( $<0$) for
    $\mathbf{p^-}$

-   $\mathbf{v_1}, \mathbf{v_2}$: Velocities of the two particles:

    -   $\mathbf{v_1}$ is an emitting particle. $\mathbf{v_2}$ the
        absorbing particle

-   $\hat{o_1}$: Unit vector of the corpuscle orientation at the time of
    its origin by the emitting particle.

-   $r$ is the distance traveled by the corpuscle form the emission
    origin to it absorption. i.e.,
    $$ r=\| \mathbf{x_1}(t_o)-\mathbf{x_2}(t) \| $$

This force generalizes Coulomb’s law by incorporating velocity-dependent
effects.

## Mass and Inertia

Inertia is implicit in the model, and the mass of particles in a direct
consequence of the finite-sized particle model. The inertial mass of a
finite-sized particle is proportional to its radius and charge
distribution. For charged particles, the mass is related to the
**electrostatic stored energy** divided by ( c\^2 ):

$$ m \propto \frac{\text{Electrostatic energy}}{c^2}.$$

The Mass folder derives the E=mc\^2 formula.

## Magnetic Force

The **magnetic force** and the **vacuum magnetic permeability** emerge
from interactions between **neutral currents** (e.g., moving neutral
particles) and moving charged particles , driven by the exchange of
vector corpuscles. Evidence is presented in the MagneticForce folder for
a current in a loop.

## Gravitational Force

**Gravity** arises as a net attractive force between **neutral composite
particles** (e.g., atoms with a positively charged nucleus and a
negatively charged shell). The slight difference in stocastic velocity
distributions between positive and negative charges results in an
Expected small attraction. The equations of interaction are derived in
the GravityForce/SteadySate folder

## Galactic Rotation Curves

The velocity-dependent nature of the gravitational force implies that
**"hot" ionized particles** in neutral plasma (with high velocities)
behave differently from **"cold" neutral particles**. This difference
accounts for the observed **flattening of galactic rotation curves**,
potentially eliminating the need for dark matter in this model.
Simulations in this repository demonstrate how these effects align with
observed galactic dynamics. The results of fitting observed rotation
curves of 22 SPARC galaxies and the local group are showcased in the
Galaxies folder.

![M33 Rotation Curve](images/clipboard-1366545713.png){width="392"}

![M33 Surface mass
densities](images/clipboard-149743154.png){width="383"}

![Rotation Curve](Galaxies/SPARC/DensityProfiles/NGC7331.jpg)

The evaluation of the 22 galaxies implies that there is no need for dark
matter for the explanation of the rotation curves:

![](images/clipboard-4116685001.png)

Collaboration is welcome. Please send me comments, corrections and
suggestions.
