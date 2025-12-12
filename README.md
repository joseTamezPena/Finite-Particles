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

# Corpuscular Action Model for Particle–Particle Interaction

I develop a mechanistic action model describing the interaction between
a receiver particle and corpuscles emitted by an independent emitting
particle. The model is grounded in a corpuscular paradigm, wherein the
universe is permeated by discrete corpuscles propagated at the speed of
light. Particles emit and absorb these corpuscles at characteristic
rates, mediating momentum transfer analogous to a radiative exchange
process.

## Model Foundations

The interaction framework posits that emitting particles radiate
corpuscles isotropically at a rate of $q$ corpuscles per second, with
propagation occurring at the vacuum speed of light $c$. Receiver
particles absorb incident corpuscles at a rate parameterised by $\mu$
(dimensions of inverse time), contingent upon the corpuscular flux
intersecting the receiver.

Observations are conducted in a coordinate system centred on a
stationary observer, from whose frame the velocity evolution of the
receiver particle is described.

The system under consideration comprises a single emitting particle
$\mathbf{p_e}$ and a single receiver particle $\mathbf{P_r}$, both
potentially in motion. At emission time $t_o$, the emitter occupies
position $\mathbf{x_e}(t_o)$ with velocity $\mathbf{v_e}(t_o)$,
radiating at rate $q_e$. At absorption time $t$, the receiver is at
$\mathbf{x_r}(t)$ with velocity $\mathbf{v_r}(t)$.

## Corpuscular Density and Flux

The number density of corpuscles at the receiver position is

$$
\rho(\mathbf{x_r}, t, t_o) = \frac{q_e}{4\pi (c \Delta t)^2 \|\mathbf{c}(t_o) + \mathbf{v_e}(t_o)\|},
$$

where $\Delta t$ is the propagation duration, and $\mathbf{c}(t_o)$ is
the light-speed directional vector from emitter to receiver in the
observer frame. The separation distance is
$r(t, t_o) = \|\mathbf{x_r}(t) - \mathbf{x_e}(t_o)\|$, yielding an
effective propagation speed $\|\mathbf{c}(t_o) + \mathbf{v_e}(t_o)\|$.
Thus,

$$
\Delta t = \frac{r(t, t_o)}{\|\mathbf{c}(t_o) + \mathbf{v_e}(t_o)\|}.
$$

Substitution provides the density

$$
\rho(\mathbf{x_r}, t, t_o) = \frac{q_e \|\mathbf{c}(t_o) + \mathbf{v_e}(t_o)\|}{4\pi c^2 r^2(t, t_o)}.
$$

The vectorial flux relative to the moving receiver is

$$
\boldsymbol{\Phi}(\mathbf{x_r}, t, t_o) = \rho \, [\mathbf{c}(t_o) + \mathbf{v_e}(t_o) - \mathbf{v_r}(t)] = \frac{q_e \|\mathbf{c}(t_o) + \mathbf{v_e}(t_o)\|}{4\pi c^2 r^2(t, t_o)} \, [\mathbf{c}(t_o) + \mathbf{v_e}(t_o) - \mathbf{v_r}(t)].
$$

## Absorption Rate

Absorption is modelled as proportional to the incident flux, the
receiver's effective cross-sectional area
$\mathbf{A} = A \hat{\mathbf{o}}_e$ (where $\hat{\mathbf{o}}_e$ defines
the corpuscle orientation at emission), and absorption efficiency. The
absorbed corpuscle count per event is

$$
\text{Abs}(\mathbf{x_r}, t, t_o) = [\boldsymbol{\Phi}(\mathbf{x_r}, t, t_o) \cdot \mathbf{A}] \, (1 - e^{-\mu \delta t}),
$$

with transit time across the receiver

$$
\delta t = \frac{r_p}{\|\mathbf{c}(t_o) + \mathbf{v_e}(t_o) - \mathbf{v_r}(t)\|},
$$

and $r_p$ the receiver radius. For $r_p \to 0$ ($\mu \delta t \ll 1$),
this approximates to

$$
\text{Abs}(\mathbf{x_r}, t, t_o) = \frac{[\boldsymbol{\Phi}(\mathbf{x_r}, t, t_o) \cdot \mathbf{A}] \, \mu r_p}{\|\mathbf{c}(t_o) + \mathbf{v_e}(t_o) - \mathbf{v_r}(t)\|}.
$$

## Momentum Transfer and Action

Corpuscle absorption imparts momentum along the incident orientation.
The rate of change of receiver momentum satisfies

$$
\frac{d (m \mathbf{v_r}(t))}{dt} \propto \text{Abs}(\mathbf{x_r}, t, t_o) \left( \frac{[\mathbf{c}(t_o) + \mathbf{v_e}(t_o) - \mathbf{v_r}(t)] \cdot \hat{\mathbf{o}}_e(t_o)}{\|\mathbf{c}(t_o) + \mathbf{v_e}(t_o) - \mathbf{v_r}(t)\|} \right) \hat{\mathbf{o}}_e(t_o).
$$

Consolidating terms yields the action expression

$$
\frac{d (m \mathbf{v_r}(t))}{dt} \propto A r_p \frac{q_e \mu_r}{4 \pi r^2} \frac{\|\mathbf{c} + \mathbf{v_e}\|}{c^2} \left( \frac{\hat{\mathbf{o}}_e \cdot (\mathbf{c} + \mathbf{v_e} - \mathbf{v_r})}{\|\mathbf{c} + \mathbf{v_e} - \mathbf{v_r}\|} \right)^2 \hat{\mathbf{o}}_e,
$$

where arguments $(t, t_o)$ are implied for brevity.

Assuming geometric scaling $A \propto r_p^2$, the action becomes

$$
\frac{d (m \mathbf{v_r}(t))}{dt} \propto r_p^3 \frac{q_e \mu_r}{4 \pi r^2} \frac{\|\mathbf{c} + \mathbf{v_e}\|}{c^2} \left( \frac{\hat{\mathbf{o}}_e \cdot (\mathbf{c} + \mathbf{v_e} - \mathbf{v_r})}{\|\mathbf{c} + \mathbf{v_e} - \mathbf{v_r}\|} \right)^2 \hat{\mathbf{o}}_e.
$$

## Dimensional Consistency

-   Distance terms: $r^2$ [m²], $r_p^3$ [m³]
-   Rates: $q_e$ [corpuscles s⁻¹], $\mu_r$ [s⁻¹]
-   Velocities: $c, v$ [m s⁻¹]

The rate of corpuscle absorption (corpuscles s⁻¹) conforms dimensionally
to

$$
\text{corpuscles s}^{-1} \sim \text{m}^3 \cdot \frac{(\text{corpuscles s}^{-1}) \cdot (\text{s}^{-1})}{\text{m}^2} \cdot \frac{\text{m s}^{-1}}{(\text{m s}^{-1})^2},
$$

confirming that the derived action is proportional to the absorption
rate per unit time.

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
