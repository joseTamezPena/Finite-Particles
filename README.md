# Finite-Sized Particle Interaction Model

This repository is a work in progress and contains the MATLAB scripts of a novel model of particle-particle interactions, demonstrating how **electrostatic**, **magnetic**, and **gravitational forces** emerge from a Galilean framework of finite-sized particles emitting and absorbing **vector corpuscles**. The model challenges conventional physics by proposing that these forces arise from corpuscle exchanges, potentially eliminating the need for dark matter in explaining galactic rotation curves.

## Model Overview

The model assumes a universe filled with finite-sized **vector corpuscles**—hypothetical entities carrying position, velocity, and orientation—that mediate interactions between finite-sized particles. Each corpuscle is characterized by:

-   **Size**: $l$

-   **Position**: $\mathbf{x_c}$,

-   **Velocity**: $\mathbf{v_c}$,

-   **Orientation**: $\mathbf{o_c}$ a unit vector.

The corpuscle velocity follows **Galilean relativity**:

$$ \mathbf{v_c} = \mathbf{v_p} + \mathbf{c},$$

where $\mathbf{v_p}$ is the velocity of the emitting particle, and $\mathbf{c}$ is a vector with magnitude equal to the speed of light $c$. The orientation, $\mathbf{o_c}$, is defined at the moment the corpuscle is emitted by the particle.

## Finite-Sized Particles

The model defines two types of particles, each emitting $q$ corpuscles per second:

-   **Positive particles** $\mathbf{p^+}$ : Emit corpuscles with orientation **parallel** to their velocity:

$$ \mathbf{o_c} \cdot \mathbf{c} = c.$$

-   **Negative particles** $\mathbf{p^-}$ : Emit corpuscles with orientation **anti-parallel** to their velocity:

$$ \mathbf{o_c} \cdot \mathbf{c} = -c.$$

The emission of the corpuscles is isotropic, and at the speed of light $c$.

Particles **absorb** the corpuscles from the surrounding medium. The absorption process generates an **action** that alters the particle's velocity based on the corpuscles' density, velocity, and orientation:

-   For $\mathbf{p^+}$ , the action aligns with the absorbed corpuscle’s orientation.

-   For $\mathbf{p^-}$ , the action opposes the absorbed corpuscle’s orientation.

# Corpuscular Action Model for Particle–Particle Interaction

This section describes the mechanistic action model describing the interaction between a receiver particle and corpuscles emitted by an independent emitting particle. The model is grounded in a corpuscular paradigm, wherein the universe is permeated by discrete corpuscles propagated at the speed of light. Particles emit and absorb these corpuscles at characteristic rates, mediating momentum transfer analogous to a radiative exchange process.

## Model Foundations

The interaction framework posits that emitting particles radiate corpuscles isotropically at a rate of $q$ corpuscles per second, with propagation occurring at the vacuum speed of light $c$. Receiver particles absorb incident corpuscles at a rate parameterised by $\mu$ (dimensions of inverse time), contingent upon the corpuscular flux intersecting the receiver.

Observations are conducted in a coordinate system centered on a stationary observer, from whose frame the velocity evolution of the receiver particle is described.

The system under consideration comprises a single emitting particle $\mathbf{p_e}$ and a single receiver particle $\mathbf{P_r}$, both potentially in motion. At emission time $t_o$, the emitter occupies position $\mathbf{x_e}(t_o)$ with velocity $\mathbf{v_e}(t_o)$, radiating at rate:

$$
q_e = \frac{4}{3} \pi r^3_p\rho_u\mu_e;
$$

where $\rho_u$ is the universal vacuum corpuscle density, $\mu_e$ the emission rate, and $r_p$ the particle radius. At absorption time $t$, the receiver is at $\mathbf{x_r}(t)$ with velocity $\mathbf{v_r}(t)$.

## Corpuscular Density

The relative observer density of corpuscles at the receiver position is

$$
\rho(\mathbf{x_r}, t, t_o) = \frac{q_e}{4\pi (c \Delta t)^2 [\mathbf{c}(t_o) + \mathbf{v_e}(t_o)] \cdot [\mathbf{c}(t_o)/c]},
$$

where $\Delta t$ is the courpuscle travel time, and $\mathbf{c}(t_o)$ is the light-speed vector from emitter to receiver in the observer frame. The emission-reception separation distance is $r(t, t_o) = \|\mathbf{x_r}(t) - \mathbf{x_e}(t_o)\|$, and the relative propagation speed is $\|\mathbf{c}(t_o) + \mathbf{v_e}(t_o)\|$. Thus,

$$
\Delta t = \frac{r(t, t_o)}{\|\mathbf{c}(t_o) + \mathbf{v_e}(t_o)\|}.
$$

Substitution provides the corpuscle density at the reciever

$$
\rho(\mathbf{x_r}, t, t_o) = \frac{q_e \|\mathbf{c}(t_o) + \mathbf{v_e}(t_o)\|^2}{4\pi c^2 r^2(t, t_o)[\mathbf{c}(t_o) + \mathbf{v_e}(t_o)] \cdot [\mathbf{c}(t_o)/c]}.
$$

## Absorption Rate

The corpuscle absorption rate, $\mathbf{Ar}$, is modeled as proportional to the corpuscle density, the receiver's absorption coefficient $\mu$ in a small volume:

$$
\mathbf{Ar}(\mathbf{x_r}, t)=\mu \rho(\mathbf{x_r}, t, t_o) \Delta_v
$$

The flux of corpuscles travel a small distance in a given $dt$. Hence:

$$
d\mathbf{Ar}(\mathbf{x_r}, t)=\mu \rho(\mathbf{x_r}, t, t_o) dA [v dt]
$$

The total number of absorbed corpuscles per second at a specific point of the receiver approximates to:

$$
\frac{d\mathbf{Ar}(\mathbf{x_r}, t)}{dt}=\mu \rho(\mathbf{x_r}, t, t_o) v dA ;
$$

where $v$ is the magnitude of relative velocity of the corpuscles at that point.

## Particle Acceleration

The absorption rate equation computes the number of absorbed corpuscles per unit of time. The model assumes that each corpuscle has a small length and orientation. Therefore, each time an corpuscle is absorbed there is a change in the receiver particle position. The rate of change is:

$$
\frac{\Delta x(\mathbf{x_r}, t)}{dt} = \frac{d \mathbf{Ar}(\mathbf{x_r}, t)}{dt}\hat{l}= \mu \rho(\mathbf{x_r}, t , t_o) v dA \hat{l}.
$$

Hence:

$$
\frac{\Delta x(\mathbf{x_r}, t)}{dt} =  \frac{ \mu q_e \|\mathbf{c}(t_o) + \mathbf{v_e}(t_o)\|^2}{4\pi c^2 r^2(t, t_o)[\mathbf{c}(t_o) + \mathbf{v_e}(t_o)] \cdot [\mathbf{c}(t_o)/c]} v dA \hat{l}.
$$

The dimensional analysis yields that the rate of change of a specific location of the the receiver partícle is: m/s\^2. In other words: A corpuscle absortion causes a small acceleration of the receiver particle.

## Newton's First Law

In the absence of a external corpuscles, $q_e=0$, the absorption is zero, hence:

$$
\frac{\Delta x(\mathbf{x_r}, t)}{dt} = 0
$$

In other words there is no change in the particle observed velocity.

"An object at rest stays at rest, and an object in motion stays in motion with the same speed and in the same direction, unless acted upon by an unbalanced external force."

## Newton's Second Law 

The absorbed corpuscles from external sources and the internal absorbed should be in equilibrium.

$$ 
\frac{\Delta x_{ext}(\mathbf{x_r}, t)}{dt}=\frac{\Delta x_{int}(\mathbf{x_r}, t)}{dt}
$$

Here we derive the Newton's second law. First I'll assume that the absorbing particle is finite and composed by two absorbing/emitting elements separated by a distance $L$ aligned to the external flow of corpuscles. I'll assume that the internal elements have half the particle charge. The external force causes an small acceleration of the receiver particle. Therefore the front particle and the back particle at the time of action summation have a small velocity: $v_e=a \Delta t$; where $a$ is the instant acceleration due to the external action and $\Delta t=L/c$ is the amount of time required to travel the distance between the two particle elements.

The internal absorption of the dual element particle is:

$$
\frac{\Delta x_{int}(\mathbf{x_r}, t)}{dt} = \frac{q_e \mu c \Delta V  \hat{o_e}}{8 \pi} \left( \frac{1}{(c+a\Delta t)(L-1/2 a \Delta t^2)^2} - \frac{1}{(c-a\Delta t)(L+1/2 a \Delta t^2)^2} \right).
$$

The external absorption is:

$$
\frac{\Delta x_{ext}(\mathbf{x_r}, t)}{dt}=\frac{ \mu Q \|\mathbf{c}(t_o) + \mathbf{v_e}(t_o)\|^2}{4\pi c^2 r^2(t, t_o)[\mathbf{c}(t_o) + \mathbf{v_e}(t_o)] \cdot [\mathbf{c}(t_o)/c]} c \Delta V \hat{o_e}
$$

Substituting we get:

$$
\frac{\mu Q  \Delta V \|c+v_e\|}{4 \pi c R^2 cos(\theta)} = \frac{q^2}{4 \pi L  c^2}a \hat{\mathbf{o}}_e,
$$

where the left hand term is the external force induced by the external charge, and the right terms originates by the internal pull of the accelerated particle.

A close inspection of the equation it reveals that it is Newton's Second law: $F=ma$, for charged finite particles.

A second observation is that the first terms of the internal absorption is equivalent to the electrostatic energy:

$$
E_{in} = \frac{q^2}{4 \pi \epsilon_o L}.
$$

Hence, the inertial mass is:

$$
 m = \frac{E_{in}}{c^2}.
$$

The inertia and particle mass are a direct consequence of the corpuscular model on finite particles.

## Rate of emission-absorption

The internal energy equation $E_{in}$ can be used to estimate the rate of corpuscles absorbed-emitted at any given time by real particles. I'll assume the electron charge $e = 1.602176634 \times 10^{-19}$ C, and radius $r_e=10^{-20}$ m

$$
\rho_u r_p^6 \mu^2 = \frac{q^2}{\epsilon_o} \therefore \rho_u \mu^2 = \frac{e^2}{r_e^6 \epsilon_o} = 2.307 \times 10^{32} \text{ m}^{-3}\text{s}^{-2}.
$$

This implies that an electron on average has $2.307 \times 10^{32}$ corpuscular interactions per second squared per cubic meter.

## Electrostatic Force

The derived model can be described as a net force of two charged particles

$$
\mathbf{F}=\int \frac{\Delta x(\mathbf{x_r}, t)}{dt} dv
$$

of the particle absorbing corpuscles $q_2$ (Absorbing) of an emitting charge $q_1$ (Emitting) and (e.g., electron charges $e$ separated by distance $r$ (at emitting-absorbing) is:

$$ \mathbf{F} = \frac{k q_1 q_2 c}{4 \pi r^2 (\mathbf{c} + \mathbf{v_e}) \cdot \mathbf{\hat{o}} } \left( (1 - \frac{v_e \cdot v_r}{c^2} )\hat{\mathbf{o}} + \frac{v_r \cdot \hat{o}}{c^2} v_e \right),
$$

where:

-   $k$: Coulomb constant,

-   $q$ is positive ($>0$) for $\mathbf{p^+}$ and negative ( $<0$) for $\mathbf{p^-}$

-   $\mathbf{v_e}, \mathbf{v_r}$: Velocities of the two particles:

    -   $\mathbf{v_e}$ is an emitting particle. $\mathbf{v_r}$ the absorbing particle

-   $\hat{o_1}$: Unit vector of the corpuscle orientation at the time of its origin by the emitting particle.

-   $r$ is the distance traveled by the corpuscle form the emission origin to it absorption:

$$
    r=\| \mathbf{x_2}(t)-\mathbf{x_1}(t_o) \|.
$$

This force is equivalent to the Total Lorentz force of two moving charged particles.

When the two particles are stationary:

$$ \mathbf{F} = \frac{k q_1 q_2}{4 \pi r^2} \hat{o_1}, $$

The is the standard Coulomb's law.

## Magnetic Force

The **magnetic force** and the **vacuum magnetic permeability** emerge from interactions between **currents** (e.g., moving charged particles in a neutral cables) and moving charged particles.

The presented corpuscular model correctly predict that the **vacuum magnetic permeability** $\mu_0$ and the electric constant are associated by:

$\epsilon_0 =\frac{1}{c^2 \mu_0 }$

The evidence is presented in the Magnetic Force folder for a current in a loop yields the exact equation predicted by the Lorentz force of a moving charged particle $Q$ in the center of a circular loop with current $I$ and radius $R$:

$$
F_m=\left(\begin{array}{ccc}
0 & -\frac{\textrm{I} Q \mu_0 v }{2 R} & 0
\end{array}\right)
$$

## Gravitational Force

**Gravity** arises as a net attractive force between **neutral composite particles** (e.g., atoms with a positively charged nucleus and a negatively charged shell). The slight difference in stochastic velocity distributions between positive and negative charges results in an Expected small attraction. The equations of interaction are derived in the GravityForce/Steady-sate folder

## Galactic Rotation Curves

The velocity-dependent nature of the gravitational force implies that **"hot" ionized particles** in neutral plasma (with high velocities) behave differently from **"cold" neutral particles**. This difference accounts for the observed **flattening of galactic rotation curves**, potentially eliminating the need for dark matter in this model. Simulations in this repository demonstrate how these effects align with observed galactic dynamics. The results of fitting observed rotation curves of 22 SPARC galaxies and the local group are showcased in the Galaxies folder.

![M33 Rotation Curve](images/clipboard-1366545713.png){width="392"}

![M33 Surface mass densities](images/clipboard-149743154.png){width="383"}

![Rotation Curve](Galaxies/SPARC/DensityProfiles/NGC7331.jpg)

The evaluation of the 22 galaxies implies that there is no need for dark matter for the explanation of the rotation curves:

![](images/clipboard-4116685001.png)

Collaboration is welcome. Please send me comments, corrections and suggestions.
