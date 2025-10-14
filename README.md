# Finite-Sized Particle Interaction Model

This repository presents a novel model of particle-particle interactions, demonstrating how **electrostatic**, **magnetic**, and **gravitational forces** emerge from a Galilean framework of finite-sized particles emitting and absorbing **vector corpuscles**. The model challenges conventional physics by proposing that these forces arise from corpuscle exchanges, potentially eliminating the need for dark matter in explaining galactic rotation curves.

## Model Overview

The model assumes a universe filled with **vector corpuscles**—hypothetical entities carrying position, velocity, and orientation—that mediate interactions between finite-sized particles. Each corpuscle is characterized by:

-   **Position**: $\mathbf{x_c}$,

-   **Velocity**: $\mathbf{v_c}$,

-   **Orientation**: $\mathbf{o_c}$ a unit vector.

The corpuscle velocity follows **Galilean relativity**:\

$$ \mathbf{v_c} = \mathbf{v_p} + \mathbf{c},$$

where$\mathbf{v_p}$ is the velocity of the emitting particle, and $\mathbf{c}$ is a vector with magnitude equal to the speed of light $c$.

## Finite-Sized Particles

The model defines two types of particles, each emitting $q$ corpuscles per second:

-   **Positive particles** $\mathbf{p^+}$ : Emit corpuscles with orientation **parallel** to their velocity:

$$ \mathbf{o_c} \cdot \mathbf{c} = c.$$

-   **Negative particles** $\mathbf{p^-}$ : Emit corpuscles with orientation **anti-parallel** to their velocity:

$$ \mathbf{o_c} \cdot \mathbf{c} = -c.$$

Particles also **absorb** corpuscles from the surrounding medium. The absorption process generates an **action** that alters the particle's velocity based on the corpuscles' density, velocity, and orientation:

-   For $\mathbf{p^+}$ , the action aligns with the absorbed corpuscle’s orientation.

-   For $\mathbf{p^-}$ , the action opposes the absorbed corpuscle’s orientation.

## Electrostatic Force

The net force between two particles with charges $q_1$ and $q_2$ (electron charges $e$ separated by distance $r$) is:\

$$ \mathbf{f_2} = \frac{k q_1 q_2}{4 \pi r^2} \frac{\|\mathbf{c} + \mathbf{v_1}\|}{c}   \frac{\mathbf{c} \cdot ( \mathbf{c} + \mathbf{v_1} - \mathbf{v_2} )}{\|\mathbf{c} + \mathbf{v_1} - \mathbf{v_2}\|} \right)^2 \hat{o_1},$$

where:

-   $k$: Coulomb constant,

-   $\mathbf{v_1}, \mathbf{v_2}$: Velocities of the two particles,

-   $\hat{o_1}$: Unit vector of the corpuscle orientation

This force generalizes Coulomb’s law by incorporating velocity-dependent effects.

## Mass and Inertia

The **inertia** (mass) of a finite-sized particle is proportional to its radius. For charged particles, the mass is related to the **electrostatic stored energy** divided by ( c\^2 ):\

$$ m \propto \frac{\text{Electrostatic energy}}{c^2}.$$

## Magnetic Force

The **magnetic force** emerges from interactions between **neutral currents** (e.g., moving neutral particles) and charged particles, driven by the exchange of corpuscles with specific orientations.

## Gravitational Force

**Gravity** arises as a net attractive force between **neutral composite particles** (e.g., atoms with a positively charged nucleus and a negatively charged shell). The slight difference in velocity distributions between positive and negative charges results in a small but consistent attraction.

## Galactic Rotation Curves

The velocity-dependent nature of the gravitational force implies that **"hot" ionized particles** (with high velocities) behave differently from **"cold" neutral particles**. This difference accounts for the observed **flattening of galactic rotation curves**, potentially eliminating the need for dark matter in this model. Simulations in this repository demonstrate how these effects align with observed galactic dynamics.
