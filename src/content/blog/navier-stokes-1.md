---
title: 'Navier-Stokes (1): Definition and Millenium Problem'
pubDate: 2026-09-25
description: 'TODO'
author: 'Arthur Feeney'
image:
    url: 'https://docs.astro.build/assets/rose.webp'
    alt: 'TODO'
tags: []
category: 'Navier-Stokes'
---

Since OpenAI announced they had a solution to the Millenium problem, I thought it would be fun
to "learn-in-public" a bit more about the (incompressible) Navier-Stokes equations. This is
hopefully a fun intro and going to be really informal!

$$
\begin{aligned}
\text{(1) Navier-Stokes:} & \quad \partial_t \mathbf{u} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\nabla p + \nu \Delta \mathbf{u} + \mathbf{f}\\
\text{(2) Continuity:} & \quad \nabla \cdot u = 0
\end{aligned}
$$

The Navier-stokes equation (1) describes the evolution of a velocity field $\mathbf{u}(\mathbf{x}, t)$, with a pressure $p(\mathbf{x}, t)$, viscosity $\nu$, 
and external forces $\mathbf{f}(\mathbf{x}, t)$. The n-dimensional vector $\mathbf{x} \in R^{n}$ is a spatial coordinate and $t \in R_{\geq 0} = [0, \infty)$ is the time. 
The external forces can be, for example, the force of gravity. (The OpenAI solution requires using the extremely weird forcing terms to cause "blowup".)
Generally, people are interested in starting from some initial velocity $u_0(\mathbf{x})$, plugging it into the Navier-Stokes equation, and looking at how the velocity
evolves over time. The "solution" is the pair functions ($\mathbf{u}(\mathbf{x}, t), p(\mathbf{x}, t)$) that arise from the choice of initial condition $u_0(\mathbf{x})$.
For example, you could choose $u_0(\mathbf{x}) := [\text{sin}(x_1), \text{sin}(x_2), \dots, \text{sin}(x_n)]$, where the initial velocity at $\mathbf{x}$ is just the $\text{sin}(\mathbf{x})$.

The continuity equation (2) is what makes it the "incompressible" Navier-Stokes. $\nabla \cdot \mathbf{u} = \sum_{i \in [n]} \partial_{\mathbf{x}_i} \mathbf{u}$ is the "divergence".
Since we have $\nabla \cdot \mathbf{u} = 0$, the velocity is "divergence-free" which essentially means that the amount of 
liquid flowing _out_ of a point must always be equal to the amount of liquid flowing _into_ a point.

I think a slightly more clear way to write this is just to isolate the time-derivative of the left-hand-side:

$$
\partial_t \mathbf{u} = -(\mathbf{u} \cdot \nabla)\mathbf{u} -\nabla p + \nu \Delta \mathbf{u} + \mathbf{f}
$$

$\partial_t \mathbf{u}$ is the derivative of the velocity with respect to time. So, this is just saying that the velocity changes over time based on the stuff on the 
right-hand-side of the equals sign. (This is also the starting point for how to solve the Navier-Stokes on a computer. The solution function is
the time-integral of both sides: $u(\mathbf{x}, t) = \int_{s\in[0, t]} \partial_s \mathbf{u}(\mathbf{x}, s) \text{d}s$, 
and one can evaluate the integral forward in time by starting from $\mathbf{u}_0$. (The integral can be very hard to compute!)) 

## What's the right-hand side?

The terms on the right-hand-side of the equals sign have very not-at-all-obvious interpretations in my opinion.

The "convective term" $-(\mathbf{u} \cdot \nabla)\mathbf{u}$ uses a slightly funky operator: 
$\mathbf{u} \cdot \nabla := \mathbf{u}_1 \partial_{\mathbf{x}_1} + \mathbf{u}_2 \partial_{\mathbf{x}_2} + \dots + \mathbf{u}_n \partial_{\mathbf{x}_n}$.
Applying this operator to $\mathbf{u}$ is then $\sum_{i\in[n]} \mathbf{u}_i \partial_{\mathbf{x}_i} \mathbf{u}$.
This is the directional derivative of the velocity with respect to the velocity.
One way to intrepret this is that it describes how the momentum of fluid carries the itself... (I think this is hard to put into English.)

- As a slightly related aside, $\nabla \mathbf{u} = (\partial_{\mathbf{x}_i}\mathbf{u}(\mathbf{x}_j))_{ij}$ is the Jacobian (a matrix with rows and columns indexed by $i$ and $j$). Since $\nabla \cdot \mathbf{u} = 0$, the trace of the Jacobian $Tr(\nabla \mathbf{u}) = 0$. 
- This is the only non-linear part of the Navier-Stokes equations and is the main reason why analysis has been so difficult.

In the "viscous term" $\nu \Delta \mathbf{u}$, the positive scalar $\nu \in R_{>0}$ is the viscosity. 
The operator $\Delta := \nabla \cdot \nabla$ is the Laplacian. The Laplacian compares a point $\mathbf{x}$ with neighboring points.
If the velocity at the point $\mathbf{x}$ is generally lower than its neighbors, the Laplacian will be positive and so the liquid will speedup. 
Otherwise, the Laplacian will be negative and the liquid is slowed down. The Laplacian will essentially push the velocity at a point towards a local average. 
The viscosity $\nu$ essentially controls the rate of averaging. So, if the viscosity $\nu$ is very low,
there will be a very slight averaging. Generally, low viscosity fluids will have a more complex flows.

The pressure gradient $\nabla p$ is used to capture how neighboring parts of the fluid push on each other.
In incompressible flow, the pressure is basically just whatever is necessary to keep $\nabla \cdot \mathbf{u} = 0$.
The pressure $p(\mathbf{x}, t)$ is actually derivable from the velocity. Take the divergence of the Navier-Stokes equation:
- $\nabla \cdot \partial_t \mathbf{u} = \partial_t (\nabla \cdot \mathbf{u}) = 0$ since derivatices commute and incompressibility.
- $\nabla \cdot (\nu \Delta \mathbf{u}) = \nu \Delta(\nabla \cdot \mathbf{u}) = 0$ since derivatives commute and incomressibility. 
- $\nabla \cdot \nabla p = \Delta p$, the Laplacian of the pressure.
- $\nabla \cdot (\mathbf{u} \cdot \nabla) \mathbf{u}$ is _something_! This is more complex since it is nonlinear, and it does not reduce to zero.

So, $\Delta p = \nabla \cdot (\mathbf{u} \cdot \nabla) \mathbf{u}$, which is a Poisson problem that can be solved for $p$.
It's possible to reduce the right-hand side to something nicer looking, but the crux is that the pressure at a point $\mathbf{x}$ and time $t$
can be computed using only spatial derivatives of $\mathbf{u}$ at the time $t$. I.e., at a given timestep, you can use the velocity to compute the pressure that
will make the system remain incompressible. $\nabla p$ is basically chosen to cancel out any divergence that would be added by the convective term!


## Millenium Prize Problem Description

The problem defined by the Clay Mathematics Institute is essentially the same as the above, but adds some additional restrictions on
what counts as a solution: 
$$
\begin{aligned}
\text{(3) smoothness:}& \quad \mathbf{u}, p \in C^{\infty}(R^n \times [0, \infty)) \\
\text{(4) Bounded Energy:}& \quad ||\mathbf{u}(\cdot, t)||_{L^2} = \int_{R^n} |u(\mathbf{x}, t)|^2 \text{d}\mathbf{x} < C \text{ for all } t > 0
\end{aligned}
$$


In (3) $C^{\infty}$ is the set of smooth functions, defined on the domain $R^n \times [0, \infty)$.
Statement (4) says that the kinetic energy must remain finite for all timesteps.

A pair ($\mathbf{u}, p)$ is a valid solution if it satisfies all of equations (1), (2), (3), (4), and 
the initial $u(\cdot, 0) = u_0$. The millenium problem is asking for a proof (or disproof) that for _all_ valid 
initial conditions $u_0$ and forces $f$, there is exists a solution that satisfies (1), (2), (3), and (4).

The clay mathematics institute wanted a proof for one of four statements:
- (A) Existence and smoothness of Navier-Stokes on $R^3$.
- (B) Existence and smoothness of Navier-Stokes on $R^3 / Z^3$.
- (C) Breakdown of Navier-Stokes on $R^3$.
- (D) Breakdown of Navier-Stokes on $R^3 / Z^3$.

($R^3 / Z^3$ is a periodic domain, which I am not going to define or discuss.)

OpenAI has shown (C) and (D) are true, disproving (A) and (B) by counter-example. 
They start from the initial condition $\mathbf{u}_0 = 0$, but construct really specific forcing functions $f$. 
Their paper is written by an LLM, so I don't really want to study it in detail. I'm not aware of a version that
is tidied up, yet.

It has been known for a while that if $[0, \infty)$ is replaced with $[0, T_{\text{max}}]$, then a solution always exists ($T_{\text{max}}$ is problem dependent). 
If there is a valid solution for time $[0, T_{\text{max}}]$, but it then the solution breaks down (i.e., becomes unsmooth or unbounded), then $T$ is called the "blowup time."
OpenAI showed that for all positive viscocities and $\mathbf{u}_0 = 0$, they can construct a forcing $f$ such that the $\mathbf{u}$ produced via
equation (1) and (2) has finite time blow up.

I think it's sort of interesting that this result is probably not interesting or useful at all to physicists. It essentially relies on
constructing a (presumably) physically meaningless forcing function. It is definitely interesting to mathematicians, but one question I have is if
this is really useful to mathematicians? I.e., will a mathematician be able take this approach and extend it to other problems?

It seems like one of the next questions for Navier-Stokes is if (C) and (D) are still true for $f = 0$.
