---
title: 'Navier-Stokes (2): Local and Non-Local Dynamics'
pubDate: 2026-10-16
description: 'Which terms of the Navier-Stokes equations only look at a neighborhood of a point, and which ones look at the whole fluid.'
author: 'Arthur Feeney'
image:
    url: 'https://docs.astro.build/assets/rose.webp'
    alt: 'TODO'
tags: []
category: 'Navier-Stokes'
---

In the [last post](/posts/navier-stokes-1), I went through the terms in the (incompressible) Navier-Stokes equations.
Some of the terms only care about what is happening right around a point $\mathbf{x}$, and some of the terms care about what is happening
everywhere. I think it is extremely unintuitive from looking at the equation, and is also kind of unphysical since it implies that some
phenomena propagate across the fluid with infinite speed.

This post assumes the forcing $\mathbf{f} = 0$ and space is $R^3$:

$$
\begin{aligned}
\text{(1) Navier-Stokes:} & \quad \partial_t \mathbf{u} = -(\mathbf{u} \cdot \nabla)\mathbf{u} - \nabla p + \nu \Delta \mathbf{u}\\
\text{(2) Continuity:} & \quad \nabla \cdot \mathbf{u} = 0
\end{aligned}
$$

## What does local mean?

An operator $A : \mathcal{U} \rightarrow \mathcal{V}$ maps a function $u$ to some other function $v$. I.e., $A[u] = v$. In these examples,
$u$ can be the velocity (and so is a function of space and time: $u(\mathbf{x}, t)$). $A$ can be any of the operators used in the navier-stokes
equations, like the gradient or divergence.

An $A$ is _local_ if the value $A[\mathbf{u}](\mathbf{x}) = v(\mathbf{x})$ only depends on $\mathbf{u}$ in an (arbitrarily small) neighborhood around $\mathbf{x}$.
If $\mathbf{u}$ were changed somewhere far away from $\mathbf{x}$, $A[\mathbf{u}](\mathbf{x})$ would not change at all.
Anything built out of derivatives at a point is local: derivatives are defined by a limit that only looks at points closer and closer to $\mathbf{x}$.

An operator is _non-local_ if $A[\mathbf{u}](\mathbf{x})$ depends on $\mathbf{u}$ at points that are far away from $\mathbf{x}$. The typical example is a kernel
integral operator:

$$
\mathbf{v}(\mathbf{x}) = A[\mathbf{u}](\mathbf{x}) = \int_{R^3} K(\mathbf{x} - \mathbf{y}) \mathbf{u}(\mathbf{y}) \text{d}\mathbf{y}
$$

If the kernel $K$ is non-zero far from the origin, then every point $\mathbf{y} \in R^3$ contributes something to the value at $v(\mathbf{x})$.
How much "far away" matters is basically controlled by how fast $K$ decays.

## Going through the terms

**1. The convective term** $-(\mathbf{u} \cdot \nabla)\mathbf{u} = -\sum_{i} \mathbf{u}_i \partial_{\mathbf{x}_i} \mathbf{u}$ is local.
It only uses the value of $\mathbf{u}$ and its first derivatives at $\mathbf{x}$. It's nonlinear, but it's local.

**2. The viscous term** $\nu \Delta \mathbf{u}$ is also local, since it's just second derivatives at $\mathbf{x}$.
- This is actually much more complicated and technically does introduce, what I consider, global effects...! This is more complicated
and explained at the end of the blog. 

**3. The pressure term** $-\nabla p$ is global. As mentioned in the previous post, the pressure at any time $t$ can be found by taking the
divergence of every term in the navier-stokes equation. This results in the Laplacian of the pressure, equal to the negative convective term.

$$
\Delta p = -\nabla \cdot (\mathbf{u} \cdot \nabla)\mathbf{u}
$$

To actually get $p$, you have to invert the Laplacian $\Delta$. In $R^3$, the inverse Laplacian is an integral against the kernel $-\frac{1}{4\pi |\mathbf{x}|}$, so

$$
p(\mathbf{x}) = \frac{1}{4\pi} \int_{R^3} \frac{1}{|\mathbf{x} - \mathbf{y}|} \sum_{i,j} \partial_{\mathbf{y}_i} \partial_{\mathbf{y}_j} (\mathbf{u}_i \mathbf{u}_j)(\mathbf{y}) \text{d}\mathbf{y}
$$

This is extremely similar to the non-local operator, with non-zero kernel, shown above. The pressure at $\mathbf{x}$ depends on the velocity _everywhere_, and the kernel
$K(\mathbf{x} - \mathbf{y}) = \frac{1}{|\mathbf{x} - \mathbf{y}|}$ actually decays quite slowly.

Another way to think about it: if you modify the speed of the fluid at one point, in order for the whole fluid to stay divergence free, any divergence 
introduced by the modification has to be cancelled out. To do this cancellation, the fluid has to adjust instantly at every other point.
Real fluids are slightly compressible and information travels across the fluid at the speed of sound. Incompressible flow acts like the speed of 
sound is infinite, so the adjustment from the pressure is propagated with infinite speed.

## Hiding the pressure (Leray projection)

Since the pressure is determined by $\mathbf{u}$, you can get rid of it entirely. The Leray projection $\mathbb{L}$ takes
a vector field and removes its "divergence part", leaving only a divergence-free field:

$$
\mathbb{L} := I - \nabla \Delta^{-1} \nabla \cdot
$$

If we apply the Leray projection to every term in the Navier-Stokes equation, we end up with
- $\mathbb{L}[\partial_t \mathbf{u}] = \partial_t \mathbf{u}$, since $\mathbf{u}$ is already divergence free.
- $\mathbb{L}[\nu \Delta \mathbf{u}] = \nu \Delta \mathbf{u}$, since it is also already divergence free.
- $\mathbb{L}[(\mathbf{u} \cdot \nabla) \mathbf{u}]$, is something else and non-zero, since the convective term has non-zero divergence.
- $\mathbb{L}[\nabla p] = 0$.

So, after applying the Leray projection, we can write the incompressible Navier-Stokes like this

$$
\partial_t \mathbf{u} = \nu \Delta \mathbf{u} - \mathbb{L}[(\mathbf{u} \cdot \nabla)\mathbf{u}]
$$

Now there's no pressure term, but the nonlinear term is non-local since the Leray projection uses the very non-local inverse Laplacian. 
You also do not need the extra continuity equation, since divergence-free-ness is imposed with the Leray projection.
This form is like a heat equation ($\partial_t \mathbf{u} = \nu \Delta \mathbf{u}$), but with an extra complicated forcing term.
It seems like a lot of analysis uses this form of the equation.

## Non-locality(?!) of the viscous term

I think this is extremely not obvious. The viscous term also has non-local effects. If you think about a computer simulation,
it just uses second derivatives, which can be implemented using a small, local finite-difference stencil..., 
so how does it have non-local effects...?! Again, this is going to be pretty informal.

Integrating the Navier-Stokes equation over time, we get

$$
\mathbf{u}(\cdot, t) = \mathbf{u}_0 + \int_{[0, t]} \nu \Delta \mathbf{u}(\cdot, s) ds - \int_{[0, t]} \mathbb{L}[(\mathbf{u} \cdot \nabla) \mathbf{u}](\cdot, s) ds
$$

We'll replace the integral of the non-linear part with $I_1(t) = \int_{[0, t]}\mathbb{L}[(\mathbf{u} \cdot \nabla) \mathbf{u}](\cdot, s) ds$.
This gives the solution of the form

<a id="eq-u-integral"></a>
$$
\mathbf{u}(\cdot, t) = \mathbf{u}_0 + \int_{[0, t]} \nu \Delta \mathbf{u}(\cdot, s) ds - I_1(t) \quad \tag{Soln. 1}
$$

The integral and $\Delta$ commute, so can rewrite the integral of the viscous term:

$$
\nu \Delta \int_{[0, t]} \mathbf{u}(\cdot, s) ds
$$

This is an integral with $\mathbf{u}(\cdot, s)$, which can be replaced with equation [(Soln. 1)](#eq-u-integral). 
After doing this replacement, we get 

$$
\begin{aligned}
\nu \Delta \int_{[0, t]} \mathbf{u}(\cdot, s) ds &= \nu \Delta \int_{[0, t]} [\mathbf{u}_0 + \int_{[0, s]} \nu\Delta\mathbf{u}(\cdot, r) dr - I_1(s)] ds \\
&= t \nu \Delta \mathbf{u}_0 + \nu \Delta \int_{[0, t]} \int_{[0, s]} \nu\Delta\mathbf{u}(\cdot, r) dr ds - \nu \Delta \int_{[0, t]} I_1(s) ds
\end{aligned}
$$

We can simplify the term $\nu \Delta \int_{[0, t]} I_1(s) ds$. We just define a new $I_2(t) = \int_{[0, t]}\nu\Delta I_1(s) ds$.

$$
\begin{aligned}
& t \nu \Delta \mathbf{u}_0 + \nu \Delta \int_{[0, t]} \int_{[0, s]} \nu\Delta\mathbf{u}(\cdot, r) dr ds - I_2(t) \\
= \; & t \nu \Delta \mathbf{u}_0 + \nu^2 \Delta^2 \int_{[0, t]} \int_{[0, s]} \mathbf{u}(\cdot, r) dr ds - I_2(t)
\end{aligned}
$$

$\Delta^n$ means $n$ applications of the Laplace operator: $\Delta^n[f] = \Delta(\Delta (\dots \Delta f))$.

You will notice that this is quite similar to the original form, where there is a $\mathbf{u}(\cdot, r)$ inside the integrals.
We can essentially repeat this process a second time:

$$
t \nu \Delta \mathbf{u}_0 + \frac{t^2}{2} \nu^2 \Delta^2 \mathbf{u}_0 + \nu^3 \Delta^3 \int_{[0,t]}\int_{[0, s]} \int_{[0, r]} \mathbf{u}(\cdot, q) dq dr ds - I_2(t) - I_3(t)
$$

Repeating this process $n-1$ times, we find

$$
\mathbf{u}(\cdot, t) = \sum_{i=0}^{n-1} \frac{(t\nu\Delta)^i}{i!} \mathbf{u}_0 + (\nu \Delta)^n \int_{[0, t]} \int_{[0, s_1]} \dots \int_{[0, s_{n-1}]} \mathbf{u}(\cdot, s_n) ds_n \dots ds_1 - \hat{I}_n(t)
$$

Where we use $\hat{I}_n := \sum_{i=1}^n I_i$ to remove the individual $I_i$ functions.

Using [Cauchy's formula](https://en.wikipedia.org/wiki/Cauchy_formula_for_repeated_integration) the nested integrals can be simplified to

$$
\int_{[0, t]} \int_{[0, s_1]} \dots \int_{[0, s_{n-1}]} \mathbf{u}(\cdot, s_n) ds_n \dots ds_1  = \frac{1}{(n-1)!} \int_{[0, t]} (t - s)^{n-1}\mathbf{u}(\cdot, s) ds
$$

This gives a solution $\mathbf{u}$ that looks like 

$$
\mathbf{u}(\cdot, t) = \sum_{i=0}^{n-1} [\frac{(t\nu\Delta)^i}{i!} \mathbf{u}_0] + \frac{(\nu \Delta)^n}{(n-1)!} \int_{[0, t]} (t-s)^{n-1} \mathbf{u}(\cdot, s) ds - \hat{I}_n(t)
$$

You may notice that as $n$ goes to infinity (i.e., we apply this process infinitely)
- (1) $\sum_{i=0}^\infty \frac{(t\nu\Delta)^i}{i!}\mathbf{u}_0 = e^{t\nu\Delta} \mathbf{u}_0$. (Note that this is the _operator_ $e^{t\nu\Delta}$ applied to $\mathbf{u}_0$.)
- (2) $\frac{(\nu\Delta)^n}{(n-1)!} \int_{[0, t]}(t-s)^{n-1} \mathbf{u}(\cdot, s) ds$, the division by $(n-1)!$ can make the whole thing go to zero as $n$ increases. 

So we further simplify it to 

$$
\mathbf{u}(\cdot, t) = e^{t\nu\Delta} \mathbf{u}_0 - \hat{I}_\infty(t)
$$

If we can figure out what $\hat{I}_n(t)$ is as $n$ goes to infinity, we can end up with a really nice form...

Recall, the $I_n$ are essentially defined inductively: 
- (1) $I_1(t) = \int_{[0, t]} F(s) ds$, where $F$ is the non-linear term.
- (2) $I_{n+1}(t) = \int_{[0, t]} \nu \Delta I_n(s) ds$

All of the $\nu \Delta$ can be factored out the integrals, so we get another repeated integral $I_{n+1}(t) = (\nu\Delta)^n \int_{[0, t]} \int_{[0, s_1]} \dots \int_{[0, s_n]} F(s_{n+1}) ds_{n+1} \dots ds_1$.
We can apply Cauchy's formula here as well. This gives

$$
I_{n+1}(t) = \frac{(\nu\Delta)^n}{n!} \int_{[0, t]} (t - s)^{n} F(s) ds
$$

So, $\hat{I}_n(t) = \sum_{i=1}^n \frac{(\nu\Delta)^{i-1}}{(i-1)!} \int_{[0, t]} (t-s)^{i-1} F(s) ds$. We can then move the integral out.

$$
\hat{I}_n(t) = \int_{[0, t]} \Big[ \sum_{i=1}^n \frac{(\nu\Delta)^{i-1}}{(i-1)!} (t-s)^{i-1} F(s) \Big] ds
$$

If we let $n$ go to infinity (and assume that the series converges),

$$
\begin{aligned}
\hat{I}_\infty(t) &= \int_{[0, t]} \Big[ \sum_{i=1}^\infty \frac{(\nu\Delta(t-s))^{i-1}}{(i-1)!} F(s) \Big] ds \\
                  &= \int_{[0, t]} \Big[ e^{(t-s)\nu\Delta} F(s) \Big] ds
\end{aligned}
$$

We end up with the solution

$$
\begin{aligned}
\mathbf{u}(\cdot, t) &= e^{t\nu\Delta} \mathbf{u}_0 - \int_{[0, t]} e^{(t-s)\nu\Delta} F(s) ds \\
                     &= e^{t\nu\Delta} \mathbf{u}_0 - \int_{[0, t]} e^{(t-s)\nu\Delta} \mathbb{L}[(\mathbf{u} \cdot \nabla) \mathbf{u}](\cdot, s) ds
\end{aligned}
$$

The main weirdness with this is that the operator $e^{t\nu\Delta}$ is not local for times $t > 0$!
It has a Gaussian kernel. (I don't want to derive it! See the heat equation solution, which is the same.) So,
$e^{t\nu\Delta}[f](\mathbf{x}) = \int_{R^3} G_t(\mathbf{x} - \mathbf{y}) f(\mathbf{y}) d\mathbf{y}$. The Gaussian kernel is
$G_t(\mathbf{x} - \mathbf{y}) := (4\pi \nu t)^{-3/2}e^{-\frac{||\mathbf{x}-\mathbf{y}||^2}{4\nu t}}$. This is global, but decays extremely quickly as $||\mathbf{x} - \mathbf{y}||^2$ increases.

So, for any $t > 0$, the operator $e^{t\nu\Delta}$ is global. This means any solution depends globally on the initial condition (even for very small $t$).
Since this operator is produced from the viscosity, it implies that the viscosity has a global effect (over non-infinitesmally small time).

## Conclusion

There are two parts of the Navier-Stokes equation that have global effects:
- (1) The pressure term is global, and instantaneous. It is a global adjustment applied to every point in time.
- (2) The viscosity is global, but not instantaneous. Over any time $t > 0$, the viscosity will have a global (but rapidly decaying with distance)
  influence on the fluid's motion.

Repeatedly inserting the solution into the integrals only works for the viscous term because it is linear. I think this trick cannot be used
with the convective term, since the integrals and non-linear part cannot commute.

Of course, all of this is known, and this is just my attempt to understand it. It is apparently more common to do this type
of derivation with fourier-stuff. (I think it makes the "Gaussian kernel" part more apparent to people who know a lot about Fourier-stuff.)
At least a lot of the stuff I found online goes to Fourier space.
