---
layout: post
title: Traveling Waves for the Regularized Long-Wave Equation
---

In this notebook we'll use the Julia language to investigate special solutions of the Regularized Long-Wave (RLW) Equation. This equation is used to model small-amplitude, long-wavelength water waves - in particular, when wave amplitudes are much smaller than the fluid depth, which in turn is much smaller than the width of the body of water. This is a reasonable assumption in, for example, ocean waves. 

Let's take a look at the equation:

\begin{eqnarray}
\label{eq: rlw}
\large
u_t + u u_x - u_{xxt} = 0 
\end{eqnarray}

Here, $u = u(x, t)$ is the wave height at position $x$ and time $t$. Here $x$ is real-valued, denoting the dimension of interest. Since the width of the body of water is large, $x$ can be any real value: $- \infty < x < \infty$. 

This is a nonlinear Partial Differential Equation (PDE).The term $uu_x$ indicates that the wave accelerates itself: active transport, while the term $-u_{xxt}$ indicates dispersion.

## Traveling Wave solutions
We'll be looking for special solutions of the RLW equation. That is, we'll be looking for solutions of the form 

\begin{eqnarray}
\large
u(x, t) = \phi(x-ct)
\end{eqnarray}

where $c$ and $\phi$ are the wave-speed and wave-profile to be determined. This assumption will simply our work by turning the RLW PDE into an ODE. Letting $y = x-ct$, the chain rule gives us $\frac{\partial}{\partial x} = \frac{d}{dy}$ and $\frac{\partial}{\partial t} = -c \frac{d}{dy}$.

Further denoting $' = d/dy$, we have the ODE

\begin{eqnarray}
\large
-c \phi' + \phi \phi' + c \phi''' = 0.
\end{eqnarray}

Finally, integrating once with respect to $y$ gives us the second-order nonlinear ODE

\begin{eqnarray}
\large
-c \phi + \frac{1}{2} \phi^2 + c \phi'' = 0
\end{eqnarray}
which we will refer to as (TW).

We will be looking for solutions to (TW) satisfying zero-boundary conditions at infinity:

\begin{eqnarray}
\phi(\pm \infty) = 0.
\end{eqnarray}

Physically, this corresponds to a "swell" in the ocean: a localized pulse which would carry your boat up and then back down to sea-level after it passes. 


```julia
using Plots
x = -10:0.1:10;
plot(sech.(x))
```




    
![svg](output_7_0.svg)
    



Let's proceed by transforming (TW) into a first-order system of ODEs by setting

\begin{eqnarray}
w = \phi'.
\end{eqnarray}

This gives us the first-order system (S)
\begin{eqnarray}
\phi' & = & w \\
w' & = & \phi + \frac{1}{2c} \phi^2.
\end{eqnarray}

Let's proceed by solving system (S) numerically in Julia. We'll do so using the `DifferentialEquations` package to define our system.


```julia
using DifferentialEquations
```


```julia
function g!(du, u, p, t)
    du[1] = u[2]
    du[2] = u[1] + (1/(2p[1]))u[1]^2
end
```




    g! (generic function with 1 method)



Notice we defined our system using the above function `g!`, with system components $\phi$ and $w$ denoted by the Array `u`, and system parameter (the wave-speed) $c$ denoted by the parameter Array `p`, and position variable $y$ denoted by $t$.

Before we solve, we should translate our boundary conditions at infinity to some initial conditions, because `ODEProblem` tool we'll be using is expecting that. 

Let's think about what the center of the wave-profile might look like. The wave will arrive, reach its maximum, then vanish. At the maximum, it must be the case the case that $w = 0$ (by definition of $w$). Therefore we will try out initial conditions of the form $(\phi_0, 0)$ to discover our solution.

We'll also need to set "time-spans" and a parameter value for $c$, which we'll pick to be $(0, 20)$ and $-1$ respectively.


```julia
tspan = (0,20.0);
p = [-1.1];
```


```julia
#u0 = [3.59;0.0];
u0 = [3.3, 0.0]
  
  prob = ODEProblem(g!,u0,tspan,p);
  sol = solve(prob,Tsit5());
  plot(sol,label=["phi" "phi dot"], reltol=1e-8, abstol=1e-8)
```




    
![svg](output_15_0.svg)
    



We can see that we indeed get a localized, pulse-like solution with wave speed $c = -1.1$ and amplitude $\phi_0 = 3.3$. One thing to note is this solution seems to be periodic, with a somewhat long period, rather than satisying boundary conditions "at infinity".

## Parameter Estimation

One thing to try next would be to choose the amplitude parameter to maximize the period length. This could be attempted using Julia's Black Box Optimizer `BBOptim`.

# Comparison with Exact Solution

One nice thing about system (S) is that it is possible to find exact solutions for it. This is possible because it turns out to be a Hamiltonian system; that is, the system possesses a conserved quantity, unchanging across trajectories, that generates the equations of the system. 

A general Hamiltonian system of ODEs satisfies 

\begin{eqnarray}
\frac{dx}{dt} & = & -\frac{\partial H} {\partial y} \\
\frac{dy}{dt} & = & \frac{\partial H} {\partial x}
\end{eqnarray}

where $x, y \in \mathbb{R}^n$ and $H(x,y)$ is a smooth scalar function called the Hamiltonian. One can check that $\frac{dH}{dt}=0$ along trajectories of this system and hence $H$ is a conserved quantity.

We'll proceed by writing down the Hamiltonian for system (S):

\begin{eqnarray}
H(\phi, w) = -\frac{1}{2} w^2 + \frac{1}{2} \phi^2 + \frac{1}{6c} \phi^3.
\end{eqnarray}

Observe that $-H_w = w = \phi'$ and $H_{\phi} = \phi + \frac{1}{2c} \phi^2 = w'$, and hence system (S) is Hamiltonian.

Recall our boundary conditions at infinity; in these variables they imply that $\phi = 0$ and $w = 0$ at infinity (so that the pulse has horizontal asymptotes of zero). 

However, because the system is Hamiltonian, we can use the conserved quantity for $H(0,0)$ and solve for a solution curve of our system.

By plugging in $\phi = 0$ and $w = 0$ into $H$, we get 

\begin{eqnarray}
H(0,0) = 0.
\end{eqnarray}

Then we can solve for $w$ in terms of $\phi$, giving the solution curve

\begin{eqnarray}
w = \pm \phi \sqrt{1 + \frac{1}{3c} \phi}.
\end{eqnarray}

Let's go ahead and plot this solution curve in the $(\phi, w)$ phase plane.


```julia
c=-1.1;
ϕ = 0.0:0.01:3.28;
sqrt.(1 .+(1/(3c))*ϕ);
```


```julia
c=-1.1;
ϕ = 0.0:0.01:3.30;
plot(ϕ.*sqrt.(1 .+(1/3c)*ϕ), label="w(ϕ)")
```




    
![svg](output_23_0.svg)
    



We can see that this agrees closely with the analogous plot from our computed solution for $w$ vs $\phi$:


```julia
plot(sol[1, 21:29], sol[2, 21:29])
```




    
![svg](output_25_0.svg)
    



In fact, we can do even better: we can substitute the relation gained from using the Hamiltonian back into the first equation in (S), giving us a separable ODE for $\phi$:

\begin{eqnarray}
\int \frac{d\phi}{\phi \sqrt{1 + \frac{1}{3c} \phi}} = \int dy.
\end{eqnarray}

Carrying out the integration and solving for $\phi$ yields 
\begin{eqnarray}
\DeclareMathOperator{\sech}{sech}
\large
\phi(y) = -\frac{3}{c} \sech^2\left(\frac{y}{2}\right)
\end{eqnarray}


```julia
#circshift(sol[1, :], 15);
#x = -10:0.1:10;
#plot(sech.(x))
plot(circshift(sol[1, 1:30], 15))
```




    
![svg](output_28_0.svg)
    




```julia
x = -10:0.1:10;
plot((-3/-1.1)sech.((x/2).^2))
```




    
![svg](output_29_0.svg)
    


