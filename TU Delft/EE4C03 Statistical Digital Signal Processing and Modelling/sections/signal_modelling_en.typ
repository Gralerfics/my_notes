#import "../generic.typ": *

#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Signal Modelling

Actually, modelling can be seen as a process of signal compression. A parameterized model uses a smaller number of parameters than the number of signal samples to represent complex signals, achieving more efficient storage and transmission.

The parameters obtained from compression can also be regarded as descriptions of the essential characteristics of the signal, such as the physical laws underlying it. This allows us to use the model to perform prediction or extrapolation on the unknown parts of the signal.

Now, our goal is to model a given digital signal $x[n]$, which means finding a model $H(z)$ such that its output signal $hat(x)[n]$ is as close as possible to the target signal $x[n]$.

== Autoregressive and Moving Average (ARMA) Model <sec:signal_modelling_arma_model>

We can use different types of models to represent signals according to practical situations. Here we take the Autoregressive and Moving Average (ARMA) model, specifically $"ARMA"(p, q)$, which is commonly used in time series analysis, as an example. Its transfer function is defined as follows:

$ Y(z) / X(z) = H(z) = (sum_(k=0)^q b[k] z^(-k)) / (1 + sum_(k=1)^p a[k] z^(-k)) = B(z) / A(z) $ <equ:signal_modelling_arma_tf>

Let the input and output signals be $x[n]$ and $y[n]$, with $X(z)$ and $Y(z)$ as their corresponding Z-transform functions. Although the index of $a[n]$ starts from 1 here formally, we can actually take $a[0] = 1$ to obtain a more unified form. Thus, by definition, we have $Y(z) A(z) = X(z) B(z)$, which transforms into the time domain as:

$ a[n] * y[n] = b[n] * x[n] $

Expanding this yields the classic form of the Linear Constant Coefficient Difference Equation (LCCDE):

$ y[n] + sum_(k=1)^p a[k] y[n-k] = sum_(k=0)^q b[k] x[n-k] $

Written in detail:

$ y[n] + a[1] y[n-1] + ... + a[p] y[n-p] \ = b[0] x[n] + b[1] x[n-1] + ... + b[q] x[n-q] $

Obviously, this system is a typical Linear Shift-Invariant (LSI) system with favorable properties.

=== Autoregressive (AR) Model

If there is only the autoregressive part, i.e., $"AR"(p)="ARMA"(p, 0)$:

$ H(z) = Y(z) / X(z) = b[0] / (1 + sum_(k=1)^p a[k] z^(-k)) $ <equ:ar_model_transfer_function>

In the time domain:

$ y[n] + sum_(k=1)^p a[k] y[n-k] = b[0] x[n] $

Written in detail:

$ a[0] y[n] + a[1] y[n-1] + ... + a[p] y[n-p] = b[0] x[n] $

This model considers the output at the current time $y[n]$ to be a linear combination of the previous $p$ outputs $y[n-1], ..., y[n-p]$ and the current input $x[n]$, hence it is called an autoregressive model. Since this model has no zeros, it is also known as an All-Pole Model.

=== Moving Average (MA) Model

If there is only the moving average part, i.e., $"MA"(q)="ARMA"(0, q)$:

$ H(z) = Y(z) / X(z) = sum_(k=0)^q b[k] z^(-k) $

In the time domain:

$ y[n] = sum_(k=0)^q b[k] x[n-k] $

Written in detail:

$ y[n] = b[0] x[n] + b[1] x[n-1] + ... + b[q] x[n-q] $

This model considers the output at the current time $y[n]$ to be a linear combination of the previous $q$ inputs $x[n], x[n-1], ..., x[n-q]$, hence it is called a moving average model. Since this model has no poles, it is also known as an All-Zero Model.

== Signal Models

As a discrete-time system, $H(z)$ cannot directly represent a signal; it needs an input to generate an output. We treat the output signal as the model's estimate of the target signal $hat(x)[n]$, and encapsulate the input signal $x[n]$ as part of the model.

=== Deterministic Modelling

#figure(
    caption: [Signal model with deterministic input]
)[
    #diagram(
        spacing: (10mm, 8mm),
        node-stroke: 0.8pt,
        edge-stroke: 0.8pt,
        node((0, 0), [$H(z)$], inset: 8pt),
        edge((-3, 0), "r,r,r", "-|>", [Impulse $delta[n]$], label-pos: 0),
        edge((0, 0), (3, 0), "-|>", [$hat(x)[n] = h[n]$], label-pos: 1),
    )
] <fig:signal_model_deterministic>

We can choose a known, deterministic signal and fix it as the system input to stably obtain the desired output signal $hat(x)[n]$, making its values as close as possible to the target signal $x[n]$. This is used for modeling deterministic signals and is called Deterministic Modelling.

This input signal can be chosen according to the practical situation; an input signal that matches the characteristics of the target signal can sometimes reduce the burden of model fitting. Here we can choose to use the simplest unit impulse signal $delta[n]$ as the input signal, which makes the system output $hat(x)[n]$ the unit impulse response $h[n]$ of the system.

=== Stochastic Modelling

#figure(
    caption: [Signal model with stochastic input]
)[
    #diagram(
        spacing: (10mm, 8mm),
        node-stroke: 0.8pt,
        edge-stroke: 0.8pt,
        node((0, 0), [$H(z)$], inset: 8pt),
        edge((-3, 0), (0, 0), "-|>", [White noise $v[n]$], label-pos: 0),
        edge((0, 0), (3, 0), "-|>", [$hat(x)[n]$], label-pos: 1),
    )
] <fig:signal_model_stochastic>

We can also choose to use random noise with a known distribution as the input, obtaining an output signal $hat(x)[n]$ whose statistical characteristics (such as mean and autocorrelation function) match those of the target signal $x[n]$. This is used for modeling stochastic processes and is called Stochastic Modelling.

We can choose to use white noise $v[n]$ with zero mean and variance $sigma_v^2$ as the input signal. The basis for this is that its autocorrelation function is $r_v [k] = sigma_v^2 delta[k]$, and its Fourier transform yields a constant power spectral density $P_v (omega) = sigma_v^2$, meaning it has a uniform energy distribution across all frequencies.

This characteristic ensures that we can obtain an output signal $hat(x)[n]$ with any frequency characteristic by filtering it. Meanwhile, the uniform energy across all frequency components ensures that the statistical characteristics of the output signal are independent of the input signal.
