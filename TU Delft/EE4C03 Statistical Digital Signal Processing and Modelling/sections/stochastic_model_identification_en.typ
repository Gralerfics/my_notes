#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Stochastic Modelling Identification

There are two main differences between modelling stochastic processes and modelling deterministic signals. First, in deterministic modelling, since the specific sample values of $x[n]$ are known, the error is defined based on those values; however, in stochastic modelling, we only possess the statistical characteristics of $x[n]$, making the previous definition of $e[n]$ unsuitable. Second is the difference in input signals: since we are modelling a stochastic process, the input is no longer a unit impulse signal but rather white noise with unit variance, as shown in @fig:signal_model_stochastic.

Given these differences, we must also assume stationarity for the stochastic process being modelled, specifically that the process is Wide-Sense Stationary (WSS).

Similarly, for stochastic processes, we can replace the least squares error in @sec:dmi_ls_method with the mean square error $cal(E)_"MS" = E{abs(x[n]-hat(x)[n])^2}$ for optimization. However, this leads to the same non-linear problems that are difficult to handle, requiring alternative solutions.

#blockquote[
    Note that from here on, $x[n]$, $v[n]$, etc., no longer represent specific discrete signals but rather stochastic processes. In practical applications, we may not directly know their statistical characteristics, in which case we must estimate them from specific signals (realizations or samples). If only one realization is available for estimation, the corresponding stationarity and ergodicity assumptions must hold for the results to be meaningful.
    
    For convenience, I will not change all subsequent square brackets to parentheses.
]

== Autoregressive Moving Average (ARMA) Processes

We first define a type of stochastic process called an ARMA process. Consider filtering white noise $v[n]$ with variance $sigma_v^2$ using the transfer function of an ARMA model (@equ:signal_modelling_arma_tf) to obtain the output $x[n]$.

#figure(
    caption: [Diagram of ARMA process generation]
)[
    #diagram(
        spacing: (10mm, 8mm),
        node-stroke: 0.8pt,
        edge-stroke: 0.8pt,
        node((0, 0), [$H(z) = B(z) / A(z) = (sum_(k=0)^q b[k] z^(-k)) / (1 + sum_(k=1)^p a[k] z^(-k))$], inset: 8pt),
        edge((-2, 0), (0, 0), "-|>", [White noise $v[n]$], label-pos: 0),
        edge((0, 0), (2, 0), "-|>", [$x[n]$], label-pos: 1),
    )
] <fig:smi_arma_proc_diagram>

Assuming $H(z)$ is stable, the output stochastic process $x[n]$ will be WSS (proof omitted). Since the power spectrum of white noise is $P_v (z) = sigma_v^2$, the power spectrum of $x[n]$ is:

$
P_x (z) = sigma_v^2 (B(z) B^*(1\/z^*))/(A(z) A^*(1\/z^*))
$

In the frequency domain, this is:

$
P_x (e^(j omega)) = sigma_v^2 abs(B(e^(j omega)))^2 / abs(A(e^(j omega)))^2
$

We define a process with a power spectrum of this form as an $"ARMA"(p, q)$ process. Note that due to symmetry, its power spectrum has $2p$ poles and $2q$ zeros.

#blockquote[
    To clarify once more, @sec:signal_modelling_arma_model mentioned ARMA models, whereas here we are discussing ARMA processes. The latter refers to the stochastic process that the output signal satisfies when white noise is passed through an ARMA model.
]

=== Yule-Walker Equations <sec:signal_modelling_yule_walker>

In stochastic modelling, we want the output of the constructed model to have the same statistical characteristics as the target process, such as ensuring the output's autocorrelation $r_x(k)$ matches that of the target process. Therefore, we need to *establish the statistical relationship between the model output's autocorrelation $r_x (k)$ and the system parameters $a[dot]$, $b[dot]$, and the unit impulse response $h[n]$*.

By definition, for an ARMA process $x[n]$ derived from $v[n]$, the following equation is satisfied:

$
x[n] + sum_(l=1)^p a[l] x[n-l] = sum_(l=0)^q b[l] v[n-l]
$

We can derive a similar relationship between the autocorrelation of $x[n]$ and the cross-correlation of $x[n]$ with $v[n]$ by multiplying both sides of the equation by $x^*[n-k]$ and taking the expectation:

$
E{x[n] x^*[n-k] + sum_(l=1)^p a[l] x[n-l] x^*[n-k]} = E{sum_(l=0)^q b[l] v[n-l] x^*[n-k]}
$

Which is:

$
E{x[n] x^*[n-k]} + sum_(l=1)^p a[l] E{x[n-l] x^*[n-k]} = sum_(l=0)^q b[l] #text(fill: blue, $E{v[n-l] x^*[n-k]}$)
$

Under the assumption of stationarity, substituting the definitions of autocorrelation and cross-correlation (see @sec:fun_rp_statistic) yields:

$
r_x (k) + sum_(l=1)^p a[l] r_x (k-l) = sum_(l=0)^(q) b[l] #text(fill: blue, $r_(v x) (k-l)$)
$

The presence of the cross-correlation term $r_(v x) (k-l)$ means the equation still contains $v$. We can replace it using the unit impulse response $h[n]$, which represents system properties, by substituting $x[n] = v[n] * h[n] = sum_(m=-infinity)^infinity v[m] h[n-m]$:

$
r_(v x) (k-l) &= E{v[k] x^*[l]} \
&= E{v[k] (sum_(m=-infinity)^infinity v[m] h[l-m])^*} \
&= sum_(m=-infinity)^infinity #text(fill: purple, $E{v[k] v^*[m]}$) h^*[l-m]
$

#blockquote[
    The use of $r_(v x) (k-l) = E{v[k] x^*[l]}$ here is consistent with the previously defined $r_(v x) (k-l) = E{v[n-l] x^*[n-k]}$ due to the stationarity assumption, because $(n-l)-(n-k)=k-l$. This substitution makes the derivation more concise.
]

Since $v[n]$ is white noise with independent and identical distribution and variance $sigma_v^2$:

$
#text(fill: purple, $E{v[k] v^*[m]}$) = cases(
    sigma_v^2\, quad &m = k,
    0\, &"otherwise"
)
$

This term is zero whenever $m != k$. Substituting this into the previous expression gives:

$
#text(fill: blue, $r_(v x) (k-l)$) = sigma_v^2 h^*[l-k]
$

Thus, we obtain an expression that does not contain $v$:

$
r_x (k) + sum_(l=1)^p a[l] r_x (k-l) = sigma_v^2 sum_(l=0)^(q) b[l] h^*[l-k]
$

Finally, considering practical constraints, we *assume the system is causal*, meaning $h[n] = 0$ for $n < 0$. Thus, $h^*[l-k]$ is zero when $l < k$. We can modify the upper and lower limits of the summation on the right side of the equation and denote it as $c[k]$:

#emphasis_equbox([
$
c[k] :&= sum_(l=0)^(q) b[l] h^*[l-k] = sum_(l=k)^(q) b[l] h^*[l-k] = sum_(l=0)^(q-k) b[l+k] h^*[l] \ &= b[k] * h^*[-k]
$
])

Incidentally, this term is zero when $k > q$. Finally, we obtain the *Yule-Walker Equations*:

#emphasis_equbox([
$
r_x (k) + sum_(l=1)^p a[l] r_x (k-l) = cases(
    sigma_v^2 c[k]\, quad &0<=k<=q,
    0\, &k>q
)
$

Note that from now on, we will default to the unit variance assumption, i.e., $sigma_v^2 = 1$.
])

This achieves our goal: establishing the statistical relationship between the model output's autocorrelation $r_x (k)$, the system parameters $a[dot]$ and $b[dot]$, and the unit impulse response $h[n]$.

#blockquote[
    As a side note, for $k > q$:

    $
    r_x (k) = - sum_(l=1)^p a[l] r_x (k-l)
    $

    The values of the autocorrelation function can be extrapolated using filter parameters and known autocorrelation values.
]

For clarity, we also present its matrix form:

$
mat(
    delim: "[",
    column-gap: #1.0em,
    row-gap: #0.5em,
    augment: #(hline: 4, stroke: (dash: (2pt, 2pt))),
    r_x (0), r_x (-1), dots, r_x (-p);
    r_x (1), r_x (0), dots, r_x (-p+1);
    dots.v, dots.v, dots.down, dots.v;
    r_x (q), r_x (q-1), dots, r_x (q-p);
    r_x (q+1), r_x (q), dots, r_x (q-p+1);
    dots.v, dots.v, dots.down, dots.v;
    r_x (q+p), r_x (q+p-1), dots, r_x (q);
    dots.v, dots.v, dots.down, dots.v;
)
mat(
    delim: "[",
    row-gap: #0.5em,
    1;
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
sigma_v^2
mat(
    delim: "[",
    row-gap: #0.5em,
    augment: #(hline: 4, stroke: (dash: (2pt, 2pt))),
    c[0];
    c[1];
    dots.v;
    c[q];
    0;
    dots.v;
    0;
    dots.v;
)
=
mat(
    delim: "[",
    row-gap: #0.5em,
    augment: #(hline: 4, stroke: (dash: (2pt, 2pt))),
    c[0];
    c[1];
    dots.v;
    c[q];
    0;
    dots.v;
    0;
    dots.v;
)
$

=== Modified Yule-Walker Equation (MYWE) Method

The Yule-Walker equations can be used to solve for filter parameters from the autocorrelation function, but due to the presence of $h^*[l]$, it remains a difficult non-linear problem.

To clarify again, in this problem we are modelling the stochastic process rather than a specific signal. Thus, the statistical characteristics of the target process (such as the autocorrelation function) are considered known. If the value of the autocorrelation function $r_v (k)$ is unknown, it must be estimated as $hat(r)_v (k)$ from some realizations (samples) using statistical methods.

Returning to the problem of parameter identification, we can follow the approach of the Padé method by solving in steps to *approximate* the optimal result. First, we use the portion where $q < k <= q+p$ to estimate $a[dot]$, with the corresponding equations:

$
mat(
    delim: "[",
    r_x (q), r_x (q-1), dots, r_x (q-p+1);
    r_x (q+1), r_x (q), dots, r_x (q-p+2);
    dots.v, dots.v, dots.down, dots.v;
    r_x (q+p-1), r_x (q+p-2), dots, r_x (q);
)
mat(
    delim: "[",
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
-
mat(
    delim: "[",
    r_x (q+1);
    r_x (q+2);
    dots.v;
    r_x (q+p);
)
$

This system of equations is called the *Modified Yule-Walker equations* (MYWE), and the method is thus called the MYWE method. Note that the form of this system is identical to that of @equ:deterministic_model_identification_pade_a_bar in the Padé method, except the values of $x[n]$ are replaced by the autocorrelation function. This matrix is also Toeplitz, allowing for accelerated solutions using algorithms like the Trench algorithm.

After obtaining $a[dot]$, the *second step is to solve for $b[dot]$*. Substituting $a[dot]$ back into the Yule-Walker equations gives the values of $c[dot]$. However, since $c[k] := b[k] * h^*[-k]$ and $h[k]$ even depends on $b[k]$, solving for $b[dot]$ is extremely difficult. The lecture slides state "We skip this," seemingly not intending to address this part. Related content in the reference book starts around page 190, mentioning several methods.

First, knowing $a[dot]$, we can construct an AR filter $A(z)$ to filter $x[n]$ and obtain a new process $y[n]$:

$
P_x (z) = (B(z) B^*(1\/z^*))/(A(z) A^*(1\/z^*)) quad attach(-->, t: A(z)) quad P_y (z) = B(z) B^*(1\/z^*)
$

This process is an MA process, which can then be handled using the methods in @sec:smi_ma_processes to estimate $b[dot]$.

Second, we can avoid explicit filtering, though the essence is likely the same. After solving for $c[dot]$ through the upper part of the Yule-Walker equations, the Laplace transform of the positive axis is obtained (since only the positive axis values can be solved through the Yule-Walker equations):

$
[C(z)]_+ = sum_(k=0)^infinity c[k] z^(-k)
$

Correspondingly, though unknown, the Laplace transform of the negative axis is:

$
[C(z)]_- = sum_(k=-infinity)^(-1) c[k] z^(-k) = sum_(k=1)^infinity c[-k] z^k
$

From the definition $c[k] := b[k] * h^*[-k]$, the power spectrum of the MA process is:

$
C(z) = B(z) H^*(1\/z^*) = B(z) (B^*(1\/z^*))/(A^*(1\/z^*)) \
quad => quad P_y (z) equiv C(z) A^*(1\/z^*) = B(z) B^*(1\/z^*)
$

Expanding this:

$
P_y (z) = C(z) A^*(1\/z^*) = [C(z)]_+ A^*(1\/z^*) + [C(z)]_- A^*(1\/z^*)
$

Since the negative axis values of $a[k]$ are 0, $A^*(1\/z^*)$ contains only positive powers of $z$, as does $[C(z)]_+$ #text(fill: red, "(TODO, did the book use a minus sign here?)"). Thus, the causal part of $P_y (z)$ is:

$
[P_y (z)]_+ = [C(z)]_+ A^*(1\/z^*)]_+
$

Therefore, despite not knowing the values of $c[k]$ on the negative axis, we can use this equation and the known positive axis values of $c[dot]$ along with $a[dot]$ to solve for $[P_y (z)]_+$, then obtain the full $P_y (z)$ through conjugate symmetry. Finally, spectral factorization is performed to obtain the coefficients $b[dot]$:

$
P_y (z) = B(z) B^*(1\/z^*)
$

#text(fill: red, "(TODO, should I copy it over?)") There is a clear example on page 192 of the reference book.

=== Extended Yule-Walker Equation Method

Correspondingly, in the first step, we can use an approach similar to Prony's method by including all equations where $k > q$. The resulting overdetermined system is called the *Extended Yule-Walker equations*.

A least squares solution is then formally sought, using a process similar to Prony's method, the details of which will not be repeated.

== Autoregressive (AR) Processes <sec:smi_eyem_ar_yule_walker_method>

We again consider the all-pole case. Since only $b[0]$ remains, the equations simplify significantly:

#emphasis_equbox([
$
r_x (k) + sum_(l=1)^p a[l] r_x (k-l) = abs(b[0])^2 delta(k), quad k>=0
$
])

Since the complex issue of $c[k]$ is absent, Prony's method can be directly applied. First, solve for $a[dot]$ using all equations except the first one; note that the system is identical to @equ:dmi_prony_all_pole_matrix. Then, use the first equation to derive $b[0]$. This is called the Yule-Walker method (a confusing naming convention).

#text(fill: red, "(TODO, what is the 1/N at the front? Although it seems to cancel out when solving for a[·], and b[·] is unaffected due to direct selection. Also, what about the asterisk? See page 194, eq 4.153 of the book.)") Once again, if the autocorrelation of the target process is unknown, it must be estimated from samples (satisfying the ergodicity assumption), for example:

$
hat(r)_x (k) = 1/N sum_(n=0)^(N-1) x[n] x^*[n-k]
$

This makes the approach equivalent to the previously mentioned Auto-correlation method (@sec:dmi_finite_data_autocorrelation_method), which aligns with intuition.

== Moving Average Processes <sec:smi_ma_processes>

For MA processes, substituting into the Yule-Walker equations yields:

$
r_x (k) = sum_(l=0)^q b[l] b^*[l-k] = b[k] * b^*[-k]
$

Thus:

$
P_x (z) = B(z) B^*(1\/z^*)
$

In summary, the results are obtained by taking the z-transform of the autocorrelation function $r_x (k)$ to get the power spectrum $P_x (z)$, followed by spectral factorization. For example:

$
r_x (k) = 17 delta(k) + 4[delta(k-1) + delta(k+1)]
$

Taking the z-transform:

$
P_x (z) = 17 + 4z^(-1) + 4z = (4 + z^(-1))(4+z)
$

Thus:

$
B(z) = 4 + z^(-1) quad "or" quad B(z) = 1 + 4z^(-1)
$

Additionally, there are methods such as Durbin's method, which are not recorded here.

// == Application: Spectrum Estimation

// Using estimated r_x (k) directly for z-transform vs calculating power spectrum after estimating model parameters

// Issues: 1. Model type may be unknown; 2. Model order may be unknown
