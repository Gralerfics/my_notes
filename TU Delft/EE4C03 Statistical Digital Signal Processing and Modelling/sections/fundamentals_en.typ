#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge
#import "@preview/suiji:0.4.0": *

= Linear Algebra

#Cre("(TODO)") Solutions, rank, inverses, eigenvalues, quadratic forms and positive definiteness, matrix calculus, etc., of systems of linear equations.

= Probability Theory

== Random Variable

Certain phenomena (such as rolling dice) are difficult to predict due to their random nature; we refer to these as *random phenomena*. Random phenomena occur continuously. Each time we observe and record their outcome, it is called a *random test*, and the recorded result is called a *sample*. A single random test may yield various different outcomes (samples). The collection of all such samples is called the sample space of that random phenomenon, encompassing all possible outcomes that may occur in that random event. A *random event* is a subset of the sample space that contains a portion of the possible outcomes (samples). If the outcome (sample) $omega$ obtained from a single random trial belongs to the sample subspace included in a random event $A$, then the random event $A$ is said to have occurred; otherwise, it is said not to have occurred.

A random event can be any statement, such as "It will rain tomorrow," which corresponds to the sample space "Whether it rains tomorrow or not." For this random phenomenon, the sample space for whether it rains tomorrow is ${"rain", "no rain"}$, and "it rains tomorrow" is its subset ${"rain"}$. Similarly, "the result of rolling a fair 20-sided die is greater than 18" is a random event, and the sample space for rolling a fair 20-sided die is ${1, 2, ..., 20}$. This random event can be denoted as ${19, 20}$ as a random event of this random phenomenon. The sample space of the random phenomenon is ${1, 2, ..., 20}$, and this random event can be denoted as ${19, 20}$. A *random variable* $x$ is used to describe the outcome of a random phenomenon whose results can be quantified. For example, the outcome of rolling a twenty-sided die yields one of twenty integers from $1$ to $20$. Even sample spaces that do not explicitly involve quantities, such as "whether it will rain tomorrow," can be quantified by arbitrarily defining "rain" as 1 and "no rain" as 0.

#blockquote[
    In most lecture notes, random variables are denoted by uppercase letters, but this is merely a difference in notation. To maintain consistency with the notation used in subsequent discussions of stochastic processes, lowercase letters will be used here to represent random variables, with explicit notation indicating their randomness when necessary. Correspondingly, parameters and other independent variables (such as the parameters in the probability density function discussed later) will be denoted by Greek letters like $alpha$ to avoid confusion.
]

The probability of a random event $A$ occurring is denoted as $Pr(A)$, representing the probability of the event happening. The probability of a random event encompassing the entire sample space is $1$, since the outcome of the random experiment is necessarily included within it, signifying that the event must occur. Further details regarding probability, conditional probability, total probability, and the sample space set will not be elaborated upon here.

== Probability Density Function (PDF)

The different values of a random variable have different probabilities. For a discrete random variable (one with a finite number of possible values), a mapping that assigns each possible value to its corresponding probability is called the probability mass function (PMF). This function contains all the information about the distribution.

If we consider continuous random variables, where the number of possible values is infinite, the probability of obtaining any specific value is approximately $0$. Thus, we can no longer construct a function by assigning probability values to each specific outcome as we do with a PMF. We first define the Cumulative Distribution Function (CDF) for the random variable $x$:

$ F_x (alpha) = "Pr"{x <= alpha} $

Whether curly braces or parentheses are used is just notation and does not change the meaning. In short, it represents the probability that $x$ is less than or equal to the parameter $alpha$. Clearly, for a continuous random variable $x$ taking real values, the CDF tends toward $0$ as $alpha$ approaches negative infinity, and toward $1$ as $alpha$ approaches positive infinity. We define the Probability Density Function (PDF) via the CDF:

$ f_x (alpha) = dif / (dif alpha) F_x (alpha) $

That is, the PDF is the derivative of the CDF. According to the Fundamental Theorem of Calculus, the integral of the PDF from $alpha_1$ to $alpha_2$ (i.e., the area under the curve over this interval) equals the probability that the random variable $x$ falls within this interval:

$
integral_(alpha_1)^(alpha_2) f_x (alpha) dif alpha &= F_x (alpha_2) - F_x (alpha_1) \ &= "Pr"{x <= alpha_2} - "Pr"{x <= alpha_1} \ &= "Pr"{alpha_1 < x <= alpha_2}
$

Similarly, the PDF contains all the information regarding the distribution of a random variable.

== Joint Distribution

A joint distribution describes the synergistic distributional relationship between multiple random variables. For example, for two random variables $x$ and $y$, the joint cumulative distribution function is defined as:

$ F_(x, y) (alpha, beta) = "Pr"{x <= alpha, y <= beta} $

This is straightforward: it simply places constraints on both random variables simultaneously. Its joint probability density function is:

$ f_(x, y) (alpha, beta) = partial^2 / (partial alpha partial beta) F_(x, y) (alpha, beta) $

== Mathematical Expectation

Mathematical expectation can be thought of as the average value obtained after performing an infinite number of random trials. By definition, it is the weighted average of all possible outcomes using their respective probabilities. For example, the expectation of a continuous random variable $x$ is:

$ E(x) = integral_(-infinity)^infinity alpha f_x (alpha) dif alpha $

If it is a discrete random variable, the definition is even clearer ($alpha$ represents any possible value):

$ E(x) = sum_(alpha) alpha f_x (alpha) $

This is simple, but let's look at an example that is not just the expectation of a single random variable. For instance, the variance of a random variable $x$ is defined as $E{(x - E(x))^2}$, and its integral form is:

$ E{(x - E(x))^2} = integral_(-infinity)^infinity (alpha - E(x))^2 f_x (alpha) dif alpha $

Note that the distribution function used here is still $f_x (alpha)$ rather than something like $f_x ((alpha - E(x))^2)$ or $f_((x - E(x))^2) (alpha)$. What I want to emphasize is that we should add a subscript to $E$ to indicate which random variable's distribution the expectation is being calculated over:

$ E_x (x) = integral_(-infinity)^infinity alpha f_x (alpha) dif alpha $

This is because calculating an expectation requires two elements: first, a distribution (the averaging process occurs over this distribution); and second, an expression (the content inside the parentheses). In the defining formula for expectation, the $x$ beneath the distribution function $f_x (alpha)$ comes from the $x$ in $E_x$, representing the first element (the distribution), rather than the $x$ in the expression within $E(x)$. The $alpha$ multiplied before the distribution function is the result of substituting the parameter $alpha$ into the $x$ of the expression, as shown in example @equ:joint_expectation_example.

This is actually a simple matter, but it is easy to get stuck here if one isn't careful, as the description of the distribution in the expectation symbol is often omitted, defaulting to the distribution of the random variables present in the expression. This is mentioned here merely for clarification.

In cases involving multiple random variables, it is even more important to clarify which distribution the expectation is being taken over. For example, for the expectation of an expression involving two random variables, the integral must be changed to a double integral to traverse all combinations of values, and the distribution function must be changed to the joint probability density function:

$ E_(x, y) (x y^*) = integral_(alpha=-infinity)^infinity integral_(beta=-infinity)^infinity alpha beta^* f_(x, y) (alpha, beta) dif alpha dif beta $ <equ:joint_expectation_example>

== Moments, Mean and Variance

Next are some statistical quantities based purely on definitions.

First, we define the *$k$-th order raw moment* of a random variable $x$ as $E(x^k)$, where the first-order raw moment is the Mean of $x$, denoted as $m_x = E(x)$.

Next, we define the *$k$-th order central moment* as $E{(x - E(x))^k}$. It can be seen that the second-order central moment is the variance of $x$, denoted as $sigma^2_x = "Var"(x)$, where $sigma_x$ is called the standard deviation (the square root of the variance).

Other statistical quantities and their corresponding physical meanings will not be detailed here.

== Independence

Independence describes whether two random variables are mutually unaffected. If two random variables $x$ and $y$ are independent, it means the probability of $x$ taking a value $alpha$ multiplied by the probability of $y$ taking a value $beta$ directly yields the probability of them taking those two values simultaneously. This is because the probability of independent events occurring together equals the product of their individual probabilities (refer to the multiplication principle).

Satisfying this condition for all possible values implies that the product of their individual distributions yields the joint distribution, which is the necessary and sufficient condition for independence:

$ f_(x, y) (alpha, beta) = f_x (alpha) f_y (beta) $

== Covariance Function

Covariance measures the degree of correlation between random variables. The covariance of two random variables $x$ and $y$ is defined as:

$ c_(x y) = "Cov"(x, y) = E{(x - m_x)(y - m_y)^*} $

Incidentally, note that $c_(x x) = E{(x - m_x)^2}$ is the variance of $x$.

Covariance can be thought of as a measure of the significance of the linear relationship between random variables. Specifically, based on the definition: when $x$ is above its mean, to what extent does $y$ tend to be simultaneously above (or below) its own mean? The greater this extent, the more it indicates that when $x$ is larger, $y$ also tends to be larger (or smaller), suggesting a stronger linear relationship. By subtracting the means, the influence of the absolute values is removed, focusing primarily on the relative trends of the variables.

== Correlation Function

The Correlation Function is defined as the expectation of the inner product of two random variables:

$ r_(x y) = E{x y^*} $

It can be viewed as the covariance function without the means removed. Their relationship can be derived as:

$
c_(x y) &= E{(x - m_x)(y - m_y)^*} \
&= E{x y^*} + E{m_x m_y^*} - m_x E{y^*} - E{x} m_y^* \
&= E{x y^*} + m_x m_y^* - m_x m_y^* - m_x m_y^* \
&= r_(x y) - m_x m_y^*
$

They differ only by a constant $m_x m_y^*$, so the physical meaning they express can be considered similar.

== Correlation Coefficient and Orthogonality

Of course, covariance does not remove the units (dimension), so it lacks universality across different combinations of random variables. We typically use normalization to define the correlation coefficient to measure Correlation, such as the widely used Pearson correlation coefficient:

$ rho_(x y) = c_(x y) / (sigma_x sigma_y) $

It can be proven that its value ranges within $[-1, 1]$, providing a unified measure of the significance of the linear relationship between random variables. If this correlation coefficient is zero, the two random distributions are said to be Orthogonal.

When the correlation coefficient is $0$ (i.e., $c_(x y) = 0$), the two random variables are uncorrelated. A positive correlation coefficient represents positive correlation, and a negative one represents negative correlation. The necessary and sufficient condition for being uncorrelated is:

$
c_(x y) = 0  quad => quad r_(x y) - m_x m_y^* = 0 quad => quad E{x y^*} = E{x} E^*{y}
$

That is, the expectation of the product of the random variables equals the product of their individual expectations.

#blockquote[
   Note that correlation and independence are not the same thing.

   First, being uncorrelated does not necessarily mean independence. For example, if two random variables satisfy $y = x^2$, they are completely dependent, yet the correlation coefficient is still $0$. As mentioned, correlation mainly measures the significance of a linear relationship. In this case, when $x$ is negative, $y$ decreases as $x$ increases; when $x$ is positive, $y$ increases as $x$ increases. These two parts of the relationship cancel each other out symmetrically, resulting in no linear component, even though a quadratic relationship exists.

   However, independence does imply that variables are uncorrelated. Since the distributions are entirely unrelated, the joint density is simply the product of the individual densities, and by definition, $E{x y^*} = E{x} E^*{y}$ follows.

   In summary: $"Independence" => "Uncorrelated"$, but the reverse is not necessarily true.
]

= Random Process

A random process ${x(n)}$ is essentially a sequence composed of random variables.

These random variables can be independent and identically distributed (i.i.d.), such as white noise, but they often are not; in the real world, most are not. The index of this sequence can represent time, or it can represent something else. For convenience, we will assume we are studying time series, where the index $n$ represents different points in time.

As before, the focus of studying randomness is on finding the "invariants" within the change, such as statistical characteristics. Therefore, we begin by defining a series of statistical features.

== Statistic <sec:fun_rp_statistic>

=== Mean

A random process does not necessarily have a single uniform mean; the mean is still defined for each individual random variable, resulting in a sequence:

$ m_x (n) = E{x(n)} $

=== Auto-correlation <sec:fun_rp_autocorrelation>

For a random process $x(n)$, the auto-correlation function is the correlation function of two random variables at specified indices $k$ and $l$:

$ r_x (k, l) = r_(x(k), x(l)) = E{x(k) x^*(l)} $ <equ:fun_rp_autocor_def>

From this formula, each $r_x$ is only related to two specific random variables and has no immediate connection to the process as a whole. However, this changes when we introduce Wide-Sense Stationarity (WSS).

=== Auto-covariance

Analogous to covariance, we define auto-covariance:

$ c_x (k, l) = E{[x(k) - m_x (k)] [x(l) - m_x (l)]^*} = r_x (k, l) - m_x (k) m_x^* (l) $

=== Cross-correlation

Cross-correlation involves two random processes. For example, the cross-correlation of ${x(n)}$ and ${y(n)}$ is defined as:

$ r_(x y) (k, l) = r_(x(k), y(l)) = E{x(k) y^*(l)} $

=== Cross-covariance

Similarly, cross-covariance is defined as:

$ c_(x y) (k, l) = E{[x(k) - m_x (k)] [y(l) - m_y (l)]^*} = r_(x y) (k, l) - m_x (k) m_y^* (l) $

These definitions are listed here for completeness, though they will not be used extensively later.

== Structural Invariance

Notice that the statistics defined previously focus on discrete random variables and do not strongly relate to the overall random process. Next, we will focus on some structural invariants of the entire random process. These properties provide a rational basis for estimating distributions when the true distribution is unknown and only finite samples are available.

Specifically, while random variables are idealizations described by a PDF, we usually study actual signals without knowing their true distribution. We must rely on samples to estimate the distribution. We can view a specific signal $x[n]$ as a single Realization of a random process. Each sample at a given point in time is a single observation of the corresponding random variable in the process. @fig:fun_rp_realizations_of_rp shows a series of possible realizations obtained after many trials of a random process, with one highlighted.

#figure(
    caption: [Different realizations of a random process]
)[
    #cetz.canvas({
        import cetz.draw: *

        let samples = 100
        let len = 30

        let (x, y) = (0, 0)
        let (x-new, y-new) = (0, 0)
        let v = ()
        let col = black
        
        let rng = gen-rng-f(42)

        for idx in range(samples) {
            (rng, v) = uniform-f(rng, low: -1.0, high: 1.0, size: 1)
            (x, y) = (0, v.at(0))
            
            for k in range(len) {
                (rng, v) = uniform-f(rng, low: -0.4, high: 0.15, size: 1)
                (x-new, y-new) = (x + 0.4, y + v.at(0))


                if (idx == 0) {
                    col = color.mix(
                        (blue.transparentize(0%), 1 - k / len),
                        (green.transparentize(0%), k / len)
                    )
                } else {
                    col = color.mix(
                        (blue.transparentize(90%), 1 - k / len),
                        (green.transparentize(90%), k / len)
                    )
                }
                
                line(stroke: (paint: col, cap: "round", thickness: 1pt),
                  (x + 0.025, y), (x-new, y-new)
                )

                (x, y) = (x-new, y-new)
            }
        }
    })
] <fig:fun_rp_realizations_of_rp>

=== Stationarity

Our *first question* is: In a random process, do certain statistical features change over time?

We introduce the concept of Stationarity. The main idea is to measure whether the statistical characteristics of any sub-sequence remain unchanged after being delayed by an arbitrary amount of time. To have identical statistical features, the simplest way is to require the distributions to be identical, which defines Strict-Sense Stationarity (SSS):

$ forall k, {n_1, n_2, dots, n_m}, f_(x(n_1), x(n_2), dots, x(n_m)) (dot) = f_(x(n_1 + k), x(n_2 + k), dots, x(n_m + k)) (dot) $

The requirements for SSS are too high and usually unnecessary. Therefore, we define Wide-Sense Stationarity (WSS), focusing only on the consistency of first and second-order statistics:

#align(center)[
    + $ m_x (n) = m_x $
    + $ r_x (k, l) = r_x (k - l) $
    + $ c_x (0) < infinity $
]

That is, for a WSS process: the mean is independent of time, and the auto-correlation function depends only on the time difference, not on absolute time.

#blockquote[
    The third condition—finite variance—is often omitted in many texts.

    Finite variance is equivalent to a finite mean-square value (differing only by the square of the mean). Physically, this represents finite power. Engineers often ignore this condition because most physical processes have finite power.
    
    Mathematically, this ensures the existence of the second moment. However, since the second condition requires $r_x (0)$ to exist, and the existence of an expectation implies convergence to a finite value, the third condition is effectively implied.
]

This allows us to take sub-sequences from a single realization at different times and claim their average properties are consistent. For a stationary signal, we can split a sufficiently long sample into multiple segments and assert that these segments share the same underlying statistical features.

=== Ergodicity

The stationarity derived from the first question allows us to use segments of a long sample as different samples. Our *second question* is: Can a single, sufficiently long realization represent the entire ensemble of infinite possible realizations (the sample space)? If I split this long sample to estimate the distribution, is the result the true distribution of the random process? Can one realization represent all aspects of the process?

We use Ergodicity to describe this property. It is built upon stationarity; we must ensure statistical features don't change over time, otherwise, different segments would have different properties, making a long sample useless.

First, we need a typical example of "stationary but non-ergodic" to show why this concept is necessary. Let the random process $x(n)$ be such that $x(0)$ is a random variable $z$ drawn once, and thereafter $x(n) = x(n-1)$. This is called a "random constant" because every realization is a constant signal, but the constant value differs between realizations depending on the initial draw.

Analyzing this process: first, it is definitely stationary because a shift in time clearly does not affect the distribution. We can also prove it satisfies WSS: the mean $m_x = E{x(0)}$ is independent of $n$; the auto-correlation $r_x (k, l) = E{x(k) x^*(l)} = E{x(0) x^*(0)} = sigma_x^2$ is independent of time.

Now consider ergodicity. The ensemble mean is $E{z}$, but for a long realization, the time average is whatever value $x(0)$ happened to be for that specific trial. No matter how long the sample is, the average remains that value, which may not equal $E{z}$. In other words, we cannot estimate the properties of $z$ using only one realization, even if it is infinitely long.

*In summary*, stationarity ensures that a sufficiently long sample yields a meaningful time average, while ergodicity ensures that this time average equals the ensemble average (mean of all possible samples). Thus, estimates derived this way are correct.

#blockquote[
    In physics and statistical mechanics, this concept is often translated as "各态历经性" (all-state traversing), which is more descriptive. It suggests that as time progresses, the system will eventually experience and exhibit all its possible macroscopic states. Only when a single realization has the potential to manifest all features of the distribution can we estimate the macroscopic distribution from a single sample.

    We can also look at this from the perspective of biased and unbiased estimation. For a stationary signal, we can estimate a meaningful statistic, but we cannot guarantee if it is unbiased; ergodicity guarantees that it is.
]

The strict definition of ergodicity is abstract and complex. Like WSS, we only consider ergodicity in terms of the mean and correlation for our needs. Of course, the following discussion assumes at least WSS.

First is *Mean-Ergodicity*. We define the sample mean, which is the temporal average of a single realization:

$ hat(m)_x^((N)) = 1/N sum_(n=0)^(N-1) x[n] $

Considering the original definition, the time average (sample mean) equals the ensemble average:

$ lim_(N -> infinity) E{abs(hat(m)_x^((N)) - m_x)^2} = 0 $

Since we cannot truly verify this from the real distribution, we provide a *necessary and sufficient* condition for mean ergodicity (proof omitted):

$ lim_(N -> infinity) 1/N^2 sum_(n=0)^(N-1) sum_(m=0)^(N-1) c_x (n - m) = 0 $

That is, the auto-covariance function decays fast enough. An intuitive interpretation is that the correlation between random variables does not persist too long in time; a long-term average can "dilute" the dependencies.

There is also a *sufficient* condition for quick judgment in simple cases:

$ c_x (0) < infinity quad "and" quad lim_(k -> infinity) c_x (k) = 0 $

Or an even stronger condition: the absolute summability of the auto-covariance function is a *sufficient* condition for mean ergodicity:

$ sum_(k=-infinity)^infinity abs(c_x (k)) < infinity $ <equ:fun_rp_ergo_autocov_abssum>

Next is *Correlation-Ergodicity*. We similarly define the sample auto-correlation function. Since the auto-correlation of a WSS process only depends on the time difference, all pairs in the sample with a time difference $k$ can be used to estimate the $k$-th term of the auto-correlation function. we sum them and take the average:

$ hat(r)_x^((N)) (k) = 1/N sum_(n=0)^(N-1-k) x[n+k] x^*[n] $ <equ:fun_rp_ergo_autocor_sample>

#blockquote[
    You might notice that in this definition, a total of $N-k$ pairs were used, yet we divide by $N$ rather than $N-k$.
    
    This actually makes the estimate biased. However, if we only consider this for the definition of correlation ergodicity where $N$ tends to infinity, the effect of $k$ becomes negligible. That is, for this definition, we only require asymptotic unbiasedness.

    In @sec:se_nonparam_periodogram, we will revisit the use of this sample auto-correlation function for actual estimation and the explanation for dividing by $N$ instead of $N-k$.
]

Similarly, the definition is that the time average equals the ensemble average:

$ lim_(N -> infinity) E{abs(hat(r)_x^((N)) (k) - r_x (k))^2} = 0 $

We provide the practical *necessary and sufficient* condition without proof:

$ forall k, lim_(N -> infinity) 1/N^2 sum_(n=0)^(N-1) sum_(m=0)^(N-1) abs(r_x (n - m) - r_x (n - m + k))^2 = 0 $

Additionally, the absolute summability condition @equ:fun_rp_ergo_autocov_abssum is also a *sufficient* condition for correlation ergodicity.

== Power Spectrum

#text(fill: red, "(TODO)") Add content regarding energy signals, power signals, and Power Spectral Density (PSD).

According to the Wiener–Khinchin Theorem, the PSD of a WSS process is the Fourier transform of its auto-correlation function.

== Filtering Random Processes

#text(fill: red, "(TODO)")

If the filter coefficients $h[n]$ are finite-length and zero outside $[0, N-1]$, the power of the output process can be expressed using the auto-correlation matrix of the input $x(n)$ and the filter coefficient vector:

$
sigma_y^2 = E{abs(y(n))^2} = bold(h)^H bold(R)_x bold(h)
$ <equ:fun_rp_filtering_power_in_y_repr_in_h_and_Rx>

== Random Process and Digital Signals

Our random process $x(n)$ is a distribution, while our digital signal $x[n]$ can be seen as a realization sampled from it. We have defined characteristics like mean and auto-correlation for random processes, and we can perform operations like mean and auto-correlation on digital signals. We now want to examine the differences and connections between the two.

Clarifying this allows us to understand how to use samples to estimate distribution features, which has been our consistent goal. We will explore this from several common perspectives.

=== Auto-correlation Estimation with Multiple Realizations <sec:fun_rp_auto_cor_estimation_with_multiple_realizations>

If we have a large number of sampling results (realizations) of a random process, we can directly estimate its distribution features, such as the auto-correlation function, using the definition (@equ:fun_rp_autocor_def).

Suppose we have $N$ realizations of the process ${x(n)}$, each of length $M$. The result of the $i$-th realization is denoted as $x_i [n]$, or as a vector (where $M$ can be finite or infinite; here it is finite for demonstration):

$ bold(x)_i = vec(x_i [0], x_i [1], dots.v, x_i [M-1]) $

We arrange the data from multiple samplings into a matrix:

$
X = [bold(x)_0 quad bold(x)_1 quad dots quad bold(x)_(N-1)]
=
mat(
    delim: "[",
    x_0 [0], x_1 [0], dots, x_(N-1) [0];
    x_0 [1], x_1 [1], dots, x_(N-1) [1];
    dots.v, dots.v, dots.down, dots.v;
    x_0 [M-1], x_1 [M-1], dots, x_(N-1) [M-1];
)
$

By definition, $r_x (k, l) = E{x(k) x^*(l)}$. Since $E$ is the expectation over the ensemble, we simply use $x_i [k]$ and $x_i [l]$ from different realizations to estimate it:

$ hat(r)_x (k, l) = 1/N sum_(i=0)^(N-1) x_i [k] x_i^* [l] $

We write the estimated auto-correlation values as a matrix, which can be computed from $X$, leading to a concise matrix expression:

$
hat(bold(R))_x &:=
mat(
    delim: "[",
    column-gap: #1.0em,
    row-gap: #0.5em,
    hat(r)_x (0, 0), hat(r)_x (0, 1), dots, hat(r)_x (0, M-1);
    hat(r)_x (1, 0), hat(r)_x (1, 1), dots, hat(r)_x (1, M-1);
    dots.v, dots.v, dots.down, dots.v;
    hat(r)_x (M-1, 0), hat(r)_x (M-1, 1), dots, hat(r)_x (M-1, M-1);
) \
&=
display(1/N sum_(i=0)^(N-1)) mat(
    delim: "[",
    column-gap: #1.0em,
    row-gap: #0.5em,
    x_i [0] x^*_i [0], x_i [0] x^*_i [1], dots, x_i [0] x^*_i [M-1];
    x_i [1] x^*_i [0], x_i [1] x^*_i [1], dots, x_i [1] x^*_i [M-1];
    dots.v, dots.v, dots.down, dots.v;
    x_i [M-1] x^*_i [0], x_i [M-1] x^*_i [1], dots, x_i [M-1] x^*_i [M-1];
) \
&=
inline(
    mat(
        delim: "[",
        x_0 [0], x_1 [0], dots, x_(N-1) [0];
        x_0 [1], x_1 [1], dots, x_(N-1) [1];
        dots.v, dots.v, dots.down, dots.v;
        x_0 [M-1], x_1 [M-1], dots, x_(N-1) [M-1];
    )
    mat(
        delim: "[",
        x_0 [0], x_0 [1], dots, x_0 [M-1];
        x_1 [0], x_1 [1], dots, x_1 [M-1];
        dots.v, dots.v, dots.down, dots.v;
        x_(N-1) [0], x_(N-1) [1], dots, x_(N-1) [M-1];
    )^*
) \
&=
1/N bold(X) bold(X)^H
$ <equ:fun_rp_autocor_matrix>

#text(fill: red, "(TODO)") In cases where there are not many independent realizations, define $bold(x)_i$ as $[x[i], dots, x[i + L - 1]]^T$. Refer to Slides Lec11 P6.

=== Auto-correlation Estimation with Correlation-Ergodicity

As seen in the previous section, the larger the number of realizations $N$, the more accurate the estimate. However, if we have only one sample ($N = 1$), the estimate using that method will be extremely imprecise. In this case, if the random process is ergodic, we are allowed to use information from different times within that single sample to perform the estimation.

We define the auto-correlation function of a finite-length digital signal $x[n]$ (length $N$, note that $N$ here is not the number of realizations above) as:

$ R_(x x) (k) = sum_(n=0)^(N-1-k) x[n + k] x^*[n] $

From a formulaic perspective, this can be understood as a measure of the similarity of a signal to itself at different time delays.

In practice, any definition expressing this meaning can be called auto-correlation; there are many variations with minor differences, such as whether to divide by the number of samples (taking an average) or whether to use delay $k$ or lead $k$ (which just flips the result). Use whichever is convenient.

We choose this definition here because the MATLAB cross-correlation function `xcorr` is defined this way (from the documentation):

$
hat(R)_(x y) (m) = cases(
    sum_(n=0)^(N-m-1) x_(n+m) y_n^*", " &m >= 0",",
    hat(R)_(y x)^* (-m)", " &m < 0".",
)
$

Recall the sample auto-correlation function $hat(r)_x^((N)) (k)$ defined in @equ:fun_rp_ergo_autocor_sample; its form is almost identical to the signal auto-correlation function $R_(x x) (k)$ we defined, differing only by a coefficient.

According to previous definitions, if the process is correlation-ergodic, then this sample auto-correlation function $hat(r)_x^((N)) (k)$ can be used to correctly estimate the auto-correlation function $r_x (k)$ of the random distribution. Consequently, the signal cross-correlation function $hat(R)_(x y) (m)$, which has a nearly identical form, carries the same physical meaning and can also be used to estimate $r_x (k)$ (with appropriate coefficient adjustments).

*In summary*, the existence of correlation ergodicity allows us to use the signal auto-correlation of a sample to estimate the auto-correlation characteristics of the random process when the number of realizations is insufficient.

=== Auto-correlation Estimation (Comprehensive)

If multiple realizations (samples) are available and ergodicity holds, we can combine both advantages for estimation.

Specifically, in @sec:fun_rp_auto_cor_estimation_with_multiple_realizations, we estimated $r_x (k, l)$. If ergodicity holds (which implies stationarity), then the auto-correlation function depends only on the time difference. We can then use all $hat(r)_x (n, n+k)$ values (those located on the same diagonal in $hat(bold(R))_x$) to estimate $r_x (k)$.

= Digital Signal Processing

#text(fill: red, "(TODO)") Main content: DTFT, z-transform, frequency domain characteristics, stability, power, energy, etc.

== Spectral Analysis <sec:fun_dsp_sa>

#text(fill: red, "(TODO)") Spectra of non-periodic signals are continuous, while those of discrete signals are periodic; Dirichlet kernel, windowing, DFT; zero-padding for increased density; zero-crossing and resolution; harmonic height and resolution, other window types; WSS random signals direct transform, averaged periodogram, BPSK example converging to the Dirichlet kernel.

For a long signal that does not satisfy the stationarity assumption, analyzing its spectrum over all time is not very meaningful because it changes over time. Therefore, we usually cut the signal into small segments for analysis. The method of cutting is direct truncation, assuming function values outside the segment are $0$. We define the Dirichlet kernel function:

$
w_R [n] = cases(
    1\, quad &0<=n<=N-1,
    0\, &"otherwise"
)
$ <equ:fun_dsp_sa_dirichlet_kernel>

The truncated signal is:

$ x_N [n] = x[n] w_R [n] $

Multiplication in the time domain corresponds to convolution in the frequency domain. Specifically:

$ X_N (omega) = 1/(2 pi) {X * W_R}(omega) $

Examining this operation from a graphical and intuitive perspective:

#text(fill: red, "(TODO)") Spectrum of W_R, frequency shifting, the intuition of replacing an impulse function with a peak of a certain width.

= Optimization <sec:fun_optimization>

#text(fill: red, "(TODO)") Mainly regarding Least Squares, Lagrange Multipliers, etc.

#text(fill: red, "(TODO)") Need to write about splitting complex variables into themselves and their conjugates, and the reason for taking partial derivatives with respect to the conjugate during solving.