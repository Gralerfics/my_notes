#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Deterministic Modelling Identification

After the model is established, we still need to perform parameter identification, which means determining the values of the parameter sequences $a[k]$ and $b[k]$. The goal of selecting parameters is to make the model output $hat(x)[n]$ as close as possible to the target signal $x[n]$.

== Least Squares (LS) Method <sec:dmi_ls_method>

First, let's discuss Deterministic Modelling. We hope that the model output $hat(x)[n]$ can accurately reproduce the target signal $x[n]$, meaning the value of each sampling point should be as close as possible. Define the error signal $e'[n] = x[n] - hat(x)[n]$. The model parameters can be determined by minimizing the mean square error $cal(E)_"LS" = sum_(n=0)^infinity abs(e'[n])^2$, as shown in @fig:deterministic_model_identification_diagram_ls.

#figure(
    caption: [Diagram for deterministic model identification (intractable)]
)[
    #diagram(
        spacing: (10mm, 8mm),
        node-stroke: 0.8pt,
        edge-stroke: 0.8pt,
        node((2, 0), [$plus$], inset: 4pt),
        node((0, 0), [$H(z) = frac(B(z), A(z))$], inset: 8pt),
        edge((-2, -1), "rrrr,d", "-|>", [$x[n]$], label-pos: 0),
        edge((-2, 0), "rr", "-|>", [$delta[n]$], label-pos: 0),
        edge((0, 0), "rr", "-|>", [$h[n]$], label-pos: 0.35),
        edge((0, 0), "rr", "-|>", [$-$], label-pos: 0.94, label-side: right),
        edge((2, 0), "rr", "-|>", [$e'[n] = x[n] - h[n]$], label-pos: 1.2),
    )
] <fig:deterministic_model_identification_diagram_ls>

This optimization problem can be solved by solving the following system of equations (for the reason of taking partial derivatives with respect to the conjugate of the variables, see @sec:fun_optimization):

$
cases(
    frac(partial cal(E)_"LS", partial a^*[k]) = 0\, quad quad & k = 1\, 2\, ...\, p,
    frac(partial cal(E)_"LS", partial b^*[k]) = 0\, & k = 0\, 1\, ...\, q
)
$

However, this system of equations is non-linear and very complex to solve.

// #blockquote[
//     Specifically, from Parseval's theorem, we know:

//     $ cal(E)_"LS" = 1/(2 pi) integral_(-pi)^pi abs(E'(e^(j omega)))^2 dif omega $

//     #text(fill: red, "(TODO)")
// ]

== Padé Approximation

Notice that the ARMA model we use has $p+q+1$ parameters, which means it has $p+q+1$ degrees of freedom. Therefore, theoretically, we can use it to perfectly fit the first $p+q+1$ samples of the signal. Let's consider this task first.

We perform a transformation on the form of the transfer function:

$ H(z) = B(z) / A(z) quad => quad H(z) A(z) = B(z) quad => quad h[n] * a[n] = b[n] $

Expanding the convolution gives:

$ h[n] + sum_(k=1)^p a[k] h[n-k] = b[n] $

For the unit impulse input $delta[n]$, the system output $h[n]$ is our estimated result $hat(x)[n]$. To completely fit the first $p+q+1$ samples, directly substitute $h[n] = x[n], 0<=n<=p+q$, to get:

$
x[n] + sum_(k=1)^p a[k] x[n-k] = cases(
    b[n]\, quad &n=0\, 1\, dots\, q,
    0\, &n=q+1\, q+2\, dots\, q+p
)
$ <equ:deterministic_model_identification_equation_pade>

This is the linear system of equations we want, containing $a[dot]$, $b[dot]$, and the known constant $x[n]$.

#blockquote[
    The above transformation process, reflected in the system block diagram, is multiplying both paths by the denominator $A(z)$ of $H(z)$. That is, let the new target error be $E(z) = A(z) E'(z) = A(z) X(z) - B(z)$, as shown in @fig:deterministic_model_identification_diagram_pade.

    #figure(
        caption: [Diagram for deterministic model identification]
    )[
        #diagram(
            spacing: (10mm, 8mm),
            node-stroke: 0.8pt,
            edge-stroke: 0.8pt,
            node((2, 0), [$plus$], inset: 4pt),
            node((0, 0), [$B(z)$], inset: 8pt),
            node((0, -1), [$A(z)$], inset: 8pt),
            edge((-2, -1), "rr", "-|>", [$x[n]$], label-pos: 0),
            edge((0, -1), "rr,d", "-|>", [$hat(b)[n]$], label-pos: 0.155),
            edge((-2, 0), "rr", "-|>", [$delta[n]$], label-pos: 0),
            edge((0, 0), "rr", "-|>", [$b[n]$], label-pos: 0.35),
            edge((0, 0), "rr", "-|>", [$-$], label-pos: 0.94, label-side: right),
            edge((2, 0), "rr", "-|>", [$e[n] = hat(b)[n] - b[n]$], label-pos: 1.2),
        )
    ] <fig:deterministic_model_identification_diagram_pade>

    The output of the unit impulse signal passing through $B(z)$ is $b[n]$. By subtracting it from $hat(b)[n]$ (which can be viewed as an estimate of the coefficient sequence $b[n]$) obtained by passing the target signal $x[n]$ through $A(z)$, we get a new error.
    
    The equations listed after such operations are linear.
]

Next is the solution to @equ:deterministic_model_identification_equation_pade. This system of equations contains the same number of equations and unknowns. If it is non-singular, it can be solved for a unique solution. To be clearer, let's draw the matrix:

#let pade_matrix = [$
mat(
    delim: "[",
    column-gap: #0.7em,
    row-gap: #0.6em,
    x[0], 0, 0, dots, 0;
    x[1], x[0], 0, dots, 0;
    x[2], x[1], x[0], dots, 0;
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[p], x[p-1], x[p-2], dots, x[0];
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[q], x[q-1], x[q-2], dots, x[q-p];
    x[q+1], x[q], x[q-1], dots, x[q-p+1];
    x[q+2], x[q+1], x[q], dots, x[q-p+2];
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[q+p], x[q+p-1], x[q+p-2], dots, x[q];
)
mat(
    delim: "[",
    row-gap: #0.6em,
    1;
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
mat(
    delim: "[",
    row-gap: #0.6em,
    b[0];
    b[1];
    b[2];
    dots.v;
    b[p];
    dots.v;
    b[q];
    0;
    0;
    dots.v;
    0;
)
$ <equ:deterministic_model_identification_matrix_pade>]

#align(center)[
    #cetz.canvas({
        import cetz.draw: *

        content(
            (0mm, 0mm),
            [#pade_matrix]
        )

        set-style(
            stroke: 0.8pt,
            mark: (transform-shape: false)
        )

        let sub_box = (x, y, w, h, color, label, lx, ly, lc) => {
            set-style(stroke: (paint: color, thickness: 0.9pt, dash: (3pt, 2pt)))
            rect((x, y), (x + w, y + h))
            if (label != none) {
                content(
                    (x + w / 2 + lx, y + h / 2 + ly),
                    text(fill: lc, label)
                )
            }
        }

        sub_box(-6.25, -0.97, 9.485, 4.65, blue, $bold(X)_0$, 0, 2.6, blue)
        sub_box(-6.25, -1.06, 1.55, -2.61, purple, $bold(x)_(q+1)$, 0, -1.68, purple)
        sub_box(-4.6, -1.06, 7.835, -2.61, purple, $bold(X)_q$, 0, -1.68, purple)
        
        sub_box(3.7, -1.64, 0.79, 2.62, green, $dash(bold(a))$, 0, -1.6, green)
        
        sub_box(5.5, -0.97, 0.78, 4.65, gray, $bold(b)$, 0, 2.6, gray)
    })
]

Here, we first use the bottom half (the last $p$ rows) to solve for $dash(bold(a))$ (i.e., $a[dot]$):

$
mat(delim: "[", bold(x)_(q+1), bold(X)_q) bold(a) = bold(0)
quad &<=> quad
mat(delim: "[", bold(x)_(q+1), bold(X)_q) mat(delim: "[", 1; dash(bold(a))) = bold(0) \
quad => quad
bold(X)_q dash(bold(a)) = - bold(x)_(q+1)
quad &=> quad
dash(bold(a)) = - bold(X)^(-1)_q bold(x)_(q+1)
$ <equ:deterministic_model_identification_pade_a_bar>

Note that $bold(X)_q$ is a non-symmetric Toeplitz matrix, for which efficient specialized methods like the Trench algorithm exist for solving its inverse. Next, substitute into the upper half (the first $q+1$ rows) to obtain $b[dot]$:

$ bold(b) = bold(X)_0 mat(delim: "[", 1; dash(bold(a))) $ <equ:dmi_b_X0a>

The Padé method is straightforward, but it clearly presents several issues:
+ It does not guarantee that the resulting system is stable;
+ It only constrains the first $p + q + 1$ samples of the model output $hat(x)[n]$ and the target signal $x[n]$ to be identical, and the matching performance beyond that may be poor;
+ $bold(X)_q$ might be singular and thus unsolvable.

#blockquote[
    Regarding the case where $bold(X)_q$ is singular and unsolvable, it can be considered that there is an issue with the default assumption of $a[0] = 1$ in the model. If the model is modified to let $a[0] = 0$, then the equations, although singular, are not unsolvable; rather, the solution is not unique.

    If it is not an all-pole model, the new transfer function obtained this way will have a factor $z$ in both the numerator and denominator, which can be canceled out. This is essentially a pole-zero cancellation occurring at zero. In terms of results, it is equivalent to a reduction in the model order, meaning there is redundancy in the model order.
]

== Prony's Method

The Padé method uses all degrees of freedom on the first $p+q+1$ terms of the sequence, whereas the idea of Prony's method is simple: reduce the fitting requirements for this initial segment of the sequence to obtain a better fit from the perspective of the overall signal.

=== Prony Normal Equations <equ:dmi_prony_normal_equations>

We first derive Prony's method from a more formal perspective. Specifically, we follow the idea of the Padé method to transform the problem into a linear one, as shown in @fig:deterministic_model_identification_diagram_pade. We write the expression for the entire signal error, not just the first $p+q+1$ terms:

$
e[n] = cases(
    x[n] + sum_(k=1)^p a[k] x[n-k] - b[n]\, quad &n=0\, 1\, dots\, q,
    x[n] + sum_(k=1)^p a[k] x[n-k]\, &n>q
)
$ <equ:deterministic_model_identification_error_prony>

In Prony's method, we first solve for $a[dot]$ by minimizing the mean square error:

$
epsilon_(p, q) = sum_(n=q+1)^infinity abs(e[n])^2 = sum_(n=q+1)^infinity abs(x[n] + sum_(k=1)^p a[k] x[n-k])^2
$ <equ:dmi_epsilon_pq>

Considering only the error for the $n>q$ part is to make this part depend only on $a[dot]$. This is based on the need to solve for $a[dot]$ and $b[dot]$ in steps. It may sacrifice some accuracy in terms of definition, but its impact relative to an infinitely long $x[n]$ is not significant. Next, we formally take the partial derivatives and set them to zero to calculate the optimal values:

$
(partial epsilon_(p, q))/(partial a^*[k]) = sum_(n=q+1)^infinity (partial [e[n] e^*[n]])/(partial a^*[k]) = sum_(n=q+1)^infinity e[n] (partial e^*[n])/(partial a^*[k]) = 0, quad k=1, 2, dots, p
$ <equ:dmi_prony_normal_equ_derivation_startpoint>

From the definition @equ:deterministic_model_identification_error_prony, we know that $(partial e^*[n])/(partial a^*[k]) = x^*[n-k]$. Substituting this, we get:

$
sum_(n=q+1)^infinity e[n] x^*[n-k] = 0, quad k=1, 2, dots, p
$ <equ:dmi_prony_orthogonality_principle>

This equation expresses the orthogonal relationship between the minimum error and the signal, known as the Orthogonality principle. We continue by substituting the definition @equ:deterministic_model_identification_error_prony (note that the letter $k$ is used, so we use $l$ instead to avoid confusion):

$
sum_(n=q+1)^infinity (x[n] + sum_(l=1)^p a[l] x[n-l]) x^*[n-k] = 0, quad k=1, 2, dots, p
$

Moving terms and rearranging the order of the summation symbols, we can obtain:

$
sum_(l=1)^p a[l] 

(sum_(n=q+1)^infinity x^*[n-k] x[n-l]) = -sum_(n=q+1)^infinity x^*[n-k] x[n], quad k=1, 2, dots, p
$

To simplify the expression, we denote:

#emphasis_equbox([
$
r_x (k, l) := sum_(n=q+1)^infinity x^*[n-k] x[n-l]
$ <equ:deterministic_model_identification_prony_ne_rx>
])

We can incidentally observe that $r_x (k, l) = r_x^* (l, k)$. Substituting this into the original equation gives:

#emphasis_equbox([
$
sum_(l=1)^p a[l] r_x (k, l) = -r_x (k, 0), quad k = 1, 2, dots, p
$ <equ:dmi_prony_normal_equ>
])

These are called the *Prony normal equations*. Written in matrix form:

#emphasis_equbox([
$
mat(
    delim: "[",
    r_x (1, 1), r_x (1, 2), dots, r_x (1, p);
    r_x (2, 1), r_x (2, 2), dots, r_x (2, p);
    dots.v, dots.v, dots.down, dots.v;
    r_x (p, 1), r_x (p, 2), dots, r_x (p, p);
)
mat(
    delim: "[",
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
-mat(
    delim: "[",
    r_x (1, 0);
    r_x (2, 0);
    dots.v;
    r_x (p, 0);
)
$
])

Denoted as:

#emphasis_equbox([
$
bold(R)_x dash(bold(a)) = - bold(r)_x
$ <equ:dmi_prony_normal_equ_matrix>
])

It can be found that $bold(R)_x$ is a Hermitian matrix. After obtaining $a[dot]$, we can substitute it back into @equ:deterministic_model_identification_error_prony, setting the error to $0$ for $n=1, 2, dots, q$, to obtain $b[k]$.

#blockquote[
    However, we need to note one thing: the original problem of minimizing $e[n]$ is still a joint non-linear least squares problem, where $a[dot]$ and $b[dot]$ are coupled. Therefore, our method of solving for $a[dot]$ and $b[dot]$ sequentially in two steps is actually a simplification and *cannot guarantee global minimization of the original error*.
]

// #blockquote[
//     *Note*, although the definition of $r_x (k, l)$ is very similar to the sample autocorrelation function of a signal, here it is essentially just a notation substitution—do not confuse the two.

//     We are not estimating the autocorrelation function, nor are we making any assumptions about the stationarity or ergodicity of the stochastic process generating the signal. Our logic is this: we naturally derived the equation using Prony's method, which happens to take a form that can be substituted with an autocorrelation matrix, rather than requiring autocorrelation calculation in principle. Therefore, we do not obsess over the detailed differences in formulaic form from @equ:fun_rp_autocor_matrix, as it holds little significance.
    
//     We only need to recognize that the formula indeed shares a physical meaning similar to autocorrelation by definition—that is, the degree of similarity between a signal and its delayed version—with the primary purpose of broadening our understanding.

//     #text(fill: red, "（TODO）") // If the signal is truly stationary, it means this can estimate the distribution behind the signal; if the signal itself is not, it can only fit this specific signal.
// ]

=== An Equivalent Perspective from Pseudoinverse <sec:dmi_equivalent_perspective_from_pseudoinverse>

In the derivation of the previous section, we naturally applied the least squares method to treat the problem as an optimization problem. In fact, we can also directly let all $hat(x)[n]=x[n]$, resulting in an overdetermined system of equations:

#let prony_matrix = [$
mat(
    delim: "[",
    column-gap: #0.7em,
    row-gap: #0.6em,
    x[0], 0, 0, dots, 0;
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[q], x[q-1], x[q-2], dots, x[q-p];
    x[q+1], x[q], x[q-1], dots, x[q-p+1];
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[q+p], x[q+p-1], x[q+p-2], dots, x[q];
    dots.v, dots.v, dots.v, dots.down, dots.v;
)
mat(
    delim: "[",
    row-gap: #0.6em,
    1;
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
mat(
    delim: "[",
    row-gap: #0.6em,
    b[0];
    dots.v;
    b[q];
    0;
    dots.v;
    0;
    dots.v;
)
$ <equ:deterministic_model_identification_matrix_prony>]

#align(center)[
    #cetz.canvas({
        import cetz.draw: *

        content(
            (0mm, 0mm),
            [#prony_matrix]
        )

        set-style(
            stroke: 0.8pt,
            mark: (transform-shape: false)
        )

        let sub_box = (x, y, w, h, color, label, lx, ly, lc) => {
            set-style(stroke: (paint: color, thickness: 0.9pt, dash: (3pt, 2pt)))
            rect((x, y), (x + w, y + h))
            if (label != none) {
                content(
                    (x + w / 2 + lx, y + h / 2 + ly),
                    text(fill: lc, label)
                )
            }
        }

        sub_box(-6.25, 0.4, 9.485, 1.95, blue, $bold(X)_0$, 0, 1.25, blue)
        sub_box(-6.25, 0.3, 1.55, -2.65, purple, $bold(x)_(q+1)$, 0, -1.7, purple)
        sub_box(-4.6, 0.3, 7.835, -2.65, purple, $bold(X)_q$, 0, -1.7, purple)
        
        sub_box(3.7, -1.64, 0.79, 2.62, green, $dash(bold(a))$, 0, -1.6, green)
    })
]

Next, we can "internalize" the process of minimizing the mean square error into the process of solving for the least squares solution of this system using the pseudoinverse. The two approaches are essentially equivalent; the former has a smoother logic, while the latter helps in understanding the problem from the perspective of linear space. Solving using the pseudoinverse yields:

$ dash(bold(a)) = - text(fill: #purple, bold(X)^+_q) bold(x)_(q+1) = - text(fill: #purple, (bold(X)_q^H bold(X)_q)^(-1) bold(X)_q^H) bold(x)_(q+1) $

That is, the optimal coefficient $dash(bold(a))$ will be the solution to the following system of equations:

$ (bold(X)_q^H bold(X)_q) dash(bold(a)) = - bold(X)_q^H bold(x)_(q+1) $

After making the following substitutions, we again obtain the Prony normal equations as in @equ:dmi_prony_normal_equ_matrix:

$ bold(R)_x = bold(X)_q^H bold(X)_q, quad bold(r)_x = bold(X)_q^H bold(x)_(q+1) $

It can be verified by calculation that $bold(R)_x$ is consistent with the definition in the previous section:

$
bold(R)_x
&=
bold(X)_q^H bold(X)_q \
&=
inline(
    mat(
        delim: "[",
        x^*[q], x^*[q+1], x^*[q+2], dots;
        x^*[q-1], x^*[q], x^*[q+1], dots;
        dots.v, dots.down, dots.v, dots.down;
        x^*[q-p+1], x^*[q-p+2], x^*[q-p+3], dots;
    )
    mat(
        delim: "[",
        x[q], x[q-1], dots, x[q-p+1];
        x[q+1], x[q], dots, x[q-p+2];
        x[q+2], x[q+1], dots, x[q-p+3];
        dots.v, dots.v, dots.down, dots.v;
    )
) \
&=
display(sum_(n=q+1)^infinity) inline(mat(
    delim: "[",
    column-gap: #1.0em,
    row-gap: #0.5em,
    x^*[n-1] x[n-1], x^*[n-1] x[n-2], dots, x^*[n-1] x[n-p];
    x^*[n-2] x[n-1], x^*[n-2] x[n-2], dots, x^*[n-2] x[n-p];
    dots.v, dots.v, dots.down, dots.v;
    x^*[n-p] x[n-1], x^*[n-p] x[n-2], dots, x^*[n-p] x[n-p];
)) \
&=
inline(
    mat(
        delim: "[",
        column-gap: #1.0em,
        row-gap: #0.5em,
        r_x (1, 1), r_x (1, 2), dots, r_x (1, p);
        r_x (2, 1), r_x (2, 2), dots, r_x (2, p);
        dots.v, dots.v, dots.down, dots.v;
        r_x (p, 1), r_x (p, 2), dots, r_x (p, p);
    )
)
$

This perspective also provides other useful information: regarding matrices of the form $A^H A$, for any vector $bold(a)$, we have $bold(a)^H (A^H A) bold(a) = (bold(A a))^H (bold(A a)) = norm(bold(A a))^2 >= 0$. This indicates that $A^H A$ is a positive semi-definite matrix.

Consequently, the aforementioned Hermitian matrix $bold(R)_x = bold(X)_q^H bold(X)_q$ is also a (positive semi-)definite matrix. This property determines that $A(z)$ is (marginally) stable, thereby addressing a drawback of the Padé method #text(fill: red, "(TODO, is it? Does it also need to be Toeplitz?)").

Furthermore, if $bold(R)_x$ is a positive definite matrix, its eigenvalues are all positive, meaning the determinant is non-zero, the matrix is invertible, and a solution exists. If $bold(R)_x$ is a positive semi-definite matrix containing zero eigenvalues, it is singular, but this actually indicates redundancy in the model order, which can be tried again after reduction.

=== The Minimum Error and Augmented Normal Equations

Since we are seeking the least squares solution, a minimum error will still exist between the final fitted signal and the true signal. Continuing the derivation from the definitions of $e[n]$ (@equ:deterministic_model_identification_error_prony}) and $epsilon_(p, q)$ (@equ:dmi_epsilon_pq):

$
epsilon_(p, q)
=
sum_(n=q+1)^infinity abs(e[n])^2
&=
sum_(n=q+1)^infinity e[n] (x[n] + sum_(k=1)^p a[k] x[n-k])^* \
&=
sum_(n=q+1)^infinity e[n] x^*[n] + sum_(n=q+1)^infinity e[n] (sum_(k=1)^p a[k] x[n-k])^* \
&=
sum_(n=q+1)^infinity e[n] x^*[n] + sum_(k=1)^p a^*[k] (sum_(n=q+1)^infinity e[n] x^*[n-k])
$

Substituting the Orthogonality principle (@equ:dmi_prony_orthogonality_principle), which holds at the optimal solution, and the definition of $e[n]$ (@equ:deterministic_model_identification_error_prony), we get:

$
epsilon_(p, q)
=
sum_(n=q+1)^infinity e[n] x^*[n]
&=
sum_(n=q+1)^infinity (x[n] + sum_(k=1)^p a[k] x[n-k]) x^*[n] \
&=
(sum_(n=q+1)^infinity x[n] x^*[n]) + sum_(k=1)^p a[k] (sum_(n=q+1)^infinity x[n-k] x^*[n])
$

Using the autocorrelation sequence $r_x (k, l)$ (@equ:deterministic_model_identification_prony_ne_rx), this simplifies to:

$
epsilon_(p, q) = r_x (0, 0) + sum_(k=1)^p a[k] r_x (0, k)
$ <equ:dmi_prony_mini_error>

Once in this form, we can unify $epsilon_(p, q)$ into the equation $bold(R)_x dash(bold(a)) = - bold(r)_x$ (by moving constants to the leftmost column of the matrix):

#emphasis_equbox([
$
mat(
    delim: "[",
    column-gap: #1.0em,
    row-gap: #0.5em,
    augment: #(hline: 1, vline: 1, stroke: (dash: (2pt, 2pt))),
    text(fill: #blue, r_x (0, 0)), text(fill: #blue, r_x (0, 1)), text(fill: #blue, r_x (0, 2)), dots, text(fill: #blue, r_x (0, p));
    text(fill: #red, r_x (1, 0)), r_x (1, 1), r_x (1, 2), dots, r_x (1, p);
    text(fill: #red, r_x (2, 0)), r_x (2, 1), r_x (2, 2), dots, r_x (2, p);
    dots.v, dots.v, dots.v, dots.down, dots.v;
    text(fill: #red, r_x (p, 0)), r_x (p, 1), r_x (p, 2), dots, r_x (p, p);
)
mat(
    delim: "[",
    row-gap: #0.5em,
    augment: #(hline: 1, stroke: (dash: (2pt, 2pt))),
    text(fill: #red, 1);
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
mat(
    delim: "[",
    row-gap: #0.5em,
    augment: #(hline: 1, stroke: (dash: (2pt, 2pt))),
    text(fill: #blue, epsilon_(p, q));
    0;
    0;
    dots.v;
    0;
)
$
])

Or (where $bold(u)_1$ is a unit vector with the first element as 1 and others as 0):

#emphasis_equbox([
$
dash(bold(R))_x bold(a) = epsilon_(p, q) bold(u)_1
$
])

This form is known as the *Augmented normal equations*. #text(fill: red, "(TODO, the book used the same R_x here; for now, an overline is added to distinguish them)")

#blockquote[
    We can also derive this using matrix form, which is more concise. The vector composed of the error sequence is:

    $ bold(e) = bold(X)_q dash(bold(a)) + bold(x)_(q+1) $

    According to the properties of the least squares solution, when we take the optimal solution that minimizes this error, it must be orthogonal to all columns of $bold(X)_q$. That is:

    $ bold(X)_q^H bold(e) = 0 quad <=> quad bold(X)_q^H (bold(X)_q dash(bold(a)) + bold(x)_(q+1)) = 0 quad <=> quad bold(R)_x dash(bold(a)) = - bold(r)_x $

    For a filter of order $p, q$, the minimum error we are concerned with is:

    $ epsilon_(p, q) = norm(bold(e))^2 = bold(e)^H bold(e) = (bold(X)_q dash(bold(a)) + bold(x)_(q+1))^H bold(e) = bold(x)^H_(q+1) bold(e) $

    The last step is because $bold(X)_q^H bold(e) = 0$. Continuing to substitute $bold(e)$, we eventually get the same equation:

    $
    epsilon_(p, q) &= bold(x)^H_(q+1) bold(e) = bold(x)^H_(q+1) (bold(X)_q dash(bold(a)) + bold(x)_(q+1)) \
    &= bold(x)^H_(q+1) bold(x)_(q+1) + (bold(x)^H_(q+1) bold(X)_q) dash(bold(a)) \
    &= r_x (0, 0) + [r_x (0, 1), r_x (0, 2), dots, r_x (0, p)] dash(bold(a)) \
    &= r_x (0, 0) + sum_(k=1)^p a[k] r_x (0, k)
    $
]

== Special Case: All-pole Modelling

We shall study the all-pole model (@equ:ar_model_transfer_function), which is common in many physical processes.

=== All-pole Normal Equations

First, as usual, we treat the process of solving for $a[dot]$ as an optimization problem. Referring to the error definition in @equ:dmi_epsilon_pq, for $q = 0$ we have:

$
epsilon_(p, 0) = sum_(n=1)^infinity abs(e[n])^2
$

The definition of $e[n]$ still comes from @equ:deterministic_model_identification_error_prony, but since $n-k<0$ for $n=0$, the value of the coefficient $x[n-k]$ for $a[dot]$ is $0$, leaving only $x[0] - b[0]$:

$
e[n] = cases(
    x[0] - b[0]\, quad &n=0,
    x[n] + sum_(k=1)^p a[k] x[n-k]\, &n>0
)
$

Now we apply a slight modification. Note that $e[0] = x[0] - b[0]$ can be viewed as a constant with respect to $a[dot]$. Therefore, when solving for $a[dot]$, minimizing $epsilon_(p, 0)$ is equivalent to minimizing a new error we define, $epsilon_(p)$:

$
epsilon_(p) = sum_(n=0)^infinity abs(e[n])^2
$

The difference is simply the inclusion of $e[0]$. Fortunately, after this change to the error function, we still derive the form of the Prony normal equations (see @equ:dmi_prony_normal_equ), *but the definition of $r_x (k, l)$ changes* to:

$
r_x (k, l) := sum_(n=0)^infinity x^*[n-k] x[n-l]
$ <equ:dmi_allpole_rxkl_new>

The difference is that the lower limit of the summation has changed from $q + 1 = 1$ to $0$. *This is the effect of changing the error function*. The specific reason requires re-derivation using the new error definition starting from @equ:dmi_prony_normal_equ_derivation_startpoint.

#blockquote[
    Actually, it's quite simple. After changing to $epsilon_(p)$, we set its partial derivative with respect to $a[dot]$ to zero:

    $
    (partial #text(fill: red, $epsilon_(p)$))/(partial a^*[k]) = sum_(n=#text(fill: red, $0$))^infinity (partial [e[n] e^*[n]])/(partial a^*[k]) = sum_(n=#text(fill: red, $0$))^infinity e[n] (partial e^*[n])/(partial a^*[k]) = 0, quad k=1, 2, dots, p
    $

    From the definition of $e[n]$ above, we still have $(partial e^*[n])/(partial a^*[k]) = x^*[n-k]$. Since $x^*[n-k]=0$ when $n=0$, it happens to be consistent with $(partial (x[0] - b[0]))/(partial a^*[k]) = 0$. Thus, it can be unified and substituted to get:

    $
    sum_(n=#text(fill: red, $0$))^infinity e[n] x^*[n-k] = 0, quad k=1, 2, dots, p
    $

    Continuing to substitute the definition of $e[n]$, and again due to the presence of the $x^*[n-k]$ term, the case for $n=0$ can be merged, directly yielding:

    $
    sum_(n=#text(fill: red, $0$))^infinity (x[n] + sum_(l=1)^p a[l] x[n-l]) #text(fill: blue, $x^*[n-k]$) = 0, quad k=1, 2, dots, p
    $

    Moving terms and rearranging the order of the summation symbols, we can obtain:

    $
    sum_(l=1)^p a[l] 

    (sum_(n=#text(fill: red, $0$))^infinity x^*[n-k] x[n-l]) = -sum_(n=q+1)^infinity x^*[n-k] x[n], quad k=dots
    $

    This equation is consistent with the form of the Prony normal equations, with the only difference being that the definition of $r_x (k, l)$ needs to be as shown in @equ:dmi_allpole_rxkl_new, changing to a summation starting from $0$.
]

At this point, we have obtained the "normal equations" for solving $a[dot]$ in the case of an all-pole model. However, observing further, since $x[n] = 0$ for $n < 0$, substituting into @equ:dmi_allpole_rxkl_new yields:

$
r_x (k+1, l+1) &= sum_(n=0)^infinity x^*[n-(k+1)] x[n-(l+1)] \
&= sum_(n=0)^infinity x^*[n-1-k] x[n-1-l] \
&= sum_(n=-1)^infinity x^*[n-k] x[n-l] \
&= x^*[-1-k] x[-1-l] + sum_(n=0)^infinity x^*[n-k] x[n-l] \
&= 0 + sum_(n=0)^infinity x^*[n-k] x[n-l] \
&= r_x (k, l), quad (forall k, l>=0)
$

Thus, we can let:

$
r_x (k-l) := r_x (k, l) = sum_(n=0)^infinity x^*[n-k] x[n-l]
$

Which is:

#emphasis_equbox([
$
r_x (k) = sum_(n=0)^infinity x^*[n-k] x[n]
$
])

// #blockquote[
//     *Note again*, although this looks very much like the definition of stationarity, and the notations for matrices and vectors follow those of the autocorrelation definition, it must be emphasized that this is merely a coincidental similarity in form.
    
//     Specifically, we find that $bold(R)_x$ carries the meaning of the signal $x[n]$'s autocorrelation matrix, and in the All-pole model, this matrix becomes Toeplitz. We would rather say that this matrix no longer represents the signal's autocorrelation matrix than say that the All-pole model imposes constraints on the signal's properties, causing the signal distribution to become stationary (how could a model affect the signal?).

//     However, this does indicate one thing: although parameter identification methods like Prony can be applied to any signal, if their physical meaning is to align well, it should be assumed that the signal is ...

//     #text(fill: red, "(TODO)") Not quite right either—where does stationarity come from for a deterministic signal? It's still just pure notation.
// ]

Observation shows that $r_x (k)$ is conjugate symmetric, i.e., $r_x (k) = r_x^* (-k)$. Substituting this yields a more concise system of equations applicable to the All-pole model:

#emphasis_equbox([
$
sum_(l=1)^p a[l] r_x (k - l) = -r_x (k), quad k=1, 2, dots, p
$ <equ:dmi_all_pole_norm_equ>
])

Or:

#emphasis_equbox([
$
mat(
    delim: "[",
    r_x (0), r_x^* (1), dots, r_x^* (p-1);
    r_x (1), r_x (0), dots, r_x^* (p-2);
    dots.v, dots.v, dots.down, dots.v;
    r_x (p-1), r_x (p-2), dots, r_x (0);
)
mat(
    delim: "[",
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
-mat(
    delim: "[",
    r_x (1);
    r_x (2);
    dots.v;
    r_x (p);
)
$ <equ:dmi_prony_all_pole_matrix>
])

This is known as the *All-pole normal equations*. Since the matrix $bold(R)_x$ is conjugate symmetric and Toeplitz, it allows us to use the Levinson-Durbin algorithm for efficient computation.

For the calculation of the minimum error value, similarly, we have:

#emphasis_equbox([
$
epsilon_(p) = r_x (0) + sum_(k=1)^p a[k] r_x^*(k)
$ <equ:dmi_all_pole_mini_error>
])

It can also be written together with the system of equations in a form similar to the Augmented normal equations, which will not be repeated here.

=== Issues on the Numerator Selection

By conventional methods, after obtaining $a[dot]$, we would use @equ:dmi_b_X0a to get $b[0] = x[0]$. However, at the end of @equ:dmi_prony_normal_equations, we mentioned that the step-by-step solution method in the Prony method does not guarantee global optimality. Modifying the value of $b[0]$ here does not necessarily destroy the optimality of the result, because the result wasn't optimal to begin with; conversely, we might even achieve better results by changing the way $b[dot]$ is selected.

The so-called better result is not necessarily a reduction in the mean square error value, but might also integrate considerations of other factors. The all-pole model here is an example: if the value of $x[0]$ in the original signal is not so credible due to noise or other interference, to prevent the entire set of model parameters from being biased by this single $b[0] = x[0]$, we prefer to make the energy of the fitted signal $hat(x)[n]$ (which equals the unit impulse response $h[n]$ in our model) equal to the energy of the target signal $x[n]$:

$
r_(hat(x)) (0) = r_h (0) = r_x (0)
$

Derivation shows that one should take $b[0] = sqrt(epsilon_(p))$.

#text(fill: red, "(TODO)") I haven't quite figured out how to derive this value for $b[0]$; the reference book says it will be covered in its section 5.2.3.

== Finite Data Records for All-pole Cases

The previous analysis of the Prony method was based on the assumption that $x[n]$ is defined over the entire positive time domain, from $0$ to $infinity$. Now we need to consider the case where we only have *$N+1$ samples* on $[0, N]$. As for why it is $N+1$ samples instead of $N$, I don't know; perhaps the author felt the subsequent indexing would be more concise. To avoid errors, I have written it this way as well, although it is somewhat inconvenient.

The following two approaches are commonly used for all-pole models, so we will only discuss all-pole models by default.

=== Auto-correlation Method <sec:dmi_finite_data_autocorrelation_method>

The first method is called the autocorrelation method. We consider applying a rectangular window to $x[n]$, or in other words, treating the parts of $x[n]$ outside $[0, N]$ as having a value of $0$:

$
x_N [n] = cases(
    x[n]\, quad &0<=n<=N,
    0\, &"otherwise"
)
$

Then we directly apply the Prony method for the solution. Note that the $r_x (k)$ estimated using the windowed $x_N [n]$ will become:

#emphasis_equbox([
$
r_x (k) = sum_(n=0)^infinity x^*_N [n-k] x_N [n] = sum_(n=k)^N x^*[n-k] x[n], quad k = 0, 1, dots, p
$ <equ:dmi_allpole_finite_autocor_rx>
])

For ease of presentation, we write out the overdetermined system:

$
bold(X)_p dash(bold(a))_p = -bold(x)_1 \
=> mat(
    delim: "[",
    // column-gap: #1.0em,
    row-gap: #0.5em,
    augment: #(hline: (4, 9), stroke: (dash: (2pt, 2pt))),
    x[0], 0, 0, dots, 0;
    x[1], x[0], 0, dots, 0;
    x[2], x[1], x[0], dots, 0;
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[p-1], x[p-2], x[p-3], dots, x[0];
    x[p], x[p-1], x[p-2], dots, x[1];
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[N-2], x[N-3], x[N-4], dots, x[N-p-1];
    x[N-1], x[N-2], x[N-3], dots, x[N-p];
    x[N], x[N-1], x[N-2], dots, x[N-p+1];
    0, x[N], x[N-1], dots, x[N-p+2];
    dots.v, dots.v, dots.v, dots.down, dots.v;
    0, 0, 0, dots, x[N];
)
mat(
    delim: "[",
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
-mat(
    delim: "[",
    row-gap: #0.5em,
    augment: #(hline: (4, 9), stroke: (dash: (2pt, 2pt))),
    x[1];
    x[2];
    x[3];
    dots.v;
    x[p];
    x[p+1];
    dots.v;
    x[N-1];
    x[N];
    0;
    0;
    dots.v;
    0;
)
$ <equ:dmi_allpole_finite_autocor_overdeter_equ>

The normal equations, except for the definition of $r_x (k)$ being modified as described in @equ:dmi_allpole_finite_autocor_rx, remain formally the same as the previous All-pole normal equations (see @equ:dmi_all_pole_norm_equ):

$
mat(
    delim: "[",
    r_x (0), r_x^* (1), dots, r_x^* (p-1);
    r_x (1), r_x (0), dots, r_x^* (p-2);
    dots.v, dots.v, dots.down, dots.v;
    r_x (p-1), r_x (p-2), dots, r_x (0);
)
mat(
    delim: "[",
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
-mat(
    delim: "[",
    r_x (1);
    r_x (2);
    dots.v;
    r_x (p);
)
$

The form of the minimum error for this method is also consistent with @equ:dmi_all_pole_mini_error in the All-pole analysis, requiring only the modification of the autocorrelation function.

The autocorrelation method truncates the signal directly, even if the signal values outside the interval are non-zero, so the results provided may be *biased* compared to the actual solution.

#text(fill: red, "(TODO)") However, this method has an important property: it guarantees that the resulting model is *stable*, which is very useful for cases requiring significant extrapolation or analysis. The book mentions that the proof is covered in Chapter 5 and is omitted here for now.

=== Covariance Method

The Auto-correlation Method sets parts outside the domain to $0$, which essentially changes the form of $x[n]$, as $0$ is also a normal signal value. In some cases, this does not achieve the best results.

The second method, the Covariance Method, typically yields more accurate results. It makes no assumptions about the signal itself but instead ignores samples outside the domain during the optimization process.

Following the standard procedure of defining error and solving the optimization problem, if we cannot consider samples outside the domain, the error can only be defined over the valid interval. Based on the previous error definition, calculating $e[n]$ requires $x[n], x[n-1], dots, x[n-p]$, so we can only define the error on $[p, N]$:

$
cal(E)_p^C = sum_(n=p)^N abs(e[n])^2
$

We can then take partial derivatives of this error with respect to the coefficients and repeat the derivation process to obtain the normal equations, but we will not expand on that here.

Looking at it from the perspective of an overdetermined system of equations is simpler: it essentially involves deleting the expressions in the autocorrelation method that involve samples outside the domain (such as $x[N+1]$, etc.). Referring to the overdetermined system of the autocorrelation method in @equ:dmi_allpole_finite_autocor_overdeter_equ, extracting only the part between the dashed lines gives the overdetermined system for the covariance method:

$
mat(
    delim: "[",
    // column-gap: #1.0em,
    // row-gap: #0.5em,
    x[p-1], x[p-2], x[p-3], dots, x[0];
    x[p], x[p-1], x[p-2], dots, x[1];
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[N-2], x[N-3], x[N-4], dots, x[N-p-1];
    x[N-1], x[N-2], x[N-3], dots, x[N-p];
)
mat(
    delim: "[",
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
-mat(
    delim: "[",
    // row-gap: #0.5em,
    x[p];
    x[p+1];
    dots.v;
    x[N-1];
    x[N];
)
$

Its normal equations are the same as the original Prony normal equations (see @equ:dmi_prony_normal_equ):

$
mat(
    delim: "[",
    r_x (1, 1), r_x (1, 2), dots, r_x (1, p);
    r_x (2, 1), r_x (2, 2), dots, r_x (2, p);
    dots.v, dots.v, dots.down, dots.v;
    r_x (p, 1), r_x (p, 2), dots, r_x (p, p);
)
mat(
    delim: "[",
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
-mat(
    delim: "[",
    r_x (1, 0);
    r_x (2, 0);
    dots.v;
    r_x (p, 0);
)
$

Doing this actually discards the Toeplitz property of the matrix in the All-pole analysis. Similarly, the autocorrelation function needs to be changed to the finite-data version:

$
r_x (k, l) := sum_(n=p)^N x^*[n-k] x[n-l]
$

The minimum error value also formally follows @equ:dmi_prony_mini_error, needing only the modification of the autocorrelation function.

// === Generalizations to the Covariance Method

// #text(fill: red, "(TODO)")

== Example: Channel Inversion

#text(fill: red, "(TODO)")
