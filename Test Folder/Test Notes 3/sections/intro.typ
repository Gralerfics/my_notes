#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/cetz-plot:0.1.3"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Introduction <sec:intro>

TODO.

== Test Subsection 1

TODO.

=== Haaaa

Haaaa.

#figure(
    caption: [The Fourier transforms of Bartlett windows]
)[
    #cetz.canvas({
        import cetz.draw: *
        import cetz-plot: plot

        let samples = 2048
        
        let func(N) = (w) => {
            let sden = calc.sin(w / 2)
            if (calc.abs(sden) < 1e-5) {
                N
            } else {
                (1.0 / N) * (calc.pow(calc.sin(N * w / 2), 2)) / (calc.pow(sden, 2))
            }
        }

        set-style(
            axes: (stroke: .5pt, tick: (stroke: .5pt)),
            legend: (stroke: none, orientation: ttb, item: (spacing: .1), scale: 80%)
        )

        plot.plot(
            size: (12, 6),
            x-tick-step: calc.pi / 6,
            x-format: plot.formats.multiple-of,
            x-label: $omega$,
            y-tick-step: 20,
            y-min: -0.2,
            y-max: 100.0,
            y-label: $W_B (e^(j omega))$,
            legend: "inner-north-east",
            {
                let domain = (-calc.pi / 3, calc.pi / 3)
                plot.add(
                    func(32),
                    domain: domain,
                    samples: samples,
                    label: $N = 32$,
                    style: (stroke: black)
                )
                plot.add(
                    func(64),
                    domain: domain,
                    samples: samples,
                    label: $N = 64$,
                    style: (stroke: blue)
                )
                plot.add(
                    func(256),
                    domain: domain,
                    samples: samples,
                    label: $N = 256$,
                    style: (stroke: purple)
                )
            }
        )
    })
] <fig:se_bartlett_window_freq>

Test @fig:se_bartlett_window_freq.

=== OHHHH

OHHHH.

$
mat(
    delim: "[",
    #Cpu($y_n$) ;
    y_(n-1);
    y_(n-2);
    dots.v;
    y_p;
)
=
mat(
    delim: "[",
    #Cre($x_n$), x_(n-1), x_(n-2), dots, #Cbl($x_(n-p+1)$) ;
    x_(n-1), x_(n-2), x_(n-3), dots, x_(n-p);
    x_(n-2), x_(n-3), x_(n-4), dots, x_(n-p-1);
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x_p, x_(p-1), x_(p-2), dots, x_1;
)
mat(
    delim: "[",
    #Cre($theta_1$) ;
    theta_2;
    theta_3;
    dots.v;
    #Cbl($theta_p$) ;
)
:= F theta
$

== Test Subsection $theta_1$

TODO.
