#import "@preview/ilm:1.4.1": *

#import "@preview/cuti:0.2.1": show-cn-fakebold
#show: show-cn-fakebold

// #set text(lang: "en")
#set text(
    font: (
        (name: "libertinus serif", covers: "latin-in-cjk"),
        // (name: "Times New Roman", covers: "latin-in-cjk"),
        // "SimSun"
    ),
    lang: "en"
)

#show: ilm.with(
    title: [Notes of Estimation and Detection],
    author: "Gralerfics",
    // date: datetime(year: 2025, month: 10, day: 12),
    date: datetime.today(),
    // abstract: [
    //     Abstract.
    // ],
    // preface: [
    //     #align(center + horizon)[
    //         Preface.
    //     ]
    // ],
    // bibliography: bibliography("refs.bib"),
    // figure-index: (enabled: true),
    // table-index: (enabled: true),
    // listing-index: (enabled: true),
)

#include "sections/introduction.typ"
#include "sections/mvue.typ"
#include "sections/crb.typ"
#include "sections/blue_mle.typ"
#include "sections/ls.typ"
#include "sections/bayesian.typ"
#include "sections/wiener_kalman.typ"
#include "sections/detection_intro.typ"
#include "sections/detection_deterministic.typ"
#include "sections/detection_stochastic.typ"
#include "sections/detection_cht.typ"
