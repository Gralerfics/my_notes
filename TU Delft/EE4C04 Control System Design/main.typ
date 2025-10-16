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
    title: [Notes of Control System Design],
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

#include "sections/fundamentals.typ"
#include "sections/state_feedback.typ"
#include "sections/output_feedback_and_gcs.typ"
#include "sections/discrete_time_lsc.typ"
#include "sections/mpc.typ"
