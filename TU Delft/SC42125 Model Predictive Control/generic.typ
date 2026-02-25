#import "@preview/ilm:1.4.1": *

#let emphasis_equbox(x) = [
    #box(
        inset: 8pt,
        fill: rgb("#f3f3f3")
    )[#x]
]

#let resize_box(body) = layout(
    container => {
        let size = measure(body)
        let ratio = calc.min(container.width / size.width, container.height / size.height) * 100%
        scale(ratio, body, reflow: true)
    }
)
