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
    title: [TODO],
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

= 一些在不同地方出现的玩意

== Kalman Filtering

+ 贝叶斯滤波的一种实现（相较粒子滤波）
+ 从观测器角度
+ SDSP 中的推导
+ F&I 中的推导（最小二乘角度、……）
+ ……

== RLS

+ SDSP 中的推导
+ F&I 中的推导

= 一些在不同地方出现的概念

== 最小方差无偏估计 (MVUE)
