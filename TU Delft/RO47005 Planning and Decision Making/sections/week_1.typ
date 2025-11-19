#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Introduction

autonomous robots
    perception & localization (not covered)
    cognition
    motion control

planning
    route planning
    behavioral layer
    motion planning
    feedback control

= Representations and Transformations

P6, A. Geometric modelling (map, obstacles)

    world, obstacles notation

    representations:
        pointcloud,
        map (with obstacles): bitmaps/gridmaps, volumetric, octomap, occupancy map
        obstacles: (non-/convex) polygon, polyhedron, semi-algebraic models, 3d triangles, NURBS
    summary & others:
        signed distance fields
        Gaussian splats
        world models

P31, B. Rigid-body transformations

    coordinate frames (world frame, body frame)
    
    R, T (2D)
    3D translation + rotation: euler angles, quaternions, rotation matrices
    homogeneous transformation (chains)

P37, assumptions

= Configuration Space and Planning Fundamentals

== Introduction

P4, motion planning, intro

== Configuration Space

=== PDM Definitions

P9, definitions, W, O, C
    examples

=== Topological Concepts

P20, topological spaces

P21, homeomorphism

P25, manifolds
    (cartesian product)
    examples

=== Groups

P33, groups
    matrix groups, GL, O, SO, E, SE

=== TODO

distance metric
    3D rotations
    R3 x SO(3)

Notes & typical configuration spaces

== Obstacles in Configuration Space

P55, definitions, Cobs, Cfree, Wfree
    free/collision configurations

representations of the (obstacle-free) C space

== Motion Planning Formulation

P70, definitions, Path

problem formulation (nonholonomic)

== Time-varying Problems

modified definitions/solutions

== State Space

planning problem formulation in state spaces

== Properties

P83

complete, optimal, anytime, probabilistically complete, asymptotic optimal

// == Summary

= Robot Models

== Manipulators

P4, Joints

configurations

FK (forward kinematics)

IK (inverse kinematics)

== Wheeled Robots

P14, (non-)holonomic
    examples

mobile robot maneuverability（机动性）
    ICR, instantaneous center of rotation

examples & dynamics

== Aerial Vehicles

TODO

= Local Collision Avoidance

(Reactive Planning)

bug, potential fields, dynamic window

dynamic environments?

    velocity obstacles, VO
        ORCA, NH-ORCA
        extensions
    
higher dimensions?

    geometric fabrics
