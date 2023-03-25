# Finite Element Method

This repository contains my first project entirely written in Rust for Differential Equation course at AGH University of Science and Technology.

## Overview

Main aim was to solve equation below with Finite Element Method 

$$-k(x) u''(x) = 100x$$

for boundary conditions:

$$ u(2) = 0 $$

$$ u'(0) + u(0) = 20$$

where 

$$ k(x) = 
    \begin{cases}
        x + 1 & \text{dla } x \in [0,1]\\
        2x & \text{dla } x \in (1,2]
    \end{cases}$$

and for  a sought-after function $u$ such that:

$$[0,2] \ni x \mapsto  u(x) \in \mathbb{R}$$

## Technology

This project was written in Rust (edition 2021). To compile it, you need to have the Rust compiler installed.

## Running

To run this project use:

```
cargo run
```

## Build

To build this project, use:

```
cargo build --release
```

## Elaboration

[`/latex`](latex) directory contains the task elaboration along with relevant calculations in Polish.