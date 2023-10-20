# Murphy et al. (2022) Formation and growth of co-culture tumour spheroids: new mathematical models and experiments https://doi.org/10.1101/2022.12.21.521515

Preprint available on bioRxiv: [https://doi.org/10.1101/2022.12.21.521515](https://doi.org/10.1101/2022.12.21.521515)

This repository holds key Julia code and all experimental data used to generate figures in the the manuscript (figure numbering corresponds to version 3 on bioRxiv).

Please contact Ryan Murphy for any queries or questions.

Code developed and run in December 2022 using:

- Julia Version  1.7.2 (see https://julialang.org/downloads/ )
- Julia packages: Plots, LinearAlgebra, NLopt, .Threads, Interpolations, Distributions, Roots, LaTeXStrings, DifferentialEquations, CSV, DataFrames, Parsers

## Guide to using the code
The script InstallPackages.jl can be used to install packages (by uncommenting the relevant lines). There are five scripts and four data files summarised in the table below. Each script loads experimental data, plots the experimental data, simulates the mathematical model, computes the MLE (Fig 5,6,8 only), profile likelihoods (Fig 5,6,8 only), and bounds for approximate confidence intervals (Fig 5,6,8 only).



| | Script        | Model           | Data  | Data Source |
| :---:   | :---: | :---: |:---: |:---: |
|1| Fig005.jl    |biphasic model | 1205Lu  spheroid size  | Fig005data.csv |
|2| Fig006.jl      | linear model      |   WM983b spheroid size   |  Fig006data.csv|
|3|  Fig008.jl |  monoculture reduced Greenspan model     |   WM983b spheroid size and structure  | Fig008data.csv|
|4| Fig010AB.jl  | two compartment two population reduced Greenspan model (heterogeneous s)    |    WM983b spheroid size and structure |  Fig010data.csv|
|5| Fig010GH.jl | three compartment two population reduced Greenspan model (heterogeneous Rd and additional cell migration)       |   WM983b spheroid size and structure | Fig010data.csv|****
