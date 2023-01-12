# Murphy et al. (2022) Formation and growth characteristics of co-culture tumour spheroids https://doi.org/10.1101/2022.07.27.501797 

Preprint available on bioRxiv: https://doi.org/10.1101/2022.07.27.501797 


This repository holds key Julia code and all experimental data used to generate figures in the manuscript.

Please contact Ryan Murphy for any queries or questions.

Code developed and run in December 2022 using:

- Julia Version  1.7.2 (see https://julialang.org/downloads/ )
- Julia packages: Plots, LinearAlgebra, NLopt, .Threads, Interpolations, Distributions, Roots, LaTeXStrings, DifferentialEquations, CSV, DataFrames, Parsers


Preprint available on bioRxiv: https://doi.org/10.1101/2022.07.27.501797 

## Guide to using the code
The script InstallPackages.jl can be used to install packages (by uncommenting the relevant lines). There are five scripts and four data files summarised in the table below. Each script loads experimental data, plots the experimental data, simulates the mathematical model, computes the MLE (Fig 2,3,5 only), profile likelihoods (Fig 2,3,5 only), and bounds for approximate confidence intervals (Fig 2,3,5 only).



| | Script        | Model           | Data  | Data Source |
| :---:   | :---: | :---: |:---: |:---: |
|1| Fig2.jl    |biphasic model | 1205Lu  spheroid size  | Fig2data.csv |
|2| Fig3.jl      | linear model      |   WM983b spheroid size   |  Fig3data.csv|
|3|  Fig5.jl |  monoculture reduced Greenspanb model     |   WM983b spheroid size and structure  | Fig5data.csv|
|4| Fig7AB.jl  | two compartment two population reduced Greenspan model (heterogeneous s)    |    WM983b spheroid size and structure |  Fig7data.csv|
|5| FigGH.jl | three compartment two population reduced Greenspan model (heterogeneous Rd and additional cell migration)       |   WM983b spheroid size and structure | Fig7data.csv|****

| Attempt | #1    | #2    |
| :---:   | :---: | :---: |
| Seconds | 301   | 283   |****
