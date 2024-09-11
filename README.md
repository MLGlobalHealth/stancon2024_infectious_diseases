# Infection Disease Modelling using Stan
Here is the code for the this introductory tutorial at StanCon 2024. It is split into two parts with a Rmd file associated with each part:
1) Compartmental model fitting.
2) Hawkes Process model fitting.

We will start by briefly covering the disease modelling theory associated with each model and then there will be plenty of time to work through the two examples.

## Packages
The following packages are necessary for this tutorial:

```install.packages(c("outbreaks", "tidyverse", "rstan", "gridExtra", "bayesplot", "loo"))```

## Getting involved
I am part of a Turing funded initative to develop a community of probabilistic programme users. Our core areas of focus are Health, Climate & Environment, and Spatial analyses, which align with the [Grand Challenges of the Turing Institute](https://www.turing.ac.uk/research). Please click [here](https://turingppl.github.io) for more information about our forecasting competitions and training material.