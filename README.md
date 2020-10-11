# Bayes-INLA-presence-only
Files for analysing spatio-temporal plant species presence trends using presence-only data.

The idea of this project is to analyse the spatial and temporal trends of rare and scarce plant species observations gathered by Thames Valley Environmental Records Centre (TVERC).
The data are presence-only records, meaning that there is no direct information on species absence locations. This forbids the direct inference species observations probabilites
using e.g. logistic regression. Instead, we build a Bayesian additive generative model that allows us to extract relative changes in the local rare species observation probabilities
in comparison to relative changes for a selection of widespread `baseline' species.

The TVERC data are not freely available as per licence agreement.
