# psweight: IPW- and CBPS-type propensity score reweighting, with various extensions

# Description

`psweight()` is a Mata class that computes inverse-probability weighting (IPW)
weights for average treatment effect, average treatment effect on the treated,
and average treatment effect on the untreated estimators for observational data.
IPW estimators use estimated probability weights to correct for the missing data on the
potential outcomes. Probabilities of treatment--propensity scores--are
computed for each observation with one of variety of methods, including
logistic regression (traditional IPW), covariate-balancing propensity scores
(CBPS), penalized covariate-balancing propensity scores (PCBPS), prognostic score-balancing
propensity scores, and other methods.  It also constructs balance tables and
assesses the distribution of the IPW weights.

`psweight` is a Stata command that offers Stata users easy access to the class.
However, the class offers more flexibility and can conduct some analyses
unavailable with the Stata command.

# The model

`psweight::solve()` and `psweight subcmd` solve for propensity score model
coefficients, propensity scores, and IPW weights as follows:

The first step involves computing coefficients for the propensity
score model, `b`.  The propensity score model takes the form of a logit
regression model.  Specifically, the propensity score for each row in
the data is defined as

```
                    p = invlogit(X * b')
```

where X is the vector of matching variables (tmvarlist) for the
respective row.

You specify a `subcmd` to control how the vector `b` is computed in the
internal numerical optimization problem.  As discussed in Kranker,
Blue, and Vollmer Forrow (2019), we can set up optimization problems
to solve for the `b` that produces the best fit in the propensity score
model, the `b` that produces the best balance on matching variables, the `b`
that produces the best balance on prognostic scores, or something
else.  The `subcmd` also determines how the term "best balance" is
defined in the previous sentence.  That is, for a given `subcmd`, we
can generically define `b` as the vector that solves the problem:

```
                    b = argmin L(X,T,W)
```

where `L(X,T,W)` is a "loss function" that corresponds to the specified
subcmd (e.g., logit regression or CBPS), given the data (`X,T`) and a
vector of weights (`W`).  (The weights are computed using the propensity
scores, as we describe below.  The propensity scores are calculated
using b, the data, and the formula given above.) The available `subcmd`s
are listed in the documentation and include logit regression and
CBPS (Imai and Ratkovic 2014).

In Kranker, Blue, and Vollmer Forrow (2019), we proposed adding a
"penalty" to the loss function that lets you effectively
prespecify the variance (or higher-order moments) of the IPW weight
distribution.  By constraining the distribution of the weights, you
can choose among alternative sets of matching weights, some of which
produce better balance and others of which yield higher statistical
power.  The penalized method solves for `b` in:

```
                b = argmin L(X,T,W) + f(W)
```

where `f(W)` is a smooth, flexible function that increases as the vector
of observation weights (`W`) becomes more variable.  The penalty
options control the functional form of `f(W)`; see details below.

Once the `b` is estimated, we can compute propensity scores (`p`) for
each observation with the formula given above and the observation's
matching variables (`tmvarlist`).  The propensity scores are returned
in a variable named `_pscore`.

Once propensity scores are computed for each observation, we can
compute IPW "matching weights" for each observation.  The formulas
for the IPW weights depend on whether you request weights for
estimating the average treatment effect (`ate`), the average treatment
effect on the treated (`atet`), or the average treatment effect on the
untreated (`ateu`).

Next, the weights are normalized to have mean equal to 1 in each
group, and returned in the variable named `_weight_mtch`.

Finally, the final weights (a variable named _weight) are set equal
to: `_weight = W :* _weight_mtch`, where W are the sample weights.


# Author

Keith Kranker

The code for implementing the CBPS method is based on work by Fong et al.
(2018), namely the CBPS package for R.  I also reviewed the Stata CBPS
implementation by Filip Premik.


# Suggested Citation

* Kranker, K. "psweight: IPW- and CBPS-type propensity score reweighting, with various extensions," Statistical Software Components S458657, Boston College Department of Economics, 2019. Available at https://ideas.repec.org/c/boc/bocode/s458657.html.

or

* Kranker, K., L. Blue, and L. Vollmer Forrow.  "Improving Effect Estimates by Limiting the Variability in Inverse Propensity Score Weights." Manuscript under review, 2019.

Source code is available at https://github.com/kkranker/psweight.
Please report issues at  https://github.com/kkranker/psweight/issues.

# Installation

To install official releases from SSC, type this from your Stata command line:

```stata
. net describe psweight, from(http://fmwww.bc.edu/RePEc/bocode/p)
```

To install the latest version from Github, type this from your Stata command line:

```stata
. net from https://raw.githubusercontent.com/kkranker/psweight/master/
```


# References

* Fong, C., M. Ratkovic, K. Imai, C. Hazlett, X. Yang, and S. Peng.  2018. CBPS: Covariate Balancing Propensity Score, Package for the Rprogramming langauage, The Comprehensive R Archive Network.Available at: https://CRAN.R-project.org/package=CBPS

* Imai, K. and M. Ratkovic.  2014.  "Covariate Balancing Propensity Score."Journal of the Royal Statistical Society: Series B (StatisticalMethodology), 76(1): 243–263, doi:10.1111/rssb.12027.

* Kranker, K., L. Blue, and L. Vollmer Forrow.  2019.  "Improving EffectEstimates by Limiting the Variability in Inverse PropensityScore Weights." Manuscript under review.

# Potential improvements or extensions:
[ ] using constraints to specify CV instead of a penalty
[ ] compute standard errors, etc. in optimization
[ ] run impacts (single- or multi-equation model) and/or pomeans option to estimate potential-outcome means inside optimization (to get SEs)
[ ] noconstant option
[ ] predict command
