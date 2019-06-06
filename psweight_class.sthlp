Title

psweight()
IPW- and CBPS-type propensity score reweighting, with various extensions{p_end}


Syntax

P = psweight()



====================================================

where

    tvar, tmvarlist, tousevar, wgtvar, depvarnames, varnames, and are string scalars (they must contain the names of Stata variables)

    t, w, and p are real column vectors
    x and y0 are real matrices

    denominator, a real (integer) scalar, determines how standardized differences are calculated:

        denominator  variance used
        0            the control groups' variances
        1            the treatment groups' variances; the default
        2            the pooled variances
        3            (control groups' variances + treatment groups' variances)/2  (the definition from Stata's {help tbalance} command)

    stat, a string scalar, specifies whether the estimand of interest for computing IPW weights:

        stat          action
        ate           computes weights for the average treatment effect; the default
        atet          computes weights for the average treatment effect on the   treated
        ateu          computes weights for the average treatment effect on the untreated

        When summarizing the distribution of weights (e.g., the wgt_cv() function), the stat specifies whether to summarize weights for the whole sample (ate, the default), the control group (atet), or the treatment group (ateu).

    subcmd, a string scalar, is defined below.

    cvopt, a real rowvector, is defined below.

    src is a (scalar) instance of the class psweight.

====================================================


Setup:

void P.st_set(treatvar, tmvarlist, tousevar, | wgtvar)

    Loads Stata data for the treatment model into the Mata class, using views wherever possible.

    tvar, tmvarlist, tousevar, and wgtvar contain the names of variables in the Stata data.
    tvar is a variable that must contain values 0 or 1, representing the treatment (1) and comparison (0) groups.
    tmvarlist specifies one or more variables that predict treatment assignment in the treatment model.
    tousevar is a variable that must contain values 0 or 1, representing the rows to include (1) or exclude (0).
    wgtvar (optional) is a variable that specifies sample weights.

void P.st_set_depvar(depvarnames, tousevar)

    Loads Stata data for the dependent variable (control group only) into the Mata class, using views wherever possible.

    depvarnames and tousevar contain the names of variables in the Stata data.
    depvarnames are the variable(s) containing the dependent variable(s).
    tousevar is a variable that must contain values 0 or 1, representing the rows to include (1) or exclude (0).

void P.set(t, x, | w)

    Loads Mata data into the Mata class.

    t must contain values 0 or >0, representing the treatment (>0) and comparison (0) groups.
    x specifies the data that predict treatment assignment in the treatment model.
    w (optional) specifies sample weights; these are treated as iweights.

void P.set_depvar(y0)

    Loads Mata data for the dependent variable (control group only) into the Mata class.

    y0 specifies the dependent variable(s) data for the control group.


void P.st_set(treatvar, tmvarlist, tousevar, | wgtvar)
void P.set(t, x, |  w)


Functions to estimate IPW weights:

real rowvector P.solve(| stat, subcmd, denominator, cvopt)

    This is the psweight class's premier solver function. See the remarks in {help psweight##remarks:help psweight} for a description of the methods.

    The function returns the vector of coefficients for the propensity score model ({it:b}).

    In addition, the function updates (private) variables containing the propensity scores (predicted values) and IPW weights. The type of IPW weights is specified by {it:stat}, and could be weights for estimating average treatment effects {ate}, average treatment effects on the treated {atet}, or average treatment effects on the untreated {ateu}.

    subcmd, a string scalar, specifies which method is used to compute coefficients, {it:b}, for the propensity score model.  In some cases, the method requires defining how standardized differences are calculated ({it:denominator}).

        subcmd       Method used to compute {it:b}
        ipw          logit regression; the default
        cbps         just-identified covariate-balancing propensity score
        cbpsoid      over-identified covariate-balancing propensity score
        mean_sd_sq   minimize mean(stddiff())^2
        sd_sq        minimize sum(stddiff()^2)
        stdprogdiff  minimize sum(stdprogdiff()^2)

    cvopt is analogue to the cvtarget(# # #), skewtarget(# # #), and kurttarget(# # #)  options in {help psweight}.

        Let loss_0 be the objective function corresponding to the specified {it:subcmd}.

        Size of cvopt   Penalty
        1 x 3           If cvopt=(a, b, c),                             the loss function is modified as: {break} loss = (loss_0 + a * abs((wgt_cv() - b)^c))
        1 x 6           If cvopt=(a, b, c, d, e, f),                    the loss function is modified as: {break} loss = (loss_0 + a * abs((wgt_cv() - b)^c) + e * abs((wgt_skewness() - e)^f))
        1 x 9           If cvopt=(a, b, c, d, e, f, g, h, i),           the loss function is modified as: {break} loss = (loss_0 + a * abs((wgt_cv() - b)^c) + e * abs((wgt_skewness() - e)^f) + g * abs((wgt_kurtosis() - h)^i))
        1 x 12          If cvopt=(a, b, c, d, e, f, g, h, i, j, k, l),  the loss function is modified as: {break} loss = (loss_0 + a * abs((wgt_cv() - b)^c) + e * abs((wgt_skewness() - e)^f) + g * abs((wgt_kurtosis() - h)^i)  + j * abs((wgt_max() - k)^l))

        The default is cvopt=(0, 0, 2) (the loss function is unmodified)

real rowvector P.ipw(| stat)

    Fits a logit regression model and then computes predicted values (propensity scores).
    Propensity scores are then used to compute IPW weights, and the weights are applied to the class instance.

    This function is a synonym for
        : P.psweight(stat, "ipw")

real rowvector P.cbps(| stat, denominator)

    Computes the just-identified covariate-balancing propensity scores from Imai and Ratkovic (2014).
    Propensity scores are then used to compute IPW weights, and the weights are applied to the class instance.

    This function is a synonym for
        : P.psweight(stat, "cbps", denominator)

real rowvector P.cbpsoid(| stat, denominator)

    Computes the over-identified covariate-balancing propensity scores from Imai and Ratkovic (2014).
    Propensity scores are then used to compute IPW weights, and the weights are applied to the class instance.

    This function is a synonym for
        : P.psweight(stat, "cbpsoid", denominator)

Post-estimation functions:

void P.reweight(| w,  p)

    Updates the class instance's (private) member variables with the supplied matching weights (w) and propensity scores (p).

    After reweight(), downstream functions (e.g., to construct the balance table) will use the reweighted sample.

    w is treated as a set of matching weights; they will be multiplied by the sample weights (if sample weights exist).

    If no arguments are provided, the matching weights are reset (that is, the sample is no longer weighted with IPW weights are set to one and the propensity scores are set to missing).

    There is no need to call reweight() when estimating IPW weights. Newly calculated IPW weights will automatically be applied to the class instance.

real colvector get_pscore()

    Returns the propensity scores (predicted values from the propensity score model).

real colvector get_weight_mtch()

    Returns the matching weights.

real colvector get_weight()

    Returns the final weights (the matching weights times the sample weights).

void P.get_scores(varnames, tousevar)

    Copies the class instance's IPW weights and propensity scores into Stata data.

    varnames contains the names of four new variables to be updated in the Stata data.
        The first  variable will receive the final weight        (typically _weight)
        The second variable will receive the matching weight     (typically _weight_mtch)
        The third  variable will receive the propensity scores   (typically _pscore)
        The fourth variable will receive the treatment indicator (typically _treated)
    tousevar is a variable that must contain values 0 or 1, representing the rows to include (1) or exclude (0).

real rowvector P.pomean()

   Returns the (weighted) mean of the dependent variable(s) for the control group.

   This function requires that a dependent variable exists.


Functions to assess balance: comparing the treatment and control groups:

real matrix P.balancetable(| denominator)

    Returns the balance table, a k x 6 matrix:
            first column:   mean for the treatment group
            second column:  mean for the control group
            third column:   difference in means, diff()
            fourth column:  standardized differences, stddiff(denominator)
            fifth column:   the standard deviation (denominator) used to compute the standardized difference
            sixth column:   the ratio of variances (treatment variance/control variance), varratio(denominator)

    The table is also returned to Stata in r(bal).

real rowvector diff()

    Returns the difference in means between the treatment and control groups for each variable in tmvarlist.

    The vector is also returned to Stata in r(diff).

real rowvector P.stddiff(| denominator)

    Returns the standardized difference in means between the treatment and control groups for each variable in tmvarlist.

    The vector is also returned to Stata in r(stddiff).

real rowvector varratio()

    Returns the ratio of variances (treatment variance :/ control variance).

    The vector is also returned to Stata in r(varratio).

real scalar P.mean_sd(| denominator)

    Returns the average of the standardized differences, mean(stddiff(denominator)').

    The value is also returned to Stata in r(mean_sd).

    Sometimes you may see r(mean_sd_sq). That value is just r(mean_sd)^2.

real scalar P.mean_asd(| denominator)

    Returns the average of the absolute standardized differences, mean(abs(stddiff(denominator))').

    The value is also returned to Stata in r(mean_asd).

real scalar max_asd(| denominator)

    Returns the maximum absolute standardized difference, max(abs(stddiff(denominator))').

    The value is also returned to Stata in r(max_asd).

real rowvector P.progdiff(| denominator)

    Returns the prognostic score balance table, a L x 5 matrix:
            first column:   mean prognostic score for the treatment group
            second column:  mean prognostic score for the control group
            third column:   difference in mean prognostic scores
            fourth column:  standardized differences, stddiff(denominator)
            fifth column:   the actual mean of the dependent variables in teh control group

    Prognostic scores are generated by: regressing depvar on the tmvarlist using OLS and the control group's data, then computing predicted values (prognostic scores) for all observations.

    The table is also returned to Stata in r(progdiff).

    This function requires that a dependent variable exists.


Functions to summarize the distribution of the IPW weights:

real scalar wgt_cv(stat)

    Returns the coefficient of variation of the IPW weights.

    The value is also returned in Stata in r(wgt_cv).

real scalar wgt_sd(stat)

    Returns the standard deviation of the IPW weights.

    The value is also returned in Stata in r(wgt_sd).

real scalar wgt_skewness(stat)

    Returns the skewness of the IPW weights.

    The value is also returned in Stata in r(wgt_skewness).

real scalar wgt_kurtosis(stat)

    Returns the excess kurtosis of the IPW weights.

    The value is also returned in Stata in r(wgt_kurtosis).

real scalar wgt_max(stat)

    Returns the maximum value of the IPW weights.

    The value is also returned in Stata in r(wgt_max).

Miscellaneous functions:

void P.balanceresults(| stat, denominator)

    This function is a 1-stop shop to call a selection of the functions defined above. The balance table (balancetable()) is always computed. The weight distribution is summarized if any of the current IPW weights do not equal 1. The prognostic scores are compared (progdiff()) if a dependent variable exists.

    Depending on which functions were called, results will be returned to Stata in r().

void P.clone(src)

    Clones variables in src into the current instance of the class, but nothing related to weighting/analyses. (Matching weights are reset to 1 and sample sizes are recalculated.) Any views will be turned into regular Mata variables.

    The last line of code in this example will produce an unweighted balance table:
        : P1 = P2 = psweight()
        : P1.st_set(...)
        : P1.cbps(...)
        : P1.balancetable()
        : P2.clone(P1)
        : P2.balancetable()


====================================================

Conformability

    t                 : n x 1   real
    x                 : n x k   real
    w                 : n x 1   real
    y0                : n0 x l  real, where n0 is the number of rows with t==0
    newweight         : n x 1   real
    colvector         : n x 1   real

    denominator       : 1 x 1   real
    stat              : 1 x 1   string
    subcmd            : 1 x 1   string
    cvopt             : 1 x 3, 1 x 6, 1 x 9, or 1 x 12 real

    src               : scalar, instance of class psweight

    solve()           : 1 x k   real
    ipw()             : 1 x k   real
    cbps()            : 1 x k   real
    cbpsoid()         : 1 x k   real

    get_pscore()      : n x 1  real
    get_weight_mtch() : n x 1  real
    get_weight()      : n x 1  real

    balancetable()    : k x 6   real
    diff()            : 1 x k   real
    stddiff()         : 1 x k   real
    varratio()        : 1 x k   real

    progdiff()        : 1 x l   real
    stdprogdiff()     : 1 x l   real
    pomean()          : 1 x l   real

    mean_sd()         : 1 x 1   real
    mean_asd()        : 1 x 1   real
    max_sd()          : 1 x 1   real

    wgt_cv()          : 1 x 1   real
    wgt_sd()          : 1 x 1   real
    wgt_skewness()    : 1 x 1   real
    wgt_kurtosis()    : 1 x 1   real
    wgt_max()         : 1 x 1   real
