{smcl}
{* $Id$}{...}
{* Copyright (C) Mathematica This code cannot be copied, distributed or used without the express written permission of Mathematica Policy Research, Inc.}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[TE] teffects intro" "mansection TE teffectsintro"}{...}
{vieweralsosee "[TE] teffects ipw" "mansection TE teffectsipw"}{...}
{vieweralsosee "[R] logit" "mansection R logit"}{...}
{viewerjumpto "Title" "psweight##title"}{...}
{viewerjumpto "Syntax" "psweight##syntax"}{...}
{viewerjumpto "Description" "psweight##description"}{...}
{viewerjumpto "Options" "psweight##options"}{...}
{viewerjumpto "Remarks" "psweight##remarks"}{...}
{viewerjumpto "Author" "psweight##author"}{...}
{viewerjumpto "Examples" "psweight##examples"}{...}
{viewerjumpto "Stored results" "psweight##results"}{...}
{viewerjumpto "References" "psweight##references"}{...}
{marker title}{...}
{title:Title}

{p2colset 1 13 15 2}{...}

{p2col:{bf:psweight} {hline 2}}IPW- and CBPS-type propensity score reweighting, with various extentions{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 12 2}
. {cmd:psweight} {it:{help psweight##subcommand:subcommand}}
                      {it:{help varname:tvar}} {it:{help varlist:tmvarlist}}
                      {ifin}
                      [{it:{help psweight##weight:weight}}]
                      [{cmd:,}
                      {it:{help psweight##stat:stat}}
                      {it:{help psweight##penalty:penalty}}
                      {it:{help psweight##variance:variance}}
                      {it:{help psweight##options_table:options}}]

{p 8 12 2}
. {cmd:psweight} {opt balance:only}
                      {it:{help varname:tvar}} {it:{help varlist:tmvarlist}}
                      {ifin}
                      [{it:{help psweight##weight:weight}}]
                      [{cmd:,}
                      {opth mw:eight(varname)}
                      {it:{help psweight##variance:variance}}
                      {it:{help psweight##options_table:options}}]

{p 8 12 2}
. {cmd:psweight} {opt call} {it:classfunction()}


{phang}
{it:tvar} must contain values 0 or 1, representing the treatment (1) and comparison (0) groups.

{phang}
{it:tmvarlist} specifies the variables that predict treatment assignment in
the treatment model.

{phang}
{it:classfunction()} is a function of the psweight Mata class. Functions may take arguments.
{error:  << link to Mata docs >>}


{marker subcommand}{...}
{synoptset 16}{...}
{synopthdr:subcommand}
{synoptline}
{synopt :{opt ipw}}logit regression{p_end}
{synopt :{opt cbps}}just identified covariate-balancing propensity score{p_end}
{synopt :{opt cbpsoid}}over identified covariate-balancing propensity score {p_end}
{synopt :{opt pcbps}}penalized covariate-balancing propensity score{p_end}
{synopt :{opt mean_sd_sq}}minimize mean(stddiff())^2{p_end}
{synopt :{opt sd_sq}}minimize sum(stddiff()^2){p_end}
{synopt :{opt stdprogdiff}}minimize sum(stdprogdiff()^2){p_end}
{synopt :{opt mean_sd}}synonym for {it:mean_sd_sq}{p_end}
{synopt :{opt sd}}synonym for {it:sd_sq}{p_end}
{synoptline}


{marker stat}{...}
{synopthdr:stat}
{synoptline}
{synopt :{opt ate}}estimate average treatment effect in population; the default{p_end}
{synopt :{opt atet}}estimate average treatment effect on the treated{p_end}
{synopt :{opt ateu}}estimate average treatment effect on the untreated{p_end}
{synoptline}


{marker penalty}{...}
{synopthdr:penalty}
{synoptline}
{synopt :{opt cvtarget(# # #)}}applies penalty using the coefficient of variation of the weight distribution; default is no penalty: cvtarget(0 0 2){p_end}
{synopt :{opt skewtarget(# # #)}}applies penalty using the skewness of the weight distribution; default no penalty: skewtarget (0 0 2){p_end}
{synopt :{opt kurttarget(# # #)}}applies penalty using the excess kurtosis of the weight distribution; default no penalty: kurttarget(0 0 2){p_end}
{synoptline}
{p 4 6 2}One or more penalty options may be specified.{p_end}


{marker variance}{...}
{synopthdr:variance}
{synoptline}
{synopt :{opt pool:edvariance}}uses the pooled (treatment plus control) sample's variances to calculate standardized differences; the default{p_end}
{synopt :{opt con:trolvariance}}uses the control groups' variances to calculate standardized differences{p_end}
{synopt :{opt tre:atvariance}}uses the treatment groups' variances to calculate standardized difference{p_end}
{synopt :{opt a:veragevariance}}uses (the control groups' variances + treatment groups' variances)/2 to calculate standardized differences{p_end}
{synoptline}


{marker options_table}{...}
{synopthdr}
{synoptline}
{synopt :{opth dep:vars(varlist)}}outcome variables{p_end}
{synopt :{it:{help psweight##display_options:display_options}}}control
INCLUDE help shortdes-displayoptall
{synopt :{it:{help psweight##maximize_options:maximize_options}}}control
the maximization process; seldom used {* includes from()}{p_end}
INCLUDE help shortdes-coeflegend
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{it:tmvarlist} may contain factor variables; see {help fvvarlists}.{p_end}
{p 4 6 2}
{marker weight}{...}
{opt fweight}s and {opt iweight}s are allowed; see {help weight}.{p_end}
{p 4 6 2}


{marker description}{...}
{title:Description}

{pstd}
{cmd:psweight} computes inverse-probability weighting (IPW) weights for average treatment effect,
average treatment effect on the treated, and average treatment effect estimators for observational data.
IPW estimators use estimated probability weights to correct for the missing data on the potential outcomes).
Probabilties of treatment--propensity scores--are computed for each observation with one of variety of methods, including
logistic regresson (traditional IPW),
covariate balance propensity scores (CBPS),
penalized balance propensity scores (PCBPS)
prognostic score balancing propensity scores, and
other methods.

{pstd}
{cmd:psweight} constructs several variables:{p_end}
{phang2}(1) The propensity scores are returned in {it:_pscore}.{p_end}
{phang2}(2) The treatment indicator is returned in {it:_treated}.}.{p_end}
{phang2}(3) The IPW weights, computed from the propensity scores, are returned in {it:_weight_mtch}.{p_end}
{phang2}(4) The final weights--the product of _weight_mtch and the sample weights--are returned in {it:_weight}.{p_end}
{pstd}If these variables exist before running {cmd:psweight}, they will be replaced.{p_end}

{pstd}
You can also use {cmd:psweight balance} to produce balanced tables to compare the treated and untreated observations
after (or before) reweighting.

{pstd}
This Stata command is a "wrapper" around the psweight Mata class.
{error:  << link to Mata docs >>}
After running the command you can access apply any of the class functions
to your data or results through {cmd:psweight call}.


{marker options}{...}
{title:Options}

{dlgtab:Subcommand}

{pstd}
The {it:subcommand} specifies which method is used to compute coefficients, {it:b}, for the
propensity score model.
In all cases, propensity score for each row in the data is defined as {p_end}

{center:p = {help mf_logit:invlogit}({it:X} {help [M-2] op_arith:*} {it:b}')}

{pstd} where {it:X} is the vector of matching variables ({it:tmvarlist}) for the respective row.
The {it:subcommand} controls how the vector {it:b} is computed.
The seven available estimation methods are:

{pmore}
The {opt ipw} subcommand fits a {help logit:logit regression model} by maximim likelihood.

{pmore2}
The logit regression solves for {it:b} in the model:

{center:Prob({it:tvar} = 1 | {it:X}) = {help mf_logit:invlogit}({it:X} {help [M-2] op_arith:*} {it:b})}

{pmore}
The {opt cbps} subcommand computes just identified covariate-balancing propensity scores
(Imai and Ratkovic 2014).

{pmore}
The {opt cbpsoid} subcommand computes over identified covariate-balancing propensity scores .
(Imai and Ratkovic 2014).

{pmore}
The {opt pcbps} subcommand impliments penalized covariate-balancing propensity score (Kranker, Blue, and Vollmer Forrow 2019).
The {opt pcbps} subcommand requires at least one
{it:{help psweight##penalty:penalty}} option be specified.
This is a synonym for the {it:cbps} subcommand when one or
more of the {it:{help psweight##penalty:penalty}} option is included.

{pmore}
The {opt mean_sd_sq} or {opt mean_sd} subcommands find the model coefficients that minimize the quantity:
mean(stddiff())^2.
{error:  << link to Mata docs >>}
That is, the weights mimizize the mean standardized difference, squared.

{pmore}
The {opt sd_sq} or {opt sd} subcommands find the model coefficients that minimize the quantity:
sum(stddiff()^2).
{error:  << link to Mata docs >>}
That is, the weights mimizize the squared (standardized) differences
in means of {it:tmvarlist} between the treatment and control groups.
If you have more than one variable in {it:tmvarlist},
the squared (standardized) differences in prognositic scores are summed.

{pmore}
The {opt stdprogdiff} subcommand finds the model coefficients that minimize the quantity:
sum(stdprogdiff()^2).
{error:  << link to Mata docs >>}
That is, the weights mimizize the squared (standardized) differences
in mean prognositic scores between the treatment and control groups.
If you have more than one outcome variable,
the squared (standardized) differences in prognositic scores are summed.

{pmore2}
The {it:stdprogdiff} subcommand requires that dependent
variables be specified (through the {opt depvar}() option).

{pmore2}
Prognoistic scores are computed by fitting a linear regression model
of the {it:tmvarlist} on the dependent variable(s) by ordinary least squares using only
control group observationsy, and then computing predicted values (for the whole sample).
This method of computing prognositic scores follows Hansen (2008) and Stuart, Lee, and Leacy (2013).


{dlgtab:Stat}

{pstd}
{it:stat} is one of three statistics: {opt ate}, {opt atet}, or {opt ateu}.
{opt ate} is the default.
The {it:stat} dictates how the command uses propensity scores ({it:p}) to compute
IPW "matching weights" (the variable named {it:_weight_mtch}).

{pmore}
{opt ate} specifies that the average treatment effect be estimated.

{pmore2}
The unnormalized {opt ate} weights are
1/{it:p} for treatment group observations and
1/(1-{it:p}) for control group observations.
The weights are then normalzied to have mean=1 in each group.

{pmore}
{opt atet} specifies that the average treatment effect on the treated be estimated.

{pmore2}
The unnormalized {opt atet} weights are
1 for treatment group observations and
{it:p}/(1-{it:p}) for control group observations.
The weights are then normalzied to have mean=1 in the control group.

{pmore}
{opt ateu} specifies that the average treatment effect on the untreated be estimated.

{pmore2}
The unnormalized {opt ateu} weights are
(1-{it:p})/{it:p}  for treatment group observations and
1 for control group observations.
The weights are then normalzied to have mean=1 in the treatment group.

{pstd}
Once the matching weights ({it:_weight_mtch}) are available,
the final weights ({it:_weight}) are set equalt to:{p_end}
{center:{it:_weight} = {it:W} {help [M-2] op_colon::*} {it:_weight_mtch}}
{pstd} where {it:W} are the {it:{help psweight##weight:sample weights}}.
The variable {it:_weight} equals {it:_weight_mtch} if no sample weights are provided.


{dlgtab:Penalty}

{pstd}
The {it:penalty} options modify the "loss function"
used to solve for the propensity score model coefficents ({it:b}).
Let the nonpenalized method for a given {it:subcommand} be written as: {p_end}
{center:{it:b} = argmin {it:L(X,T,W)}}
{pstd} where {it:L(X,T,W)} is the loss function that corresponds to the specified {it:subcommand}
(e.g., logit regression or CBPS),
given the data ({it:(X,T)} and vector of weights {it:W}.
(The weights are computed using the propensity scores,and the propensity scores
are calcullated usign the ceofficents and data using the formulas given above.)
Then the penalized method would be written as:{p_end}
{center:{it:b} = argmin {it:L(X,T,W)} + {it:f(W)}}
{pstd} where {it:f(W)} is smooth, flexible function that increases as the vector of observation weights (W) becomes more variable. Specifically:

{phang2}{opt cvtarget(# # #)} applies a penalty using the coefficient of variation of the weight distribution.
If {opt cvopt(a, b, c)} is specified, then the loss function is modified as:{p_end}
{center:{it:L'(X,T,W) = L(X,T,W) + a * abs((wgt_cv() - b)^c)}}
{error:  << link to Mata docs >>}
{phang3}The default is no penalty: cvtarget(0 0 2).

{phang2}{opt skewtarget(# # #)} applies a penalty using the skewness of the weight distribution.
If {opt skewtarget(d, e, f)} is specified, then the loss function is modified as:{p_end}
{center:{it:L'(X,T,W) = L(X,T,W) + e * abs((wgt_skewness() - e)^f)}}
{error:  << link to Mata docs >>}
{phang3}The default is no penalty: skewtarget(0 0 2).{p_end}

{phang2}{opt kurttarget(# # #)} applies a penalty using the excess kurtosis of the weight distribution.
If {opt kurttarget(g, h, i)} is specified, then the loss function is modified as:{p_end}
{center:{it:L'(X,T,W) = L(X,T,W) + g* abs((wgt_kurtosis() - h)^i)}}
{error:  << link to Mata docs >>}
{phang3}The default is no penalty: kurttarget(0 0 2).{p_end}

{* maxtarget option is undocumented}{...}


{dlgtab:Variance}

{pstd}
{it:variance} is one of three statistics: {opt pooledvariance}, {opt controlvariance}, {opt treatvariance}, or {opt averagevariance}.
{opt pooledvariance} is the default.
The {it:variance} dictates how the command standardizes the difference in means between the treatment and control groups.

{phang2}
{opt pooledvariance} uses the pooled (treatment plus control) sample's variances to calculate standardized differences; the default

{phang2}
{opt controlvariance} uses the control groups' variances to calculate standardized differences

{phang2}
{opt treatvariance} uses the treatment groups' variances to calculate standardized difference

{phang2}
{opt averagevariance} uses (the control groups' variances + treatment groups' variances)/2 to calculate standardized differences (as in {help tebalance})


{dlgtab:Other options}

{phang}
{opth depvars(varlist)} specificies dependent variables (outcome variables).
The data for the treatment group observations are ignored, but the data for the
the control group are used to compute prognositic scores (see above).

INCLUDE help displayopts_list

{phang2}
The following options control the display of the balance tables:
{opt formats(%fmt)},
{opt noomit:ted},
{opt vsquish},
{opt noempty:cells},
{opt base:levels},
{opt allbase:levels},
{opt nofvlab:el},
{opt fvwrap(#)},
{opt fvwrapon(style)}, and
{opt nolstretch};
see {help _matrix_table##display_options: _matrix_table}.

{marker maximize_options}{...}
{phang}
{it:maximize_options} control the maximization process:
{opt from(init_specs)},
{opt tr:ace},
{opt grad:ient},
{opt hess:ian},
{cmd:showstep},
{opt tech:nique(algorithm_specs)},
{opt iter:ate(#)},
{opt tol:erance(#)},
{opt ltol:erance(#)},
{opt nrtol:erance(#)},
{opt qtol:erance(#)},
{opt nonrtol:erance},
{opt showtol:erance}, and
{cmdab:dif:ficult}.
For a description of these options, see {manhelp maximize R} and
{help mf_moptimize_init_mlopts:moptimize_init_mlopts()}.

{phang}
{cmd:coeflegend};
{helpb estimation options##coeflegend:[R] estimation options}.

{* [opt noconstant]; see [helpb estimation options:[R] estimation options].}{...}

{dlgtab:psweight balance}

{pstd}
{cmd:psweight} {cmd:balanceonly} constructs a balance table instead of computing IPW weights.
If {opth mweight(varname)} is not specified, the balance table is constructed with "unweighted" data.
(Only the {it:{help psweight##weight:sample weights}} are applied.)
If the option is specified, then teh reweighted data are used to construct the balance table.

{phang2}
{opth mweight(varname)} is a variable containing user-specified "matching" weights.
This variable is the analogue to the _weight_mtch variable;
it should not be multiplied with the {it:{help psweight##weight:sample weights}}.

{pstd}
As explained in the {help psweight##syntax:syntax}, {cmd:psweight} {cmd:balanceonly} also allows many of the options listed above.


{dlgtab:pweight call}

{pstd}
The {cmd:psweight} {it:subcommand} program is a "wrapper" around the psweight Mata {mansection M-2 class:class}.
{error:  << link to Mata docs >>}
When the command is run, it constructs an instance of the class, and this instance
remains accessible afterward the command has finished.
Your can access any of the class functions through {cmd:psweight call}.

{pmore}
Note that any default options that were overridden when {cmd:psweight} {it:subcommand} was called
will continue to be applied with {cmd:psweight call}.

{pmore}For exmaple, supposed you wrote:{p_end}
{phang2}{cmd:. psweight ipw mbsmoke mmarried mage fbaby medu, treatvariance}{p_end}
{phang2}{cmd:. psweight call balanceresults()}{p_end}
{pmore}Then the balanced table will use the treatment group's variance to calculate standardized differences.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse cattaneo2}{p_end}

{pstd}
Balance before reweighting{p_end}
{phang2}{cmd:. psweight balanceonly mbsmoke mmarried mage fbaby medu}{p_end}

{pstd}
Estimate the average treatment effect of smoking on birthweight, using a
logit model to predict treatment status{p_end}
{phang2}{cmd:. psweight ipw mbsmoke mmarried mage fbaby medu}{p_end}
{phang2}{cmd:. psweight call balanceresults()}{p_end}

{pstd}Estimate the average treatment effect on the treated with CBPS{p_end}
{phang2}{cmd:. psweight cbps mbsmoke mmarried mage fbaby medu, atet}{p_end}
{phang2}{cmd:. psweight call balanceresults()}{p_end}

{pstd}Estimate the average treatment effect on the treated with Penalized CBPS{p_end}
{phang2}{cmd:. psweight pcbps mbsmoke mmarried mage fbaby medu, atet cvtarget(1 .5 6)}{p_end}
{phang2}{cmd:. psweight call balanceresults()}{p_end}

{phang}For more examples, see psweight_example.do{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:psweight} stores the following in {cmd:e()}:

{phang}lorem ipsum


{marker author}{...}
{title:Author}

{pstd}By Keith Kranker{break}
Mathematica Policy Research{p_end}

{pstd}This help file last updated $Date${p_end}

{pstd}My coauthors, Laura Blue and Lauren Vollmer Forrow, were closely involved with the
developement of the Penalized CBPS methodology,
and we also received many helpful suggestions from our colleages at Mathematica.
I thank Liz Potamites for testing early versions of the program and providing helpful feedback.{p_end}


{marker references}{...}
{title:References}

{psee}
Fong, C., M. Ratkovic, K. Imai, C. Hazlett, X. Yang, and S. Peng.
2018.
CBPS: Covariate Balancing Propensity Score,
Package for the R programming langauage,
{it:The Comprehensive R Archive Network}.
Available at: https://CRAN.R-project.org/package=CBPS

{psee}
Hansen, B. B.
2008.
"The Prognostic Analogue of the Propensity Score."
{it:Biometrika},
95(2): 481–488, doi:10.1093/biomet/asn004.

{psee}
Imai, K. and M. Ratkovic.
2014.
"Covariate Balancing Propensity Score."
{it:Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
76(1): 243–263, doi:10.1111/rssb.12027.

{psee}
Kranker, K., L. Blue, and L. Vollmer Forrow.
2019.
"Improving Effect Estimates by Limiting the Variability in Inverse Propensity Score Weights."
Manuscript under review.

{psee}
Stuart, E. A., B. K. Lee, and F. P. Leacy.
2013.
"Prognostic Score–based Balance Measures Can Be a Useful Diagnostic for Propensity Score Methods in Comparative Effectiveness Research."
{it:Journal of Clinical Epidemiology},
66(8): S84–S90.e1, doi:10.1016/j.jclinepi.2013.01.013.
