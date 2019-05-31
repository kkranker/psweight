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
. {cmd:psweight} {it: {help psweight##subcommand:subcommand}}
                      {it:{help varname:tvar}} {it:{help varlist:tmvarlist}}
                      {ifin}
                      [{it:{help psweight##weight:weight}}]
                      [{cmd:,}
                      {it:{help psweight##stat:stat}}
                      {it:{help psweight##variance:variance}}
                      {it:{help psweight##penalty:penalty}}
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


{marker subcommand}{...}
{synoptset 16}{...}
{synopthdr:subcommand}
{synoptline}
{synopt :{opt balanceonly}}computes a balance test (without computing weights){p_end}
{synopt :{opt ipw}}computes weights with loss function: logit regression{p_end}
{synopt :{opt cbps}}computes weights with loss function: CBPS (just identified){p_end}
{synopt :{opt cbpsoid}}computes weights with loss function: CBPS (over identified){p_end}
{synopt :{opt pcbps}}computes weights with loss function: penalized CBPS{p_end}
{synopt :{opt mean_sd_sq}}computes weights with loss function: mean(stddiff())^2{p_end}
{synopt :{opt sd_sq}}computes weights with loss function: sum(stddiff()^2){p_end}
{synopt :{opt stdprogdiff}}computes weights with loss function: sum(stdprogdiff()^2){p_end}
{synopt :{opt mean_sd}}synonym for mean_sd_sq{p_end}
{synopt :{opt sd}}synonym for sd_sq{p_end}
{synoptline}


{marker stat}{...}
{synopthdr:stat}
{synoptline}
{synopt :{opt ate}}estimate average treatment effect in population; the default{p_end}
{synopt :{opt atet}}estimate average treatment effect on the treated{p_end}
{synopt :{opt ateu}}estimate average treatment effect on the untreated{p_end}
{synoptline}


{marker variance}{...}
{synopthdr:variance}
{synoptline}
{synopt :{opt pool:edvariance}}uses the pooled (treatment plus control) sample's variances to calculate standardized differences; the default{p_end}
{synopt :{opt con:trolvariance}}uses the control groups' variances to calculate standardized differences{p_end}
{synopt :{opt tre:atvariance}}uses the treatment groups' variances to calculate standardized difference{p_end}
{synopt :{opt a:veragevariance}}uses (the control groups' variances + treatment groups' variances)/2 to calculate standardized differences (as in {help tebalance}){p_end}
{synoptline}


{marker penalty}{...}
{synopthdr:penalty}
{synoptline}
{synopt :{opt cvtarget(# # #)}}applies penalty using the coefficient of variation of the weight distribution; default is no penalty: cvtarget(0 0 2){p_end}
{synopt :{opt skewtarget(# # #)}}applies penalty using the skewness of the weight distribution; default no penalty: skewtarget (0 0 2){p_end}
{synopt :{opt kurttarget(# # #)}}applies penalty using the excess kurtosis of the weight distribution; default no penalty: kurttarget(0 0 2){p_end}
{synopt :{opt maxtarget(# # #)}}applies penalty using hte maximum in the weight distribution{p_end}
{synoptline}
{p 4 6 2}One or more penalty options may be specified.{p_end}


{marker options_table}{...}
{synopthdr}
{synoptline}
{synopt :{opth dep:vars(varlist)}} outcome variables{p_end}
{synopt :{opth mw:eight(varname)}} variable containing weights for observations;
only used for the {opt balanceonly} {help psweight##subcommand:subcommand}.{p_end}
{synopt :{it:{help psweight##display_options:display_options}}}control
INCLUDE help shortdes-displayoptall
{* incliudes options for _matrix_table: formats(passthru) NOOMITted vsquish NOEMPTYcells BASElevels ALLBASElevels NOFVLABel fvwrap(passthru) fvwrapon(passthru) nolstretch *}{...}
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
Probabilties of treatment--propensity scores--are computed for each observation with one of variety of methods, including{p_end}
{phang}- logistic regresson (traditional IPW),{p_end}
{phang}- covariate balance propensity scores (CBPS),{p_end}
{phang}- penalized balance propensity scores (PCBPS){p_end}
{phang}- prognostic score balancing propensity scores, and{p_end}
{phang}- other methods{p_end}

{pstd}
You can also use {cmd:psweight} to produce balanced tables to compare the treated and untreated observations
after (or before) reweighting.

{pstd}
This Stata command is a "wrapper" around the psweight Mata class.
After running the command you can access apply any of the class functions
to your data or results through {cmd:psweight call}.


{marker options}{...}
{title:Options}

{phang}lorem ipsum


{marker remarks}{...}
{title:Remarks}

{phang}lorem ipsum


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
and we received many helpful suggestions from our colleages at Mathematica.
I thank Liz Potamites for testing early versions of the program and
providing helpful feedback.{p_end}


{marker references}{...}
{title:References}

{psee}CBPS paper{p_end}
{psee}R CRAN package{p_end}
{psee}what else? {p_end}
