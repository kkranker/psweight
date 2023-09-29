*! $Id$
*! One replication of the Simulations with DGP_KSIR.ADO
// Here the estimatand is the treament group rather than the ATE (as in onrep_ksir.ado)
// User specifies
//   - one or more sample sizes of dataset(s)
//   - models to run (on each dataset)
// The command then loops over DGPs, Ns, and Impacts.
// In each loop it
//   - generates the data (calls dgp_ssbgc.ado)
//   - estimates one or more models
//   - computes various statistics
// All the results are returned in one wide matrix (e(b))
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.

program define onerep_ksir_means, eclass
  version 15.1

  local cmd `"onerep_ksir_means `0'"'
  syntax,  ///
         n(numlist min=1 >0 integer sort) /// loop over one or more sample sizes
         ESTimators(namelist) /// pick wich estimators to include (IPW, CBPS, etc.) see below
       [ CVTargets(numlist >=10 <=100) /// passed to gmatch
         /// ate atet ateu /// atet is the default; passed to reweighting commands (teffects, gmatch, etc)
         TRUEpscore /// use Z, rather than X in p-score model
         TRUEoutcome /// use Z, rather than X in outcome model
         /// AUGmented /// run OLS models to estimate impacts
         HISTogram irversion /// passed to DGP
         /// vce(passthru) /// e.g., add robust standard errors in outcome models
         elasticopts(string) rfopts(string) /// passed to elasticregress and randomforest, respectively
         ITERate(integer 0) /// number of itations
         QUIETly /// supress a lot of the output
         *] // display options passed everywhere
  _get_diopts diopts options, `options'

  // process options
  // throw error if invalid options
  local valid_est raw ipw_true_ps ipw ipw_te stdprogdiff cbps ipwcbps elastic rf
  // if ("`estimators'"=="") local estimators : copy local valid_est
  if !`:list estimators in valid_est' {
    di as error "estimators(`: list estimators-valid_est') invalid"
    error 198
  }
  if ("`ate'`atet'`ateu'"=="") {
    local atet atet
  }
  if ("`truepscore'"=="truepscore") {
    local pscorevarlist "z1-z4"
  }
  else {
    local pscorevarlist "x1-x4"
  }
  if ("`trueoutcome'"=="trueoutcome") {
    local omvarlist "z1-z4"
    local augmented "augmented" // trueoutcome implies augmented
    local aug "aug"
  }
  else if ("`augmented'"=="augmented") {
    local omvarlist "x1-x4"
    local aug "aug"
  }
  if ("`aug'"=="aug") local omopt omvarlist(\`omvarlist')
  if (`iterate'!=0) {
    set maxiter `iterate'
    local options `options' iterate(`iterate')
  }
  if ("`quietly'"=="") local quietly noisily
  set matsize `c(max_matsize)'

  tempname _b_ add from
  local c = 0
  local scenario I
      local L : word count `n'
      dgp_ksir, n(`: word `L' of `n'') `histogram' `irversion'

      scalar _g1 = 210 // per K-S (page 529) and I-R (page 251), this is the quantity of interest to estimate

      local l = `L'
      while (`l' > 0) {
        local thisN : word `l' of `n'
        local prefix = strtoname(trim("`scenario'_`++c'"))

        keep in 1/`thisN'
        di _n(3) as txt _dup(20) "-" _n ///
                 as txt "`prefix'" _n ///
                 as txt _dup(20) "-" _n ///
                 as txt "Scenario: " as res "`scenario'" _n ///
                 as txt "True mu:  " as res _g1  _n ///
                 as txt "N:        " as res "`thisN'" _n ///
                 as txt _dup(20) "-"  _n(2)

        qui count if a
        matrix          `add' = (`thisN', r(N), _g1)
        matrix colnames `add' = N Nt true
        matrix coleq    `add' = `prefix'
        matrix          `_b_' = (nullmat(`_b_'), `add')
        ereturn clear
        return  clear
        cap mata: mata drop gmatch_ado_most_recent

        // Difference in means ("raw")
        local e "raw"
        if (`: list e in estimators') cap `quietly' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `truepscore' `trueoutcome'" _n(2)
          summ a, mean
          gen double rawps = r(mean)
          ipw_pop rawps, `diopts'
          addstats `_b_' 1.a `prefix'_`e'
        }
        cap mata: mata drop gmatch_ado_most_recent

        // Difference in means, weighted using true propensity scores ("true")
        local e "ipw_true_ps"
        if (`: list e in estimators') cap `quietly' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `truepscore' `trueoutcome'" _n(2)
          ipw_pop ps, `diopts'
          addstats `_b_' 1.a `prefix'_`e'
        }
        cap mata: mata drop gmatch_ado_most_recent

        // IPW model (with teffects)
        local e "ipw_te"
        if (`: list e in estimators') cap `quietly' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `truepscore' `trueoutcome'" _n(2)
          teffects ipw (y) (a `pscorevarlist'), `ate'`atet'`ateu' aeq `diopts'
          matrix `from' = e(b)
          matrix `from' = `from'[1, "TME1:"]
          local fromopt from(`from')
          tempvar  ipw_te_ps
          predict `ipw_te_ps', ps tlevel(1)
          ipw_pop `ipw_te_ps', `diopts'
          addstats `_b_' 1.a `prefix'_`e'
        }

        // IPW model with elastic net used to estimate P-scores
        local e "elastic"
        if (`: list e in estimators') cap `quietly' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `truepscore' `trueoutcome'" _n(2)
          elasticregress a c.(`pscorevarlist')##c.(`pscorevarlist'), `elasticopts'
          tempvar elasticW elasticPS
          predict `elasticPS'
          ipw_pop `elasticPS', `diopts'
          addstats `_b_' 1.a `prefix'_`e'
        }
        cap mata: mata drop gmatch_ado_most_recent

        // IPW model with randomforest model used to estimate P-scores
        local e "rf"
        if (`: list e in estimators') cap `quietly' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `truepscore' `trueoutcome'" _n(2)
          randomforest a `pscorevarlist', type(class) `rfopts'
          tempvar rfW rfPS rfPS0
          predict `rfPS0' `rfPS', pr
          ipw_pop `rfPS', `diopts'
          addstats `_b_' 1.a `prefix'_`e'
        }
        cap mata: mata drop gmatch_ado_most_recent

        // IPW model (with gmatch.ado)
        local e "ipw"
        if (`: list e in estimators') cap `quietly' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `truepscore' `trueoutcome'" _n(2)
          gmatch a `pscorevarlist', ipw `ate'`atet'`ateu' `fromopt' `options' `diopts'
          matrix `from' = e(b)
          local fromopt from(`from')
          ipw_pop _ps, `diopts'
          addstats `_b_' 1.a `prefix'_`e'
        }
        cap mata: mata drop gmatch_ado_most_recent

        // Minimize difference in prognostic scores model (with gmatch.ado)
        local e "stdprogdiff"
        if (`: list e in estimators') cap `quietly' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `truepscore' `trueoutcome'" _n(2)
          gmatch a `pscorevarlist', stdprogdiff depvar(y) `ate'`atet'`ateu' `fromopt' `options' `diopts'
          ipw_pop _ps, `diopts'
          addstats `_b_' 1.a `prefix'_`e'
        }
        cap mata: mata drop gmatch_ado_most_recent

        // CBPS overidentified model (with gmatch.ado)
        local e "ipwcbps"
        if (`: list e in estimators') cap `quietly' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `truepscore' `trueoutcome'" _n(2)
          gmatch a `pscorevarlist', cbps ipw `ate'`atet'`ateu' `fromopt' `options' `diopts'
          ipw_pop _ps, `diopts'
          addstats `_b_' 1.a `prefix'_`e'
        }
        cap mata: mata drop gmatch_ado_most_recent

        // CBPS model (with gmatch.ado)
        local e "cbps"
        if (`: list e in estimators') cap `quietly' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `truepscore' `trueoutcome'" _n(2)
          gmatch a `pscorevarlist', cbps `ate'`atet'`ateu' `fromopt' `options' `diopts'
          matrix `from' = e(b)
          local fromopt from(`from')
          gmatchcall wgt_cv("`ate'`atet'`ateu'")
          local cvcbps = r(wgt_cv)
          ipw_pop _ps, `diopts'
          addstats `_b_' 1.a `prefix'_`e'
        }
        cap mata: mata drop gmatch_ado_most_recent

        // CBPS model (with gmatch.ado), with CV at X%
        foreach cut of local cvtargets {
          // for some weird reason cap `quietly' doesn't work" see https://www.statalist.org/forums/forum/general-stata-discussion/general/1482533-macro-expansion-in-a-while-loop
          cap nois {
          local e = strtoname("cbps`cut'")
          local confirm cbps
          if (!`: list confirm in estimators') continue
          if (`cut'==100) continue
          local cvtarget = (`cvcbps'*`cut'/100)
          `quietly' di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `truepscore' `trueoutcome' `augmented'" _n ///
                   as txt "CV target: cvtarget(20 " as res %7.4f `cvtarget' as txt " 6)" _n(2)
          `quietly' gmatch a `pscorevarlist', cbps cvtarget(20 `cvtarget' 6) `ate'`atet'`ateu' `fromopt' `options' `diopts'
          ipw_pop _ps, `diopts'
          addstats `_b_' 1.a `prefix'_`e'
          }
          cap mata: mata drop gmatch_ado_most_recent
        }

        local --l

      }

  ereturn clear
  cap mata: mata drop gmatch_ado_most_recent
  // mat list `_b_'
  ereturn post `_b_'
  ereturn local cmd "onerep_ksir_means"
  ereturn local cmdline "`cmd'"
  ereturn local title "Simulation results"
  ereturn display, `diopts'
end

// Imai & Ratkovic's IPW estimator (page 251 of the article)
program define ipw_pop, eclass
  ereturn clear
  syntax varname [,*] // p-score and display options
  tempvar numer denom
  tempname b

  egen `numer' = total(y * a / `varlist')
  egen `denom' = total(    a / `varlist')
  matrix `b'   = (`numer'[1] / `denom'[1])
  matrix colnames `b' = 1.a
  ereturn post `b'
  ereturn display, `options'
end

// helper program to add impact esimate (named "coef" to matrix
program define addstats
  version 15.1
  args matname coef eqname
  tempname add

  // impact estimate, standard error, bias, MSE
  cap nois {
    mataddscalar `add' impact_est = _b[`coef']              // beta
    /*
    mataddscalar `add' sd_error   = _se[`coef']             // standard error
    */
    mataddscalar `add' bias       =  _b[`coef'] - _g1       // impact estimate - true
    mataddscalar `add' error_sqr  = (_b[`coef'] - _g1)^2    // (impact estimate - true)^2

    /*
    // power to detect true effect given SE
    cap nois {
      mata: power_zstat(st_matrix("`add'")[2], st_numscalar("_g1"), .05, 2)
      mataddscalar `add' power_zstat_0 = r(beta)
    }
    */
  }

  /*
  // Ho: null of effect = 0
  // lincom is just a conventient way to get p-value and CIs.
  cap nois {
    return clear
    qui lincom _b[`coef']
    mataddscalar `add' p_0      = r(p)                              // grab p-value
    mataddscalar `add' reject_0 = (r(p) <= (1 - c(clevel) / 100))   // reject null?
    mataddscalar `add' covered  = ((r(lb) <= _g1) & (_g1 <= r(ub))) // coverage
  }

  cap nois {
    gmatchcall balanceresults()
    mataddscalar `add' bal_max_asd  = r(max_asd)       // Maximum absolute standardized diff.
    mataddscalar `add' bal_mean_asd = r(mean_asd)      // Mean absolute standardized diff.
    mataddscalar `add' wgt_sd       = r(wgt_sd)        // S.D. of matching weights
    mataddscalar `add' wgt_cv       = r(wgt_cv)        // C.V. of matching weights
    mataddscalar `add' wgt_skewness = r(wgt_skewness)  // Skewness of matching weights
    mataddscalar `add' wgt_kurtosis = r(wgt_kurtosis)  // Kurtosis of matching weights
    mataddscalar `add' wgt_max      = r(wgt_max)       // Maximum matching weight
  }
  */
  matrix coleq `add' = `eqname'
  matrix `matname' = (nullmat(`matname'), `add')
  return clear
  ereturn clear

end

program define mataddscalar
  syntax namelist(min=1 max=2) =/exp
  if (missing((`exp'))) exit
  gettoken matrix colname: namelist
  tempname cell
  mat `cell' = (`exp')
  if ("`colname'"!="") mat colnames `cell' = `colname'
  mat `matrix' = (nullmat(`matrix'), `cell')
end

// converts p-scores into IPW weights
program define mkwgt

  syntax newvarname =/exp [, ate atet ateu trim]

  tempvar ps
  qui gen double `ps' = `exp'

  if ("`trim'"=="trim") {
    di "Trimming p-scores <.0001 or >.9999 when contructing weights"
    cap nois replace `ps' = .0001 if `ps' < .0001
    cap nois replace `ps' = .9999 if `ps' > .9999 & !missing(`varlist')
  }
  if inlist("`ate'`atet'`ateu'", "ate", "") {
    gen `varlist' = cond(a, 1/`ps', 1/(1-`ps'))
    summ `varlist' if a, mean
    qui replace `varlist' = `varlist' / r(mean) if a
    summ `varlist' if !a, mean
    qui replace `varlist' = `varlist' / r(mean) if !a
  }
  else if "`ate'`atet'`ateu'"=="atet" {
    gen `varlist' = cond(a, 1, `ps'/(1-`ps'))
    summ `varlist' if !a, mean
    qui replace `varlist' = `varlist' / r(mean) if !a
  }
  else if "`ate'`atet'`ateu'"=="ateu" {
    gen `varlist' = cond(a, (1-`ps')/`ps', 1)
    summ `varlist' if a, mean
    qui replace `varlist' = `varlist' / r(mean) if a
  }
  else {
    di as error "`ate' `atet' `ateu' invalid"
    error 198
  }

end

// calculates power from a SE and effect size, assuming normal distribution
mata:
  mata set matastrict on
  void power_zstat(real scalar se, real scalar effect, real scalar alpha , real scalar sides) {
    if (alpha<=0 | alpha>=1)  _error("The third argument (alpha) must be between 0 and 1")
    if (sides!=2 & sides!=1)  _error("The fourth argument (sides) must be 1 or 2")
    numeric scalar beta
    beta = normal(abs(effect/se) + invnormal(alpha/sides))
    st_numscalar("r(beta)", beta)
    strofreal(round(beta#100,.01)) + " percent power to detect an effect of " + strofreal(effect) + ", assuming a standard error of " + strofreal(se) + " and a " + strofreal(sides) + "-sided test and alpha = " + strofreal(alpha)
  }
end