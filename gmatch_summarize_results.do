assert "`: environment computername'"=="M116" | "`: environment computername'"=="NJ1STATA1"
clear all
cap log close gmatch_test_postestimation
cls
set linesize 180
set varabbrev off
set cformat %9.2f

// Program: gmatch_summarize_results.do
// Project: CPC+ 50139
// Programmer: Richard Chapman, Keith Kranker
// Researcher: Laura Blue, Keith Kranker
// Specs: Lauren Vollmer
// Created: 2/16/18

cd              "G:/Secure_Data/50319_CPC_ Eval/TASKS/SBCG/Pgms/"
adopath ++      "G:/Secure_Data/50319_CPC_ Eval/TASKS/SBCG/Pgms/gmatch_source"
local dtafolder "G:/Secure_Data/50319_CPC_ Eval/TASKS/SBCG/Data"

log using "gmatch_log_files/gmatch_summarize_results.log", replace name(gmatch_test_postestimation)

mata: mata mlib index
which gmatch

cd "`dtafolder'/Output/"

dir "gmatch_*.dta"
dir "`dtafolder'\Input\*.dta"

local fn=0
foreach track in 1 2 {

  foreach drop in 0 15 30 50 {

    if `track'==1 & `drop'==50 continue

    di _n(3) "Begin loop with track `track' and drop=`drop'" _n(2)

    use "`dtafolder'\Input\stataInput_track`track'drop`drop'.dta",clear
    qui include "G:\Secure_Data\50319_CPC_ Eval\TASKS\SBCG\Pgms\variable_lists.do"

    di "Unweighted" _n(2)
    gen byte all=1
    _rmcoll treatment `track`track'_vars' if all [iw = weight_cbps_practicesize], expand logit touse(all)
    local track`track'_vars `r(varlist)'
    gettoken trash track`track'_vars : track`track'_vars
    isid cpc_plus_practice_id
    sort cpc_plus_practice_id
    mata: D = gmatch()
    mata: D.set("treatment", "`track`track'_vars'", "all", "weight_cbps_practicesize")
    mata: D.balanceresults("atet",1)

    matrix thisrow = (r(mean_sd_sq), r(mean_asd), r(max_asd), r(wgt_cv), r(wgt_sd), r(wgt_skewness), r(wgt_kurtosis), r(wgt_max), `track', `drop')
    matrix rownames thisrow = unmatched
    matrix roweq    thisrow = track`track'_drop`drop'
    count  if treatment
    matrix thisrow = (thisrow, r(N))
    count  if !treatment
    matrix thisrow = (thisrow, r(N))
    matrix summary = (nullmat(summary) \ thisrow)

    preserve
    dir     "gmatch_track`track'_drop`drop'_*.dta"
    qui fs  "gmatch_track`track'_drop`drop'_*.dta"
    foreach f in `r(files)' {

      local fabbrev: subinstr local f "gmatch_track`track'_drop`drop'_" ""
      local fabbrev : subinstr local fabbrev "  " " ", all
	  local fabbrev : subinstr local fabbrev "  " " ", all
	  local fabbrev : subinstr local fabbrev "  " " ", all
	  local fabbrev =strtoname("`fabbrev'")

      di _n(2) "Begin inner loop with track `track' and drop=`drop'" _n "File #`++fn': `f'" _n "(`fabbrev' in the table below)" _n(2)

      merge 1:1 cpc_plus_practice_id using "`f'", assert(3) noreport

      tabstat _weight _weight_mtch _pscore, by(treatment) c(s) s(N mean sd min p1 p10 p25 p50 p75 p90 p99 max) format

      // hist _pscore  , kdens by(treatment, col(1) title("Propensity scores" "`file'")) ysize(10) name(yhat_rf`run')
      // hist _weight_mtch if !treatment, kdens name(C_wt_rf`run') title("Comparison group weights" "`file'")

      di _n(2) "Reweighted"
      isid cpc_plus_practice_id
      sort cpc_plus_practice_id
      mata: D.reweight(st_data(.,"_weight_mtch", "all"))
      mata: D.balanceresults("atet",1)

      matrix thisrow = (r(mean_sd_sq), r(mean_asd), r(max_asd), r(wgt_cv), r(wgt_sd), r(wgt_skewness), r(wgt_kurtosis), r(wgt_max), `track', `drop')
      matrix rownames thisrow = `fabbrev'
      matrix roweq    thisrow = track`track'_drop`drop'
      count  if treatment & _weight_mtch
      matrix thisrow = (thisrow, r(N))
      count  if !treatment & _weight_mtch
      matrix thisrow = (thisrow, r(N))
      matrix summary = (summary \ thisrow)
      restore, preserve
	  
	  // di "Stored estimates"
	  // estimates use `"`: subinstr local f ".dta" ".ster"'"'
	  // ereturn list
    }
    restore, not
  }
}

// unsorted
matrix colnames summary = mean_sd_sq mean_asd max_asd wgt_cv wgt_sd wgt_skewness wgt_kurtosis wgt_max track drop N_t N_c 
_matrix_table summary

// sort the matrix
mata: 
	summ = st_matrix("summary")
	rown = st_matrixrowstripe("summary")
	coln = st_matrixcolstripe("summary")
	s1 = order(summ, (9, 10, 2, 3, 4, 5, 6, 7, 1)) 
	s2 = order(summ, (9, 10, 4, 5, 6, 7, 2, 3, 1))
	st_matrix("summary_sort1",summ[s1,1::8])
	st_matrix("summary_sort2",summ[s2,1::8])
	st_matrixrowstripe("summary_sort1",rown[s1,.])
	st_matrixrowstripe("summary_sort2",rown[s2,.])
	st_matrixcolstripe("summary_sort1",coln[1::8,.])
	st_matrixcolstripe("summary_sort2",coln[1::8,.])
end
_matrix_table summary_sort1
_matrix_table summary_sort2

log close gmatch_test_postestimation
