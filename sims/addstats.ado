program define addstats
  args matname eqname coef
  tempname add cell

  mat `add' = (_b[`coef'], _se[`coef'], (_b[`coef'] - _g1), (_b[`coef'] - _g1)^2)
  mat colnames  `add' = impact_est sd_error bias error_sqr

  cap nois {
    test _b[`coef']=0
    mat `cell' = (r(p), (r(p) <= (1-c(clevel)/100)))
    mat colnames `cell' = p_null_0 reject_null_0
    mat `add' = (`add', `cell')
  }

  cap nois {
    test _b[`coef']=_g1
    mat `cell' = (r(p), (r(p) <= (1-c(clevel)/100)))
    mat colnames `cell' = p_null_true reject_null_true
    mat `add' = (`add', `cell')
  }

  cap nois {
    gmatchcall balanceresults()
    mat `cell' = (r(max_asd), /// Maximum absolute standardized diff.
                  r(mean_asd), /// Mean absolute standardized diff.
                  r(wgt_sd), /// S.D. of matching weights:
                  r(wgt_cv), /// C.V. of matching weights:
                  r(wgt_skewness), /// Skewness of matching weights:
                  r(wgt_kurtosis), /// Kurtosis of matching weights:
                  r(wgt_max)) // Maximum matching weight:
    mat colnames `cell' = bal_max_asd bal_mean_asd wgt_sd wgt_cv wgt_skewness wgt_kurtosis wgt_max
    mat `add' = (`add', `cell')
  }

  matrix coleq `add' = `eqname'
  matrix `matname' = (nullmat(`matname'), `add')
  mat list `add'

end
