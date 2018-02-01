 void gmatch::cbpseval2( real    scalar        todo,
                        real    rowvector     beta,
                        string  scalar        est,
                        string  scalar        fctn,
                        real   scalar         lnf,
                        real    rowvector     g,
                        real    matrix        H)
{

}  


 
    // Designate sample size, number of treated and control observations,
    // theta, which are used to generate probabilities.
    // Trim probabilities, and generate weights.
    theta = (this.X*beta)
    probs = (1+exp(-theta))^-1
    probs = rowmax(J(rows(this.X),1,1e-6),rowmin(1:-J(rows(this.X),1,1e-6),probs))	
    if(ATT) {
      w = ATT.wt.func(beta)
    }
    else {
      w = (probs-1+treat)^-1
    }
    
    // Generate the vector of mean imbalance by weights.
    w.del = 1/(n)*t(this.W*X) * (w)
    w.del = (w.del)
    w = (w)
    
    // Generate g-bar, as in the paper.
    gbar = c( 1/n*t(this.W*X) * (treat-probs),w.del)
    
    // Generate the covariance matrix used in the GMM estimate.
    // Was for the initial version that calculates the analytic variances.
    if(is.null(invV))
    {
      if(ATT){
        X.1   = this.W^.5*X*((1-probs)*probs)^.5
        X.2   = this.W^.5*X*(probs/(1-probs))^.5
        X.1.1 = this.W^.5*X*(probs)^.5
      }
      else{
        X.1   = this.W^.5*X*((1-probs)*probs)^.5
        X.2   = this.W^.5*X*(probs*(1-probs))^-.5		
        X.1.1 =  this.W^.5*X
      }
      if (ATT){
        V = rbind(1/n*cbind(t(X.1) * X.1,t(X.1.1) * X.1.1)*n/sum(treat),
                 1/n*cbind(t(X.1.1) * X.1.1*n/sum(treat),t(X.2) * X.2*n^2/sum(treat)^2))
      }
      else{
        V = rbind(1/n*cbind(t(X.1) * X.1,t(X.1.1) * X.1.1),
                 1/n*cbind(t(X.1.1) * X.1.1,t(X.2) * X.2))
      }
      invV = ginv(V)
    }			   
    
    // Calculate the GMM loss.
    loss1 = (t(gbar) * invV * (gbar))		
    out1 = list("loss"=loss1, "invV"=invV)
    out1
  }