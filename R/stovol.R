#' @title stovol
#'
#' @description Fits two types of stochasitc volatility models
#'
#' @param data A data set object, which is one dimension time series
#' @param burn_in The sample size of the burn in before estimating the parameters
#' @param N Sample size of the totally sampled time series
#' @param model_type Specify the type of the model that will fitted.(BSV, ASV)
#'
#' @return (1) The estimated parameters and their confidence intervals,
#' (2) the sampled time series after a burn_in period for parameter estimation,
#' (3) the probability integral transform of the fitted model based on the training data,
#' (4) the standardized residuals , and (5) the AIC an BIC.
#'
#' @import nimble
#' @import HDInterval
#' @examples
#' stovol(data, burn_in, N, model_type)
#' @export



stovol = function(data, burn_in, N, model_type){
  y = data
  if (sum(is.na(y)) != 0){
    return("There are some missing values in the data set!")
  }


  if (!((is.numeric(burn_in)) & (burn_in%%1 ==0)& (burn_in>0))){
    return("burn_in has to be a positive integer!")
  }


  if (!((is.numeric(N)) & (N%%1 ==0)& (N>0))){
    return("The sample size of MCMC time sereis has to be a positive integer!")
  }


  if(burn_in >= N){
    return ("Burin_in can not be larger than the simulated sample size")
  }

  if (!(toupper(model_type) %in% c("BSV", "ASV"))){
    return("model_type has to be of of the values c('BSV', 'ASV')!")
  }

  if (toupper(model_type) %in% c("BSV")){

    T = length(y)


    mu = -10
    phi = 0.7
    sigma1 = 0.19

    mu_est = 0
    phi_est = 0
    sigma1_est = 0
    likelihood_est = 0
    z_est = rep(0, T)


    n = burn_in
    N = N
    z = rep(0, T) # the latent state vector

    est = matrix(rep(0, (N-n)*3), nrow = N-n)


    # initialized the states before the estimation process
    temp = rnorm(1, 0, 1)
    z[1] = mu + sigma1/sqrt(1.0 - phi*phi)*temp

    for (i in 2:T){
      temp = rnorm(1, 0, 1)
      z[i] = mu + phi* (z[i-1]-mu)+ sigma1*temp
    }

    sigma1 =sigma1^2

    for (k in 1:N){
      #print(k)
      u1 = runif(1,0,1)
      a1 = u1*exp(-z[1]/2.0)
      a1_right = -2.0*log(a1)

      u2 = runif(1,0,1)
      b1 = u2*exp(-y[1]*y[1]/2.0*exp(-z[1] ) )

      if (y[1] != 0.0){
        # b1_left = -log(-2.0/y[1]/y[1]*log(b1))
        b1_left = -log(-2.0/y[1]/y[1]*( log(u2) -y[1]*y[1]/2.0*exp(-z[1] )  ))
      }

      u3 = runif(1,0,1)
      d = mu + phi*(z[2]- mu)
      c1 = u3*exp(-(z[1]-d)^2/2.0/sigma1)
      c1_left = d - sqrt(-2.0*sigma1*(log(u3)-(z[1]-d)^2/2.0/sigma1))
      c1_right = d + sqrt(-2.0*sigma1*(log(u3)-(z[1]-d)^2/2.0/sigma1))

      if (y[1] != 0.0){
        left = max(b1_left, c1_left)
        right = min(c1_right, a1_right)
      } else{
        left=c1_left
        right=min(c1_right,a1_right)
      }

      u = runif(1, 0,1)
      z[1] = left + (right-left)*u


      for( t in 2:(T-1)){

        u1 = runif(1, 0, 1)
        a1 = u1*exp(-z[t]/2.0)
        a1_right = -2.0*log(a1)

        if (y[t]  != 0.0){
          u2 = runif(1, 0,1)
          b1 = u2*exp(-y[t]*y[t]/2.0*exp(-z[t]))
          b1_left = -log(-2.0/y[t]/y[t]*( log(u2) -y[t]*y[t]/2.0*exp(-z[t] )  ))

          #b1_left = -log(-2.0/y[t]/y[t]*log(b1))
        }

        u3 = runif(1, 0,1)
        d = mu + phi* ( (z[t-1]-mu) +(z[t+1]-mu))/(1+phi^2)
        beta1 = sigma1/(1 +phi^2)

        c1 = u3*exp(-(z[t]-d)^2/2.0/beta1)
        c1_left = d - sqrt(-2.0*beta1*log(c1))
        c1_right = d + sqrt(-2.0*beta1*log(c1))

        if (y[t] != 0.0){
          left = max(b1_left, c1_left)
          right = min(c1_right, a1_right)
        }else{
          left = c1_left
          right = min(c1_right, a1_right)
        }

        u = runif(1,0,1)
        z[t] = left + (right-left)*u

      }

      u1 = runif(1, 0,1)
      a1 = u1*exp(-z[T]/2.0)
      a1_right = -2.0*log(a1)

      if (y[T]  != 0.0){
        u2 = runif(1, 0,1)
        b1 = u2*exp(-y[T]*y[T]/2.0*exp(-z[T]))
        b1_left = -log(-2.0/y[T]/y[T]*( log(u2) -y[T]*y[T]/2.0*exp(-z[T] )  ))
        #b1_left = -log(-2.0/y[T]/y[T]*log(b1))
      }

      u3 = runif(1, 0,1)
      d = mu + phi*(z[T-1]-mu)

      c1 = u3*exp(-(z[T]-d)^2/2.0/sigma1)
      c1_left = d-sqrt(-2.0*sigma1*log(c1))
      c1_right = d+sqrt(-2.0*sigma1*log(c1))

      if (y[T] != 0.0){
        left = max(b1_left,c1_left)
        right = min(c1_right,a1_right)
      }else{
        left = c1_left
        right = min(c1_right,a1_right)
      }


      u = runif(1,0,1)
      z[T] = left + (right-left)*u

      # Now we simulate the parameters


      a1 = (T-1)*(1.0-phi)^2 +(1-phi^2)
      a1 = a1/sigma1

      b1 = 0.0
      for (h in 1:(T-1)){
        b1 = b1 + (1.0-phi)*(z[h+1]-phi*z[h])
      }

      b1 = b1 +z[1]*(1-phi^2)
      b1 = b1/sigma1


      temp1 = rnorm(1,0,1)
      mu = b1/a1+1.0/sqrt(a1)*temp1



      # u = runif(1, 0,1)
      # latent = sqrt(1.0-phi^2)*u
      # bounder = min(1.0, sqrt(1.0-latent*latent))
      bounder = 1

      a1 = 0.0
      for(h in 1:(T-1)){
        a1 = a1+(z[h]-mu)^2
      }

      a1 = a1-(z[1]-mu)^2
      a1 = a1/sigma1

      b1 = 0.0
      for (h in 1:(T-1)){
        b1 = b1+(z[h]-mu)*( z[h+1]- mu )
      }

      b1 = b1/sigma1

      mu1 = b1/a1
      sig = sqrt(1/a1)


      u = runif(1, 0,1)

      p = pnorm((-bounder-mu1)/sig)+(pnorm((bounder-mu1)/sig)-pnorm((-bounder-mu1)/sig))*u;
      phi = mu1 + sig* qnorm(p, 0, 1);

      a1 = 0.0
      for (h in 1:(T-1)){
        a1 = a1 + ( z[h+1]-mu - phi*(z[h]-mu) )^2
      }
      a1 = a1+(z[1]-mu)^2 * (1.0-phi*phi )
      a1 = a1/2.0

      sigma1 = nimble::rinvgamma(1,(T/2-1), a1)

      if (k > n){

        mu_est = mu_est + mu
        phi_est = phi_est + phi
        sigma1_est = sigma1_est + sqrt(sigma1)
        est[k-n,]=c(mu, phi, sqrt(sigma1))
        z_est = z_est +z

        tem = 0
        for (i in 1:T){
          tem = tem  -2*log(pnorm(y[i], 0, exp(z[i]/2)))
        }
        likelihood_est = likelihood_est + tem
      }
    }


    pk = 3
    likehood1 = likelihood_est/(N-n)
    BIC1 = likehood1+pk*log(T)
    AIC1 = likehood1+2*pk


    mu = mu_est/(N-n)
    phi = phi_est/(N-n)
    sigma1 = sigma1_est/(N-n)
    z_est  = z_est/(N-n)


    mean(est[,1])
    mean(est[,2])
    mean(est[,3])

    sd(est[,1])
    sd(est[,2])
    sd(est[,3])


    pit = rep(0, T-1)

    for ( t in 2:T){
      temp = mu +phi*(z_est[t-1]-mu)
      temp = exp(temp/2)
      pit[t-1] = pnorm(y[t],0,  temp)
    }


    std_error =  rep(0, T)
    for ( t in 1:T){
      std_error[t] = y[t]*exp(-z_est[t]/2)
    }



    estimate = c(mu, phi, sigma1)
    sd1 = apply(est,2, sd)
    cf = t(hdi(est))
    df= data.frame(estimate, sd1, cf)
    colnames(df) = c("Est.", "Std.", "HPD CI(95%)_lower", "HPD CI(95%)_upper")

    df_simu_time_series = data.frame(est)
    colnames(df_simu_time_series) = c("mu", "phi", "simga")


    lt = list("estimate" = df, "standard_residuals"=std_error , "mcmc_series" = df_simu_time_series, "Aic" = AIC1, "Bic"=BIC1, "pit"=pit, "est_volatilities"= z_est)
    return(lt)
  }


  # now we fit ASV model

  if (toupper(model_type) %in% c("ASV")){

    T = length(y)

    mu = -10
    phi = 0.7
    rho = -0.41
    sigma1 = 0.19

    mu_est = 0
    phi_est = 0
    rho_est = 0
    sigma1_est = 0
    likelihood_est = 0
    likelihood = rep(0, N)
    z_est = rep(0, T)

    n = burn_in
    N = N
    z = rep(0, T) # the latent state vector

    est = matrix(rep(0, (N-n)*4), nrow = N-n)


    # initialized the states before the estimation process
    temp = rnorm(1, 0, 1)
    z[1] = mu + sigma1/sqrt(1.0 - phi*phi)*temp

    for (i in 2:T){
      temp = rnorm(1, 0, 1)
      z[i] = mu + phi* (z[i-1]-mu)+ sigma1*temp
    }

    tau = sigma1*sigma1*(1.0-rho^2)
    psi=rho*sigma1


    for (k in 1:N){
      #print(k)
      u1 = runif(1,0,1)
      a1 = u1*exp(-z[1]/2.0)
      a1_right = -2.0*log(a1)

      u2 = runif(1,0,1)
      b1 = u2*exp(-y[1]*y[1]/2.0*exp(-z[1] ) )

      if (y[1] != 0.0){
        # b1_left = -log(-2.0/y[1]/y[1]*log(b1))
        b1_left = -log(-2.0/y[1]/y[1]*( log(u2) -y[1]*y[1]/2.0*exp(-z[1] )  ))
      }

      u3 = runif(1,0,1)
      d = mu
      c1 = u3*exp(-(z[1]-d)^2/2.0/tau*(1-phi^2))
      c1_left = d - sqrt(-2.0*tau/(1-phi^2)*(log(u3)-(z[1]-d)^2/2.0/tau*(1-phi^2)))
      c1_right = d + sqrt(-2.0*tau/(1-phi^2)*(log(u3)-(z[1]-d)^2/2.0/tau*(1-phi^2)))

      if (y[1] != 0.0){
        left = max(b1_left, c1_left)
        right = min(c1_right, a1_right)
      } else{
        left=c1_left
        right=min(c1_right,a1_right)
      }

      u = runif(1, 0,1)
      rnd = left + (right-left)*u

      #a1=exp(-(z[2] -mu -phi*(rnd-mu)-psi*exp(-rnd/2.0)*y[1])^2/2.0/tau)/5
      a1 = min(1, exp(-(z[2] -mu -phi*(rnd-mu)-psi*exp(-rnd/2.0)*y[1])^2/2.0/tau+
                        (z[2] -mu -phi*(z[1]-mu)-psi*exp(-z[1]/2.0)*y[1])^2/2.0/tau ) )

      u = runif(1,0,1)
      if (u <= a1){
        z[1] = rnd
      }



      for( t in 2:(T-1)){

        u1 = runif(1, 0, 1)
        a1 = u1*exp(-z[t]/2.0)
        a1_right = -2.0*log(a1)

        if (y[t]  != 0.0){
          u2 = runif(1, 0,1)
          b1 = u2*exp(-y[t]*y[t]/2.0*exp(-z[t]))
          b1_left = -log(-2.0/y[t]/y[t]*( log(u2) -y[t]*y[t]/2.0*exp(-z[t] )  ))

          #b1_left = -log(-2.0/y[t]/y[t]*log(b1))
        }

        u3 = runif(1, 0,1)
        d = mu + phi*(z[t-1]-mu) + psi*exp(-z[t-1]/2.0)*y[t-1]

        c1 = u3*exp(-(z[t]-d)^2/2.0/tau)
        c1_left = d - sqrt(-2.0*tau*log(c1))
        c1_right = d + sqrt(-2.0*tau*log(c1))

        if (y[t] != 0.0){
          left = max(b1_left, c1_left)
          right = min(c1_right, a1_right)
        }else{
          left = c1_left
          right = min(c1_right, a1_right)
        }

        u = runif(1,0,1)
        rnd = left + (right-left)*u

        #a1=exp(-(z[t+1] -mu -phi*(rnd-mu)-psi*exp(-rnd/2.0)*y[t])^2/2.0/tau)/2
        a1 = min(1, exp(-(z[t+1] -mu -phi*(rnd-mu)-psi*exp(-rnd/2.0)*y[t])^2/2.0/tau+
                          (z[t+1] -mu -phi*(z[t]-mu)-psi*exp(-z[t]/2.0)*y[t])^2/2.0/tau ) )

        u = runif(1,0,1)
        if (u <= a1){
          z[t] = rnd
        }
      }

      u1 = runif(1, 0,1)
      a1 = u1*exp(-z[T]/2.0)
      a1_right = -2.0*log(a1)

      if (y[T]  != 0.0){
        u2 = runif(1, 0,1)
        b1 = u2*exp(-y[T]*y[T]/2.0*exp(-z[T]))
        b1_left = -log(-2.0/y[T]/y[T]*( log(u2) -y[T]*y[T]/2.0*exp(-z[T] )  ))
        #b1_left = -log(-2.0/y[T]/y[T]*log(b1))
      }

      u3 = runif(1, 0,1)
      d = mu + phi*(z[T-1]-mu) + psi*exp(-z[T-1]/2.0)*y[T-1]

      c1 = u3*exp(-(z[T]-d)^2/2.0/tau)
      c1_left = d-sqrt(-2.0*tau*log(c1))
      c1_right = d+sqrt(-2.0*tau*log(c1))

      if (y[T] != 0.0){
        left = max(b1_left,c1_left)
        right = min(c1_right,a1_right)
      }else{
        left = c1_left
        right = min(c1_right,a1_right)
      }


      u = runif(1,0,1)
      z[T] = left + (right-left)*u


      # Now we simulate the parameters


      a1 = (T-1)*(1.0-phi)^2 +(1-phi^2)
      a1 = a1/tau

      b1 = 0.0
      for (h in 1:(T-1)){
        b1 = b1 + (1.0-phi)*(z[h+1]-phi*z[h]-psi*exp(-z[h]/2.0)*y[h])
      }

      b1 = b1 +z[1]*(1-phi^2)
      b1 = b1/tau


      temp1 = rnorm(1,0,1)
      mu = b1/a1+1.0/sqrt(a1)*temp1



      # u = runif(1, 0,1)
      # latent = sqrt(1.0-phi^2)*u
      # bounder = min(1.0, sqrt(1.0-latent*latent))
      bounder = 1

      a1 = 0.0
      for(h in 1:(T-1)){
        a1 = a1+(z[h]-mu)^2
      }

      a1 = a1-(z[1]-mu)^2
      a1 = a1/tau

      b1 = 0.0
      for (h in 1:(T-1)){
        b1 = b1+(z[h]-mu)*( z[h+1]- mu -psi*exp( -z[h]/2.0 )*y[h] )
      }

      b1 = b1/tau

      mu1 = b1/a1
      sig = sqrt(1/a1)


      u = runif(1, 0,1)

      p = pnorm((-bounder-mu1)/sig)+(pnorm((bounder-mu1)/sig)-pnorm((-bounder-mu1)/sig))*u
      phi = mu1 + sig* qnorm(p, 0, 1)




      a1=0.0
      for ( h in 1:(T-1)){
        a1 = a1 + (exp(-z[h]/2.0)*y[h])^2
      }
      a1=a1/tau



      b1 = 0.0
      for ( h in 1:(T-1)){
        b1 = b1+exp(-z[h]/2.0)*y[h]*(z[h+1]-mu-phi*(z[h]-mu))
      }
      b1 = b1/tau


      mu1=b1/a1
      sig=1.0/sqrt(a1)

      temp = rnorm(1,0,1)
      psi = mu1 +sig*temp



      a1 = 0.0
      for (h in 1:(T-1)){
        a1 = a1 + ( z[h+1]-mu - phi*(z[h]-mu)-psi*exp(-z[h]/2.0)*y[h] )^2
      }
      a1 = a1+(z[1]-mu)^2 * (1.0-phi*phi )
      a1 = a1/2.0

      tau = nimble::rinvgamma(1,(T/2-1), a1)

      sigma1 = sqrt(psi*psi +tau)
      rho = psi/sigma1


      if (k > n){

        mu_est = mu_est + mu
        phi_est = phi_est + phi
        rho_est = rho_est +rho
        sigma1_est = sigma1_est + sigma1
        est[k-n,]=c(mu, phi, rho, sigma1)
        z_est = z_est +z

        tem = 0
        for (i in 1:T){
          tem = tem  -2*log(pnorm(y[i], 0, exp(z[i]/2)))
        }
        likelihood_est = likelihood_est + tem
      }
    }

    pk = 4
    likehood1 = likelihood_est/(N-n)
    BIC1 = likehood1+pk*log(T)
    AIC1 = likehood1+2*pk

    mu = mean(est[,1])
    phi = mean(est[,2])
    rho = mean(est[,3])
    sigma1 = mean(est[,4])
    z_est = z_est/(N-n)


    # pit = rep(0, T-1)
    #
    # for ( t in 2:T){
    #   temp = mu +phi(z_est[t-1]-mu)
    #   temp_v = exp(temp/2)
    #   pit[t] = pnorm(y[t],0,  temp_v)
    # }


    pit = rep(0, T-1)

    for ( t in 2:T){
      temp2 = mu + phi*(z_est[t-1]-mu)+ sigma1*rho*y[t-1]*exp(-z_est[t-1]/2)
      temp2 = exp(temp2/2)
      pit[t-1] = pnorm(y[t],0,  temp2)
    }


    std_error =  rep(0, T)
    for ( t in 1:T){
      std_error[t] = y[t]*exp(-z_est[t]/2)
    }



    estimate = c(mu, phi, rho, sigma1)
    sd1 = apply(est,2, sd)
    cf = t(hdi(est))
    df= data.frame(estimate, sd1, cf)
    colnames(df) = c("Est.", "Std.", "HPD CI(95%)_lower", "HPD CI(95%)_upper")

    df_simu_time_series = data.frame(est)
    colnames(df_simu_time_series) = c("mu", "phi","rho" ,"simga")


    lt = list("estimate" = df, "standard_residuals"=std_error , "mcmc_series" = df_simu_time_series, "Aic" = AIC1, "Bic"=BIC1, "pit"=pit, "est_volatilities"= z_est)
    return(lt)
  }
}


