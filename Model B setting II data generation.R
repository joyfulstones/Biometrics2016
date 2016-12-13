data.gen=function(alpha.0, beta.0, delta.0, theta.0, gamma.0, shape.R, shape.D, n.covar, n.id, r.0, h.0, p)
# To generate the recurrent event data with Weibull hazard, parameter r0 and k
# Hazard = k [r0 * exp(X beta) ^ 1/k] ^ k * t ^(k-1) 
{
  max.obs=n.id*16
  id=rep(0,max.obs)
  starttime=rep(0,max.obs)
  stoptime=rep(0,max.obs)
  trt=rep(0,max.obs)
  age=rep(0,max.obs)
  X1=rep(0, max.obs)
  cure=rep(0,max.obs)
  event=rep(0, max.obs)
  death.time=rep(0, max.obs)
  enum=event
  covar=matrix(0, max.obs, n.covar)
  j=1

    for (i in 1:n.id)
    {

      #censor.time=runif(1, 0, 1)
       censor.time=1 + runif(1, 0, 4)
      # censor.time=runif(1, 0, 2)
      X1=ifelse(runif(1) < p , 1, 0)
      covar.i=c(1, X1)
      event.time=0
      a=rnorm(1, 0, sqrt(theta.0))    # normal frailty

      #event.rate=r0 * exp(beta.0 %*% covar.i + a)
      haz.recur=r.0 * exp(beta.0 * X1 +a) ^ (1/shape.R)

    haz.death=h.0 * exp(delta.0 * X1 + gamma.0 * a) ^ (1/shape.D)
    death.time=rweibull(1, shape.D, 1/haz.death)

    follow.up=ifelse(censor.time>=death.time, death.time, censor.time)
    status=ifelse(censor.time>=death.time, 2 , 0)


      enum[j]=1

      prob=1/(1+exp(-alpha.0 %*% covar.i ))

      runif.1=runif(1)

      if (runif.1<prob)
      {
        starttime[j]=0
        stoptime[j]=follow.up
        event[j]=status
        covar[j, ]=X1
        id[j]=i
        cure[j]=1

        j=j+1
        # enum[j]=enum[j-1]+1
 
      }
      if (runif.1>=prob)
      {
          repeat
          {
            id[j]=i
            #gap.time=rexp(1, event.rate)
            increment= - log(runif(1)) / haz.recur ^shape.R
            starttime[j]=event.time
            event.time=(event.time^shape.R + increment) ^ (1/shape.R)
            cure[j]=0
            covar[j, ]=X1

            if (event.time> follow.up)
            {
              stoptime[j]=follow.up
              event[j]=status
       #       event[j]=0
              j=j+1
              # enum[j]=enum[j-1]+1
              break
            }
            else
            {
              stoptime[j]=event.time
              event[j]=1
              enum[j+1]=enum[j]+1
              j=j+1
            } # End of else


          } # End of Repeat

        } # End of if (runif.1>=prob)


    } # End of for (i in 1:n.id)

  x=data.frame(id, covar, starttime, stoptime, event, enum, cure)
  x[id>0,]
}



##############################################################################################

library(MASS)

# True values of parameters
  n.id=500
  alpha.0=c(0, -1)
  beta.0=1.2
  delta.0=1
  gamma.0=1
  theta.0=.4
  shape.R=1.25
  shape.D=1.25

  n.covar=1
  r.0=.25
  h.0=.1
  p=.5

  covar.names="covar"


#######

rand.num=1001:1600
n.replicate=length(rand.num)
recur=rep(0, n.replicate)
pct.zero=recur
pct.cure=recur
pct.censor=recur



for (ii in 1:n.replicate)
{
# Generate the data

  set.seed(rand.num[ii])


  data.all=data.gen(alpha.0, beta.0, delta.0, theta.0, gamma.0, shape.R, shape.D, n.covar , n.id , r.0, h.0,  p)

#  attach(data.all)
  max.obs=length(data.all$id)
  recur[ii]=sum(data.all$event==1)/n.id
  pct.zero[ii]=1-sum(data.all$enum==2)/n.id
  pct.cure[ii]=sum(data.all$cure)/n.id
  pct.censor[ii]=sum(data.all$event==0)/n.id



    filedir="c:\\liulei\\paper\\"
    filename=paste(filedir, ii, ".csv", sep="") # Generate comma delimited file

    write.matrix(as.matrix(data.all), filename, sep=",")


}

 mean(recur)     

 mean(pct.zero)     

 mean(pct.cure)     

 mean(pct.censor)  

>  mean(recur)     
[1] 0.8091033
> 
>  mean(pct.zero)     
[1] 0.6441333
> 
>  mean(pct.cure)     
[1] 0.3835233
> 
>  mean(pct.censor)  
[1] 0.6622367
