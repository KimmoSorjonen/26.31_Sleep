

############################# BIENVENUE #############################

############# I SOLEMNLY SWEAR THAT I AM UP TO NO GOOD ##############

## Loading package

library(MASS)
library(metafor)
library(lavaan)

#################
## Data
## Correlations reported in Peng et al. (2026)
## order: igd1, sleep1, igd2, sleep2, igd3, sleep3

rm <- matrix(c(
  
  1.000, 0.402, 0.455, 0.210, 0.394, 0.183,
  0.402, 1.000, 0.219, 0.457, 0.194, 0.379,
  0.455, 0.219, 1.000, 0.400, 0.530, 0.272,
  0.210, 0.457, 0.400, 1.000, 0.265, 0.497,
  0.394, 0.194, 0.530, 0.265, 1.000, 0.385,
  0.183, 0.379, 0.272, 0.497, 0.385, 1.000), nrow=6)

n <- 20137  ## sample size

df <- data.frame(mvrnorm(n=n, ## generating data frame with required size and corr.
            mu=rep(0,6), Sigma=rm, empirical=T))

## The six models, with dy1 = Y2-Y1 and dy2 = Y3-Y2

m1 <- "dy1 ~ x1 + y1"
m2 <- "dy2 ~ x2 + y2"
m3 <- "dy1 ~ x1 + y2"
m4 <- "dy2 ~ x2 + y3"
m5 <- "dy1 ~ x1"
m6 <- "dy2 ~ x2"

modlist <- list(m1,m2,m3,m4,m5,m6) ## list with all six models

#################
## Figure

panlab <- c("Gaming → Sleep", 
            "Sleep → Gaming")

slab <- c("7.RMA","",
          "6.y3-y2.no.adj", "5.y2-y1.no.adj", ## y-labels
          "4.y3-y2.adj.y3", "3.y2-y1.adj.y2", 
          "2.y3-y2.adj.y2", "1.y2-y1.adj.y1")

cx <- 0.7 ## sizing factor
r.low <- -0.4 ## range, lower
r.upp <- 0.2 ## range, upper
f.upp <- r.upp+1*(r.upp-r.low) ## room for text
tic <- 0.2 ## distance, tics

if(dev.cur()==2) dev.off() ## removing earlier plots
par(mar=c(1,1,1.2,0), oma=c(1.5,4,0,1), mfrow=c(2,1)) ## setting margins and layout

for(i in 1:2){ ## 1. ST -> LO; 2. LO -> ST
  
  if(i==1) names(df) <- c("x1","y1","x2","y2","x3","y3") ## igd as X and sleep as Y
  if(i==2) names(df) <- c("y1","x1","y2","x2","y3","x3") ## sleep as X and igd as Y
  
  dy1 <- df$y2 - df$y1 ## difference score 1
  dy2 <- df$y3 - df$y2 ## difference score 2
  
  numeff <- 8 ## number of rows per panel (one is empty)
  
  plot(c(r.low,f.upp),c(0.5,8.5), type="n",xaxt="n",yaxt="n",xlab="",ylab="") ## empty plot
  
  mtext(panlab[i],3, line=0.2, cex=cx) ## panel label
  
  axis(1,at=seq(r.low,r.upp,tic),labels=F, cex.axis=cx) ## x-labels
  if(i==2) axis(1,at=seq(r.low,r.upp,tic),labels=seq(r.low,r.upp,tic), cex.axis=cx) ## x-labels
  #axis(2,at=1:8,labels=F, las=1, cex.axis=cx) ## y-tics
  axis(2,at=1:8,labels=slab, las=1, cex.axis=cx) ## y-labels
  
  lines(c(0,0),c(-1,15),col="gray") ## vertical gray line at x=0
  
  ball <- vector() ## to be filled with effects below
  seall <- vector() ## to be filled with standard errors below
  
  for(j in 1:6){ ## for the six models
    
    fit <- lm(modlist[[j]], data=df) ## fitting the model
    b <- fit$coefficients[2] ## reg. coefficient
    low <- confint(fit)[2,1] ## CI, low
    upp <- confint(fit)[2,2] ## CI, upp
    se <- summary(fit)$coefficients[2,2] ## standard error
    tx <- paste(round(b,2)," [", round(low,2),"; ", 
                round(upp,2),"]", sep="") ## string with values
    
    ball <- c(ball,b) ## adding effect to object
    seall <- c(seall,se) ## adding se to object
    
    lines(c(-2,r.upp),c(numeff,numeff), col="gray") ## vertical gray line
    arrows(low,numeff,upp,numeff,angle=90,code=3,length=0.05,lwd=2) ## CI
    points(b,numeff,pch=21,col="black",bg="black") ## point for effect
    text(r.upp,numeff,tx,pos=4,cex=cx) ## adding string with values
    
    numeff <- numeff-1 ## next row, please
  }
  
  ## Meta-analysis of the six effects
  
  b.fish <- 0.5*log((1+ball)/(1-ball)) ## Fisher's transformation of effects
  se.fish <- 0.5*log((1+seall)/(1-seall)) ## Fisher's transformation of SE
  
  ###### Following Bartos et al.
  
  we <- rep(1/6,6) ## equal weight to all effects
  
  random1 <- rma(yi=b.fish, vi=se.fish^2)
  res.ma <- rma(yi=b.fish, vi=(se.fish^2)/we, tau=random1$tau2)
  
  ######
  
  pred <- predict(res.ma, transf=transf.ztor) ## transforms back from Fisher's
  metb <- pred$pred ## estimate
  metlow <- pred$ci.lb ## lower CI
  metupp <- pred$ci.ub ## upper CI
  
  ##
  
  lines(c(-2,4),c(2,2),lty=2) ## dashed line
  lines(c(-2,r.upp),c(1,1),col="gray") ## gray vertical line
  
  polygon(c(metlow,metb,metupp,metb), c(1,1.5,1,0.5), ## diamond for RMA 
          col = "black", border = "black", lwd = 1)
  
  tx.ma <- paste(round(metb,2)," [", round(metlow,2),"; ", 
                 round(metupp,2),"]", sep="") ## string with values
  text(r.upp,1,tx.ma,pos=4,cex=cx) ## adding string of values
  legend("topleft", title=LETTERS[i], legend="", bty="n", inset=0, cex=1.5*cx) ## A-B legend
}


####################
## Alternative model

mosla <- "

## Trait

tx =~ 1*x1+1*x2+1*x3
ty =~ 1*y1+1*y2+1*y3

tx ~~ ty

## State

st1 =~ sx*x1+sy*y1
st2 =~ sx*x2+sy*y2
st3 =~ sx*x3+sy*y3

st2 ~ as*st1
st3 ~ as*st2

## (Error) variances

x1 ~~ x1
x2 ~~ x2
x3 ~~ x3

y1 ~~ y1
y2 ~~ y2
y3 ~~ y3

tx ~~ tx
ty ~~ ty

st1 ~~ 1*st1
st2 ~~ 1*st2
st3 ~~ 1*st3

"

names(df) <- c("x1","y1","x2","y2","x3","y3") ## igd as X and sleep as Y

fit <- lavaan(mosla, meanstructure=F, data=df) ## fitting mosla
summary(fit, fit.measures=T, standardized=T, ci=T) ## let's have a look



########################## MISCHIEF MANAGED #########################

############################# AU REVOIR #############################

