
library(PracTools) #📦
####################################################################### Simple Random Sampling
################   Sample size for a target CV ###################

## the population CV of some variable is CV = 2.0
## If the population is extremely large and CV0 (the target CV ) is set to 0.05

(nCont(CV0 = 0.05, CVpop = 2))

# If the population size is N = 500

(nCont(CV0 = 0.05,CVpop = 2, N = 500))  # 381 after roundin

######################## Sample sizes for a vector of target CVs ####################
## Often it will be useful to show a client the sample sizes for a series of precision target
# ample sizes for a vector of  values of CV0 from 0.01 to 0.21 in increments of 0.02

ceiling(nCont(CV0 = seq(0.01, 0.21, 0.02), CVpop=2))

######## Sample sizes for proportions based on an MOE #############################

# Suppose that we want to estimate a proportion for a characteristic where an advance estimate is Pu = 0.5
# The MOE is to be e when  α = 0.05

ceiling(nPropMoe(moe.sw=1, e=seq(0.01,0.08,0.01), alpha=0.05, pU=0.5))

####################################################################################  Stratified Simple Random Sampling


## n.tot = fixed total sample size
## Nh = vector of pop stratum sizes or pop stratum proportions (required parameter)
## Sh = stratum unit standard deviations, required unless alloc = "prop"
## cost = total variable cost
## ch = vector of costs per unit in strata
## V0 = fixed variance target for estimated mean
## CV0 = fixed CV target for estimated mean
## ybarU = pop mean of y
## alloc = type of allocation, must be one of "prop", "neyman", "totcost", "totvar"

# Neyman allocation
Nh <- c(215, 65, 252, 50, 149, 144)
Sh <- c(26787207, 10645109, 6909676, 11085034, 9817762, 44553355)
strAlloc(n.tot = 100, Nh = Nh, Sh = Sh, alloc = "neyman")

# cost constrained allocation
ch <- c(1400, 200, 300, 600, 450, 1000)
strAlloc(Nh = Nh, Sh = Sh, cost = 100000, ch = ch, alloc = "totcost")

# allocation with CV target of 0.05
 strAlloc(Nh = Nh, Sh = Sh, CV0 = 0.05, ch = ch, ybarU = 11664181, alloc = "totvar")








################    Probability proportional to size ####################




data("smho.N874")

## Extract EXPTOTAL & BEDS

y <- smho.N874[,"EXPTOTAL"]
x <- smho.N874[, "BEDS"]
y <- y[x>0]
x <- x[x>0]
ybarU <- mean(y)

(N<- length(x))
CV0<-0.15

# calculate V1 based on pp(x) sample
pik <- x/sum(x)
T <- sum(y)
(V1 <- sum( pik*(y/pik - T)^2))


# calculate V1 based on pp(x) sample
n <- V1 / (N*ybarU*CV0)^2
(n <- ceiling(n))

# Anticipated SE for the pps sample
(cv.pps <- sqrt(V1/(N^2*n)) / ybarU)

# sample size for an srs to produce the same SE
ceiling(nCont(CV0 = cv.pps, S2 = var(y), ybarU = ybarU, N = N))

#
# Variance Component Estimation in Multistage Sampling
# link:  https://cran.r-project.org/web/packages/PracTools/vignettes/Varcomps-multistage.pdf




