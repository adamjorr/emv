#The largest plant genome is 150GBP (150e9 bases)
#The probability of having any number of incorrect basecalls for a genome of this size
#	is given by a binomial distribtion with parameters 150e9 and p,
#	where p is the probability that a base is called incorrectly.
#	To find the probability of calling each base such that there is a 95%
#	probability that there are no errors in the genome, we use the pdf of the
#	binomial to find F(X = 0).

x <- seq(from = -14, to = -10, length.out = 100)
p_0_incorrect = function(x){ dbinom(0,150e9,x) }
plot(x, sapply(10^x, FUN = p_0_incorrect))
abline(a = .95, b = 0)

#Approximately 95% is ~ x = 10^-12.466
#Analytically the solution is 1-exp(log(.95)/150e9) (appx 3.419487e-13)
#Human: 3.3e9
#Amoeba:670e9

##
## plot a 3d space of to find appropriate
##
p_incorrect = function(cov, percent, error){ sum(dbinom( floor(percent * cov ):cov, cov, error  )) }

#full plot
x <- 1:50
y<-seq(.01,1,by=.01)
z<- outer(x,y,Vectorize(p_incorrect), error=.01)
p <- list(x=x,y=y,z=z)
image2D(p,xlab="Coverage",ylab="Consensus Level",main="Probability of an Incorrect Call",contour=list(levels=c(.05,.001,1e-7,1e-11),col="white",labcex=.8,lwd=4),lighting=TRUE,rasterImage=TRUE)

#Note that as the amount of coverage increases, the level of consensus required to call a base at a particular level of certainty decreases exponentially

#This also means that at a fixed consensus level, the probability of an incorrect call decreases exponentially with coverage
x <- 2:50
y <- sapply(x, FUN=p_incorrect, percent = .4, error = .01)
ggplot() + geom_line(aes(x,y)) + scale_y_log10()

#We expect the runtime of tools like GATK to increase exponentially with coverage as well, but this requires data collection



p_incorrect_cutoff = function(cov, cutoff, error){sum(dbinom(cutoff:cov,cov,error))}
minimum_cutoff = function(coverages, p){
	coverage = 1:coverages
	min(which(sapply(coverage, FUN=p_incorrect_cutoff, cov = max(coverage), error = .01) < p))
	} #gives first cutoff where p_incorrect < p
coverage <- 1:50
cutoffs <- sapply(coverage, FUN=minimum_cutoff, p=1e-11) #seems some kind of number pattern here, possibly associated with the binomial
#For p=1e-11, the number of times the next cutoff appears equals the number of times the current cutoff appears + 3 if the current number is even, +2 if the current number is odd
ggplot() + geom_line(aes(coverage,cutoffs))






