inFile = '../hydrogen/molecule/energies.dat' #Set input file

# Values needed to set up matrices
numsVals = 20  # Number of separate values tried in outer loop
numBetaVals = 20   # Number of separate values tried in inner loop

# Read in data
dat<-read.table(inFile)

# Shape list of values into matrix
s=matrix(dat$V1,nrow=numsVals,ncol=numBetaVals)
beta=matrix(dat$V2,nrow=numsVals,ncol=numBetaVals)
z=matrix(dat$V3,nrow=numsVals,ncol=numBetaVals)

# Temporary code to generate some distribution to make pretty plot
# Get rid of when energy properly implemented
z = exp(-(s-1.75)^2*10)*exp(-(beta-.5)^2*50)/5

# Get a column/row of s/beta accordingly for plotting
sVals=s[1,]
betaVals=beta[,1]

#Set up plot coloring, relate to z value
jet.colors<-colorRampPalette(c('blue','green','yellow'))
nbcol<-100
color<-jet.colors(nbcol)
nrz<-nrow(z)
ncz<-ncol(z)
zfacet<-z[-1,-1]+z[-1,-ncz]+z[-nrz,-1]+z[-nrz,-ncz]
facetcol<-cut(zfacet,nbcol)

# Make 3d plot and view at angle theta, phi with colors set up
persp(sVals,betaVals,z,theta=45,phi=45,col=color[facetcol],ticktype='detailed')
