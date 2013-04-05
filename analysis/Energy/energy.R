library("Hmisc")

colors<- c( rgb(189,215,231,maxColorValue=555),
            rgb(189,215,231,maxColorValue=255),
            rgb(107,174,214,maxColorValue=555),
            rgb(107,174,214,maxColorValue=255),
            rgb(186,228,179,maxColorValue=555),
            rgb(186,228,179,maxColorValue=255))

energies <- read.table("./data.dat")

pdf(file="energy.pdf",height=6.0,width=7.5)

plot(energies$V1,energies$V2,ann=FALSE, ylim = c(-1.2,-0.8))
#errbar(energies$V1, energies$V2, energies$V2 + energies$V3, energies$V2 - energies$V3, ann=FALSE)
#lines(energies$V1,energies$V2 + energies$V3,col=colors[2])
#lines(energies$V1,energies$V2 - energies$V3,col=colors[2])
curve(0.2018*(1 - exp(-1.178*(x - 1.371)))^2 - 1.165, add = TRUE, col = "blue", lty=5, ann=FALSE)
curve(0.1633*(1 - exp(-1.087*(x - 1.408)))^2 - 1.152, add = TRUE, col = "red", lty=5, ann=FALSE)
curve(1/2*0.372453129147*(x-1.4239868277)**2-1.15121538293, add = TRUE, col = "green", lty=5, ann=FALSE)
#curve(0.173096419211*x^2 - 0.493225078702*x - 0.799972863613, add = TRUE, col = "green", lty=5, ann=FALSE)
title(xlab="Distance between nuclei (Bohr radii)",
        ylab="Energy (Hartree)")
legend(2.65,-1.1, # places a legend at the appropriate place
c("Data","Wide range Morse fit","Narrow range Morse fit","Harmonic fit"), # puts text in the legend
lty=c(1,1,1,1), # gives the legend appropriate symbols (lines)

lwd=c(1,1,1,1),col=c("black","blue","red","green"))


dev.off()
