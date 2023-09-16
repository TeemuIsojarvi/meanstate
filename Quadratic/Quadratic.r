# This program calculates the ground state wave function of a harmonic 1D oscillator with an unusual method
# based on propagating an initial state along the real time axis, and finding the average wave function during
# a given time interval. This is repeated for several iterative steps where the final result of previous iteration
# is set to be the initial state in a new one. The result is then compared to the result obtained in a more usual
# way, propagating along the negative imaginary time axis.
#
# Teemu Isojärvi 9/16/2023


library(graphics)        							 #load the graphics library needed for plotting

# In this first part of the calculation, time evolution is computed along the imaginary time axis

lx <- 6.0              						       	         #length of the computational domain
lt <- 3.0		        						 #length of the simulation time interval
nx <- 120		  							 #number of discrete lattice points
nt <- 500 		  							 #number of timesteps
dx <- lx/nx		 							 #length of one discrete lattice cell
dt <- -1i*lt/nt 	       						      	 #length of timestep, set negative imaginary

V = c(1:nx)									 #potential energies at discrete points

for(j in c(1:nx)) {
V[j] = as.complex(2*(j*dx-3)*(j*dx-3))			
}

kappa1 = (1i)*dt/(2*dx*dx)							 #an element needed for the matrices	        			        
kappa2 <- c(1:nx)								 #another element

for(j in c(1:nx)) {
kappa2[j] <- as.complex(kappa1*2*dx*dx*V[j])
}

psi_ground = as.complex(c(1:nx))                                                        #array for the wave function values

for(j in c(1:nx)) {
psi_ground[j] = as.complex(exp(-1.5*(j*dx-3)*(j*dx-3)))					 #Gaussian initial wavefunction
}
	
xaxis1 <- c(1:nx)*dx								 #the x values corresponding to the discrete lattice points

A = matrix(nrow=nx,ncol=nx)								 #matrix for forward time evolution
B = matrix(nrow=nx,ncol=nx)								 #matrix for backward time evolution

for(j in c(1:nx)) {
for(k in c(1:nx)) {
A[j,k]=0
B[j,k]=0
if(j==k) { 
A[j,k] = 1 + 2*kappa1 + kappa2[j]						 #diagonal elements
B[j,k] = 1 - 2*kappa1 - kappa2[j]
}
if((j==k+1) || (j==k-1)) {
A[j,k] = -kappa1									 #off-diagonal elements
B[j,k] = kappa1
}
}
}

for (k in c(1:nt)) {                     					 #main time stepping loop

sol <- solve(A,B%*%psi_ground)								 #solve the system of equations

for (l in c(1:nx)) { 
psi_ground[l] <- sol[l] 
}
                                   			         
}

nrm <- 0

for(j in c(1:nx)) {
nrm <- nrm + abs(psi_ground[j])*abs(psi_ground[j])*dx
}

for(j in c(1:nx)) psi_ground[j] <- psi_ground[j]/sqrt(nrm)


# In this next part of calculation, the time interval and discretization parameters can be slightly different,
# and time evolution is taken along real time axis. After several iterations and time averaging of the wave function,
# the result is compared to that of the imaginary time version

              						       	         	 #length of the computational domain
lt <- 3.141592*1.5		        				         #length of the simulation time interval
nx <- 140		  							 #number of discrete lattice points
nt <- 150		  							 #number of timesteps
dx <- lx/nx		 							 #length of one discrete lattice cell
dt <- lt/nt 	       						      	         #length of timestep
nr <- 4                                                                          #number of iterations

xaxis2 <- c(1:nx)*dx

V = c(1:nx)									 #potential energies at discrete points

for(j in c(1:nx)) {
V[j] = as.complex(2*(j*dx-3)*(j*dx-3))			
}

kappa1 = (1i)*dt/(2*dx*dx)							 #an element needed for the matrices	        			        
kappa2 <- c(1:nx)								 #another element

for(j in c(1:nx)) {
kappa2[j] <- as.complex(kappa1*2*dx*dx*V[j])
}

psi = as.complex(c(1:nx))                                                        #array for the wave function values
psi_avg = as.complex(c(1:nx))                                                    #array for mean wave function values
psi_rec = matrix(nrow=nt,ncol=nx)                                                #matrix for wave function evolution

psi_avg_nrm = as.complex(c(1:nx))

for(j in c(1:nx)) {
psi[j] = as.complex(exp(-1.5*(j*dx-3)*(j*dx-3)))					 #Gaussian initial wavefunction
}

nrm <- 0

A = matrix(nrow=nx,ncol=nx)								 #matrix for forward time evolution
B = matrix(nrow=nx,ncol=nx)								 #matrix for backward time evolution

for(j in c(1:nx)) {
for(k in c(1:nx)) {
A[j,k]=0
B[j,k]=0
if(j==k) { 
A[j,k] = 1 + 2*kappa1 + kappa2[j]						 #diagonal elements
B[j,k] = 1 - 2*kappa1 - kappa2[j]
}
if((j==k+1) || (j==k-1)) {
A[j,k] = -kappa1									 #off-diagonal elements
B[j,k] = kappa1
}
}
}

for (r in c(1:nr)) {

for (j in c(1:nt)) {
for (l in c(1:nx)) {
psi_rec[j,l] <- 0
}
}

for (k in c(1:nt)) {                     					 #main time stepping loop

sol <- solve(A,B%*%psi)								 #solve the system of equations

for (l in c(1:nx)) { 
psi[l] <- sol[l] 
psi_rec[k,l] <- sol[l]					#The wave function values are recorded throughout the simulation
}

for (p in c(1:nx)) psi_avg[p] <- 0

for (q in c(1:k)) {
for (p in c(1:nx)) {
psi_avg[p] <- psi_avg[p] + psi_rec[q,p]
}
}

for (p in c(1:nx)) psi_avg[p] <- psi_avg[p]/k		#Calculate the time average wave function from recorded values

nrm <- 0

for(j in c(1:nx)) {
nrm <- nrm + abs(psi_avg[j])*abs(psi_avg[j])*dx
}

for (p in c(1:nx)) psi_avg[p] <- psi_avg[p]/sqrt(nrm)

# The time averaged WF is plotted on every timestep of each iteration, and compared in the same image to the imaginary time
# result and the shape of the potential energy function V(x)

jpeg(file = paste("iter_",r,"_plot_",k,".jpg",sep=""),width = 2.00, height = 2.00, units = "in", res = 300, pointsize = 3)
par(mar = c(6, 6, 5, 5))
plot(xaxis2,abs(psi_avg)*abs(psi_avg),xlab="position (x)", ylab=expression(abs(psi[avg])^2),type='l',col='red',ylim=c(0,1.2),cex.lab=2)
lines(xaxis1,abs(psi_ground)*abs(psi_ground),type='l',col='blue',ylim=c(0,1.2))
lines(xaxis2,0.2*V,col='black',ylim=c(0:1.2))
legend(0.1, 1.1, legend=c("Ground state from time averaging", "Ground state from ITP", "Shape of V(x)"), fill =c("red","blue","black"))
title(paste("Probability density at t = ",round(k*dt,digits=3)," on iteration ",r))
dev.off()

}

for (p in c(1:nx)) psi[p] <- psi_avg[p]

}
