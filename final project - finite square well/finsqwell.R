###############################################################################
# Author: 		Adam J Jaso, Jr
# Title:		Calculation of the Energy Levels of the Infinite Square Well in 1 Dimension
# Date:			12/05/2011
# Class:		PHYS 140 Computational Methods for Physicists
# Instructor:	Dr. Olencka Hubickyj
###############################################################################

start = Sys.time();
###############################################################################
# INPUTS
	# potential
	V0 = 40; #eV
	# well width
	a = 10; #nm
# END INPUTS
###############################################################################
# CONTROL VARIABLES
# search method
	search_method = 0;
	# bisection			0
	# newton-raphson	1
# search method  tolerance
	root_dx = 0.0001;
# outputprecision
	precision = 9;
# plot settings
	# search function
		# the cutoff point for calculated points
			fceil = 100;
		# the maximum allowed distance between two points
			# this value should be adjusted appropriately with the well depth and widths,
			# since it affects the resolution required to find smaller and smaller spaced roots
			maxdist = 4.9; 
		# the default dx interval
			defdx = .00001;
	# wave equation
		# maximum distance between plot points
			psi_maxdist = .05;
		# factor of the well width to be shown on either side of the origin
			psi_n = 2;
		# mode to be plotted
			psi_mode = 29;
# END CONTROL VARIABLES
###############################################################################
# CONSTANTS
	# m*c^2 (m = mass of the electron)
	mc2 = 0.511e6; #eV
	# hbar*c
	hbarc = 197.3; #eV * nm
	# 2*m/hbar^2
	twomhbar2 = 2 * mc2 / hbarc^2;
	K = twomhbar2 * V0;
# END CONSTANTS
###############################################################################
# FUNCTIONS
	# primary function
	f <- function(x,K,a) sqrt(K)*(sqrt(1-x)*tan(sqrt(K*(1-x))*a)-sqrt(x));

	# first derivative
	fp <- function(x,K,a) {
		z = 1-x;
		-sqrt(K)/2*(1/sqrt(z)*tan(sqrt(K*z)*a)+a*sqrt(K)*1/cos(sqrt(K*z)*a)^2+1/sqrt(x));
	}
	
	# second derivative
	fpp <- function(x,K,a) {
		A = -1/(2*z)*(1/sqrt(z)*tan(sqrt(K*z)*a)+a*sqrt(K)*1/cos(sqrt(K*z)*a)^2);
		B = -K*a^2*1/sqrt(z)*1/cos(sqrt(K*z)*a)^2*tan(sqrt(K*z)*a);
		C = -1/2*1/sqrt(x)^3;
		
		-sqrt(K)/2*(A+B+C);
	}
	
	# kappa constant calculation
	k <- function(E,twomhbarc) sqrt(twomhbarc)*sqrt(E);
	
	# L constant calculation
	el <- function(E,twomhbarc,V0) sqrt(twomhbarc)*sqrt(V0-E);
	
	# energy level constant calculation from the solution, x
	Energy <- function(V0,x) V0*x;
	
	# normalization constants
	normconst <- function(x,el,kappa,a) {
		D=(a+1/kappa*(sin(el*a)*cos(el*a)+cos(el*a)^2))^(-1/2);
		if(x < -a) {
			return(D*exp(kappa*a)*cos(el*a));
		} else if(x >= -a && x <= a) {
			return(D);
		} else if(x > a) {
			return(D*exp(kappa*a)*cos(el*a))
		} else {
			return(0);
		}
	}
	
	# dx calculator
	dx_calculator <- function(dp,yp) dp/sqrt(yp^2+1);
	
	# newton raphson dx
	nr_dx_calculator <- function(fx,fpx) {
		-fx/fpx;
	}
	
	# point distance calculator
	dist <- function(x1,y1,x2,y2) sqrt((x1-x2)^2 + (y1-y2)^2);
	
	# the wave function
	psi <- function(x,el,kappa,a) {
		const = normconst(x,el,kappa,a)
		if( x < -a) {
			return(const*exp(kappa*x))
		} else if(x >= -a && x <= a) {
			return(const*cos(el*x));
		} else if(x > a) {
			return(const*exp(kappa*-x));
		} else {
			return(0);
		}
	}
	
	# the wave function first derivative
	psip <- function(x,el,kappa,a) {
		if(x < -a) {
			return(kappa*psi(x,el,kappa,a));
		} else if(x >= -a && x <= a) {
			return(-el*normconst(x,el,kappa,a)*sin(el*x));
		} else if(x > a) {
			return(-kappa*psi(x,el,kappa,a));
		} else {
			return(0);
		}
	}
	
	# the probability distribution (psi*psi) first derivative
	psi2p <- function(x,el,kappa,a) {
		if(x < -a) {
			return(2*kappa*psi(x,el,kappa,a)^2);
		} else if(x >= -a && x <= a) {
			return(-2*el*normconst(x,el,kappa,a)^2*cos(el*x)*sin(el*x));
		} else if(x > a) {
			return(-2*kappa*psi(x,el,kappa,a)^2);
		} else {
			return(0);
		}
	}
	
	# the probability distribution (psi*psi) second derivative
	psi2pp <- function(x,el,kappa,a) {
		if(x < -a) {
			return(2*kappa*psi(x,el,kappa,a)^2);
		} else if(x >= -a && x <= a) {
			return(-2*el*normconst(x,el,kappa,a)^2*cos(el*x)*sin(el*x));
		} else if(x > a) {
			return(-2*kappa*psi(x,el,kappa,a)^2);
		} else {
			return(0);
		}
	}
	
	# midpoint calculator
	mid <- function(lo,hi) (lo+hi)/2;
# END FUNCTIONS
###############################################################################
# GENERATE PLOT
par(mfrow=c(2,2));
# initial values
	root_x=0;
	root_xs = c();
	root_fxs = c();
	intervals = c();
	root_fofx_1 = f(root_x,K,a);
	root_fpofx_1 = fp(root_x,K,a);
# end initial values

# stat variables
	# for dx stats
	root_dxsum = 0;
	root_dx2sum = 0;
	# for distance stats (distance between graphed points)
	root_distsum = 0;
	root_dist2sum = 0;
# end stat variables

# counter variables
	# iteration counter
	root_iter = 0;
	# plot point counter
	root_ptcounter = 0;
	root_droppeddist = 0;
	# domain pair counter;
	intervalcounter = 0;
# end counter variables

# analyze interval
	while(root_x < 1) {
		root_fofx = f(root_x,K,a);
		root_fpofx = fp(root_x,K,a);
		root_dxcalc = dx_calculator(maxdist,root_fpofx);
		
		if(abs(root_fofx) < fceil && abs(root_fpofx) < Inf) {
			root_xs = c(root_xs,root_x);
			root_fxs = c(root_fxs,root_fofx);
			root_ptcounter = root_ptcounter+1;
			cat("x = ",root_x,"\tdxcalc = ",root_dxcalc,"\tf(x) = ",f(root_x,K,a));
			
			if(root_ptcounter > 1) {
				distance = dist(root_xs[root_ptcounter-1],root_fxs[root_ptcounter-1],root_xs[root_ptcounter],root_fxs[root_ptcounter]);
				if(distance < maxdist+1) {
					cat("\tdist = ",distance);
					root_distsum = root_distsum + distance;
					root_dist2sum = root_dist2sum + distance^2;
					if(root_fofx_1 > 0 && root_fofx < 0) {
						intervals = c(intervals,root_xs[root_ptcounter-1],root_xs[root_ptcounter]);
					}
				} else {
					cat("\tf'(x) = ",root_fpofx);
					root_droppeddist = root_droppeddist+1;
				}
			}
			cat("\n");
			
			root_x = root_x+root_dxcalc;
			root_dxsum = root_dxsum+root_dxcalc;
			root_dx2sum = root_dx2sum = root_dx2sum+root_dxcalc^2;
			defdx = root_dxcalc;
		} else {
			root_x = root_x+defdx;
			root_dxsum = root_dxsum+defdx;
			root_dx2sum = root_dx2sum+defdx^2;
		}
		root_fofx_1 = root_fofx;
		root_iter = root_iter+1;
	}
	
	root_dxavg = root_dxsum/root_iter;
	root_dxstd = sqrt(root_dx2sum / root_iter - root_dxavg^2);
	root_distavg = root_distsum / (root_ptcounter-root_droppeddist-1);
	root_diststd = sqrt(root_dist2sum / (root_ptcounter-root_droppeddist-1) - root_distavg^2);
# end analyze interval

# plot
	plot(root_xs,root_fxs,main="Roots to be Searched",
			xlab="E/V_0", ylab="Root Value");
	grid();
	#title(main="main title", sub="sub-title",
	#		xlab="x-axis label", ylab="y-axis label") 
# end plot

# END GENERATE PLOT
###############################################################################
# ROOT SEARCH METHODS

E = c();
kappa = c();
L = c();

root_intervalsum = 0;
root_itersum = 0;
root_iter2sum = 0;
root_iterstd = 0;
root_iteravg = 0;
root_itercount = 0;
root_itermax = 0;
root_itermin = 1000000000;

# newton-raphson search method
if(search_method == 1) {
	
	cat("\nEnergy Levels of the Finite Square Well (V0 =",V0,"eV  and  a =",a,"nm)\n");
	for(i in 1:(length(intervals)/2)) {
		root_intervalsum = root_intervalsum + intervals[2*i] - intervals[2*i-1];
		root_min = intervals[2*i-1];
		root_max = intervals[2*i];
		
		root_x = root_min;
		root_iter = 0;
		
		# error management
		root_x_1 = root_x;
		root_errorcount = 0;
		root_maxerrorcount = 10;
		root_errorlist = c();
		root_erroryes = 0;
		# yes		1
		# no		0
		
		while(abs(f(root_x,K,a)) > root_dx) {
			root_nr_dx =  nr_dx_calculator(f(root_x,K,a),fp(root_x,K,a));
			if(root_x+root_nr_dx < root_max) {
				root_x = root_x + root_nr_dx;
			} else {
				root_x = root_max;
			}
			
			root_iter = root_iter+1;
			
			if(root_x_1 == root_x) {
				root_errorcount = root_errorcount+1;
				if(root_errorcount > root_maxerrorcount) {
					root_errorlist = c(root_errorlist,i);
					root_erroryes = 1;
					root_errorcount = 0;
					break;
				}
			}
			root_x_1 = root_x;
		}
		
		# calculate results
		E = c(E,Energy(root_x,V0));
		kappa = c(kappa,k(E[i],twomhbar2));
		L = c(L,el(E[i],twomhbar2,V0));
		
		# solution output
		cat("\tLevel ",length(intervals)/2-i,"\tx  = ",format(root_x,digits=precision),"\tE = ",format(-E[i],digits=precision)," +/- ",V0*root_dx," eV\titer = ",root_iter,"\tkappa = ",kappa[i],"\tL = ",L[i]);
		if(root_erroryes == 1) {
			cat("\tWARNING: THE ZERO MAY NOT HAVE BEEN FOUND EXACTLY!");
			root_erroryes = 0;
		}
		
		# iteration statistics
		root_itersum = root_itersum + root_iter;
		root_iter2sum = root_iter2sum + root_iter^2;
		root_itercount = root_itercount + 1;
		
		if(root_iter > root_itermax) {
			root_itermax = root_iter;
		}
		if(root_iter < root_itermin) {
			root_itermin = root_iter;
		}
		
		
		if(root_itercount != 0) {
			root_iteravg = root_itersum / root_itercount;
			root_iterstd = sqrt(root_iter2sum / root_itercount - root_iteravg^2);
		}
		
		cat("\n");
	}
}
# end newton-raphson search method

# bisection search method
if(search_method == 0) {
	
	cat("\nEnergy Levels of the Finite Square Well (V0 =",V0,"eV  and  a =",a,"nm)\n");
	for(i in 1:(length(intervals)/2)) {
		root_intervalsum = root_intervalsum + intervals[2*i] - intervals[2*i-1];
		root_xlo = intervals[2*i-1];
		root_xhi = intervals[2*i];
		root_xmid = mid(root_xlo,root_xhi);
		root_x = root_xmid;
		root_iter = 0;
		
		# error management
		root_x_1 = root_x;
		root_errorcount = 0;
		root_maxerrorcount = 10;
		root_errorlist = c();
		root_erroryes = 0;
			# yes		1
			# no		0
		
		while(abs(f(root_x,K,a) ) > root_dx) {
			if(f(root_xlo,K,a) * f(root_xmid,K,a) < 0) {
				root_xhi = root_xmid;
				root_xmid = mid(root_xlo,root_xhi);
				root_x = root_xmid;
			} else if(f(root_xhi,K,a) * f(root_xmid,K,a) < 0) {
				root_xlo = root_xmid;
				root_xmid = mid(root_xlo,root_xhi);
				root_x = root_xmid;
			} else if(f(root_xhi,K,a) == 0) {
				root_x = root_xhi;
				break;
			} else if(f(root_xmid,K,a) == 0) {
				root_x = root_xmid;
				break;
			} else if(f(root_xlo,K,a) == 0) {
				root_x = root_xlo;
				break;
			}
			root_iter = root_iter+1;
			
			if(root_x_1 == root_x) {
				root_errorcount = root_errorcount+1;
				if(root_errorcount > root_maxerrorcount) {
					root_errorlist = c(root_errorlist,i);
					root_erroryes = 1;
					root_errorcount = 0;
					break;
				}
			}
			root_x_1 = root_x;
		}
		
		# calculate results
		E = c(E,Energy(root_x,V0));
		kappa = c(kappa,k(E[i],twomhbar2));
		L = c(L,el(E[i],twomhbar2,V0));
		
		# solution output
		cat("\tLevel ",length(intervals)/2-i,"\tx  = ",format(root_x,digits=precision),"\tE = ",format(-E[i],digits=precision)," +/- ",V0*root_dx," eV\titer = ",root_iter,"\tkappa = ",kappa[i],"\tL = ",L[i]);
		if(root_erroryes == 1) {
			cat("\tWARNING: THE ZERO MAY NOT HAVE BEEN FOUND EXACTLY!");
			root_erroryes = 0;
		}
		
		# iteration statistics
		root_itersum = root_itersum + root_iter;
		root_iter2sum = root_iter2sum + root_iter^2;
		root_itercount = root_itercount + 1;
		
		if(root_iter > root_itermax) {
			root_itermax = root_iter;
		}
		if(root_iter < root_itermin) {
			root_itermin = root_iter;
		}
		
		
		if(root_itercount != 0) {
			root_iteravg = root_itersum / root_itercount;
			root_iterstd = sqrt(root_iter2sum / root_itercount - root_iteravg^2);
		}
		
		cat("\n");
	}
	
}
# end bisection search method
# END ROOT SEARCH METHODS
###############################################################################
# WAVE EQUATION PLOT
#psi_maxdist = .0125;
#psi_n = 3;
psi_x = -psi_n*a;
psi_xs = c();
psi_psis = c();
psi_max = c();
psi_min = c();
psi_iter = 0;
#for(i in 1:length(E)) {
#psi_i = 1;
	while(psi_x < psi_n*a) {
		psi_xs = c(psi_xs,psi_x);
		psi_psis = c(psi_psis,psi(psi_x,L[psi_mode],kappa[psi_mode],a));
		#cat(dx_calculator(psi_maxdist,psip(psi_x,L[psi_mode],kappa[psi_mode],a)),"\n");
		psi_x = psi_x + dx_calculator(psi_maxdist,psip(psi_x,L[psi_mode],kappa[psi_mode],a));
		psi_iter = psi_iter+1;
		
		if(length(psi_xs) > 1) {
			psi_psi2p = psi2p(psi_xs[psi_iter],L[psi_mode],kappa[psi_mode],a);
			psi_psi2p_1 = psi2p(psi_xs[psi_iter-1],L[psi_mode],kappa[psi_mode],a);
			if(psi_psi2p < 0 && psi_psi2p_1 > 0) {
				psi_max = c(psi_max,psi_xs[psi_iter-1],psi_xs[psi_iter]);
				#cat("psi_intervals = ",psi_xs[psi_iter-1]," -> ",psi_xs[psi_iter]," max\n");
			} else if(psi_psi2p > 0 && psi_psi2p_1 < 0) {
				psi_min = c(psi_min,psi_xs[psi_iter-1],psi_xs[psi_iter]);
				#cat("psi_intervals = ",psi_xs[psi_iter-1]," -> ",psi_xs[psi_iter]," min\n");
			}
		}
	}
	
	plot(psi_xs,psi_psis,main="The Wave Function",
			xlab="x (nm)", ylab="psi(x)");
	grid();
#}
# END WAVE EQUATION PLOT
###############################################################################
# PSI * PSI PLOT
psi_psis2 = c();
for(i in 1:length(psi_xs)) {
	psi_psis2 = c(psi_psis2,psi_psis[i]^2);
}
plot(psi_xs,psi_psis2,main="PSI * PSI\nThe Probability Distribution",
		xlab="x (nm)", ylab="|psi(x)|^2");
grid();
# END PSI * PSI PLOT
###############################################################################
# PSI * PSI: MAXIMUM AND MINIMUM PROBABILITIES
cat("Root Search of PSI * PSI ****************************************\n");
cat("\tMaxima: mode =",psi_mode,"\n");

# newton-raphson search method
if(search_method == 1) {
	for(i in 1:(length(psi_max)/2)) {
		psimm_min = psi_max[2*i-1];
		psimm_max = psi_max[2*i];
		psimm_x = psimm_min;
		psimm_psi2p_x = psi2p(psimm_x,L[psi_mode],kappa[psi_mode],a);
		
		while(abs(psi2p(psimm_x,L[psi_mode],kappa[psi_mode],a)) > root_dx) {
			psimm_nr_dx =  nr_dx_calculator(psi(psimm_x,L[psi_mode],kappa[psi_mode],a)^2,psi2p(psimm_x,L[psi_mode],kappa[psi_mode],a));
			if(psimm_x+psimm_nr_dx < psimm_max) {
				psimm_x = psimm_x + psimm_nr_dx;
			} else {
				psimm_x = psimm_max;
			}
		}
		cat("\t\tpsi * psi max = ",psi(psimm_x,L[psi_mode],kappa[psi_mode],a)^2,"\tx = ",psimm_x,"\n");
	}
	cat("\tMinima: mode =",psi_mode,"\n");
	if(length(psi_min) > 1)
	for(i in 1:(length(psi_min)/2)) {
		psimm_min = psi_min[2*i-1];
		psimm_max = psi_min[2*i];
		psimm_x = psimm_min;
		psimm_psi2p_x = psi2p(psimm_x,L[psi_mode],kappa[psi_mode],a);
		
		while(abs(psi2p(psimm_x,L[psi_mode],kappa[psi_mode],a)) > root_dx) {
			psimm_nr_dx =  nr_dx_calculator(psi(psimm_x,L[psi_mode],kappa[psi_mode],a)^2,psi2p(psimm_x,L[psi_mode],kappa[psi_mode],a));
			if(psimm_x+psimm_nr_dx < psimm_max) {
				psimm_x = psimm_x + psimm_nr_dx;
			} else {
				psimm_x = psimm_max;
			}
		}
	}
}
# end newton-raphson search method

# bisection search method
if(search_method == 0) {
	for(i in 1:(length(psi_max)/2)) {
		psimm_xlo = psi_max[2*i-1];
		psimm_xhi = psi_max[2*i];
		psimm_xmid = mid(psimm_xlo,psimm_xhi);
		psimm_x = psimm_xmid;
		psimm_psi2p_xhi = psi2p(psimm_xhi,L[psi_mode],kappa[psi_mode],a);
		
		while(abs(psi2p(psimm_xhi,L[psi_mode],kappa[psi_mode],a)) > root_dx) {
			if(psi2p(psimm_xlo,L[psi_mode],kappa[psi_mode],a) * psi2p(psimm_xmid,L[psi_mode],kappa[psi_mode],a) < 0) {
				psimm_xhi = psimm_xmid;
				psimm_xmid = mid(psimm_xlo,psimm_xhi);
				psimm_x = psimm_xmid;
			} else if(psi2p(psimm_xhi,L[psi_mode],kappa[psi_mode],a) * psi2p(psimm_xmid,L[psi_mode],kappa[psi_mode],a) < 0) {
				psimm_xlo = psimm_xmid;
				psimm_xmid = mid(psimm_xlo,psimm_xhi);
				psimm_x = psimm_xmid;
			} else if(psi2p(psimm_xhi,L[psi_mode],kappa[psi_mode],a) == 0) {
				psimm_x = psimm_xhi;
				break;
			} else if(psi2p(psimm_xmid,L[psi_mode],kappa[psi_mode],a) == 0) {
				psimm_x = psimm_xmid;
				break;
			} else if(psi2p(psimm_xlo,L[psi_mode],kappa[psi_mode],a) == 0) {
				psimm_x = psimm_xlo;
				break;
			}
		}
		cat("\t\tpsi * psi max = ",psi(psimm_x,L[psi_mode],kappa[psi_mode],a)^2,"\tx = ",psimm_x,"\n");
	}
	cat("\tMinima: mode =",psi_mode,"\n");
	if(length(psi_min) > 1)
	for(i in 1:(length(psi_min)/2)) {
		psimm_xlo = psi_min[2*i-1];
		psimm_xhi = psi_min[2*i];
		psimm_xmid = mid(psimm_xlo,psimm_xhi);
		psimm_x = psimm_xmid;
		psimm_psi2p_xhi = psi2p(psimm_xhi,L[psi_mode],kappa[psi_mode],a);
		
		while(abs(psi2p(psimm_xhi,L[psi_mode],kappa[psi_mode],a)) > root_dx) {
			if(psi2p(psimm_xlo,L[psi_mode],kappa[psi_mode],a) * psi2p(psimm_xmid,L[psi_mode],kappa[psi_mode],a) < 0) {
				psimm_xhi = psimm_xmid;
				psimm_xmid = mid(psimm_xlo,psimm_xhi);
				psimm_x = psimm_xmid;
			} else if(psi2p(psimm_xhi,L[psi_mode],kappa[psi_mode],a) * psi2p(psimm_xmid,L[psi_mode],kappa[psi_mode],a) < 0) {
				psimm_xlo = psimm_xmid;
				psimm_xmid = mid(psimm_xlo,psimm_xhi);
				psimm_x = psimm_xmid;
			} else if(psi2p(psimm_xhi,L[psi_mode],kappa[psi_mode],a) == 0) {
				psimm_x = psimm_xhi;
				break;
			} else if(psi2p(psimm_xmid,L[psi_mode],kappa[psi_mode],a) == 0) {
				psimm_x = psimm_xmid;
				break;
			} else if(psi2p(psimm_xlo,L[psi_mode],kappa[psi_mode],a) == 0) {
				psimm_x = psimm_xlo;
				break;
			}
		}
		cat("\t\tpsi * psi min = ",psi(psimm_x,L[psi_mode],kappa[psi_mode],a)^2,"\tx = ",psimm_x,"\n");
	}
}
# end bisection search method
# END PSI * PSI: MAXIMUM AND MINIMUM PROBABILITIES
###############################################################################
# DISPLAY STATS
# root plot stats
	cat("Root Search Plot Statistics ********************************************\n");
	cat("\tdx avg =\t",root_dxavg," +/- ",root_dxstd,"\n");
	cat("\tdist avg =\t",root_distavg," +/- ",root_diststd,"\n");
	cat("\tmaxdist =\t",maxdist,"\n");
	cat("\t# points =\t",root_ptcounter,"\n");
# end root plot stats

# bisection stats
	search_method_name = "Search Method";
	if(search_method == 0) {
		search_method_name = "Bisection ";
	} else if(search_method == 1) {
		search_method_name = "Newton-Raphson "
	}

	cat(search_method_name,"Statistics ***************************************\n")
	cat("\t",search_method_name,"iteration avg =\t",root_iteravg,"+/-",root_iterstd,"\n");
	cat("\t",search_method_name,"total iterations = ",root_itersum,"\n");
	cat("\t",search_method_name,"max iterations = ",root_itermax,"\n");
	cat("\t",search_method_name,"min iterations = ",root_itermin,"\n");
	cat("\t",search_method_name,"avg interval   = ",root_intervalsum/root_itercount,"\n");
	cat("\t# of energy levels       = ",length(intervals)/2,"\n");
# end bisection stats
# END DISPLAY STATS
cat(Sys.time()-start,"s\n");
###############################################################################
# TEST
#cat("test\n");
#cat("K = ",K,"\na = ",a,"\n");
#cat("f (",.1,") = ",f(.1,K,a),"\n");
#cat("f'(",.1,") = ",fp(.1,K,a),"\n");
