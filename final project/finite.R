###############################################################################
# Author: 		Adam J Jaso, Jr
# Title:		Calculation of the Energy Levels of the Infinite Square Well in 1 Dimension
# Date:			12/4/2011
# Class:		PHYS 140 Computational Methods for Physicists
# Instructor:	Dr. Olencka Hubickyj
###############################################################################
# BASIC APPROACH
# This program accepts the potential and the width of the square well, and
# 	from these calculates the energy levels from the zeros of the equation
#	shown below by the bisection method.
#		f(x) = sqrt(K) * ( sqrt(1-x) * tan( sqrt(1-x) * a ) - sqrt(x))
#		V0 = potential
#		a = well width
#		K = 2*m*V0/hbar^2
#		x = E_n / V0
#
###############################################################################
#CONSTRAINTS
#	Given values for the current situation
	V0 = 75; #eV
	a = .001;#0.05; #nm
#END CONSTRAINTS


##############################################################################
#CONTROL VARIABLES
#	ALGORITHM
		algorithm=2;
#		Do nothing 		-1
#		Newton-Raphson	0
#		Brute Force		1
#		Bisection		2
			dx_bisection = 0.000000001;
			outputprecision = 9;
#	END ALGORITHM
#	PLOT
		plot=1;
#		on		1
#		off 	0

#		the ceiling for all displayed points, that is, f(x) < fceil
		fceil = 100;
#		maximum distance allowed between two graphed points
		maxdist=fceil/100;
#		default distance between dx increments of the graphed points
#			DO NOT ADJUST THIS RATIO UNLESS ABSOLUTELY NECESSARY!
		defdx = 1/10000;
#	END PLOT
#END CONTROL VARIABLES


###############################################################################
#FUNCTION DEFINITIONS
#
#	the function whose zeros represent energy levels of the finite square well in 1-D
f <- function(x,K,a) sqrt(K)*(sqrt(1-x)*tan(sqrt(1-x)*a) - sqrt(x));

#	1st derivative of the above function
#fp <- function(x,K,a) -1/2*sqrt(K)*(sqrt(1-x)*tan(sqrt(1-x)*a) + a*1/cos(sqrt(1-x)*a)^2 + 1/sqrt(x));
fp <- function(x,K,a) -1/2*sqrt(K)*(sqrt(1-x)*tan(sqrt(1-x)*a) + a*1/cos(sqrt(1-x)*a)^2 + 1/sqrt(x));

#	2nd derivative of the above function
fpp <- function(x,K,a) 
{
	z = 1-x;
	za = sqrt(z)*a;
	A = -1/2*1/sqrt(z)*tan(za) + sqrt(z)*1/cos(za)^2*(-1/2*a*1/sqrt(z));
	B = 2*a*1/cos(za)*(1/cos(za)*tan(za))*(-1/2*a*1/sqrt(z));
	C = -1/2*1/sqrt(x)^3;
	
	-1/2*sqrt(K)*(A+B+C);
}
#E <- function(x,V0) x*V0;
#END FUNCTION DEFINITIONS

#COLLECTED CONSTANTS
#	calculation constants
mc2 = 100*0.511E6; #eV
hbarc = 197.3; #eV . nm
K=2*mc2*V0/hbarc^2;
#END COLLECTED CONSTANTS

###############################################################################
#PLOT
#	Generate the plot of the energies as percentages of x
#	Plot x vs f(x) on 0 <= x <= 1

#	PLOT CONSTRAINT VARIABLES
#		the ceiling for all displayed points, that is, f(x) < fceil
		#fceil = 100;
#		maximum distance allowed between two graphed points
		#maxdist=2;
#		default distance between dx increments of the graphed points
#			DO NOT ADJUST THIS RATIO UNLESS ABSOLUTELY NECESSARY!
		#defdx = maxdist/10000;
#	END PLOT CONSTRAINT VARIABLES


#	FUNCTIONS
#		function to find the distance between two points
dist <- function(x1,fx1,x2,fx2) sqrt((x1-x2)^2 + (fx1-fx2)^2);
#		function to calculate the approximate next increment to maintain a distance, maxdist, between all points
dx <- function(dp,yp) dp / sqrt(yp^2 + 1);
#	END FUNCTIONS


#	INITIAL VALUE VARIABLES
#		starting plot points
x=0;
xs = c();
fxs = c();
domains = c();
fofx_1 = f(x,K,a);
fpofx_1 = fp(x,K,a);
fppofx_1 = fpp(x,K,a);
#	END INITIAL VALUE VARIABLES


#	STAT VARIABLES
#		for dx stats
dxsum = 0;
dx2sum = 0;
#		for distance stats (distance between graphed points)
distsum = 0;
dist2sum = 0;
#	END STAT VARIABLES


#	COUNTER VARIABLES
#		iteration counter
iter = 0;
#		plot point counter
ptcounter = 0;
droppeddist = 0;
#		domain pair counter;
domaincounter = 0;
#	END COUNTER VARIABLES


#	GENERATE PLOT POINTS
while(x < 1-.01) {
	#	calculate the function value and derivate at the point, x
	fofx = f(x,K,a);
	fpofx = fp(x,K,a);
	fppofx = fpp(x,K,a);
	dxcalc = dx(maxdist,fpofx);
	
	#	if the absolute value of the function is less than the ceiling value, 
	#		store this point for plotting
	if(abs(fofx) < fceil) {
		#	add value to the plot points
		xs = c(xs,x);
		fxs = c(fxs,fofx);
		#	count this plot point
		ptcounter = ptcounter + 1;
		
		#	if the 1st derivative (slope) isn't infinite (or at least very large), 
		#		calculate the next dx increment and the distance between the two 
		#		most recently added points
		if(abs(fpofx) != Inf) {
			#x	dxcalculated	dist
			
			cat("x = ",x,"\tdxcalc = ",dxcalc);
			#	Calculate the distance between the current point and the previous point,
			#		excluding points greater than one unit above the set max distance,
			#		and summing the distances, counting dropped points above the max distance,
			#		finally displaying the distance.
			distance = dist(xs[ptcounter-1],f(xs[ptcounter-1],K,a),xs[ptcounter],f(xs[ptcounter],K,a));
			
			#	If the point counter is greater than 1, and the distance is less than the maxdist
			#		plus some buffer, then check the distance between the two most recent points.
			#		the purpose of the buffer is to exclude huge distances that occur at the asymptotes.
			if(ptcounter > 1 && distance < maxdist+1) {
				cat("\tdist = ",distance,"\n");
				#	include this distance in the distance sum
				distsum = distsum + distance;
				dist2sum = dist2sum + distance^2;
				#	check the 2nd derivative for concavity change (+) to (-)
				if(fofx_1 > 0 && fofx < 0) {
					domains = c(domains,xs[ptcounter-1],xs[ptcounter]);
				}
			} else {
				cat("\n");
				#	drop this distance from the average
				droppeddist = droppeddist+1;
			}
			
			#	generate next x from the calculated differential correction
			x=x+dxcalc;
			
			#	add this dx to the total of dx and dx^2
			dxsum = dxsum + dxcalc;
			dx2sum = dx2sum + dxcalc^2;
		} else {
			#	generate next x from the default differential correction specified above
			#		this is needed to get past the asymptotes
			x = x + defdx;
			
			#	add this dx to the total of dx and dx^2
			dxsum = dxsum + defdx;
			dx2sum = dx2sum + defdx^2;
		}
	} else {
		#	generate next x from the default differential correction specified above
		#		this is needed to get past the asymptotes
		x = x + defdx;
		
		#	add this dx to the total of dx and dx^2
		dxsum = dxsum + defdx;
		dx2sum = dx2sum + defdx^2;
	}
	#	store fppofx in fppofx_1
	fofx_1 = fofx;
	#	count iteration
	iter = iter+1;
}
#	END GENERATE PLOT POINTS

#	CALCULATE AND DISPLAY DX AND DISTANCE STATS
#		Stats on the dx's

dxavg = dxsum / iter;
dxstd = sqrt(dx2sum / iter - (dxsum / iter)^2);
distavg = distsum / (ptcounter-droppeddist-1);
diststd = sqrt(dist2sum / (ptcounter-droppeddist) - (distsum / (ptcounter-droppeddist))^2);
#	cat("avg dx = ",dxavg," +/- ",dxstd,"\n");
#	cat("avg dist = ",distavg,"+/-",diststd,"\n");
#	cat("dist2sum = ",dist2sum,"\tdistsum = ",distsum,"\tptcounter = ",ptcounter,"\tdroppeddist = ",droppeddist,"\n");
#	END CALCULATE AND DISPLAY DX AND DISTANCE STATS

plot(xs,fxs);
grid();

#END PLOT


##############################################################################
#ALGORITHMS

#Print the number of pi's inside the tan
#cat("k*a / pi = ",sqrt(K*(1-x))*a / pi,"\n");

#Find the intervals on which to search


#Newton-Raphson Algorithm
if(algorithm == 0) {
	
	x = .5;
	dx = -f(x,K,a)/fp(x,K,a);
	
	while(abs(f(x,K,a)) > 0.00001) {
		x = x + dx;
		dx = -f(x,K,a)/fp(x,K,a);
	}
}

#Brute Force Algorithm
if(algorithm == 1) {
	
	x = 0.002491;
	dx = .0000000000005;
	min = 1;
	count = 0;
	
	while(x < 1-dx && count < 10) {
		if(abs(f(x,K,a)) < abs(f(min,K,a))) {
			min = x;
			cat("x = ",x,"\tf(x) = ",f(x,K,a),"\tE = ",x*V0,"\n");
		} else {
			count = count + 1;
		}
		x=x+dx;
		cat("x = ",x,"\tf(x) = ",f(x,K,a),"\tE = ",x*V0,"\n");
	}
	
}

#Bisection Algorithm
if(algorithm == 2) {
	mid <- function(lo,hi) (lo + hi) / 2;
	
	itersum = 0;
	iter2sum = 0;
	iterstd = 0;
	iteravg = 0;
	itercount = 0;
	itermax = 0;
	itermin = 1000000000;
	
	#	take each domain and bisect it
	for(i in 1:(length(domains)/2)) {
		
		#	set the range to bisect
		xlo = domains[2*i-1];
		xhi = domains[2*i];
		#cat("xlo = ",xlo,"\tf(xlo) = ",f(xlo,K,a),"\txhi = ",xhi,"\tf(xhi) = ",f(xhi,K,a),"\n");
		xmid = mid(xlo,xhi);
		x = xmid;
		#dx_bisection = 0.00001;
		iter = 0;
		
		#	error management
		x_1 = x;
		errorcount = 0;
		maxerrorcount = 10;
		errorlist = c();
		erroryes = 0;
		#	yes 	1
		#	no		0
		
		while(abs(f(x,K,a)) > dx_bisection) {
			if(f(xlo,K,a) * f(xmid,K,a) < 0) {
				xhi = xmid;
				xmid = mid(xlo,xhi);
				x = xmid;
			} else if(f(xhi,K,a) * f(xmid,K,a) < 0) {
				xlo = xmid;
				xmid = mid(xlo,xhi);
				x = xmid;
			} else if(f(xhi,K,a) == 0) {
				x = xhi;
				break;
			} else if(f(xmid,K,a) == 0) {
				x = xmid;
				break;
			} else if(f(xlo,K,a) == 0) {
				x = xlo;
				break;
			}
			iter = iter+1;
			
			if(x_1 == x) {
				errorcount = errorcount + 1;
				if(errorcount > maxerrorcount) {
					errorlist = c(errorlist,i);
					erroryes = 1;
					errorcount = 0;
					break;
				}
			}
			x_1 = x;
			cat("f(x) = ",f(x,K,a),"\n");
			
		}
		#Bisection Output
		cat("Level ",length(domains)/2-i,"\tx  = ",format(x,digits=outputprecision),"\tf(x) = ",f(x,K,a),"\tE = ",format(-x*V0,digits=outputprecision)," +/- ",V0*dx_bisection," eV\titer = ",iter,"\tdx_bisection = ",format(dx_bisection,digits=outputprecision));
		if(erroryes == 1) {
			cat("\tWARNING: THE ZERO MAY NOT HAVE BEEN FOUND EXACTLY!\n");
			erroryes = 0;
		} else {
			cat("\n");
		}
		
		#Iteration statistics
		itersum = itersum + iter;
		iter2sum = iter2sum + iter^2;
		itercount = itercount + 1;
		
		if(iter > itermax) {
			itermax = iter;
		}
		if(iter < itermin) {
			itermin = iter;
		}
	}
	
	if(itercount != 0) {
		iteravg = itersum / itercount;
		iterstd = sqrt(iter2sum / itercount - iteravg^2);
	}
}

#END ALGORITHMS

#SHOW PLOT STATS
cat("avg dx = ",dxavg," +/- ",dxstd,"\n");
cat("avg dist = ",distavg," +/- ",diststd,"\n");
cat("avg iteration = ",iteravg," +/- ",iterstd,"\n");
cat("itermin = ",itermin,"\titermax = ",itermax,"\n");

#cat("K = ",K,"\na = ",a,"\nf  = ",f(.01,K,a),"\nfp = ",fp(.01,K,a),"\nfpp=",fpp(.01,K,a),"\n");