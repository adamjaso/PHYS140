function energies=finitesquarewell(V0,a)
%finitesquarewell.m
%Author:    Adam J. Jaso, Jr.
%Date:      12/05/2011
%Funcion:   Calculation of the Energies of the finitesquarewell Square Well

%*************************************************************************
%   CONSTANTS

        %V0 = 40;
        %a = .05;
        mc2 = 0.511E6; %eV
        hbarc = 197.3; %eV . nm
        qc=2*mc2/hbarc^2;
        K=qc*V0;

%   END CONSTANTS
%*************************************************************************
%   CONTROL VARIABLES
    
    %Solution Search Method
    algorithm = 2;
		%Newton-Raphson Method  1
		%Bisection Method       2
    
    %Solution Tolerance Settings
    	%Bisection Tolerance
        dx_bis = 0.0001;
        %Newton-Raphson Tolerance
        dx_nr = 0.0001;
    %Plot Settings
        %the ceiling for all displayed points, that is, f(x) < fceil
		fceil = 100;
        %maximum distance allowed between two graphed points
		maxdist=fceil/100;
        %default distance between dx increments of the graphed points
            %DO NOT ADJUST THIS RATIO UNLESS ABSOLUTELY NECESSARY!
        defdx=1/10000;
%   END CONTROL VARIABLES
%*************************************************************************
%   GENERATE PLOT & FIND SEARCH INTERVALS

%	INITALIZE LOOP VARIABLES
		x=0;
		xs=[];
		fxs=[];
		domains=[];
		fofx_1 = f(x,K,a);
		fpofx_1 = fp(x,K,a);
%	END INITIALIZE LOOP VARIABLES

%	STAT VARIABLES
	%	for differentials between points
			dxsum = 0;
			dx2sum = 0;
	%	for distance stats (distance between graphed points)
			distsum = 0;
			dist2sum = 0;
%	END STAT VARIABLES

%	COUNTER VARIABLES
	%	iteration counter
			iter = 0;
	%	plot point counter
			ptcounter = 0;
			droppeddist = 0;
	%	domain pair counter;
			domaincounter = 0;
%	END COUNTER VARIABLES

%	GENERATE PLOT POINTS
while x < 1
	%	calculate the function value and derivate at the point, x
    fofx = f(x,K,a);
    fpofx = fp(x,K,a);
    %fppofx = fpp(x,K,a);
    dxcalc = dx_calculator(maxdist,fpofx);
	
	%	if the absolute value of the function is less than the ceiling value, 
	%		store this point for plotting
    fprintf('fofx = %f\n',fofx);
	if abs(fofx) < fceil
		%	add value to the plot points
		xs = [xs x];
		fxs = [fxs fofx];
        
		%	count this plot point
		ptcounter = ptcounter + 1;
		
		%	if the 1st derivative (slope) isn't infinitesquarewell (or at least very large), 
		%		calculate the next dx increment and the distance between the two 
		%		most recently added points
        %fprintf('fpofx = %f\n',fpofx);
		%if abs(fofx) < 100
			%	Display the current x and differential
			fprintf('x = %3.8f\tdxcalc = %1.8f',x,dxcalc);
			
			%	If the point counter is greater than 1, and the distance is less than the maxdist
			%		plus some buffer, then check the distance between the two most recent points.
			%		the purpose of the buffer is to exclude huge distances that occur at the asymptotes.
			if ptcounter > 1
				%	Calculate the distance between the current point and the previous point,
				%		excluding points greater than one unit above the set max distance,
				%		and summing the distances, counting dropped points above the max distance,
				%		finally displaying the distance.
				distance = dist([xs(ptcounter-1) f(xs(ptcounter-1),K,a)], [xs(ptcounter) f(xs(ptcounter),K,a)]);
				
                if distance < maxdist+1
					fprintf('\tdist = %1.8f\n',distance);
					%	include this distance in the distance sum
					distsum = distsum + distance;
					dist2sum = dist2sum + distance^2;
					%	check the 2nd derivative for concavity change (+) to (-)
					if fofx_1 > 0 && fofx < 0
						domains = [domains [xs(ptcounter-1) xs(ptcounter)]];
					end
				end
			else
				fprintf('\n');
				%	drop this distance from the average
				droppeddist = droppeddist+1;
            end

            %	generate next x from the calculated differential correction
            x=x+dxcalc;

            %	add this dx to the total of dx and dx^2
            dxsum = dxsum + dxcalc;
            dx2sum = dx2sum + dxcalc^2;
        %else
            %	generate next x from the default differential correction specified above
            %		this is needed to get past the asymptotes
            x = x + defdx;

            %	add this dx to the total of dx and dx^2
            dxsum = dxsum + defdx;
            dx2sum = dx2sum + defdx^2;
            
    else
        x = x + defdx;
    end
end
    
	%	store fppofx in fppofx_1
	fofx_1 = fofx;
	%	count iteration
	iter = iter+1;
end

%	CALCULATE AND DISPLAY DX AND DISTANCE STATS
%		Stats on the dx's

dxavg = dxsum / iter;
dxstd = sqrt(dx2sum / iter - (dxsum / iter)^2);
distavg = distsum / (ptcounter-droppeddist-1);
diststd = sqrt(dist2sum / (ptcounter-droppeddist) - (distsum / (ptcounter-droppeddist))^2);
%	cat("avg dx = ",dxavg," +/- ",dxstd,"\n");
%	cat("avg dist = ",distavg,"+/-",diststd,"\n");
%	cat("dist2sum = ",dist2sum,"\tdistsum = ",distsum,"\tptcounter = ",ptcounter,"\tdroppeddist = ",droppeddist,"\n");
%	END CALCULATE AND DISPLAY DX AND DISTANCE STATS

plot(xs,fxs);
grid;
%	END GENERATE PLOT POINTS

%   END GENERATE PLOT & FIND SEARCH INTERVALS
%*************************************************************************
%   EQUATION SOLTUION ALGORITHMS



%   END EQUATION SOLUTION ALGORITHMS
%*************************************************************************
%   DISPLAY SUMMARY STATS



%   END DISPLAY SUMMARY STATS
end
%   END OF MAIN FUNCTION: finitesquarewell.M





%*************************************************************************
%   FUNCTION DEFINITIONS

%   SEARCH FUNCTION
function f=f(x,K,a)
    f=sqrt(K)*(sqrt(1-x)*tan(sqrt(K*(1-x))*a) - sqrt(x));
end

%   1ST DERIVATIVE OF SEARCH FUNCTION
function fp=fp(x,K,a)
    fp=-1/2*sqrt(K)*(sqrt(1-x)*tan(sqrt(1-x)*a) + a*1/cos(sqrt(1-x)*a)^2 + 1/sqrt(x));
end

%   2ND DERIVATIVE OF SEARCH FUNCTION
function fpp=fpp(x,K,a)
    z = 1-x;
	za = sqrt(z)*a;
	A = -1/2*1/sqrt(z)*tan(za) + sqrt(z)*1/cos(za)^2*(-1/2*a*1/sqrt(z));
	B = 2*a*1/cos(za)*(1/cos(za)*tan(za))*(-1/2*a*1/sqrt(z));
	C = -1/2*1/sqrt(x)^3;
	
	fpp=-1/2*sqrt(K)*(A+B+C);
end

%   KAPPA CONSTANT CALCULATION
function kappa=k(E,qc)
    kappa=sqrt(qc)*sqrt(E);
end

%   L CONSTANT CALCULATION
function el=el(E,qc,V0)
    el=sqrt(qc)*sqrt(V0-E);
end

%   ENERGY LEVEL CONSTANT CALCULATION FROM x SOLUTION
function E=E(V0,x)
    E=V0*x;
end

%   NORMALIZATION CONSTANTS
function const=normconst(x,el,kappa,a)
    D=(a + 1/kappa * ( sin( el * a) * cos( el * a) + exp(2 * kappa * a) * cos(el * a) ^ 2))^(-1/2);
    if x < -a
        const=D*exp(-2 * kappa * a) * cos( el * a);
    elseif x >= -a && x <= a
        const=D;
    elseif x > a
        const=D*exp(-2 * kappa * a) * cos( el * a);
    else
        const=0;
    end
end

%   DX CALCULATOR
function dx=dx_calculator(dp,yp) 
    dx=dp / sqrt(yp^2 + 1);
end

%   POINT DISTANCE CALCULATOR
function dist=dist(p1,p2) 
    dist=sqrt( ( p1(1) - p2(1) )^2 + ( p1(2) - p2(2) )^2 );
end

%   THE WAVE FUNCTION
function wavfunc=wavf(x,el,kappa,a)
    if x < -a
        wavfunc=normconst(x,el,kappa,a) * exp(kappa * x);
    elseif x >= -a && x <= a
        wavfunc=normconst(x,el,kappa,a) * cos(el * x);
    elseif x > a
        wavfunc=normconst(x,el,kappa,a) * exp(-kappa * x);
    else
        wavfunc=0;
    end
end


