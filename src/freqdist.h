/*
 Ckrige, a program for universal, constrained and covariance-matching
 constrained point and block kriging.
 Copyright 2010 (C) Christoph Hofer
 
 Christoph Hofer christoph.hofer@alumni.ethz.ch
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version. As a special exception, linking 
 this program with the Qt library is permitted.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
// Headerfile Definition der Distanz-Dichtefunktionen zwischen
// einem Punkt und einem Rechteck

//
// f_int_boundaries berechent die Koordinaten des Pixels für den rechten
// oberen Quadraten sowie die Grenzen des Integrals
//

// Author: Christoph Hofer
// Date: Dezember 2009

#include <stdio.h>
#include <R.h>



//*******************************************
// f_dist_freq Berechnung der Distanzdichte
// zwischen einem Punkt und einem Rechteck
double f_dist_freq( double lag_dist, double *p_sc_pixcoords);
//*******************************************


double f_dist_freq( double lag_dist, double *p_sc_pixcoords) 
{ 
	double lag_dist_sq, a, a_sq, b, b_sq, c, c_sq, d, d_sq, K;
	double dist_freq = 0;
	//double *psc_v;
	//psc_v = (double*) modparam; // cast void pointer in double pointer
	 //// printf("a %f \n", *(psc_v - 1));
	 //// printf("b %f \n", *(psc_v - 2));
	 //// printf("c %f \n", *(psc_v - 3));
	 //// printf("d %f \n", *(psc_v - 4));
	
	lag_dist_sq = pow(lag_dist,2);
	
	a = *(p_sc_pixcoords); // -4 Adresse des Wertes von a bei 0 und grösser sind die Parameter der Kovarianzfunktion gespeichert
	a_sq= pow(a, 2);
	
	b =  *(p_sc_pixcoords + 1);
	b_sq= pow(b, 2);
	
	c =  *(p_sc_pixcoords + 2);
	c_sq= pow(c, 2);
	
	d =  *(p_sc_pixcoords + 3);
	d_sq= pow(d, 2);
	
	
	/*  printf("a %f \n", a);
	  printf("b %f \n", b);
	  printf("c %f \n", c);
	  printf("d %f \n", d);
	 printf("a2 %f \n", a_sq);
	 printf("b2 %f \n", b_sq);
	 printf("c2 %f \n", c_sq);
	 printf("d2 %f \n", d_sq);
	
	 printf("lag dist %f \n", lag_dist);
	 printf("lag dist squared %f \n", lag_dist_sq);
	*/
	
	
	// Konstante
	
	K = ( 1 / (b - a) ) * ( 1 / (d - c) );
	
	// printf(" dist %f ", lag_dist);
	
	//case I (kommt nie vor)
	
	if( ( (lag_dist_sq) > (b_sq + c_sq) ) && ( (lag_dist_sq) < (a_sq + d_sq) ) && ( (a_sq) <= (b_sq) ) )
	{
		//printf("case I");
		
		dist_freq = K * lag_dist * ( atan( b / sqrt( lag_dist_sq - b_sq ) ) - atan( a / sqrt( lag_dist_sq - a_sq ) ) );
	
	}
	
	//case II 
	else if( ( (lag_dist_sq) > (b_sq + c_sq) ) && ( (lag_dist_sq) > (a_sq + d_sq) ) && ( (lag_dist_sq) <= (b_sq + d_sq) ) )
	{
		// printf("case II");
		
		dist_freq = K * lag_dist * ( atan( b / sqrt( lag_dist_sq - b_sq ) ) - atan( sqrt( lag_dist_sq - d_sq ) / d ) );
	
	}
	
	//case III
	else if( ( (lag_dist_sq) < (b_sq + c_sq) ) && ( (lag_dist_sq) < (a_sq + d_sq) ) && ( (lag_dist_sq) >= (a_sq + c_sq) ) )
	{	
		// printf("case III");
		
		dist_freq = K * lag_dist * ( atan( sqrt( lag_dist_sq - c_sq ) / c ) - atan( a / sqrt( lag_dist_sq - a_sq ) ) );
	
	}
	
	//case IV
	else if( ( (lag_dist_sq) < (c_sq + b_sq) ) && ( lag_dist_sq > ( a_sq + d_sq ) ) && ( ( d_sq ) >= ( c_sq ) ) )
	{	
		//printf("case IV");
		
		dist_freq = K * lag_dist * ( atan(  sqrt( lag_dist_sq - c_sq ) / c )  - atan( sqrt(  lag_dist_sq - d_sq )  / d ) );
	}


	return(dist_freq);
	
} 


