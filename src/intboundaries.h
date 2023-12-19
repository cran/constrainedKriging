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
 */// Headerfile Definition der Integrationsgrenzen zwischen
// einem Punkt und einem Rechteck


// Author: Christoph Hofer
// Date: Dezember 2009
// 2023-01-24 A. Papritz elimination of unnneeded variable sc_v

#include <stdio.h>
#include<math.h>

void f_int_boundaries(
				 void *modparam, // Pointer zeigt auf Array mit den skalierten Pixelcoordinaten
				 double *plb,  // pointer auf Array[4] mit den unteren Integralgrenzen 
				 double *pub, // pointer auf Array[4] mit den unteren Integralgrenzen
				 double *ppa); // pointer auf Array[4] mit der Fläche des Pixels 
				 

// Berechung der minimalen und maximalen euklidischen Distanz (Integrationsgrenzen)
// zwischen dem Punkt c(0,0) und dem skalierten Pixel

void f_min_max_euclidian_distance(double *psc_v, double *plb, double *pub);



void f_int_boundaries(
						void *modparam,
						double *plb,   
						double *pub,
					    double *ppa)
{

	/* 
	 * int i;
	 */
	double sc_a, sc_b, sc_c, sc_d;
	
	double *psc_v;
	psc_v = (double*) modparam; // cast void pointer in double pointer
		
	sc_a = *(psc_v + 0);
	sc_b = *(psc_v + 1);
	sc_c = *(psc_v + 2);
	sc_d = *(psc_v + 3);
	

	

	
	//  printf("sc_a %f \n", sc_a);
	//  printf("sc_c %f \n", sc_b);
	//  printf("sc_d %f \n", sc_c);
	// printf("sc_d %f \n", sc_d);
	
	
	// .....................................................................
	// Der Punkt ist ausserhalb des Pixels oder liegt auf einem seiner Ecken
	// .....................................................................
	
	
	// FALL I: px <= a und py <= c
	
	if( (sc_a >= 0) && (sc_c >= 0) ){
		
		//printf("\n \n vertex I  \n");
		
		*(psc_v + 0) = sc_a;
		*(psc_v + 1) = sc_b;
		*(psc_v + 2) = sc_c;
		*(psc_v + 3) = sc_d;
		
		
		
		f_min_max_euclidian_distance(psc_v, plb, pub);
	
		*ppa = ( *(psc_v + 1) - *(psc_v + 0) ) * ( *(psc_v + 3) - *(psc_v + 2) );
	}
	
	// FALL II: px <= a und py >= d

	if( (sc_a >= 0) && (sc_d <= 0) ){
		
		//printf("\n \n vertex II  \n");
		
		*(psc_v + 0) = sc_a;
		*(psc_v +  1) = sc_b;
		*(psc_v +2) = fabs(-1* sc_d); // c new = -1* d old
		*(psc_v + 3) = fabs(-1* sc_c); // d new = -1* c old
		
	
		
		
		
		// printf("\n \n vertex II \n \n \n");
		
		*ppa = ( *(psc_v + 1) - *(psc_v + 0)) * ( *(psc_v + 3)  - *(psc_v + 2) ) ;
		
		f_min_max_euclidian_distance(psc_v, plb, pub);
	}
	
	// FALL III: b <= px und py <=  c

	if( (sc_b <= 0) && (sc_c >= 0) ){
		
		//printf("\n \n vertex III \n");
		
		*(psc_v + 0) = fabs(-1 * sc_b); // a new = -1* b old
		*(psc_v + 1) = fabs(-1 * sc_a) ; // b new = -1* a old
		*(psc_v + 2) = sc_c;
		*(psc_v + 3) = sc_d; 
		
		
	
		
		// printf(" \n \n vertex III \n \n \n");
		
		
		
		*ppa = ( *(psc_v +1) - *(psc_v +0) )* ( *(psc_v +3) - *(psc_v +2) );
		
		
		
		f_min_max_euclidian_distance(psc_v, plb, pub);
		
	}
	
	// FALL IV: b <= px und py >=  d

	if( (sc_b <=  0) && (sc_d <= 0) ){
		
		//printf("\n \n vertex IV \n");
		
		*(psc_v + 0) =  fabs(-1 * sc_b); // a new = -1* b old
		*(psc_v + 1) =  fabs(-1 * sc_a); // b new = -1* a old
		*(psc_v + 2) =  fabs(-1 * sc_d); // c new = -1* d old
		*(psc_v + 3) =  fabs(-1 * sc_c);   // d new = -1 * c old
		
		
		//printf("sc_a %f \n",*(psc_v + 0));
		//printf("sc_b %f \n",*(psc_v + 1));
		//printf("sc_c %f \n",*(psc_v + 2));
		//printf("sc_d %f \n",*(psc_v + 3));
		
		*ppa = (*(psc_v + 1) - *(psc_v + 0))*(*(psc_v + 3) - *(psc_v + 2));
		
		f_min_max_euclidian_distance(psc_v, plb, pub);
	
		
		//printf("plb %f \n",*plb);
		//printf("pub %f \n",*pub);
	}
	
	// .....................................................................
	// Der Punkt ist innerhalb des Pixels 
	// .....................................................................
	

	if( (sc_a < 0) && (sc_c < 0) && (sc_b > 0) && (sc_d > 0)  )
	{
	 
		// printf(" \n \n Punkt liegt im Pixel \n \n");
		
		
		// oberes rechtes Teilrechteck
		// ..................................................................
		*(psc_v + 0) = 0;
		*(psc_v + 1) = sc_b;
		*(psc_v + 2) = 0;
		*(psc_v + 3) = sc_d; 
		
		
		
		f_min_max_euclidian_distance(psc_v, plb, pub);
		
		*ppa = (*(psc_v + 1) - *(psc_v + 0))*(*(psc_v + 3) - *(psc_v + 2));
		
		
	
		// printf("lower %f \n", *plb);
		// printf("upper %f \n", *pub);
		// printf("lower %f \n", *(plb + 1));
		// printf("upper %f \n", *(pub + 1));
		// printf("lower %f \n", *(plb + 2));
		// printf("upper %f \n", *(pub + 2));
		// printf("lower %f \n", *(plb + 3));
		// printf("upper %f \n", *(pub + 3));
		
		// oberes linkes Teilrechteck
		// ..................................................................
		*(psc_v + 4) = 0;
		*(psc_v + 5) = -1 * sc_a;
		*(psc_v + 6) = 0;
		*(psc_v + 7) = sc_d; 
		
		
		
		f_min_max_euclidian_distance((psc_v + 4), (plb + 1), (pub + 1));
		
		*(ppa +1) = (*(psc_v + 5) - *(psc_v + 4))*(*(psc_v + 7) - *(psc_v + 6));
									  
		
		
		
		// printf("lower %f \n", *plb);
		// printf("upper %f \n", *pub);
		// printf("lower %f \n", *(plb + 1));
		// printf("upper %f \n", *(pub + 1));
		// printf("lower %f \n", *(plb + 2));
		// printf("upper %f \n", *(pub + 2));
		// printf("lower %f \n", *(plb + 3));
		// printf("upper %f \n", *(pub + 3));
		
		// unteres linkes Teilrechteck
		// ..................................................................
		*(psc_v + 8) = 0;
		*(psc_v + 9) = fabs(-1 * sc_a);
		*(psc_v + 10) = 0;
		*(psc_v + 11) = fabs(-1 *sc_c); 
		
      
		
		f_min_max_euclidian_distance((psc_v + 8), (plb + 2), (pub + 2));
		
		*(ppa +2) = (*(psc_v + 9) - *(psc_v + 8))*(*(psc_v + 11) - *(psc_v + 10));
		
		
		// printf("lower %f \n", *plb);
		// printf("upper %f \n", *pub);
		// printf("lower %f \n", *(plb + 1));
		// printf("upper %f \n", *(pub + 1));
		// printf("lower %f \n", *(plb + 2));
		// printf("upper %f \n", *(pub + 2));
		// printf("lower %f \n", *(plb + 3));
		// printf("upper %f \n", *(pub + 3));
		
		// unteres rechtes Teilrechteck
		// ..................................................................
		*(psc_v + 12) = 0;
		*(psc_v + 13) = sc_b;
		*(psc_v + 14) = 0;
		*(psc_v + 15) = fabs(-1 *sc_c); 
		

		
		f_min_max_euclidian_distance((psc_v + 12), (plb + 3), (pub + 3));
		
		*(ppa +3) = (*(psc_v + 13) - *(psc_v + 12))*(*(psc_v + 15) - *(psc_v + 14));
		
		
		// printf("lower %f \n", *plb);
		// printf("upper %f \n", *pub);
		// printf("lower %f \n", *(plb + 1));
		// printf("upper %f \n", *(pub + 1));
		// printf("lower %f \n", *(plb + 2));
		// printf("upper %f \n", *(pub + 2));
		// printf("lower %f \n", *(plb + 3));
		// printf("upper %f \n", *(pub + 3));
		
		
		
	}
	
	
	
	// .....................................................................
	// Der Punkt ist zwischen zwei Kanten des Pixels
	// .....................................................................
	
	
	
	// zwischen den Kanten a und b und auf der Kante oder unterhalb von c
	// ..................................................................
	
	if( (sc_a < 0) && (sc_b > 0) && (sc_c >= 0))
	{
	
		// printf("\n \n Punkt liegt zwischen a und b auf oder unter c \n \n");
		
			// linkes Rechteck vom Punkt aus gesehen
			*(psc_v + 0) = 0;
			*(psc_v + 1) = fabs(-1* sc_a);
			*(psc_v + 2) = sc_c;
			*(psc_v + 3) = sc_d; 
		
			
			f_min_max_euclidian_distance(psc_v, plb, pub);
			
			*ppa = (*(psc_v + 1) - *(psc_v + 0)) * (*(psc_v + 3) - *(psc_v + 2));
		    

		
			// rechtes Rechteck vom Punkt aus gesehen
			*(psc_v + 4) = 0;
			*(psc_v + 5) = sc_b;
			*(psc_v + 6) = sc_c;
			*(psc_v + 7) = sc_d; 
		
			
			f_min_max_euclidian_distance((psc_v +4), (plb + 1), (pub + 1));
			
			*(ppa + 1) = (*(psc_v + 5) - *(psc_v + 4)) * (*(psc_v + 7) - *(psc_v + 6)); 

	} // end if
	
	// zwischen den Kanten a und b und auf der Kante oder überhalb von d
	// ..................................................................
	
	if( (sc_a < 0) && (sc_b > 0) && (sc_d <= 0))
	{
		
		 //printf("\n \n Punkt liegt zwischen a und b auf oder über d \n \n");
		
		// linkes Rechteck vom Punkt aus gesehen
		
		*(psc_v + 0) = 0;
		*(psc_v + 1) = fabs(-1* sc_a);
		*(psc_v + 2) = fabs(-1 * sc_d);
		*(psc_v + 3) = fabs(-1 * sc_c); 
		

		
		f_min_max_euclidian_distance(psc_v, plb, pub);
		
		*ppa = ( *(psc_v + 1) - *(psc_v +0) ) * ( *(psc_v + 3) - *(psc_v + 2) );
		
		// rechtes Rechteck vom Punkt aus gesehen
		*(psc_v + 4) = 0;
		*(psc_v + 5) = sc_b;
		*(psc_v + 6) = fabs(-1 * sc_d);
		*(psc_v + 7) = fabs(-1 * sc_c);  
		

		
		f_min_max_euclidian_distance((psc_v + 4), (plb + 1), (pub + 1));
		
		*(ppa + 1) = (*(psc_v + 5) - *(psc_v + 4)) * (*(psc_v +7) - *(psc_v + 6)); 
		
	} // end if
	
	// zwischen den Kanten c und d und auf der Kante oder links von a
	// ..................................................................
	
	if( (sc_c < 0) && (sc_d > 0) && (sc_a >= 0))
	{
		
		// printf("\n \n Punkt liegt zwischen c und d auf oder links von  a \n \n");
		// oberes Rechteck vom Punkt aus gesehen
		
		*(psc_v + 0) = sc_a;
		*(psc_v + 1) = sc_b;
		*(psc_v + 2) = 0;
		*(psc_v + 3) = sc_d; 
		
		
		
		f_min_max_euclidian_distance(psc_v, plb, pub);
		
		*ppa = (*(psc_v + 1) - *(psc_v + 0)) * (*(psc_v + 3) - *(psc_v + 2));
		
		// unteres Rechteck vom Punkt aus gesehen
		
		*(psc_v + 4) = sc_a;
		*(psc_v + 5) = sc_b;
		*(psc_v + 6) = 0;
		*(psc_v + 7) = fabs(-1 * sc_c);  
		
		
		
		f_min_max_euclidian_distance((psc_v + 4), (plb + 1), (pub + 1));
		
		*(ppa + 1) = (*(psc_v + 5)  - *(psc_v + 4)) * (*(psc_v + 7) - *(psc_v + 6)); 
		
	
	
	} // end if
	
	// zwischen den Kanten c und d und auf der Kante oder rechts von b
	// ..................................................................
	
	if( (sc_c < 0) && (sc_d > 0) && (sc_b <= 0))
	{
		
		// printf("\n \n Punkt liegt zwischen c und d auf oder rechts von b \n \n");
		// oberes Rechteck vom Punkt aus gesehen
		
		*(psc_v + 0) = fabs(-1 * sc_b);
		*(psc_v + 1) = fabs(-1 * sc_a);
		*(psc_v + 2) = 0;
		*(psc_v + 3) = sc_d; 
		
		
		
		f_min_max_euclidian_distance(psc_v, plb, pub);
		
		*ppa = (*(psc_v + 1) - *(psc_v + 0)) * (*(psc_v + 3) - *(psc_v + 2));
		
		// unteres Rechteck vom Punkt aus gesehen
		
		*(psc_v + 4) = fabs(-1 * sc_b);
		*(psc_v + 5) = fabs(-1 * sc_a);
		*(psc_v + 6) = 0;
		*(psc_v + 7) = fabs(-1 * sc_c);  
		
		
		
		f_min_max_euclidian_distance((psc_v + 4), (plb + 1), (pub + 1));
		
		*(ppa + 1) = (*(psc_v + 5) - *(psc_v + 4)) * (*(psc_v + 7) - *(psc_v + 6)); 
		
		
		
	} // end if
	

	

	
} // end of f_int_boundaries

void f_min_max_euclidian_distance(double *psc_v, double *plb, double *pub){
	
	int i;
	double min, max;
	double dist[4];

	
	dist[0] = pow( pow( (*(psc_v + 0) - 0), 2 ) + pow( (*(psc_v + 2) - 0), 2), 0.5 );
	dist[1] = pow( pow( (*(psc_v + 0) - 0), 2 ) + pow( (*(psc_v + 3) - 0), 2), 0.5 );
	dist[2] = pow( pow( (*(psc_v + 1) - 0), 2 ) + pow( (*(psc_v + 3) - 0), 2), 0.5 );
	dist[3] = pow( pow( (*(psc_v + 1) - 0), 2 ) + pow( (*(psc_v + 3) - 0), 2), 0.5 );

	
	
	 //printf("1 %f \n", dist[0]);
	 //printf("2 %f \n", dist[1]);
	 //printf("3 %f \n", dist[2]);
	 //printf("4 %f \n", dist[3]);
		  								
	min = dist[1];
	max = dist[1];
	
	for (i = 0; i < 4; i++)
    {
		if (dist[i] >= max)
        {
			max = dist[i];
        }
		else if (dist[i] <= min)
        {
			min = dist[i];
        }
    }
	
	*(plb) = min;
	*(pub) = max;
	
	// printf("lowerbound  %f \n", *plb);
	// printf("upperbound %f \n",  *pub);
	
	

}




