/*
 Ckrige, a program for universal, constrained and covariance-matching
 constrained point and block kriging. Check the help page of the R function
 CovarianceFct in the RandomFields package to get more information about 
 a single covariance function.
 Copyright 2010 (C) Christoph Hofer
 
 Christoph Hofer christoph.hofer@alumni.ethz.ch
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
// Integrieren der Kovarianzfunktionen zwischen den Stützpunkten
// und Rechtecken
//
// Author: Christoph Hofer
// Datum: Dezember 2009
// geändert am 29-01-2010, ch
// 2023-01-24 A. Papritz elimination of unnneeded variable setlb, setub

#include <stdio.h>
#include <stdlib.h> 
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
// Definiton der Kovarianzfunktionen
#include "covariancefunctions.h"
// Definition der Distanz-Dichtefunktion
#include "intboundaries.h"
#include "scalePixelCoord.h"
#include "freqdist.h"
//
#define datasupport pow(10, -6) // Für Pixel mit kleinerere Fläche werden als Punkt behandelt
//
// pcp = pixel covariance parameter
typedef struct pcp 
{
	double sc_pixcoords[4]; // skalierte Pixelkoordinaten
	int ncp; // Anzahl Kovarianzparameter, wird nicht gebraucht
	double *cp; // Zeiger auf die Parameter der Kovarianzfunktion
	
} pcp_type;
//
// wird von R aufgerufen
void PointRectCov(
				  double *pxcoord, // x-Koordinate des (Stütz)Punktes
				  double *pycoord, // y-Koordinate des (Stütz)Punktes
				  double *crxcoord, // x Koordinate des Zentrums des Rechteckes
				  double *crycoord, // y Koordinate des Zentrums des Rechteckes
				  double *rxwidth, // x-width of the rectangle
				  double *rywidth, // y-width of the rectangle
				  double *npoints,
				  double *npixel,
				  double *epsabs,
				  double *epsrel,
				  double *result,
				  double *abserr,
				  int *limit,
				  int *neval,
				  int *ier,
				  int *lenw,
				  int *last,
				  int *iwork,
				  double *work,
				  void *param, // param: Parameter der Kovarianzfunktionen
				  int *nmod, // nmod: Anzahl Kovarianzfunktionen
				  int *nmodparam); // Anzahl Kovarianzparameter pro Kovarianzfunktion
//
// ***************************************************
// Definiton der Funktion(en) welche Integriert werden soll(en)
void f_integrate( double *x, int n, void *p_pcp );
// ***************************************************
//
// ***************************************************
// Pointer auf Kovarianzfunktion wird in der Funktion 
// PointRectCov gesetzt
// ***************************************************
double (*f_cov)( double x,  double *covpar);
// ***************************************************
//
void PointRectCov(
				  double *pxcoord, 
				  double *pycoord, 
				  double *crxcoord, 
				  double *crycoord, 
				  double *rxwidth, 
				  double *rywidth,
				  double *npoints,
				  double *npixel,
				  double *epsabs,
				  double *epsrel,
				  double *result,
				  double *abserr,
				  int *limit,
				  int *neval,
				  int *ier,
				  int *lenw,
				  int *last,
				  int *iwork,
				  double *work,
				  void *param, // param: Parameter der Kovarianzfunktionen
				  int *nmod, // nmod: Anzahl Kovarianzfunktionen
				  int *nmodparam) // rect_area:Fläche des Pixels
{
	unsigned long int ipixel, ipoint; // Laufvariable ipixel über die 
									  // Anzahl der Pixel und ipoint über die Anzahl Punkte
	int  icov, q, inmod,j, modnr, imodr; 
	// i= Laufvariable # Cov-Modelle, j = Laufvariable # Integrationsgrenzen (1,2, oder 4), 
	//modnr = Nummer der Kovarianzfunktion, imod = Laufvariable Anzahl Kovarianzmodelle
	
	double temp_cov, cov, res, rect_area, pcxdiff, pcydiff;
	//int *neval = NULL, *ier = NULL, *lenw= NULL, *last=NULL;  
	//int *iwork= NULL;
	//double *work = NULL;
	
	double sc_pixcoords[] = {-1,-1,-1,-1,
						  -1,-1,-1,-1,
						  -1,-1,-1,-1,
						  -1,-1,-1,-1};
	
	// Pointer zeigt auf Array mit den skalierten Pixelcoordinaten
	double *p_sc_pixcoords = NULL;  
	
	p_sc_pixcoords = & sc_pixcoords[0]; 

	
	double *param_new = NULL;
	double *pres=NULL;
	
	// Array für die unteren Integrationsgrenzen (max. 4 Punkt im Pixel)
	double setlowerbound[4] = {-1,-1,-1,-1}; 
	// Array für die oberen Integrationsgrenzen (max. 4 Punkt im Pixel)
	double setupperbound[4] = {-1,-1,-1,-1}; 
	// Pointer auf die Arrays der Integrationsgrenzen
	//
	double *plb, *pub; 
	plb = &setlowerbound[0]; //
	pub = &setupperbound[0]; //
	
		
	double pixelarea[4] = {-1,-1,-1,-1};
	double *ppa;
	

	
	ppa = &pixelarea[0];
	

	
    // printf("px %f \n", *pxcoord);
	// printf("py %f \n", *pycoord);
	
	
	param_new = (double*)param; // cast void pointer in double pointer


	//pcp Struktur vom Type pcp_type
	//p_pcp pointer auf die Struktur pcp
	pcp_type pcp, *p_pcp;
	p_pcp = &pcp;
	
	// Werte der skalierten Pixelkoordinaten auf 0 setzen
	(*p_pcp).sc_pixcoords[0] = 0;
	(*p_pcp).sc_pixcoords[1] = 0;
	(*p_pcp).sc_pixcoords[2] = 0;
	(*p_pcp).sc_pixcoords[3] = 0;
	
	// Initialisierung des Kovarianzparametervektors
	(*p_pcp).cp = param_new;
	
	
	icov = 0;

	for(ipixel=0; ipixel < *npixel; ipixel++)
	{
		


	
	//Schleife für einen Pixel über alle Punkte
	ipoint = 0;
		
	for(ipoint=0; ipoint < *npoints; ipoint++)
	{
	imodr = 0;
	res = 0;
    pres = &res;
	cov = 0;
	temp_cov = 0;

    rect_area = (*rxwidth) * (*rywidth);
	//*ipxcoord = *(pxcoord + ipoint);
	//*ipycoord = *(pycoord + ipoint);
	
		
	
	for(inmod = 0; inmod < *nmod; inmod++)
	{
		//printf("imodr = %d \n", imodr);
		//printf("ipoint = %d \n", ipoint);
		//printf("i = %d \n", i);
		//printf("1 param_new(0) = %f \n", *(param_new));
		modnr = (int) *(param_new + imodr);
		
		//printf("modnr %d \n", modnr);
		
		switch(modnr)
		{
			case 0: // nugget
				//printf("Nugget\n");
				// printf("Area %f \n",(*rxwidth) * (*rywidth));
		
				
				pcxdiff = fabs( *(pxcoord + ipoint) - *(crxcoord +ipixel) );
				pcydiff = fabs( *(pycoord + ipoint) - *(crycoord +ipixel) );
				
				if( pcxdiff > 0.5* (*rxwidth) || pcydiff > 0.5* (*rywidth) )
				// Punkt auserhalb des Pixels
				{
					temp_cov = 0;
				}
				else // Punkt innerhalb des Pixels
				{
					if (rect_area > datasupport) // Pixelfläche > 10^-6
					{
						temp_cov = *(param_new + imodr + 1) * datasupport  / rect_area;
					}
					else // Pixelfläche < 10^-6
					{
						temp_cov = *(param_new + imodr + 1);
					}
				}
			
				break;
			
			case 1: // exponential
				// printf("exponential \n"); 
				
				f_cov = f_cov_exponential;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
			
			case 2: // spherical
			// printf("spherical \n");
				
				f_cov = f_cov_sphercial;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 3: // matern
				// printf("matern \n");
				
				f_cov = f_cov_matern;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
			
			case 4: // bessel
				
				f_cov = f_cov_bessel;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 5: // cauchy
				
				f_cov = f_cov_cauchy;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 6: // cauchytbm
				
				f_cov = f_cov_cauchytbm;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 7: // circular
				
				f_cov = f_cov_circular;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 8: // constant
				
				f_cov = f_cov_constant;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 9: // cubic
				
				f_cov = f_cov_cubic;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 10: // dampedcosine
				
				f_cov = f_cov_dampedcosine;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 11: // gencauchy
				
				f_cov = f_cov_gencauchy;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 12: // gengneiting;  a= 1
				
				f_cov = f_cov_gengneiting1;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 13: // gengneiting;  a= 2
				
				f_cov = f_cov_gengneiting2;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 14: // gengneiting;  a= 3
				
				f_cov = f_cov_gengneiting3;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 15: // gneiting
				
				f_cov = f_cov_gneiting;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 16: // hyperbolic
				
				f_cov = f_cov_hyperbolic;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 17: // penta
				
				f_cov = f_cov_penta;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 18: // lgd1
				
				f_cov = f_cov_lgd1;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 19: // power
				
				f_cov = f_cov_power;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 20: // wave
				
				f_cov = f_cov_wave;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			case 21: // qexponential
				
				f_cov = f_cov_qexponential;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
                
      case 22: // whittle
				
				f_cov = f_cov_whittle;
				(*p_pcp).cp =  param_new + imodr + 1;
				break;
				
			default:
				// printf("Kein gültiges Kovarianzmodell.");
				temp_cov = -999;
	
		}
		
		// printf("modnr %d \n", modnr);
		if (modnr != 0) // Nugget wird nicht integriert
			{
				// skalierung der Pixelkoordinaten bezüglich eines Stützpunktes
				scale_pixel_coord( (pxcoord+ipoint), (pycoord+ipoint), (crxcoord + ipixel), (crycoord+ipixel), rxwidth, rywidth, p_sc_pixcoords);
			
				// Berechung der Integrationsgrenzen lower und upperbound
				// ppa pointer inhalt Fläche des entsprechenden Rechteckes
				// plb pointer auf untere Integrationsgrenze
				// pub pointer auf obere Integrationsgrenze
				f_int_boundaries(p_sc_pixcoords, plb, pub, ppa);
					
				j = 0;
				q = 0;
				//rect_area = 0;
				while( *(plb + j) >= 0 && *(pub + j) >= 0 && j < 4)
				{
					(*p_pcp).sc_pixcoords[0] = *(p_sc_pixcoords + q + 0);
					(*p_pcp).sc_pixcoords[1] = *(p_sc_pixcoords + q + 1);
					(*p_pcp).sc_pixcoords[2] = *(p_sc_pixcoords + q + 2);
					(*p_pcp).sc_pixcoords[3] = *(p_sc_pixcoords + q + 3);
					 q = q + 4;
				
			    // Integral Funktion (von R)
				Rdqags(f_integrate, (void*)p_pcp, (plb + j), (pub + j),
					   epsabs, epsrel, pres, abserr, neval, ier,
					   limit, lenw, last, iwork, work);
				
					
					// die benutztenn Integrationsgrenzne werden auf -1 gesezt
					// Abruchkriterium des while loops
					*(plb +j) = -1;
					*(pub +j) = -1;
					
				//	printf("\n ier %d \n", *ier);
				//	printf("\n neval %d \n", *neval);
				//	printf("\n last %d \n", *last);
				//	printf("\n iwork %d \n", *iwork);
				//	printf("\n work %f \n", *work);
					
					
					//Rprintf("Nr:  %d  ,", ipoint);
					
				
					//printf("integral %f \n", *pres);
					
					temp_cov = temp_cov + (*pres) * (*(ppa +j)) ;
					//printf("temp_cov %f \n", temp_cov);
					//printf("temp_cov %f \n", temp_cov);
					//*pres = 0;
					
				j++; 
				}
				temp_cov =  temp_cov / rect_area;
			} 
		// printf("temp_cov %f \n", temp_cov);
		//printf("temp_cov %f \n", temp_cov);
		

		cov= cov +  temp_cov ;
		temp_cov = 0;
		
				//printf("cov %f \n", cov);
		imodr = imodr + 3 + *nmodparam;
	}
	
	
	
	
		//printf("cov %f \n", cov);
		 
		*(result + icov) = cov ;// / rect_area;
	icov++;	
	} // end for ipoint Schleife
		
		
	
	} // end for ipixel Schleife 
	

	

} // end PointRectCov
	
	
	
// ***************************************************
//Die zu integrierende Funktion 
//***************************************************************

void f_integrate( double *x, int n, void *p_pcp ) 
{ 
	int i;

	// Vektor auf den Array der Pixelkoordinaten
	double *p_sc_pixcoods;
	
	//Cast void pointer auf pcp_type pointer
	pcp_type *p_pixel_cov_param;
	p_pixel_cov_param = (pcp_type*)p_pcp;
	
	p_sc_pixcoods = &(*p_pixel_cov_param).sc_pixcoords[0];
	
	for(i=0;i<n;i++)
	{ 
		
		x[i] =  f_dist_freq( x[i], p_sc_pixcoods) *  f_cov( x[i], (*p_pixel_cov_param).cp );
		
	} // end for
} // end f_integrate 
