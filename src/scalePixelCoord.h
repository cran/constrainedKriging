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
// Skalierung der Pixelkoordinaten bezüglich eines Punktes 
//
//
// Author: Christoph Hofer
// Date: Dezember 2009

#include <stdio.h>


void scale_pixel_coord(
					   double *pxcoord,   // x-Koordinate des (Stütz)Punktes
					   double *pycoord,   // y-Koordinate des (Stütz)Punktes
					   double *crxcoord,  // x Koordinate des Zentrums des Rechteckes
					   double *crycoord,  // y Koordinate des Zentrums des Rechteckes
					   double *rxwidth,   // x-width of the rectangle
					   double *rywidth,   // y-width of the rectangledouble x
					   double *psc_v);    // Point auf den Array mit den Skalierten Koordinaten des Pixels


void scale_pixel_coord(
						double *pxcoord, 
						double *pycoord, 
						double *crxcoord, 
						double *crycoord, 
						double *rxwidth,
					    double *rywidth,
						double *psc_v)
{

	//printf("px %f \n", *pxcoord);
	// printf("py %f \n", *pycoord);
	

	
	// Skalierung der Eck-Koordinaten des Pixels bezüglich des Stützpunktes
	
	*(psc_v + 0)  = (*crxcoord - 0.5 * (*rxwidth)) - (*pxcoord); // skalierte x-coord der linken untere Ecke a
	*(psc_v + 1)  = (*crxcoord + 0.5 * (*rxwidth)) - (*pxcoord); // skalierte x-coord der rechten untere Ecke b
	*(psc_v + 2)  = (*crycoord - 0.5 * (*rywidth)) - (*pycoord); // skalierte y-coord der linken unteren Ecke c
	*(psc_v + 3)  = (*crycoord + 0.5 * (*rywidth)) - (*pycoord); // skalierte y-coord der linken oberen Ecke d
	

	//printf("a %f \n", *(psc_v + 0));
	//printf("b %f \n", *(psc_v + 1));
	//printf("c %f \n", *(psc_v + 2));
	//printf("d %f \n", *(psc_v + 3));
	
	
} // end scale_pixel_coord
