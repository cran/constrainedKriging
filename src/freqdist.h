/*
 * Ckrige, a program for universal, constrained and covariance-matching
 * constrained point and block kriging.
 * Copyright 2010 (C) Christoph Hofer
 *
 * Christoph Hofer christoph.hofer@alumni.ethz.ch
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version. As a special exception, linking
 * this program with the Qt library is permitted.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

// Headerfile Definition der Distanz-Dichtefunktionen zwischen
// einem Punkt und einem Rechteck in ersten Quadranten

// Author: Christoph Hofer
// Date: Dezember 2009

#include <stdio.h>
#include <R.h>


double f_dist_freq(
  double lag_dist,         // distance between origin and a point uniformly in rectangle in first quadrant
  double *p_sc_pixcoords); // pointer to array with coordinates of lower left and upper right corner of rectangle


double f_dist_freq( double lag_dist, double *p_sc_pixcoords)
{

  // density of distance between origin and a point uniformly distributed
  // in a rectangle that lies in the first quadrant
  // cf Hofer & Papritz, 2011, constrainedKriging: An R-package for customary, constrained
  // and covariance-matching constrained point or block kriging, Computers &amp; Geosciences
  // Vol. 37, No. 10, p. 1562-1569, Appendix A

  // the function returns the density for a given distance and rectangle

  // 2009-12-xx C. Hofer
  // 2023-01-24 A. Papritz editing of comments, verification of computations
  // 2025-02-01 A. Papritz, minor changes (definition of angle_a, angle_c)

  double angle_a, angle_c;
  double ba, dc, lag_dist_sq, a, a_sq, b, b_sq, c, c_sq, d, d_sq, K;
  double dist_freq = 0;

  lag_dist_sq = pow(lag_dist,2);

  a = *(p_sc_pixcoords);       // x-coordinate of lower left corner of rectangle
  a_sq= pow(a, 2);

  b =  *(p_sc_pixcoords + 1);  // x-coordinate of upper right corner of rectangle
  b_sq= pow(b, 2);

  c =  *(p_sc_pixcoords + 2);  // y-coordinate of lower left corner of rectangle
  c_sq= pow(c, 2);

  d =  *(p_sc_pixcoords + 3);  // y-coordinate of upper right corner of rectangle
  d_sq= pow(d, 2);

  // lengths of rectangle

  ba = b - a;
  dc = d - c;

  // inverse of area of rectangle

  K = ( 1 / ba ) * ( 1 / dc );

	// angles for positions where rectangle touches axes

	if( c_sq <= 0. ){
	   // c = 0
	  angle_c = 0.5 * M_PI;
	} else {
	  angle_c = atan( sqrt( lag_dist_sq - c_sq ) / c );
	}
	if( a_sq <= 0. ){
	  // a = 0
	  angle_a = 0.;
	} else {
	  angle_a = atan( a / sqrt( lag_dist_sq - a_sq ) );
	}

  // convention: angles positive from north in clockwise direction

	// case BII (zone II if (b - a) < (d - c))
	if( ( lag_dist_sq > (b_sq + c_sq) ) && ( lag_dist_sq <= (a_sq + d_sq) ) )
	{
	  dist_freq = K * lag_dist * ( atan( b / sqrt( lag_dist_sq - b_sq ) ) - angle_a );
	}

	//cases AI and BI (zone I)
	else if( ( lag_dist_sq > (b_sq + c_sq) ) && ( lag_dist_sq > (a_sq + d_sq) ) && ( lag_dist_sq <= (b_sq + d_sq) ) )
	{
	  dist_freq = K * lag_dist * ( atan( b / sqrt( lag_dist_sq - b_sq ) ) - atan( sqrt( lag_dist_sq - d_sq ) / d ) );
	}

	//cases AIII and BIII (zone III)
	else if( ( lag_dist_sq <= (b_sq + c_sq) ) && ( lag_dist_sq <= (a_sq + d_sq) ) && ( lag_dist_sq >= (a_sq + c_sq) ) )
	{
	  dist_freq = K * lag_dist * ( angle_c - angle_a );
	}

	//case AII (zone II if (b - a) >= (d-c))
	else if( ( lag_dist_sq <= (c_sq + b_sq) ) && ( lag_dist_sq > ( a_sq + d_sq ) ) )
	{
	  dist_freq = K * lag_dist * ( angle_c  - atan( sqrt( lag_dist_sq - d_sq ) / d ) );
	}

  return(dist_freq);

}
