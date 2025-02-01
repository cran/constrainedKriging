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

// Author: Christoph Hofer
// Date: Dezember 2009

#include <stdio.h>
#include<math.h>

void f_int_boundaries(void *cpixcoords, // pointer to the array[16] with the coordinates of the pixel
                      // or its respective 2 - 4 parts, centered before by the coordinates of a support point
                      double *plb,  // pointer to the array with the lower integration limits of the pixel or its respective 2 - 4 parts
                      double *pub,  // pointer to the array with the upper integration limits of the pixel or its respective 2 - 4 parts
                      double *ppa); // pointer to the array with the area of the pixel or its respective 2 - 4 parts


void f_min_max_euclidian_distance(double *psc_v, double *plb, double *pub);



void f_int_boundaries(void *cpixcoords,
                      double *plb,
                      double *pub,
                      double *ppa)
{

  // function shifts the pixel into the first quadrant around a support point
  // if the pixel pixel lies in a single quadrant.  if the pixel intersects
  // with 2 quadrants then the pixel is split into 2 parts and the 2 parts
  // are shifted again into the first quadrant.  if the pixel intersects with
  // all 4 quadrants then the pixel is split into 4 parts and the 4 parts are
  // shifted again into the first quadrant.

  // the function then computes for the single pixel or for each of the 2
  // or 4 parts minimum and maximum distance(s) of the 4 corners from the
  // support points (these are used as lower and upper integration limits)

  // the function returns by the pointers plb, puk and ppa the lower and
  // upper integration limits for the pixel or its 2-4 parts and the areas
  // of the pixel or its 2-4 parts

  // 2009-12-xx C. Hofer
  // 2023-01-24 A. Papritz elimination of unnneeded variable sc_v
  // 2025-01-19 A. Papritz pointer *modparam renamed to *cpixcoords
  //                       Rprintf() substituted for printf()
  //                       editing of comments and code reformatting

	
	/* 
	 * int npart;  // number of parts into which pixel is subdivided
	 */
	

  double sc_a, sc_b, sc_c, sc_d;

  double *psc_v;
  psc_v = (double*) cpixcoords; // cast void pointer in double pointer

	/* 
	 * npart = 1;
	 */

  sc_a = *(psc_v + 0);
  sc_b = *(psc_v + 1);
  sc_c = *(psc_v + 2);
  sc_d = *(psc_v + 3);

	/*
	 * Rprintf("\nf_int_boundaries:centred coordinates of lower left (a,c) and upper right pixel corners (b,d)\n");
	 * Rprintf("f_int_boundaries:a %f \n", sc_a);
	 * Rprintf("f_int_boundaries:b %f \n", sc_b);
	 * Rprintf("f_int_boundaries:c %f \n", sc_c);
	 * Rprintf("f_int_boundaries:d %f \n", sc_d);
	 */

  // .....................................................................
  // CASES I-IV:
  // support points lies outside of pixel or lies on one of its 4 corners
  // pixel intersects with only one quadrant around support point
  // .....................................................................


  // CASE I: pixel intersects only with first quadrant (no reflection)
  // .....................................................................

  if( (sc_a >= 0) && (sc_c >= 0) )
  {

		/*
		 * Rprintf("\nf_int_boundaries: CASE I: pixel intersects only with first quadrant  \n");
		 */

    *(psc_v + 0) = sc_a;
    *(psc_v + 1) = sc_b;
    *(psc_v + 2) = sc_c;
    *(psc_v + 3) = sc_d;

    f_min_max_euclidian_distance(psc_v, plb, pub);

    *ppa = ( *(psc_v + 1) - *(psc_v + 0) ) * ( *(psc_v + 3) - *(psc_v + 2) );

  }


  // CASE II: pixel intersects only with second quadrant (reflection at ordinate)
  // .....................................................................

  if( (sc_b <= 0) && (sc_c >= 0) )
  {

		/*
		 * Rprintf("\nf_int_boundaries: CASE III: pixel intersects only with second quadrant  \n");
		 */

    *(psc_v + 0) = fabs(-1 * sc_b); // a new = -1* b old
    *(psc_v + 1) = fabs(-1 * sc_a) ; // b new = -1* a old
    *(psc_v + 2) = sc_c;
    *(psc_v + 3) = sc_d;

    *ppa = ( *(psc_v +1) - *(psc_v +0) )* ( *(psc_v +3) - *(psc_v +2) );

    f_min_max_euclidian_distance(psc_v, plb, pub);

  }


  // CASE III: pixel intersects only with third quadrant (reflection at origin)
  // .....................................................................

  if( (sc_b <= 0) && (sc_d <= 0) )
  {

		/*
		 * Rprintf("\nf_int_boundaries: CASE IV: pixel intersects only with third quadrant  \n");
		 */

    *(psc_v + 0) =  fabs(-1 * sc_b); // a new = -1 * b old
    *(psc_v + 1) =  fabs(-1 * sc_a); // b new = -1 * a old
    *(psc_v + 2) =  fabs(-1 * sc_d); // c new = -1 * d old
    *(psc_v + 3) =  fabs(-1 * sc_c); // d new = -1 * c old

    *ppa = (*(psc_v + 1) - *(psc_v + 0))*(*(psc_v + 3) - *(psc_v + 2));

    f_min_max_euclidian_distance(psc_v, plb, pub);

  }


  // CASE IV: pixel intersects only with fourth quadrant (reflection at abscissa)
  // .....................................................................

  if( (sc_a >= 0) && (sc_d <= 0) )
  {

		/*
		 * Rprintf("\nf_int_boundaries: CASE II: pixel intersects only with fourth quadrant  \n");
		 */

    *(psc_v + 0) = sc_a;
    *(psc_v + 1) = sc_b;
    *(psc_v + 2) = fabs(-1 * sc_d); // c new = -1* d old
    *(psc_v + 3) = fabs(-1 * sc_c); // d new = -1* c old

    f_min_max_euclidian_distance(psc_v, plb, pub);

    *ppa = ( *(psc_v + 1) - *(psc_v + 0)) * ( *(psc_v + 3)  - *(psc_v + 2) ) ;

  }


  // .....................................................................
  // CASE V: support points lies in interior of pixel
  // pixel intersects with all 4 quadrants around support point
  // the intersections of the pixel with the 4 quadrants define the 4 parts
  // which are then point- or line-symmetrically reflected into the first quadrant
  // .....................................................................

  if( (sc_a < 0) && (sc_c < 0) && (sc_b > 0) && (sc_d > 0)  )
  {

		/*
		 * Rprintf("\nf_int_boundaries: CASE V: pixel intersects with 4 quadrants \n");
		 */

		/* 
		 * npart = 4;
		 */

    // upper right part (no reflection)

    *(psc_v + 0) = 0;
    *(psc_v + 1) = sc_b;
    *(psc_v + 2) = 0;
    *(psc_v + 3) = sc_d;

    f_min_max_euclidian_distance(psc_v, plb, pub);

    *ppa = (*(psc_v + 1) - *(psc_v + 0))*(*(psc_v + 3) - *(psc_v + 2));


    // upper left part (reflection at ordinate)

    *(psc_v + 4) = 0;
    *(psc_v + 5) = fabs(-1 * sc_a);
    *(psc_v + 6) = 0;
    *(psc_v + 7) = sc_d;

    f_min_max_euclidian_distance((psc_v + 4), (plb + 1), (pub + 1));

    *(ppa +1) = (*(psc_v + 5) - *(psc_v + 4))*(*(psc_v + 7) - *(psc_v + 6));


    // lower left part (reflection at origin)

    *(psc_v + 8)  = 0;
    *(psc_v + 9)  = fabs(-1 * sc_a);
    *(psc_v + 10) = 0;
    *(psc_v + 11) = fabs(-1 * sc_c);

    f_min_max_euclidian_distance((psc_v + 8), (plb + 2), (pub + 2));

    *(ppa +2) = (*(psc_v + 9) - *(psc_v + 8))*(*(psc_v + 11) - *(psc_v + 10));


    // lower right part (reflection at abscissa)

    *(psc_v + 12) = 0;
    *(psc_v + 13) = sc_b;
    *(psc_v + 14) = 0;
    *(psc_v + 15) = fabs(-1 * sc_c);

    f_min_max_euclidian_distance((psc_v + 12), (plb + 3), (pub + 3));

    *(ppa +3) = (*(psc_v + 13) - *(psc_v + 12))*(*(psc_v + 15) - *(psc_v + 14));

  }


  // .....................................................................
  // CASES VI - IX: pixel intersects with 2 quadrants around support point
  // the intersections of the pixel with the 2 quadrants define the 2 parts
  // which are then point or line-symmetrically reflected into the first quadrant
  // .....................................................................

  // CASE VI: pixel intersects only with first and second quadrants
  // ..................................................................

  if( (sc_a < 0) && (sc_b > 0) && (sc_c >= 0))
  {

		/*
		 * Rprintf("\nf_int_boundaries: CASE VI: pixel intersects only with first and second quadrants\n");
		 */

		/* 
		 * npart = 2;
		 */

    // left part (reflection at ordinate)

    *(psc_v + 0) = 0;
    *(psc_v + 1) = fabs(-1 * sc_a);
    *(psc_v + 2) = sc_c;
    *(psc_v + 3) = sc_d;

    f_min_max_euclidian_distance(psc_v, plb, pub);

    *ppa = (*(psc_v + 1) - *(psc_v + 0)) * (*(psc_v + 3) - *(psc_v + 2));


    // right part (no reflection)

    *(psc_v + 4) = 0;
    *(psc_v + 5) = sc_b;
    *(psc_v + 6) = sc_c;
    *(psc_v + 7) = sc_d;

    f_min_max_euclidian_distance((psc_v +4), (plb + 1), (pub + 1));

    *(ppa + 1) = (*(psc_v + 5) - *(psc_v + 4)) * (*(psc_v + 7) - *(psc_v + 6));

  } // end if


  // CASE VII: pixel intersects only with third and fourth quadrant
  // ..................................................................

  if( (sc_a < 0) && (sc_b > 0) && (sc_d <= 0))
  {

		/*
		 * Rprintf("\nf_int_boundaries: CASE VII: pixel intersects only with third and fourth quadrants\n");
		 */

		/* 
		 * npart = 2;
		 */

    // left part (reflection at origin)

    *(psc_v + 0) = 0;
    *(psc_v + 1) = fabs(-1 * sc_a);
    *(psc_v + 2) = fabs(-1 * sc_d);
    *(psc_v + 3) = fabs(-1 * sc_c);

    f_min_max_euclidian_distance(psc_v, plb, pub);

    *ppa = ( *(psc_v + 1) - *(psc_v +0) ) * ( *(psc_v + 3) - *(psc_v + 2) );


    // right part (reflection at abscissa)

    *(psc_v + 4) = 0;
    *(psc_v + 5) = sc_b;
    *(psc_v + 6) = fabs(-1 * sc_d);
    *(psc_v + 7) = fabs(-1 * sc_c);

    f_min_max_euclidian_distance((psc_v + 4), (plb + 1), (pub + 1));

    *(ppa + 1) = (*(psc_v + 5) - *(psc_v + 4)) * (*(psc_v +7) - *(psc_v + 6));

  } // end if


  // CASE VIII: pixel intersects only with first and fourth quadrant
  // ..................................................................

  if( (sc_c < 0) && (sc_d > 0) && (sc_a >= 0))
  {

		/*
		 * Rprintf("\nf_int_boundaries: CASE VIII: pixel intersects only with first and fourth quadrants\n");
		 */

		/* 
		 * npart = 2;
		 */

    // upper part (no reflection)

    *(psc_v + 0) = sc_a;
    *(psc_v + 1) = sc_b;
    *(psc_v + 2) = 0;
    *(psc_v + 3) = sc_d;

    f_min_max_euclidian_distance(psc_v, plb, pub);

    *ppa = (*(psc_v + 1) - *(psc_v + 0)) * (*(psc_v + 3) - *(psc_v + 2));


    // lower part (reflection at abscissa)

    *(psc_v + 4) = sc_a;
    *(psc_v + 5) = sc_b;
    *(psc_v + 6) = 0;
    *(psc_v + 7) = fabs(-1 * sc_c);

    f_min_max_euclidian_distance((psc_v + 4), (plb + 1), (pub + 1));

    *(ppa + 1) = (*(psc_v + 5)  - *(psc_v + 4)) * (*(psc_v + 7) - *(psc_v + 6));

  } // end if


  // CASE IX: pixel intersects only with second and third quadrant
  // ..................................................................

  if( (sc_c < 0) && (sc_d > 0) && (sc_b <= 0))
  {

		/*
		 * Rprintf("\nf_int_boundaries: CASE IX: pixel intersects only with second and third quadrants\n");
		 */

		/* 
		 * npart = 2;
		 */

    // upper part (reflection at ordinate)

    *(psc_v + 0) = fabs(-1 * sc_b);
    *(psc_v + 1) = fabs(-1 * sc_a);
    *(psc_v + 2) = 0;
    *(psc_v + 3) = sc_d;

    f_min_max_euclidian_distance(psc_v, plb, pub);

    *ppa = (*(psc_v + 1) - *(psc_v + 0)) * (*(psc_v + 3) - *(psc_v + 2));


    // lower part (reflection at origin)

    *(psc_v + 4) = fabs(-1 * sc_b);
    *(psc_v + 5) = fabs(-1 * sc_a);
    *(psc_v + 6) = 0;
    *(psc_v + 7) = fabs(-1 * sc_c);

    f_min_max_euclidian_distance((psc_v + 4), (plb + 1), (pub + 1));

    *(ppa + 1) = (*(psc_v + 5) - *(psc_v + 4)) * (*(psc_v + 7) - *(psc_v + 6));

  } // end if


	/*
	 * // debugging output
	 *
	 * Rprintf("\nf_int_boundaries: transformed centred coordinates of pixel or its parts\n");
	 * Rprintf("  a_1 %f, b_1 %f, c_1 %f, d_1 %f, lwrlim_1 %f, uprlim_1 %f, area_1 %f\n",
	 *         *(psc_v + 0), *(psc_v + 1), *(psc_v + 2), *(psc_v + 3), *(plb + 0), *(pub + 0), *(ppa + 0));
	 *
	 * if( npart == 2 )
	 * {
	 *   Rprintf("  a_2 %f, b_2 %f, c_2 %f, d_2 %f, lwrlim_2 %f, uprlim_2 %f, area_2 %f\n",
	 *         *(psc_v + 4), *(psc_v + 5), *(psc_v + 6), *(psc_v + 7), *(plb + 1), *(pub + 1), *(ppa + 1));
	 * }
	 *
	 * if( npart == 4 )
	 * {
	 *   Rprintf("  a_2 %f, b_2 %f, c_2 %f, d_2 %f, lwrlim_2 %f, uprlim_2 %f, area_2 %f\n",
	 *         *(psc_v + 4), *(psc_v + 5), *(psc_v + 6), *(psc_v + 7), *(plb + 1), *(pub + 1), *(ppa + 1));
	 *   Rprintf("  a_3 %f, b_3 %f, c_3 %f, d_3 %f, lwrlim_3 %f, uprlim_3 %f, area_3 %f\n",
	 *         *(psc_v + 8), *(psc_v + 9), *(psc_v + 10), *(psc_v + 11), *(plb + 2), *(pub + 2), *(ppa + 2));
	 *   Rprintf("  a_4 %f, b_4 %f, c_4 %f, d_4 %f, lwrlim_4 %f, uprlim_4 %f, area_4 %f\n",
	 *         *(psc_v + 12), *(psc_v + 13), *(psc_v + 14), *(psc_v + 15), *(plb + 3), *(pub + 3), *(ppa + 3));
	 * }
	 */

} // end of f_int_boundaries


void f_min_max_euclidian_distance(
                                  double *psc_v, double *plb, double *pub)
{

	/*
	 * function distances from origin to lower left and upper right corner of
	 * pixel or of a respective part of the pixel after the reflection into the
	 * first quadrant aroung a support point; these distances form the lower
	 * upper integration limits
	 */

  // 2009-12-xx C. Hofer
  // 2025-01-19 A. Papritz Rprintf() substituted for printf()
  //                       editing of comments and code reformatting
  // 2025-01-21 A. Papritz simplifiction because the pixel or its parts lie always in first quadrant

	/*
   * // 2025-01-21 begin changes
	 * int i;
	 * double min, max;
	 * double dist[4];
	 *
	 * dist[0] = pow( pow( (*(psc_v + 0) - 0), 2 ) + pow( (*(psc_v + 2) - 0), 2), 0.5 ); // sqrt(a^2+c^2)
	 * dist[1] = pow( pow( (*(psc_v + 0) - 0), 2 ) + pow( (*(psc_v + 3) - 0), 2), 0.5 ); // sqrt(a^2+d^2)
	 * dist[2] = pow( pow( (*(psc_v + 1) - 0), 2 ) + pow( (*(psc_v + 3) - 0), 2), 0.5 ); // sqrt(b^2+d^2) sollte wohl sqrt(b^2+c^2) sein
	 * dist[3] = pow( pow( (*(psc_v + 1) - 0), 2 ) + pow( (*(psc_v + 3) - 0), 2), 0.5 ); // sqrt(b^2+d^2)
	 *
	 * min = dist[1];
	 * max = dist[1];
	 *
	 * for (i = 0; i < 4; i++)
	 * {
	 *   if (dist[i] >= max)
	 *   {
	 *     max = dist[i];
	 *   }
	 *   else if (dist[i] <= min)
	 *   {
	 *     min = dist[i];
	 *   }
	 * } // 2025-01-21 end change
	 */

  double min, max;

  min = pow( pow( (*(psc_v + 0) - 0), 2 ) + pow( (*(psc_v + 2) - 0), 2), 0.5 ); // sqrt(a^2+c^2)
  max = pow( pow( (*(psc_v + 1) - 0), 2 ) + pow( (*(psc_v + 3) - 0), 2), 0.5 ); // sqrt(b^2+d^2)

  *(plb) = min;
  *(pub) = max;

}
