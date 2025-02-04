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

// Integrieren der Kovarianzfunktionen zwischen den Stützpunkten
// und Rechtecken
//
// Author: Christoph Hofer
// Datum: Dezember 2009

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

#include "covariancefunctions.h"
#include "intboundaries.h"
#include "scalePixelCoord.h"
#include "freqdist.h"

#define datasupport pow(10, -6)  // pixels with areas smaller than datasupport are considered as points

// definition of structure to pass centred coordinates of pixel (or its parts) and
// extra covariance parameters to integrand function
typedef struct pcp
{
  double sc_pixcoords[4];  // coordinates of a pixel centered before by the coordinates of a support point
  int ncp;                 // number of extra covariance parameters, not used
  double *cp;              // pointer to the extra covariance parameters

} pcp_type;

// main C function called from R
void PointRectCov(
  double *pxcoord,  // pointer to x-coordinates of all support points
  double *pycoord,  // pointer to y-coordinates of all support points
  double *crxcoord, // pointer to x-coordinates of all pixel centers
  double *crycoord, // pointer to y-coordinates of all pixel centers
  double *rxwidth,  // pointer to width of pixel in x-directions
  double *rywidth,  // pointer to height of pixel in y-directions
  double *npoints,  // pointer to number of support points
  double *npixel,   // pointer to number of pixels
  double *epsabs,   // pointer to absolute error accuracy for integration
  double *epsrel,   // pointer to absolute error accuracy for integration
  double *result,   // pointer to results
  double *abserr,
  int *limit,
  int *neval,
  int *ier,
  int *lenw,
  int *last,
  int *iwork,
  double *work,
  void *param,    // pointer to array of stacked covariance model and covariance parameters
                  // {modelno, varriance, scale, extra parameters, modelno, variance, scale, extra parameters, ...}
  int *nmod,      // pointer to number of stacked covariance models
  int *nmodparam);// pointer to array with numbers of extra parameters of stacked covariance models

// integrand function
void f_integrate( double *x, int n, void *p_pcp );

// pointer to covariance function that is used in integrand
double (*f_cov)( double x,  double *covpar);

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
  void *param,
  int *nmod,
  int *nmodparam)
{

  // 2009-12-xx C. Hofer
  // 2023-01-24 A. Papritz elimination of unnneeded variable sc_v
  // 2025-01-26 A. Papritz various variables renamed, added dimension (16) in definition of
  //                       array sc_pixcoords, corrected update of index variable imodr
  //                       ("imodr = imodr + 3 + *nmodparam;" replaced by "imodr = imodr + 3 + *(nmodparam + kmod);")
  //                       minor changes and editing of comments
  // 2025-02-03 A. Papritz extended length of arrays setlowerbound and setupperbound to avoid stack-buffer-overflow
  //                       in while loop over all pixel parts
  
  unsigned long int ipixel; // index variable in loop over all pixels
  unsigned long int jpoint; // index variable in loop over all support points

  int ijcov;       // index variable for results array
  int q;           // index variable for array for centred coordinates
                   // of lower left and upper right corners of pixel parts
  int kmod;        // index variable for stacked covariance models
  int lpixpart;    // index variable for pixel parts
  int modnr;       // model number of covariance model
  int imodr;       // index variable for model number and parameters in array of stacked covariance models

  double cov_current_model;  // point-block average of covariance for current model
  double cov_all_models;     // point-block average of covariance summed up over all covariance models
  double rect_area;          // area of pixel
  double pcxdiff, pcydiff;   // differences of x- and y-coordinates between support point and pixel center

  // array with centred coordinates of pixel (or its parts)
  double sc_pixcoords[16] = {-1,-1,-1,-1,
              -1,-1,-1,-1,
              -1,-1,-1,-1,
              -1,-1,-1,-1};
  double *p_sc_pixcoords = NULL;
  p_sc_pixcoords = &sc_pixcoords[0];

  // pointer to array of stacked covariance model and covariance parameters
  double *param_new = NULL;
  param_new = (double*)param; // cast void pointer in double pointer

  double res;            // pointer to result passed from integration function
  double *pres = NULL;

  // 5 elements to avoid stack-buffer-overflow in following while statement
  double setlowerbound[5] = {-1,-1,-1,-1,-1}; // array for lower integration limit(s) of pixel, 
  double setupperbound[5] = {-1,-1,-1,-1,-1}; // array for lower integration limit(s) of pixel (or its parts)

  // pointers to array for lower and upper integration limits
  double *plb, *pub;
  plb = &setlowerbound[0];
  pub = &setupperbound[0];

  // pointer to area of pixel (or its parts)
  double pixelarea[4] = {-1,-1,-1,-1};
  double *ppa;
  ppa = &pixelarea[0];

  // pcp structure of pcp_type
  // p_pcp pointer on pcp
  pcp_type pcp, *p_pcp;
  p_pcp = &pcp;

  // initialize values of centred coordinates of pixel corners in structure
  (*p_pcp).sc_pixcoords[0] = 0;
  (*p_pcp).sc_pixcoords[1] = 0;
  (*p_pcp).sc_pixcoords[2] = 0;
  (*p_pcp).sc_pixcoords[3] = 0;

  // initialize values of covariance model and parameters in structure
  (*p_pcp).cp = param_new;

  // compute area of pixels

  rect_area = (*rxwidth) * (*rywidth);
	/*
	 * Rprintf("PointRectCov:area of pixel %f \n", rect_area);
	 */

  ijcov = 0; // initialize index variable for result array


  // loop over all pixels

  for(ipixel=0; ipixel < *npixel; ipixel++)
  {

    // loop over all support points

    for(jpoint=0; jpoint < *npoints; jpoint++)
    {

      // center cooordinates of pixel corners by current support point
      // pxcoord, pycoord pointers to coordinates of current support points
      // crxcoord, crycoord pointers to coordinates of center of current pixel
      // p_sc_pixcoords pointer to coordinates of corners of pixel centred by current support point

			/*
			 * Rprintf("\nPointRectCov:processing pixel %lu and support point %lu\n", ipixel + 1, jpoint + 1);
			 */

      scale_pixel_coord((pxcoord+jpoint), (pycoord+jpoint),
                        (crxcoord + ipixel), (crycoord+ipixel),
                        rxwidth, rywidth, p_sc_pixcoords);

      // initialization for loop over stacked covariance models

      imodr = 0;   // index variable for array of covariance model number and covariance parameters
      cov_all_models = 0;

      // initialization of result

      res = 0;
      pres = &res;

      // loop over all covariance models

      for(kmod = 0; kmod < *nmod; kmod++)
      {

        // initializations for current covariance model

        cov_current_model = 0;              // integreated point-block covariance
        modnr = (int) *(param_new + imodr); // model number

				/*
				 * Rprintf(
				 *   "\nPointRectCov:current model %d: no %d, variance %f, scale %f, extra parameter %f\n",
				 *   kmod, modnr, *(param_new + imodr + 1), *(param_new + imodr + 2), *(param_new + imodr + 3)
				 * );
				 */

        // regularize nugget effect or prepare integration of covariance function

        switch(modnr)
        {
          case 0: // nugget

          // regularization of nugget effect

          // distance between support point and pixel center
          pcxdiff = fabs( *(pxcoord + jpoint) - *(crxcoord + ipixel) );
          pcydiff = fabs( *(pycoord + jpoint) - *(crycoord + ipixel) );

          if( pcxdiff > 0.5 * (*rxwidth) || pcydiff > 0.5 * (*rywidth) )
          {
            // support point lies outside of pixel
            // cf Journel & Huijbregts, 1981, Mining geostatistids, eq. III.6c, p. 154
           cov_current_model = 0;
          }
          else
          {
            // support point lies inside pixel
            if (rect_area > datasupport)
            {
              // pixel area larger than data support: regularization of nugget
              // cf Journel & Huijbregts, 1981, Mining geostatistids, eq. III.6b, p. 154
              cov_current_model = *(param_new + imodr + 1) * datasupport / rect_area;
            }
            else
            {
              // pixel and data support identical,
              // cf Journel & Huijbregts, 1981, Mining geostatistids, eq. III.6a, p. 154
              cov_current_model = *(param_new + imodr + 1);
            }
          }

          break;

          case 1: // stable or powered exponential (including exponential and gauss
          f_cov = f_cov_exponential;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 2: // spherical
           f_cov = f_cov_sphercial;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 3: // matern
          f_cov = f_cov_matern;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 4: // bessel
          f_cov = f_cov_bessel;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 5: // cauchy
          f_cov = f_cov_cauchy;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 6: // cauchytbm
          f_cov = f_cov_cauchytbm;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 7: // circular
          f_cov = f_cov_circular;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 8: // constant
          f_cov = f_cov_constant;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 9: // cubic
          f_cov = f_cov_cubic;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 10: // dampedcosine
          f_cov = f_cov_dampedcosine;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 11: // gencauchy
          f_cov = f_cov_gencauchy;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 12: // gengneiting;  a= 1
          f_cov = f_cov_gengneiting1;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 13: // gengneiting;  a= 2
          f_cov = f_cov_gengneiting2;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 14: // gengneiting;  a= 3
          f_cov = f_cov_gengneiting3;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 15: // gneiting
          f_cov = f_cov_gneiting;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 16: // hyperbolic
          f_cov = f_cov_hyperbolic;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 17: // penta
          f_cov = f_cov_penta;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 18: // lgd1
          f_cov = f_cov_lgd1;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 19: // power
          f_cov = f_cov_power;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 20: // wave
          f_cov = f_cov_wave;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 21: // qexponential
          f_cov = f_cov_qexponential;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          case 22: // whittle
          f_cov = f_cov_whittle;
          (*p_pcp).cp = param_new + imodr + 1;
          break;

          default:
          // Rprintf("Kein gültiges Kovarianzmodell.");
          cov_current_model = -999;

        } // end switch(modnr)


        // compute point-block average of covariance

        if (modnr != 0) // point-block average of nugget covariance has already been computed
        {

            // reflect pixel (or its 2 or 4 parts) into first quadrant around support point
            // and compute for pixel (or each its parts) lower and upper integrations limits
            // and area of the part
            // plb pointer to lower integration limits of pixel or its parts
            // pub pointer to upper integration limits of pixel or its parts
            // ppa pointer to areas of pixel or its parts

            f_int_boundaries(p_sc_pixcoords, plb, pub, ppa);

            // loop over all parts into which pixel was subdivided

            lpixpart = 0;       // index variable for parts of pixel
                                // of lower left and upper right corners of pixel parts

            while( *(plb + lpixpart) >= 0 && *(pub + lpixpart) >= 0 && lpixpart < 4)
            {

              q = 4 * lpixpart;   // index variable for array for centred coordinates

              (*p_pcp).sc_pixcoords[0] = *(p_sc_pixcoords + q + 0);
              (*p_pcp).sc_pixcoords[1] = *(p_sc_pixcoords + q + 1);
              (*p_pcp).sc_pixcoords[2] = *(p_sc_pixcoords + q + 2);
              (*p_pcp).sc_pixcoords[3] = *(p_sc_pixcoords + q + 3);

							/*
							 * Rprintf( "\nPointRectCov:processing pixel part %d\n", lpixpart+1 );
							 */

              // integration

              Rdqags(f_integrate, (void*)p_pcp, (plb + lpixpart), (pub + lpixpart),
                     epsabs, epsrel,
                     pres, abserr, neval, ier,
                     limit, lenw, last,
                     iwork, work);


              // set lower and upper integrations limits for the parts that were used back
              // to initial values; this is important because plb and pub are used to exit from while loop

              *(plb + lpixpart) = -1;
              *(pub + lpixpart) = -1;

              // sum up integrated covariance for pixel parts (weighted by their area)

              cov_current_model = cov_current_model + (*pres) * (*(ppa + lpixpart));

							/*
							 * Rprintf("  average covariance = %f weighted sum = %f\n", *pres, cov_current_model);
							 */

              // increase index variable for parts of pixel by 1

              lpixpart++;

            } // end while lpixpart loop

            // compute weighted mean of the integrated covariance for the pixel parts

            cov_current_model =  cov_current_model / rect_area;

          } // end if (modnr != 0)


          // sum up mean integrated covariances over all models
          // possibly including also regularized nugget effect

          cov_all_models = cov_all_models + cov_current_model;

					/*
					 * Rprintf("  average covariance current model = %f accumulated average covariance = %f\n", cov_all_models, cov_all_models);
					 */

          // increase index variable variable for array of covariance parameters

          imodr = imodr + 3 + *(nmodparam + kmod);

        } // end for kmod loop

        // assing point-block average of covariance to result array

        *(result + ijcov) = cov_all_models;

        // increase index variable for result array by 1

        ijcov++;

      } // end for jpoint loop

   } // end for ipixel loop

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
