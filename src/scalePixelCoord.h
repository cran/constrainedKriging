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
//

#include <stdio.h>

void scale_pixel_coord(
             double *pxcoord,   // x-coordinate of one support point
             double *pycoord,   // y-coordinate of onesupport point
             double *crxcoord,  // x-coordinate of center of one pixel
             double *crycoord,  // y-coordinate of center of one pixel
             double *rxwidth,   // width of pixel in x-direction
             double *rywidth,   // height of pixel in y-direction
             double *psc_v);    // point on array with centred coordindates of pixel corners

void scale_pixel_coord(
            double *pxcoord,
            double *pycoord,
            double *crxcoord,
            double *crycoord,
            double *rxwidth,
            double *rywidth,
            double *psc_v)
{

  // Centering the coordinates of the lower left and upper right corners of a pixel
  // with the coordinates of a single support point
  // the function computes and returns (in this order)
  //   centred x-coordinate of the lower left corner,
  //   centred x-coordinate of the upper right corner
  //   centred y-coordinate of the lower left corner,
  //   centred y-coordinate of the upper right corner

  // 2009-12-xx C. Hofer
  // 2025-01-17 A. Papritz editing of comments and substituting Rprintf(...) for printf(...)

	/*
	 * Rprintf("\nscale_pixel_coord:pxcoord %f \n", *pxcoord);
	 * Rprintf("scale_pixel_coord:pycoord %f \n", *pycoord);
	 * Rprintf("scale_pixel_coord:crxcoord %f \n", *crxcoord);
	 * Rprintf("scale_pixel_coord:crycoord %f \n", *crycoord);
	 */

  // centering coordinates of corners of pixel

  *(psc_v + 0)  = (*crxcoord - 0.5 * (*rxwidth)) - (*pxcoord); // centred x-coordinate of the lower left corner (a)
  *(psc_v + 1)  = (*crxcoord + 0.5 * (*rxwidth)) - (*pxcoord); // centred x-coordinate of the upper right corner (b)
  *(psc_v + 2)  = (*crycoord - 0.5 * (*rywidth)) - (*pycoord); // centred y-coordinate of the lower left corner (c)
  *(psc_v + 3)  = (*crycoord + 0.5 * (*rywidth)) - (*pycoord); // centred y-coordinate of the upper right corner (d)


	/*
	 * Rprintf("\nscale_pixel_coord:a %f \n", *(psc_v + 0));
	 * Rprintf("scale_pixel_coord:b %f \n", *(psc_v + 1));
	 * Rprintf("scale_pixel_coord:c %f \n", *(psc_v + 2));
	 * Rprintf("scale_pixel_coord:d %f \n", *(psc_v + 3));
	 */

} // end scale_pixel_coord
