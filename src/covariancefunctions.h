/*
 Ckrige is a program for universal, constrained and covariance-matching
 constrained point and block kriging. The function name of the the covariance
 functions refer to the names used for the covariance functions in the Random Fields package.
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
//
# include<R.h>
# include<Rmath.h>
// ***************************************************
// qexponential : C(x) = (2 exp(-x)-a exp(-2x))/(2-a)
// The parameter a takes values in [0,1].
// ***************************************************
double f_cov_qexponential(double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// wave: C(x)=sin(x)/x if x>0 and C(0)=1
// This isotropic covariance function is valid only for
// dimensions less than or equal to 3. It is a special
// case of the bessel model (for a=0.5).
// ***************************************************
double f_cov_wave(double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// power: C(x) = 1-x)^a if 0<=x<=1, 0 otherwise
// This covariance function is valid for dimension d if
// a >= (d+1)/2. For a=1 we get the well-known triangle
// (or tent) model, which is valid on the real line, only.
// ***************************************************
double f_cov_power(double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// local-global distinguisher
// lgd1: C(x)= 1- b(a+b)^{-1}|x|^a for |x| <=  1 and 
//   C8x) = a(a+b)^{-1}|x|^-b for |x|> 1  
// Here b>0 and a is in (0, 1.5-d/2] for dimension d=1,2.
// The random field has fractal dimension d + 1 - a/2 and
// Hurst coefficient 1 - b/2 for b in (0,1]
// ***************************************************
double f_cov_lgd1(double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// penta: C(x) = 1 - 22/3 x^2 +33 x^4 - 77/2 x^5 + 33/2 x^7 - 11/2 x^9 + 5/6 x^11 if 0<=x<=1, 0 otherwise
// valid only for dimensions less than or equal to 3.
// This is a 4 times differentiable covariance functions with compact support
// ***************************************************
double f_cov_penta(double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// hyperbolic:
// C(x)= c^(-b) (K_b(a*c))^(-1) * (c^2 +x^2)^(0.5 b) * K_b(a sqrt(c^2 + x^2))
//K_b = bessel_k(x ,b )
//The parameters are such that
//c>=0, a>0 and b>0, or
//c>0 , a>0 and b=0, or
///c>0 , a>=0, and b<0.
// Note that this class is over-parametrised; always one of the three parameters a, c,
// and scale can be eliminated in the formula. Therefore, one of these parameters
// should be kept fixed in any simulation study.  The model contains as special cases
// the matern model and the cauchy model, for c=0 and a=0, respectively.
// ***************************************************
double f_cov_hyperbolic(double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// gneiting:
// C(x)= (1 + 8 s x + 25 s^2 x^2 + 32 s^3 x^3)*(1-s x)^8 if 0 <= x <= 1/s, 0 otherwise
// s = 0.301187465825
// ***************************************************
double f_cov_gneiting(double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// gengneiting (generalised gneiting):
// if a = 1: C(x)=[1 + (b+1) * x] * (1-x)^(b+1) if 0<=x<=1, 0 otherwise
// => f_cov_gengneiting1
// if a = 2: C(x)= [1 + (b+2) * x + ((b+2)^2-1) * x^2 / 3] * (1-x)^(b+2) if 0<=x<=1, 0 otherwise
// => f_cov_gengneiting2
// if a = 3: C(x)=[1 + (b+3) * x + (2 * (b+3)^2 - 3) * x^2 / 5 + ((b+3)^2 - 4) * (b+3) * x^3 / 15] * (1-x)^(b+3) if 0<=x<=1, 0 otherwise
// => f_cov_gengneiting3
// The parameter b is greater than or equal to (d + 2a +1)/2 where d is the dimension of the random field.
// ***************************************************
double f_cov_gengneiting1(double x, double *covpar);
double f_cov_gengneiting2(double x, double *covpar);
double f_cov_gengneiting3(double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// gencauchy (generalised cauchy): C(x)= (1+x^a)^(-b/a)
// The parameter a is in (0,2], and b is positive. 
// ***************************************************
double f_cov_gencauchy(double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// dampedcosine: C(x) = C(x)= exp(-a x) cos(x)
// This model is valid for dimension 1 iff a>=0,
//for dimension 2 iff a>=1, and for dimension 3 iff a >= sqrt(3).
// ***************************************************
double f_cov_dampedcosine(double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// cubic:
// C(x)= 1- 7 x^2 + 8.75 x^3 - 3.5 x^5 + 0.75 x^7 if 0<=x<=1,
// 0 otherwise//
// ***************************************************
double f_cov_cubic( double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// constant: C(x) = sill
// ***************************************************
double f_cov_constant( double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// circular = C(x)=1-2/pi*(x sqrt(1-x^2)+asin(x)) if 0<=x<=1, 0 otherwise
// This isotropic covariance function is valid only for
// dimensions less than or equal to 2.
// ***************************************************
double f_cov_circular( double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// cauchytbm: C(x)= (1+(1-b/c)x^a)(1+x^a)^(-b/a-1)
// The parameter a is in (0,2] and b is positive. The model
// is valid for dimensions d<=c; this has been shown for
// integer c, but the package allows real values of c. 
// It allows for simulating random fields where fractal
// dimension and Hurst coefficient can be chosen independently.
// It has negative correlations for b>c and large x
// 2023-11-20 A. Papritz parameter c set equal to 3                          
// ***************************************************
double f_cov_cauchytbm( double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// cauchy Kovarianzfunktion: C(x)=(1+x^2)^(-a)
// The parameter a is positive. The model possesses two generalisations, the gencauchy model and the hyperbolic model.
// ***************************************************
double f_cov_cauchy( double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// bessel Kovarianzfunktion: c(x) = 2^aGamma(a+1)x^(-a)J_a(x)
// The parameter a is greater than or equal to (d-2)/2, where d is the dimension of the random field.
// ***************************************************
double f_cov_bessel( double x, double *covpar);
// ***************************************************
//
//
//
// ***************************************************
// pure_nugget Kovarianzfunktion
// ***************************************************
double f_cov_nugget( double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// exponentielle Kovarianzfunktion
// ***************************************************
double f_cov_exponential( double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// sphaerisches Kovarianzfunktion
// ***************************************************
double f_cov_sphercial( double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// matérn Kovarianzfunktion
// ***************************************************
double f_cov_matern( double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// whittle Kovarianzfunktion
// ***************************************************
double f_cov_whittle( double x, double *covpar);
// ***************************************************
//
//
// ***************************************************
// exponentielle Kovarianzfunktion
// ***************************************************

// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = 0 <= Kappa <= 2
// x = lad Distanz
double f_cov_exponential( double x,  double *covpar) 
{ 
	
	double exp_cov;

	exp_cov =   (*(covpar + 0)) * exp( - (  pow( ( x / *(covpar + 1) ), *(covpar + 2) ) ) );
		
return(exp_cov);
} 


// ***************************************************
// sphearische  Kovarianzfunktion
// ***************************************************

// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// x = lad Distanz
double f_cov_sphercial( double x,  double *covpar) 
{ 
	
	double sph_cov;
	
	if( x <= ( *(covpar + 1))  )
	{
		sph_cov =   ( *(covpar + 0) ) * ( 1 - (1.5  * ( x / ( *(covpar + 1) )  ) - 0.5 * pow( ( x / ( *(covpar + 1) ) ), 3 ) ) );
	} 
	else
	{
		sph_cov = 0;
		
	}
	return(sph_cov);
} 

// ***************************************************
// matérn (bessel) Kovarianzfunktion
// ***************************************************

// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = Delta > 0
// gammafn(a) Gammafunktion in Rmath.h
// bessel_k(x,n,1) modifizierte besselfunction zweiter Art der Ordnung n
double f_cov_matern( double x, double *covpar)
{ 
	
	double mat_cov;

//mat_cov = *(covpar + 0) * ( 1 / ( pow( 2, ( *(covpar + 2)  - 1 ) ) * gammafn( *(covpar + 2) ) ) ) * pow( ( x / *(covpar + 1) ),  *(covpar + 2) ) * bessel_k( x / *(covpar + 1), *(covpar + 2), 1 );
    mat_cov = *(covpar + 0)* pow( 2, ( 1- *(covpar + 2) ) ) * pow( gammafn( *(covpar + 2) ), -1) * pow( sqrt( *(covpar + 2)  * 2 ) * (x / *(covpar + 1) ), *(covpar + 2) ) * bessel_k( x / *(covpar + 1), *(covpar + 2), 1 );
	

	return(mat_cov);
}

// ***************************************************
// whittle Kovarianzfunktion
// ***************************************************

// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = Delta > 0
// gammafn(a) Gammafunktion in Rmath.h
// bessel_k(x,n,1) modifizierte besselfunction zweiter Art der Ordnung n
double f_cov_whittle( double x, double *covpar)
{
	
	double whi_cov;
        whi_cov = *(covpar + 0)* pow( 2, ( 1- *(covpar + 2) ) ) * pow( gammafn( *(covpar + 2) ), -1) * pow( x / *(covpar + 1), *(covpar + 2) )  * bessel_k( x / *(covpar + 1), *(covpar + 2), 1 );
	return(whi_cov);
}

// ***************************************************
// pure_nugget Kovarianzfunktion
// ***************************************************
double f_cov_nugget( double x, double *covpar);
// ***************************************************
// ***************************************************
// exponentielle Kovarianzfunktion
// ***************************************************

// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = 0 <= Kappa <= 2
// x = lad Distanz
double f_cov_nugget( double x,  double *covpar) 
{ 
	
	double nugget_cov;
	
	if( x > 0)
	{
		nugget_cov = 0;
	}
	else
	{
		nugget_cov =  *(covpar + 0);
	}
	
	return(nugget_cov);
} 
// ***************************************************
// bessel Kovarianzfunktion c(x) = 2^aGamma(a+1)x^(-a)J_a(x)
// The parameter a is greater than or equal to (d-2)/2, where d is the dimension of the random field.
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = parameter a
// x = lad Distanz
// ***************************************************
double f_cov_bessel( double x, double *covpar)
{
	double bessel_cov;
	
	bessel_cov = *(covpar + 0) * pow(2, *(covpar + 2)) * gammafn( *(covpar + 2) + 1 ) * pow(x / *(covpar + 1), -1*( *(covpar + 2) ) ) * bessel_j(x / *(covpar + 1), *(covpar + 2) );
	
	return(bessel_cov );
	
}
// ***************************************************
// cauchy Kovarianzfunktion c(x) = C(x)=(1+x^2)^(-a)
// The parameter a is positive. The model possesses two generalisations, the gencauchy model and the hyperbolic model.
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = parameter a
// x = lad Distanz
// ***************************************************
double f_cov_cauchy( double x, double *covpar)
{
	double cauchy_cov;
	
	cauchy_cov =  *(covpar + 0) * pow( 1 + pow( x / *(covpar + 1) ,2 ), -1* ( *(covpar + 2) ) ) ; 
	
	return(cauchy_cov);
}
// ***************************************************
// cauchytbm = C(x)= (1+(1-b/c)x^a)(1+x^a)^(-b/a-1)
// a = (0,2] b must be positive
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = parameter a real
// *(covpar + 3) = parameter b real 
// *(covpar + 4) = parameter c real (2023-11-20, A.Papritz, parameter set equal to 3)
// x = lad Distanz
// ***************************************************
double f_cov_cauchytbm( double x, double *covpar)
{
	double cauchytbm_cov;
// 	cauchytbm_cov = *(covpar + 0) * (1 + ( 1- *(covpar + 3) / *(covpar + 4) ) * pow( x / *(covpar + 1), *(covpar + 2)) ) * pow(1 + pow(x / *(covpar + 1), *(covpar + 2)), ( -1 * (*(covpar + 3))/ (*(covpar + 2)) ) -1 );
	cauchytbm_cov = *(covpar + 0) * (1 + ( 1- *(covpar + 3) / 3 ) * pow( x / *(covpar + 1), *(covpar + 2)) ) * pow(1 + pow(x / *(covpar + 1), *(covpar + 2)), ( -1 * (*(covpar + 3))/ (*(covpar + 2)) ) -1 );
	return(cauchytbm_cov);
	
}
// ***************************************************
// circular=
// C(x)=1-2/pi*(x sqrt(1-x^2)+asin(x)) if 0<=x<=1, 0 otherwise
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// x = lad Distanz
// ***************************************************
double f_cov_circular( double x, double *covpar)
{
	double circular_cov;
	
	if( ( x / (*(covpar + 1)) ) >= 0 && ( x / (*(covpar + 1)) ) <= 1)
	{
		circular_cov = *(covpar + 0) * ( 1 - 2/PI * ( x /( *(covpar + 1)) * sqrt(1- pow( x/(*(covpar + 1)), 2 ) ) + asin( x / ( *(covpar + 1) ) ) ) );
	}
	else
	{
		circular_cov = 0;
	}
	return(circular_cov);
}
// ***************************************************
// constant C(x) = sill
// *(covpar + 0) = partial Sill
double f_cov_constant( double x, double *covpar)
{
	double constant_cov;
	constant_cov = *(covpar + 0);
	return( constant_cov );
	
}
// ***************************************************
// cubic:
// C(x)= 1- 7 x^2 + 8.75 x^3 - 3.5 x^5 + 0.75 x^7 if 0<=x<=1,
// 0 otherwise//
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// x = lad Distanz
// ***************************************************
double f_cov_cubic( double x, double *covpar)
{
	double cubic_cov;
	if(  ( x / *(covpar + 1) ) >= 0 && ( x / *(covpar + 1) )  <= 1 )
	{
	cubic_cov = *(covpar + 0) * ( 1 - 7 * pow( x / *(covpar + 1), 2 ) + 8.75 * pow(  x / *(covpar + 1), 3 ) - 3.5 * pow(  x / *(covpar + 1), 5 ) + 0.75 * pow( x / *(covpar + 1), 7 )) ;
	}
	else
	{
		cubic_cov = 0;
	}
	return( cubic_cov );
}
// ***************************************************
// dampedcosine: C(x) = C(x)= exp(-a x) cos(x)
// This model is valid for dimension 1 iff a>=0,
//for dimension 2 iff a>=1, and for dimension 3 iff a >= sqrt(3).
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = parameter a real
// x = lad Distanz
// ***************************************************
double f_cov_dampedcosine(double x, double *covpar)
{
	double dampedcosine_cov;
	
	dampedcosine_cov = *(covpar + 0) *  exp( -1 * ( *(covpar + 2) ) * ( x / (*(covpar + 1) ) ) ) * cos ( x / (*(covpar + 1) )  );
	
	return( dampedcosine_cov );
										
										
	
}
// ***************************************************
// gencauchy (generalised cauchy): C(x)= (1+x^a)^(-b/a)
// The parameter a is in (0,2], and b is positive. 
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = parameter a real; a = (0,2] 
// *(covpar + 2) = parameter b real; b > 0
// x = lad Distanz
// ***************************************************
double f_cov_gencauchy(double x, double *covpar)
{
	double gencauchy_cov;
	gencauchy_cov = *(covpar + 0) * pow( ( 1 + pow(  x / (*(covpar + 1) ), *(covpar + 2) ) )	, -1 * ( *(covpar + 3) / *(covpar + 2) ) );
	return( gencauchy_cov);
}
//
//
// ***************************************************
// gengneiting (generalised gneiting):
// if a = 1: C(x)=[1 + (b+1) * x] * (1-x)^(b+1) if 0<=x<=1, 0 otherwise
// => f_cov_gengneiting1
// if a = 2: C(x)= [1 + (b+2) * x + ((b+2)^2-1) * x^2 / 3] * (1-x)^(b+2) if 0<=x<=1, 0 otherwise
// => f_cov_gengneiting2
// if a = 3: C(x)=[1 + (b+3) * x + (2 * (b+3)^2 - 3) * x^2 / 5 + ((b+3)^2 - 4) * (b+3) * x^3 / 15] * (1-x)^(b+3) if 0<=x<=1, 0 otherwise
// => f_cov_gengneiting3
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = parameter a real; a = (0,2] 
// *(covpar + 3) = parameter b real; b > 0
// x = lad Distanz
// ***************************************************
double f_cov_gengneiting1(double x, double *covpar)
{
	double gengeneiting1_cov;
	if( ( x / *(covpar + 1) ) >= 0 && ( x / *(covpar + 1) ) <= 1 )
	{
	gengeneiting1_cov  = *(covpar + 0) * ( 1 + ( *(covpar + 3) + 1) * (  x / *(covpar + 1) ) ) * pow( (1 - (  x / *(covpar + 1) ) ), ( *(covpar + 3) + 1));
	}
	else
	{
		gengeneiting1_cov = 0;	
	}
	return( gengeneiting1_cov ) ;
}
double f_cov_gengneiting2(double x, double *covpar)
{
	double gengeneiting2_cov;
	if( ( x / *(covpar + 1) ) >= 0 && ( x / *(covpar + 1) ) <= 1 )
	{
		gengeneiting2_cov  = *(covpar + 0) * ( ( 1 + ( *(covpar + 3) + 2) * (  x / *(covpar + 1) ) )  +  
											  (pow( ( *(covpar + 3) + 2 ),2) -1 ) * (pow( ( x / *(covpar + 1) ),2) / 3 ) )
											* pow( (1 - (  x / *(covpar + 1) ) ), ( *(covpar + 3) + 2));
	}
	else
	{
		gengeneiting2_cov = 0;
	}
	return( gengeneiting2_cov ) ;
}
double f_cov_gengneiting3(double x, double *covpar)
{
	double gengeneiting3_cov;
	if( ( x / *(covpar + 1) ) >= 0 && ( x / *(covpar + 1) ) <= 1 )
	{
		gengeneiting3_cov  = *(covpar + 0) * (
											  ( 1 + ( *(covpar + 3) + 3) * (  x / *(covpar + 1) ) )  + 
											  ( 2 *  (pow( ( *(covpar + 3) + 3 ),2) - 3 ) ) * (pow( ( x / *(covpar + 1) ),2) / 5 ) +  
											  ( pow( ( *(covpar + 3) + 3 ),2) - 4 )  * ( *(covpar + 3) + 3 ) * (pow( ( x / *(covpar + 1) ),3) / 15 ))
											* pow( (1 - (  x / *(covpar + 1) ) ), ( *(covpar + 3) + 3) );
	}
	else
	{
		gengeneiting3_cov = 0;
	}
	return( gengeneiting3_cov ) ;
}
// ***************************************************
// gneiting:
// C(x)= (1 + 8 s x + 25 s^2 x^2 + 32 s^3 x^3)*(1-s x)^8 if 0 <= x <= 1/s, 0 otherwise
// s = 0.301187465825
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// ***************************************************
double f_cov_gneiting(double x, double *covpar)
{
	double s, gneiting_cov;
	s =  0.301187465825;
	
	if( ( x / *(covpar + 1) ) >= 0 && ( x / *(covpar + 1) ) <= 1 )
	{
		gneiting_cov = *(covpar + 0) * ( 1 + 8 * s * ( x / *(covpar + 1) )  + 
										25 * pow(s,2) * pow( x / *(covpar + 1) , 2 ) + 
										32 *pow( s, 3 ) *pow( x / *(covpar + 1 ) ,3 ) )  * pow( 1 - s * ( x / *(covpar + 1) ) , 8) ;
	}
	else
	{
		gneiting_cov = 0;
	}
	return( gneiting_cov );
}
// ***************************************************
// hyperbolic:
// C(x)= c^(-b) (K_b(a*c))^(-1) * (c^2 +x^2)^(0.5 b) * K_b(a sqrt(c^2 + x^2))
//K_b(x) = bessel_k(x ,b )
//The parameters are such that
//c>=0, a>0 and b>0, or
//c>0 , a>0 and b=0, or
///c>0 , a>=0, and b<0.
// Note that this class is over-parametrised; always one of the three parameters a, c,
// and scale can be eliminated in the formula. Therefore, one of these parameters
// should be kept fixed in any simulation study.  The model contains as special cases
// the matern model and the cauchy model, for c=0 and a=0, respectively.
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = parameter a real;
// *(covpar + 3) = parameter b real;
// *(covpar + 4) = parameter c real;
// ***************************************************
double f_cov_hyperbolic(double x, double *covpar)
{
	double hyperbolic_cov;
	hyperbolic_cov = *(covpar + 0) * pow( *(covpar + 4), -1 * (*(covpar + 3)) ) * pow( bessel_k( (*(covpar + 2)) * (*(covpar + 4)),*(covpar + 3),1), -1)
												* pow( pow( *(covpar + 4), 2 ) + pow(x / *(covpar + 1),2), 0.5 * (*(covpar + 3)) ) 
												*  bessel_k( *(covpar + 2) * sqrt(pow(*(covpar + 4),2) + pow( x / *(covpar + 1), 2) ), *(covpar + 3),1);
	return( hyperbolic_cov);
	
}
// ***************************************************
// penta: C(x) = 1 - 22/3 x^2 +33 x^4 - 77/2 x^5 + 33/2 x^7 - 11/2 x^9 + 5/6 x^11 if 0<=x<=1, 0 otherwise
// valid only for dimensions less than or equal to 3.
// This is a 4 times differentiable covariance functions with compact support
// ***************************************************
double f_cov_penta(double x, double *covpar)
{
	double penta_cov;
	
	if( ( x / *(covpar + 1) ) >= 0 && ( x / *(covpar + 1) ) <= 1  )
	{
		penta_cov =  *(covpar + 0) * ( 1 -  22.0/3.0 * pow( x / ( *(covpar + 1) ), 2.0 )  +( 33.0 * pow( ( x / *( covpar + 1 ) ), 4.0 ) )
									   -( (77.0/2.0) * pow( ( x / *( covpar + 1 ) ), 5.0 ) ) +( (33.0/2.0) * pow( ( x / *( covpar + 1 ) ), 7.0 ) )
									   -( (11.0/2.0) * pow( ( x / *( covpar + 1 ) ), 9.0 ) ) +( (5.0/6.0)  * pow( ( x / *( covpar + 1 ) ) ,11.0) ) );
	}
	else
	{
		penta_cov = 0;
	}
	return( penta_cov ) ;
}
// ***************************************************
// local-global distinguisher)
// lgd1: C(x)= 1- b(a+b)^{-1}|x|^a for |x| <=  1 and 
//   C8x) = a(a+b)^{-1}|x|^-b for |x|> 1  
// Here b>0 and a is in (0, 1.5-d/2] for dimension d=1,2.
// The random field has fractal dimension d + 1 - a/2 and
// Hurst coefficient 1 - b/2 for b in (0,1]
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = parameter a real;
// *(covpar + 3) = parameter b real;
// *(covpar + 4) = parameter c real;
// ***************************************************
double f_cov_lgd1(double x, double *covpar)
{
	double lgd1_cov;
	if( ( x / *(covpar + 1) ) <= 1 )  
	{
		lgd1_cov =  *(covpar + 0) * (1 - *(covpar + 3)*pow( (*(covpar + 2) + *(covpar + 3)), -1)* pow( ( x / *(covpar + 1) ), *(covpar + 2) ));
	}
	else
	{
		lgd1_cov = 	*(covpar + 2) * pow( *(covpar + 2) + *(covpar + 3), -1 )* pow( ( x / *(covpar + 1) ), -1* ( *(covpar + 3))  );
	}
	return(lgd1_cov );
}
// ***************************************************
// power: C(x) = 1-x)^a if 0<=x<=1, 0 otherwise
// This covariance function is valid for dimension d if
// a >= (d+1)/2. For a=1 we get the well-known triangle
// (or tent) model, which is valid on the real line, only.
// ***************************************************
double f_cov_power(double x, double *covpar)
{
	double power_cov;
	if( ( x / *(covpar + 1) ) >= 0 && ( x / *(covpar + 1) ) <= 1)
	{
		power_cov = *(covpar + 0) *  pow( 1 - ( ( x / *(covpar + 1) ) ), *(covpar + 2) );
	}
	else
	{
		power_cov = 0;
	}
	return( power_cov );
}
// ***************************************************
// wave: C(x)=sin(x)/x if x>0 and C(0)=1
// This isotropic covariance function is valid only for
// dimensions less than or equal to 3. It is a special
// case of the bessel model (for a=0.5).
// ***************************************************
double f_cov_wave(double x, double *covpar)
{
	double wave_cov;
	
	if( x / (*(covpar + 1)) == 0 )
	{
		wave_cov = *(covpar + 0);
	}
	else // (x / scale parameter)  > 0
	{
		wave_cov = *(covpar + 0) * sin( x / (*(covpar + 1)) ) / ( x / (*(covpar + 1)) );
	}
	return( wave_cov ) ;
}
// ***************************************************
// qexponential : C(x) = (2 exp(-x)-a exp(-2x))/(2-a)
// The parameter a takes values in [0,1].
// ***************************************************
// *(covpar + 0) = partial Sill
// *(covpar + 1) = Rangeparameter
// *(covpar + 2) = parameter a real;
// ***************************************************
double f_cov_qexponential(double x, double *covpar)
{
	double qexponential_cov;
	qexponential_cov = *(covpar + 0) *( 2.0 * exp( -x / (*(covpar + 1)) )- *(covpar + 2) * exp( -2.0 * x / (*(covpar + 1)) )) / (2.0 - *(covpar + 2));
	return( qexponential_cov );
}



