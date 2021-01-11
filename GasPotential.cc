 //-----------------------------------------------------------------------------+
//                                                                             |
// GasPotential.cc                                                             |
//                                                                             |
// Copyright (C) 2004-2006 Walter Dehnen                                       |
//               2015      Bekdaulet Shukirgaliyev			       |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Versions                                                                    |
// 0.1    09/08/2006  WD use $NEMOINC/defacc.h                                 |
//        18/11/2015  BS                                                       |
//        20/08/2018  BS  formulation of rg3() function is shortened           |
//-----------------------------------------------------------------------------+
#define POT_DEF
#include <cmath>
#include <defacc.h> // $NEMOINC/defacc.h

#define SQR(x) x*x
#define Pi 3.141592653589793

double M_sc = 1.0, a_pl = 1.0; //mass and scale radius of a star cluster
double K, RmaxGasPot;
double G_const = 1.0, e_ff, t_sf;

      
// rs() - calculates density of stars at given r (distance from the centre) for a Plummer sphere
double rs(double const r){
	double rhos = (3.0 * M_sc)/(4.0 * Pi * a_pl*a_pl*a_pl)*pow((1 + r*r/(a_pl*a_pl)),(-(5./2.)));
	return rhos;
    }

// rg3() - calculates density of residual gas at given r (distance from the centre)
//         and parameter k, which calculates as
//         K = sqrt(8.0 * G_const/(3.0 * Pi)) * e_ff * t_sf;
double rg3(double const r, double const k){

	double k4 = k*k*k*k;
	double k6 = k*k*k*k*k*k;
	double aa = k4 * rs(r) * rs(r);

	double K0 = pow(( aa*aa*aa + 36*aa*aa + 216*aa + 24*aa*sqrt(3*(aa + 27))  ),(1.0/3.0));
	double K1 = sqrt( (aa*aa + aa*(K0 + 24) + K0*(K0 + 12))/(12*k4*K0) );
	double K2 = (aa - K0 + 24)*(K0 - aa)/(3*k4*K0);

    double rhog=  1.0/(k*k) - 0.5*rs(r) + K1 - 0.5*sqrt(K2 + 8.0/(k6*K1) );
    return rhog;
}


////////////////////////////////////////////////////////////////////////////////
namespace {

  class GasPotential {

  public:
    static const char* name() { return "GasPotential"; }
    GasPotential(const double*pars,
	    int          npar,
	    const char  *file)
    {
      if(npar < 4)
	warning("%s: recognizing 4 parameters:\n"
		" e_ff       star formation efficiency per free-fall time	[0.05]\n"
		" t_sf	     duration of star formation phase [NB]		[10.0]\n"
		" a_pl	     Plummer scale radius [NB]     			[1.0]\n"
		" M_sc	     Mass of stellar cluster [NB]			[1.0]\n"
		"\n\n",name());
      if(file && file[0])
	warning("%s: file \"%s\" ignored",name(),file);
      double
	eff   = npar>0? pars[0] : 0.05,
	tsf   = npar>1? pars[1] : 10.0,
	apl   = npar>2? pars[2] : 1.0,
	mpl   = npar>3? pars[3] : 1.0;
	
      if (!((eff < 1.0))&&(eff > 0))
	error("e_ff value is out of range: %f\n"
	      "\t\t\t  correct range is : 0.0 < e_ff < 1.0",eff);
      if (tsf <=0.0)
	error("t_sf value is out of range: %f\n"
	      "\t\t\t  correct range is : 0.0 < t_sf",tsf);
    RmaxGasPot = 32.0; 
      e_ff    = eff;
      t_sf    = tsf;
      M_sc = mpl;
      a_pl = apl;
      
      if(npar > 4) warning("%s: skipped parameters beyond 4\n"
			   " parameters : e_ff  = %f\n"
		           "              t_sf  = %f\n",   
		           "              a_pl  = %f\n",   
		           "              M_sc  = %f\n",   
			   e_ff, t_sf, a_pl, M_sc, name());
      nemo_dprintf (1,
		    "initializing %s\n"
		    " parameters : e_ff  = %f\n"
		    "              t_sf  = %f\n",
		    "              a_pl  = %f\n",   
		    "              M_sc  = %f\n",   
		    name(), e_ff, t_sf, a_pl, M_sc);
    }
    template<typename scalar>
    void potacc(scalar const&Rq,
		scalar      &P,
		scalar      &T) const
    {
      register scalar R = std::sqrt(Rq);
      
      double ff_in, ff_out,
	    R_min, R_max, R_len,
	    R_len_1_2,
	    R_len_1_4,
	    R_in_min,
	    R_out_max,
	    R_in_max,
	    R_out_min;

      int i1, i2, ii, bins, indx = 512;

      double dr_in, dr_out,
	     R_in, R_out;
	    
      R_min = 0.0;
      R_max = RmaxGasPot;

      R_len = R_max-R_min;
      R_len_1_2 = R_len/2.0;
      R_len_1_4 = R_len/4.0;

      K = sqrt(8.0 * G_const/(3.0 * Pi)) * e_ff * t_sf;
// check whether we close to the center or to the outskirts in order to find most efficient binning
      if (R < R_len_1_4){
	  i1 = 1 * indx;
	  i2 = 3 * indx;
      }
      else if (R > R_len_1_2 + R_len_1_4){
	  i1 = 3 * indx;
	  i2 = 1 * indx;
      }
      else{
	  i1 = 2 * indx;
	  i2 = 2 * indx;
      }

      R_in_min = R_min;
      R_out_max = R_max;
  
      if (!(R > R_max)) {
	R_in_max  = R;
	R_out_min = R;
      }
      else{
	R_in_max  = R_max;
	ff_out = 0.0;
        nemo_dprintf (1,
		    "warning in %s\t"
		    "max_R < R\n"
		    ,name());
	goto in_edge;
      }
      
      ii = i2;
      bins = 1+ii*2;

      dr_out = 1.0*(R_out_max - R_out_min)/(bins-1.0);
      ff_out= 0.0;
      for (int i = 0; i < bins; i++){
	  R_out = dr_out * i + R_out_min;
	  if ((i > 0)&&(i < bins-1)){
	    if (i % 2 == 0)
	      ff_out += 2 * rg3(R_out,K) * (R_out);
	    else
	      ff_out += 4 * rg3(R_out,K) * (R_out);
	  }
	  else
	    ff_out += rg3(R_out,K) * (R_out);
      }  
      ff_out *= dr_out / 3.0 * 4.0 * Pi;

in_edge:      
      ii = i1;
      bins = 1+ii*2;

      dr_in = 1.0*(R_in_max - R_in_min)/(bins-1.0);

      ff_in = 0.0;
      for (int i = 0; i < bins; i++){
	  R_in = dr_in * i + R_in_min;
	  if ((i > 0)&&(i < bins-1)){
	    if (i % 2 == 0)
	      ff_in += 2 * rg3(R_in,K) * SQR(R_in);
	    else
	      ff_in += 4 * rg3(R_in,K) * SQR(R_in);
	  }
	  else
	    ff_in += rg3(R_in,K) * SQR(R_in);
      }  
      ff_in *= dr_in / 3.0 * 4.0 * Pi;

      
	    
//	ff_in - cumulative mass of gas inside R
//              that is an integral :  Int_0^R[4 Pi r^2  \rho_{gas}(r)  dr]
//	ff_out - an integral : Int_R^R_max[4 Pi r  \rho_{gas}(r)  dr]
      P = - (G_const * ff_in / R + G_const * ff_out);
      T = - G_const * ff_in / (R*R*R);
//      nemo_dprintf (1,"R = %.8E  P = %.8E   T = %.8E\n",R,P,T);
    }
  };
}
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<GasPotential>)
__DEF__POT(SphericalPot<GasPotential>)

//------------------------------------------------------------------------------
