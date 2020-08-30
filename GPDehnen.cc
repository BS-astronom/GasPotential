//-----------------------------------------------------------------------------+
//                                                                             |
// GPDehnen.cc                                                             |
//         
//  Gas Potential for the Family of Dehnen models
//                                                                  |
// Copyright (C) 2020     Bekdaulet Shukirgaliyev			       |
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
//-----------------------------------------------------------------------------+
#define POT_DEF
#include <cmath>
#include <defacc.h> // $NEMOINC/defacc.h
//#include <vector>
//#include <stdio.h>
#define SQR(x) x*x
#define Pi 3.141592653589793
//////////////////////////////////////////////////////////////////////////////////
double M_sc = 1.0, r_sc =1.0;//mass and scale radius of a star cluster
double K, RmaxGasPot;
double G_const = 1.0, e_ff, t_sf;
double eps = 1.0E-4, inner=0.0;

      
// rs() - calculates density of stars at given r (distance from the centre)
/*double ssh(double y){
	return 2.0/(exp(y) + exp(-y));
	}
double pow(double x, int p){
	double z=1;
	for(int i=0;i<p;i++)
		z*=x;
	return z;
}
*/
double rs(double const r, double gamma){
	//double x    = sqrt(r*r + r_c*r_c)/r_sc;
	//double rhos = c * ssh(r/r_t)/(pow(x,inner)*pow((x+1),(4-inner)));
	double rhos = M_sc * (3-gamma)/ (4*Pi) * r_sc/(pow(r,gamma)*pow((r+r_sc),(4-gamma)));
	return rhos;
    }

// rg3() - calculates density of residual gas at given r (distance from the centre)
//         and parameter k, which calculates as
//         K = sqrt(8.0 * G_const/(3.0 * Pi)) * e_ff * t_sf;
double rg3(double const r, double const k){

    double k4  = pow(k,4);
    double k6  = pow(k,6);
    double rs1 = rs(r,inner);
    double aa = k4 * rs1 * rs1;

    double K0 = pow(( aa*aa*aa + 36*aa*aa + 216*aa + 24*aa*sqrt(3*(aa + 27))  ),(1.0/3.0));
	double K1 = sqrt( (aa*aa + aa*(K0 + 24) + K0*(K0 + 12))/(12*k4*K0) );
	double K2 = (aa - K0 + 24)*(K0 - aa)/(3*k4*K0);


    double rhog=  1.0/(k*k) - 0.5*rs1 + K1 - 0.5*sqrt(K2 + 8.0/(k6*K1) );
    return rhog;
}
/*
double ms(double r, double gamma){
      
      double tmp_m = M_sc * pow((r/(r+r_sc)),(3-gamma));

return tmp_m;
}*/
////////////////////////////////////////////////////////////////////////////////
namespace {

  class GPDehnen {


  public:
    static const char* name() { return "GPDehnen"; }
    GPDehnen(const double*pars,
	    int          npar,
	    const char  *file)
    {
      if(npar < 5)
	warning("%s: recognizing 5 parameters:\n"
		" e_ff       star formation efficiency per free-fall time	[0.05]\n"
		" t_sf	     duration of star formation phase [NB]		[1.0]\n"
		" inner      Dehnen profile index				[0]\n"
                " r_s        Scale radius [NB]                                  [1.0]\n "
                " M_s        Cluster stellar mass [NB]                          [1.0]\n"
		"\n\n",name());
      if(file && file[0])
	warning("%s: file \"%s\" ignored",name(),file);
      double
	eff   = npar>0? pars[0] : 0.05,
	tsf   = npar>1? pars[1] : 1.0,
	innr  = npar>2? pars[2] : 0.0,
	r_s   = npar>3? pars[3] : 1.0,
	m_s   = npar>4? pars[4] : 1.0;

      if (!((eff < 1.0))&&(eff > 0))
	error("e_ff value is out of range: %f\n"
	      "\t\t\t  correct range is : 0.0 < e_ff < 1.0",eff);
      if (tsf <=0.0)
	error("t_sf value is out of range: %f\n"
	      "\t\t\t  correct range is : 0.0 < t_sf",tsf);
      if (((innr < 0.0))||(innr >= 3))
	error("inner value is out of range: %f\n"
	      "\t\t\t  correct range is : 0.0 <= e_ff < 3.0",innr);
      
      e_ff    = eff;
      t_sf    = tsf;
      inner   = innr;
      r_sc    = r_s;
      M_sc    = m_s;
      double rhalf = r_s/(pow(2,(1.0/(3.0-inner)))-1.0);
      RmaxGasPot = 500;
      eps*=rhalf;
      warning("%s: next values to parameters assigned:\n"
		" e_ff 	[%f]\n"
		" t_sf	[%f]\n"
		" inner [%f]\n"
		" r_s   [%f]\n"
		" M_s   [%f]\n"
		"\n\n",name(), e_ff, t_sf, inner, r_sc, M_sc);
      
      if(npar > 5) warning("%s: skipped parameters beyond 5",name());
      nemo_dprintf (1,
		    "initializing %s\n"
		    " parameters : e_ff  = %f\n"
		    "              t_sf  = %f\n"
		    "              inner = %f\n"
		    "              r_s   = %f\n"
		    "              M_s   = %f\n",
		    name(), e_ff, t_sf, inner, r_sc, M_sc);
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

      int i1, i2, ii, bins, indx = 1024*16;

      double dr_in, dr_out,
	     R_in, R_out;
	    
restart:
      R_min = eps;
      R_max = RmaxGasPot;

      R_len = R_max-R_min;
      R_len_1_2 = R_len/2.0;
      R_len_1_4 = R_len/4.0;

      K = sqrt(8.0 * G_const/(3.0 * Pi)) * e_ff * t_sf;

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
	R_in_max  = R;// R_max;
	ff_out = 1e-12;
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
      if (ff_in != ff_in) {
         eps *=2;
         goto restart;
         }

      
	    
//	ff_in - cumulative mass of gas inside R
//              that is an integral :  Int_0^R[4 Pi r^2  \rho_{gas}(r)  dr]
//	ff_out - an integral : Int_R^R_max[4 Pi r  \rho_{gas}(r)  dr]
      P = - (G_const * ff_in / R + G_const * ff_out);
      T = - G_const * ff_in / (R*R*R);
if (R>RmaxGasPot){
      P = - (G_const * ff_in / R );
}
      nemo_dprintf (1,"R = %.8E  P = %.8E  ff_in = %.8E   ff_out = %.8E\n",R,P,ff_in,ff_out);
//      nemo_dprintf (1,"R = %.8E  P = %.8E   T = %.8E\n",R,P,T);
    }
  };
}
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<GPDehnen>)
__DEF__POT(SphericalPot<GPDehnen>)

//------------------------------------------------------------------------------
