#ifndef HEADER_H
#define HEADER_H

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <fitsio.h>



static const std::string winpix_directory="/mnt/d/karin/PseudoPower/Healpix_3.11/data";

#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <alm.h>
#include <alm_fitsio.h>
#include <alm_powspec_tools.h>
#include <alm_healpix_tools.h>
#include <string_utils.h>
#include <datatypes.h>
#include <error_handling.h>
#include <fitshandle.h>
#include <lsconstants.h>

#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#include <healpix_data_io.h>
#include <healpix_base.h>
//#include <healpix_base2.h>

#include <healpix_data_io.h>

#include <powspec.h>
#include <powspec_fitsio.h>
#include <planck_rng.h>
//#include <wigxjpf.h>
#include <UCLWig3j.hpp>

#define PI 3.14159265
#define light 299792458
#define boltzmann 1.3806488e−23 //in J/K
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
                                                                                                                                                                                                                               #define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

extern const double Tsys, Aeff, integration, FOV, fracsense, HIline, HIfrequency;
extern const int telescopes;
extern std::string survey;


class spherical_shell{
public:
  int order, nside;
  int lmax, mmax;
  double theta_lim, phi_lim;
  Healpix_Ordering_Scheme scheme;
  int mask_area;
  double gal_number;
  double mean_gal_number;
  spherical_shell();
  spherical_shell(int l, std::string filename);
  spherical_shell(int ns);
  void set(int l, std::string inputfile);
  void set( std::string inputfile);
  void set( std::string inputfile, int nside);
  void set( std::string inputfile, int nside, int l);
  void set(int ns, int l);
  ~spherical_shell();
  Healpix_Map<double> hp_map;
  Healpix_Map<double> hp_mask;
  xcomplex<double> **I;
  double **J;
  PowSpec powspec;
  PowSpec Peebles; 
  PowSpec mix_Cl; 
  PowSpec unmix_Cl;
  Alm<xcomplex <double> > alm;
  double *partialCl;
  int *pixelpos;
  double *flux;
  double get_npix();
  void pixelised();
  void masked();
  void reset_l(int l);
  void fluxed(int b);
  void write_hp_map(std::string filename);
  void set_hp_map();
  void draw_alm(PowSpec input,  int ll);
  void draw_alm2(PowSpec input,  int ll);
  void create_map( int ll);
  void create_map();
  void calc_cl(bool MASK, std::string outputfile);
  void calc_cross(bool MASK, std::string outputfile, std::string Alm2file, std::string map2file);
  void calc_PeeblesCross(bool OVERDENSITY, std::string outputfile, std::string mapfile, bool GALAXYCOUNTS, std::string Alm2file, std::string map2file);
  void write_alm();
  void read_alm(std::string inputfile);
  void write_map(std::string outputfile);
  void calc_I_lm_J_lm(std::string IlmJlmfile );
  void calc_I_lm_J_lm_weighted(std::string IlmJlmfile ); //PCM
  void calc_I_lm_J_lm(std::string IlmJlmfile,int startl, int stopl);
  void calc_I_lm_J_lm_weighted(std::string IlmJlmfile,int startl, int stopl); //PCM
  void read_I_lm_J_lm(std::string IlmJlmfile );
  void read_map(std::string inputfile);
  void calc_mask_lognorm(std::string outputfile);
  void calc_mask_halfsky(std::string outputfile);
  void calc_mask_fg(std::string outputfile, std::string fgfile, double tlim);

  void calc_mask_accurate();
  void read_mask(std::string maskfile);

  void calc_mask(std::string file,double *theta, double *phi);
  //void calc_mask();
  void calc_gal_number();
  Healpix_Map<double> antipixelise(std::string file, Healpix_Map<double> newmap);
  void antipixelise();
  void winpix(int whichcl);

  void calc_PeeblesCl(bool OVERDENSITY, std::string outputfile, std::string mapfile, bool GALAXYCOUNTS,bool SHOTOVERDENSITY);
  void calc_mix_Cl(PowSpec fullCl);
 void read_file(std::string filename);
  void calc_alm(bool MASK);
  void write_alm(std::string outputfile);
  void calc_mix_matrix(std::string outputfile);
  void calc_mix_matrix(std::string outputfile,int lzero);
private:
  void calc_alm_galcat();
  void set_order();
  void set_scheme();

};
void  noise_map(int numb, double fwhm, std::string outputfile, int nside);
void  noise_map(double *zlim, double fwhm, std::string outputfile, int nside);
void  noise_map_SKA(double *flim, double omega, std::string outputfile, int nside);
void  noise_map_rms(double rms, double omega, std::string outputfile, int nside);

//double chi2(double** powerspectrum, double ***covariance, cosmomodel *dee,int lines);
//gsl_matrix determine_inverse_covmatrix(gsl_matrix *covariancematrix, double **powerspectrum, int multipol);


//Outputcolors



#endif
