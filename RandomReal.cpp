#include "header.h"

int main(int argc, char *argv[]){

  double zmin=0.0, zmax=0.0;
  int galaxies=0.0;
  string outputfile;
  bool HALF=false, START=false;
    int arg=1;
 int nside=512;
while (arg < argc) {
    if (argv[arg][0] == '-') {
      switch (argv[arg][1]) {
      case 'O':
 	outputfile=string(argv[arg+1]);
	arg+=2;
	break;
      case 'N':
 	nside=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'G':
	galaxies=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'H':
	HALF=true;
	arg+=1;
	break;
      case 'h':
	cout<<BOLDRED<<"Options:"<<RESET<<endl;
	cout<<endl;
	cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in fits format"<<endl;
	cout<<YELLOW<<"-N"<<RESET<<"side: "<<BLUE<<"'nside'"<<RESET<<" Healpix resolution of the output map; default 512"<<endl;	
	cout<<YELLOW<<"-G"<<RESET<<"alaxies: "<<BLUE<<"'galaxy number'"<<RESET<<" of the uniform distributed map"<<endl;
	cout<<YELLOW<<"-H"<<RESET<<"alf: upper half of the sphere will be realised, otherwise full sphere"<<endl;
	cout<<MAGENTA<<"Example:"<<endl<<"./RandomReal  -O out_map.fits -G 100000  -H"<<endl<<";"<<RESET<<endl<<endl;
	cout<<GREEN<<"Queries should be directed to wolz.laura@gmail.com"<<RESET<<endl;
	
	START=false;
	arg++;
	break;
      default:
	cout<<"bad argument: do RandomReal -h to get man page. "<<RED<<"Exit now"<<RESET<<endl;
	exit(-1);
	break;
      }
    }else{
      cout<<"bad argument: do RandomReal -h to get man page. "<<RED<<"Exit now."<<RESET<<endl;
      exit(-1);
    }
  }
 int order=log2(nside);
  Healpix_Map<double> randomreal(order, RING);
 
//Initialise random numbers
  gsl_rng *random=gsl_rng_alloc(gsl_rng_default);
  gsl_rng *random2=gsl_rng_alloc(gsl_rng_default);
  gsl_rng *random3=gsl_rng_alloc(gsl_rng_default);
  int seed=time(NULL);
  int seed2=time(NULL)+100;
  int seed3=time(NULL)+1000;
  gsl_rng_set(random, seed);
  gsl_rng_set(random2, seed2);
  gsl_rng_set(random3, seed3);
 
  
    randomreal.fill(0.0);
    pointing pointer;
    double phi, costheta;
    int pixel=0;
    for(int i=0; i<galaxies;i++){
      costheta=gsl_rng_uniform(random);
      phi=2.0*pi*gsl_rng_uniform(random2);
      if(!HALF) if(gsl_rng_uniform(random3)>0.5) costheta=-costheta;
      pointer.theta=acos(costheta);
      pointer.phi=phi;
      pixel=randomreal.ang2pix(pointer);
      randomreal[pixel]+=1.0; 
    }
  
  write_Healpix_map_to_fits<float64>(outputfile, randomreal, PLANCK_FLOAT64);

  gsl_rng_free(random);
  gsl_rng_free(random2);



  return 0;
}


