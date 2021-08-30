#include "header.h"

int main(int argc, char *argv[]){


  bool HALFSKY=false, LOGNORM=false, CONSTLIM=false, FOREGROUND=false,START=false;
  std::string inputfile, outputfile, maskfile, fgfile;
  std::string inputfileroot, outputfileroot;

  int arg=1;
  int nside=512;
  int Nmaps=1;
  double *theta, *phi, tlim=0.0;
  theta=new double[2];
  phi=new double[2];

 while (arg < argc) {
    if (argv[arg][0] == '-') {
      switch (argv[arg][1]) {
      case 'O':
 	outputfile=std::string(argv[arg+1]);
	arg+=2;
	break;
      case 'N':
 	nside=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'H':
	HALFSKY=true;
	arg++;
	break;
      case 'L':
	LOGNORM=true;
	arg++;
	break;
      case 'F':
	FOREGROUND=true;
	tlim=atof(argv[arg+1]);
	fgfile=std::string(argv[arg+2]);
	arg+=3;
	break;
      case 'C':
	CONSTLIM=true;
	theta[0]=atof(argv[arg+1])*pi/180.0;
	theta[1]=atof(argv[arg+2])*pi/180.0;
	phi[0]=atof(argv[arg+3])*pi/180.0;
	phi[1]=atof(argv[arg+4])*pi/180.0;
	arg+=5;
	break;
      case 'h':
	std::cout<<BOLDRED<<"Options:"<<RESET<<std::endl;
	std::cout<<std::endl;
	std::cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in fits format"<<std::endl;
	std::cout<<YELLOW<<"-N"<<RESET<<"side: "<<BLUE<<"'nside'"<<RESET<<" Healpix resolution of the output map; default 512"<<std::endl;
	std::cout<<YELLOW<<"-H"<<RESET<<"alfsky: A mask for the northern halfsk with theta<90 is produced"<<std::endl;
	std::cout<<YELLOW<<"-L"<<RESET<<"ognorm: A mask with 20<theta<95 and 110<phi<270 is produced"<<std::endl;
	std::cout<<YELLOW<<"-F"<<RESET<<"oreground: A mask is cut off according to the input"<<BLUE<<" 'Temperaturelimit'"<<RESET<<" out of the input "<<BLUE<<" foreground map"<<RESET<<" in fits format"<<std::endl;
	std::cout<<YELLOW<<"-C"<<RESET<<"onstlim: A mask with "<<BLUE<<"'thetamin'"<<RESET<<"<theta<"<<BLUE<<"'thetamax'"<<RESET<<" and "<<BLUE<<"'phimin'"<<RESET<<"<phi<"<<BLUE<<"'phimax'"<<RESET<<" is produced; with 0<theta<180 and 0<phi<360"<<std::endl;
	std::cout<<MAGENTA<<"Example:"<<std::endl<<"./CreateMask  -O output_mask.fits -N 256 -C 10 30 0 90 "<<std::endl
	<<" computes a mask with nside 256 with area 10<theta<30 and 0<phi<90 writes it to output file 'output_alm.fits';"<<RESET<<std::endl<<std::endl;
	std::cout<<GREEN<<"Queries should be directed to wolz.laura@gmail.com"<<RESET<<std::endl;
	START=false;
	arg++;
	break;
      default:
	std::cout<<"bad argument: do CreateMask -h to get man page. "<<RED<<"Exit now"<<RESET<<std::endl;
	exit(-1);
	break;
      }
    }else{
      std::cout<<"bad argument: do CreateMask -h to get man page. "<<RED<<"Exit now."<<RESET<<std::endl;
      exit(-1);
    }
  }
  if( !outputfile.empty()) START=true;

  //constructor for class where analysis is computed
  spherical_shell Shell;

  if(START){
    Shell.order=log2(nside);
    Shell.hp_mask.Set(Shell.order, RING);
      
    if(LOGNORM) Shell.calc_mask_lognorm(outputfile);
    if(HALFSKY) Shell.calc_mask_halfsky(outputfile);
    if(CONSTLIM) Shell.calc_mask(outputfile, theta, phi);
    if(FOREGROUND) Shell.calc_mask_fg(outputfile,fgfile, tlim);
  }else std::cout<<" Output file name empty; nothing to do;"<<std::endl;
    
 return 0;
}
