#include "header.h"

int main(int argc, char *argv[]){

bool OVERDENSITY=false, START=false,  MASK=false;
  std::string  outputfile, maskfile, IlmJlmfile;

  int arg=1;

  int lmin=0;
  int lmax=-100;
  int nside=512;
 while (arg < argc) {
    if (argv[arg][0] == '-') {
      switch (argv[arg][1]) {
      case 'l':
	lmin=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'L':
	lmax=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'O':
 	IlmJlmfile=std::string(argv[arg+1]);
	arg+=2;
	break;
      case 'm':
	MASK=true;
	maskfile=std::string(argv[arg+1]);
	std::cout<<maskfile<<std::endl;
	arg+=2;
	break;
      case 'N':
	nside=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'h':
	std::cout<<RED<<"This function can have long running times for high l and high Nside. It is advised to run it for small l ranges in parallel and copy all output files together using the 'cat' command."<<std::endl;
	std::cout<<BOLDRED<<"Options:"<<RESET<<std::endl;
	std::cout<<std::endl;

	std::cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in ascii format"<<std::endl;

	std::cout<<YELLOW<<"-m"<<RESET<<"ask: "<<BLUE<<"'inputfile'"<<RESET<<" of mask in fits format"<<std::endl;

	std::cout<<YELLOW<<"-N"<<RESET<<"side:"<<BLUE<<"'nside'"<<RESET<<" of the mask and map"<<std::endl;
	std::cout<<YELLOW<<"-l"<<RESET<<"min: "<<BLUE<<"'lmin'; "<<RESET<<"start with 0 for complete runs;"<<std::endl;
	std::cout<<YELLOW<<"-L"<<RESET<<"max: "<<BLUE<<"'lmax'"<<RESET<<std::endl;
	std::cout<<MAGENTA<<"Example:"<<std::endl<<".IlmJlm  -O ilmjlm_mask1.dat -m mask1.fits -N 128 -l 0 -L 100"<<std::endl<<"; Calculates the Ilm and Jlm for mask1 with resolution Nside=128 for the multipoles 0 to 100"<<RESET<<std::endl<<std::endl;

	START=false;
	arg++;
	break;
      default:
	std::cout<<"bad argument: do IlmJlm -h to get man page. "<<RED<<"Exit now"<<RESET<<std::endl;
	exit(-1);
	break;
      }
    }else{
      std::cout<<"bad argument: do IlmJlm -h to get man page. "<<RED<<"Exit now."<<RESET<<std::endl;
      exit(-1);
    }
  }
  if( !IlmJlmfile.empty()) START=true;

  //constructor for class where analysis is computed
  spherical_shell Shell;
  if(START){
        
    Shell.set(nside,lmax);   
    Shell.read_mask(maskfile);
    std::cout<<"here"<<std::endl;
    Shell.calc_I_lm_J_lm_weighted(IlmJlmfile,lmin,lmax);
 
 
  }else std::cout<<" Output file name empty; nothing to do;"<<std::endl;
    
 return 0;
}
