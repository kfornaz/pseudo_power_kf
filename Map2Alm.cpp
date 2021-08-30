#include "header.h"

int main(int argc, char *argv[]){


  bool OVERDENSITY=false, START=false, MASK=false;
  std::string inputfile, outputfile, maskfile;
  std::string inputfileroot, outputfileroot;

  int arg=1;
  int lmax=-100;
  int Nmaps=1;

 while (arg < argc) {
    if (argv[arg][0] == '-') {
      switch (argv[arg][1]) {
      case 'I':
	inputfile=std::string(argv[arg+1]);
	arg+=2;
	break;
      case 'O':
 	outputfile=std::string(argv[arg+1]);
	arg+=2;
	break;
      case 'L':
 	lmax=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'm':
	MASK=true;
	maskfile=std::string(argv[arg+1]);
	arg+=2;
	break;
      case 'h':
	std::cout<<BOLDRED<<"Options:"<<RESET<<std::endl;
	std::cout<<std::endl;
	std::cout<<YELLOW<<"-I"<<RESET<<"nput: name of the "<<BLUE<<"'inputfile'"<<RESET<<" in fits"<<std::endl;
	std::cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in Healpix fits format"<<std::endl;
	std::cout<<YELLOW<<"-L"<<RESET<<"max: "<<BLUE<<"'lmax'"<<RESET<<" maximum multipole of Cl; default lmax=2*Nside of input map"<<std::endl;
	std::cout<<YELLOW<<"-m"<<RESET<<"ask: "<<BLUE<<"'inputfile'"<<RESET<<" of mask in fits format; cuts and weights survey window according to mask; "<<std::endl;
	std::cout<<MAGENTA<<"Example:"<<std::endl<<"./Map2Alm -I input_map.fits -O output_alm.fits -L 100 "<<std::endl<<" computes the Alms of the input map 'input_map.fits'to lmax=100 and writes it to output file 'output_alm.fits';"<<RESET<<std::endl<<std::endl;
	std::cout<<GREEN<<"Queries should be directed to wolz.laura@gmail.com"<<RESET<<std::endl;

	START=false;
	arg++;
	break;
      default:
	std::cout<<"bad argument: do Map2Alm -h to get man page. "<<RED<<"Exit now"<<RESET<<std::endl;
	exit(-1);
	break;
      }
    }else{
      std::cout<<"bad argument: do Map2Alm -h to get man page. "<<RED<<"Exit now."<<RESET<<std::endl;
      exit(-1);
    }
  }
  if(!inputfile.empty() && !outputfile.empty()) START=true;

  //constructor for class where analysis is computed
  spherical_shell Shell;

  if(START){
 
      //for first file set parameter and set first map
	Shell.set(lmax, inputfile);
	if(MASK) Shell.read_mask(maskfile);
      
      //calculation of alm
      Shell.calc_alm(MASK);
      Shell.write_alm(outputfile);
    
  }else std::cout<<"Input or Output file name empty; nothing to do;"<<std::endl;
    
 return 0;
}
