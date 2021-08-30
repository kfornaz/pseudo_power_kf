#include "header.h"

int main(int argc, char *argv[]){


  bool OVERDENSITY=false, START=false, MULTIPLE=false, MASK=false;
  std::string inputfile, outputfile;
  std::string inputfileroot, outputfileroot;

  int arg=1;
  int lmax=-100;
  int Nmaps=1;
  double psi=0.0, phi=0.0, theta=0.0;

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
      case 'A':
	psi=atof(argv[arg+1]);
	theta=atof(argv[arg+2]);
	phi=atof(argv[arg+3]);
	arg+=4;
	break;
      case 'h':
	std::cout<<BOLDRED<<"Options:"<<RESET<<std::endl;
	std::cout<<std::endl;
	std::cout<<YELLOW<<"-I"<<RESET<<"nput: name of the "<<BLUE<<"'inputfile'"<<RESET<<" in fits"<<std::endl;
	std::cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in fits format"<<std::endl;
	std::cout<<YELLOW<<"-A"<<RESET<<"ngles: "<<BLUE<<"'psi' 'theta' 'phi'"<<RESET<<" the Euler angles of the rotation in degrees"<<std::endl;
	std::cout<<MAGENTA<<"Example:"<<std::endl<<"./RotateAlm -I input_alm.fits -O output_cl.dat -A 45.0 45.0 0.0 "<<std::endl<<" computes the Alm rotated around and writes it to output file 'output_alm.fits';"<<RESET<<std::endl<<std::endl;
	std::cout<<GREEN<<"Queries should be directed to wolz.laura@gmail.com"<<RESET<<std::endl;
	START=false;
	arg++;
	break;
	
      default:
	std::cout<<"bad argument: do RotateAlm -h to get man page. "<<RED<<"Exit now"<<RESET<<std::endl;
	exit(-1);
	break;
      }
    }else{
      std::cout<<"bad argument: do RotateAlm -h to get man page. "<<RED<<"Exit now."<<RESET<<std::endl;
      exit(-1);
    }
  }
  if(!inputfile.empty() && !outputfile.empty()) START=true;

  //constructor for class where analysis is computed
  spherical_shell Shell;

  if(START){
      Shell.set(inputfile);
      rotate_alm(Shell.alm, psi*pi/180.0, theta*pi/180.0, phi*pi/180.0);
      Shell.write_alm(outputfile);
     
  }else std::cout<<"Input or Output file name empty; nothing to do;"<<std::endl;
    
 return 0;
}
