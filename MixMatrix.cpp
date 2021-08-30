#include "header.h"

int main(int argc, char *argv[]){


  bool  START=false, MASK=false, CALCILMJLM=false, READILMJLM=false,CALCONECOLUMONLY=false;
  std::string  outputfile, maskfile,IlmJlmfile;
  std::string  outputfileroot;

  int arg=1;
  int filenumber=0;
  int lmax=-100;
  int lzero=-100;
  int nside=512;
  int Nmaps=1;
  double fwhm=0.0;

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
      case 'R':
	READILMJLM=true;
	IlmJlmfile=std::string(argv[arg+1]);
	arg+=2;
	break;
      case 'l':
  CALCONECOLUMONLY=true;
  lzero=atoi(argv[arg+1]);
  arg+=2;
  break;
      case 'L':
  lmax=atoi(argv[arg+1]);
  arg+=2;
  break;
      case 'C':
	CALCILMJLM=true;
	IlmJlmfile=std::string(argv[arg+1]);
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
	std::cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in ascii format"<<std::endl;
	std::cout<<YELLOW<<"-m"<<RESET<<"ask: "<<BLUE<<"'inputfile'"<<RESET<<" of mask in fits format"<<std::endl;
	std::cout<<YELLOW<<"-N"<<RESET<<"side: "<<BLUE<<"'nside'"<<RESET<<" Healpix resolution of the output map, should be the same as mask or the input Ilm Jlm;"<<std::endl;
	std::cout<<YELLOW<<"-L"<<RESET<<"max: "<<BLUE<<"'lmax'"<<RESET<<" maximum multipole of Cl; default lmax=2*Nside of input map"<<std::endl;
	std::cout<<YELLOW<<"-R"<<RESET<<"eadIlmJlm: "<<BLUE<<"'inputfile'"<<RESET<<" of Ilm and Jlm in ascii format"<<std::endl;
	std::cout<<YELLOW<<"-C"<<RESET<<"alcIlmJlm: Ilm and Jlm are calculated and written to "<<BLUE<<"'outputfile'"<<RESET<<";"<<std::endl;

	std::cout<<MAGENTA<<"Example:"<<std::endl<<"./MixMatrix  -O output_mixmatrix.dat -N 256 -m mask.fits -C outputIlmJlm.dat "<<std::endl<<" computes a the Mixing matix for the mask in mask.fits and writes it to output file 'output_mixmatrix.dat'; computes the Ilm and Jlm and writes them to outputIlmJlm.dat ;"<<RESET<<std::endl<<std::endl;
	START=false;
	arg++;
	break;
      default:
	std::cout<<"bad argument: do NoiseMap -h to get man page. "<<RED<<"Exit now"<<RESET<<std::endl;
	exit(-1);
	break;
      }
    }else{
      std::cout<<"bad argument: do NoiseMap -h to get man page. "<<RED<<"Exit now."<<RESET<<std::endl;
      exit(-1);
    }
  }
  if( !outputfile.empty()) START=true;
    spherical_shell Shell;
   
  if(START){
    if(!READILMJLM && !CALCILMJLM){
      std::cout<<"Neither ReadIlmJlm nor CalcIlmJlm defined."<<RED<<" Exit now."<<RESET<<std::endl;
      exit(-1);
    }
    Shell.set(nside, lmax);
    if(MASK) Shell.read_mask(maskfile);
    if(READILMJLM) Shell.read_I_lm_J_lm(IlmJlmfile);
    if(CALCILMJLM){
      if(!MASK){
	std::cout<<"For option CalcIlmJlm a mask is needed!"<<RED<<"Exit now."<<RESET<<std::endl;
	exit(-1);
      }
      Shell.calc_I_lm_J_lm_weighted(IlmJlmfile);
    }
    //int wigner=gsl_sf_coupling_3j(464, 1018, 650, 0, 0, 0);
    //std::cout << "Wigner: " << wigner << "\n" << std::endl;
    if (CALCONECOLUMONLY) Shell.calc_mix_matrix(outputfile,lzero);
    else Shell.calc_mix_matrix(outputfile);

  }else std::cout<<"Input or Output file name empty; nothing to do;"<<std::endl;
  return 0;
}
