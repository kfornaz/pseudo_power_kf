#include "header.h"

int main(int argc, char *argv[]){


  bool CROSS=false, OVERDENSITY=false, SHOTOVERDENSITY=false, START=false, MULTIPLE=false, PEEBLES=false, MASK=false, CALCILMJLM=false, READILMJLM=false, GALAXYCOUNTS=false;
  std::string map2file, inputfile,Alm2file, outputfile, maskfile, IlmJlmfile, mapfile;
  std::string inputfileroot, outputfileroot;

  int arg=1;
  int lmax=-100;
  int Nmaps=1;
  int nside=512;
  
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
      case 'P':
	PEEBLES=true;
	arg++;
	break;
      case 'o':
	OVERDENSITY=true;
	arg++;
	break;
      case 's':
	SHOTOVERDENSITY=true;
	arg++;
	break;
      case 'm':
	MASK=true;
	maskfile=std::string(argv[arg+1]);
	std::cout<<maskfile<<std::endl;
	arg+=2;
	break;
      case 'T':
	mapfile=std::string(argv[arg+1]);
	arg+=2;
	std::cout<<mapfile<<std::endl;
	break;
      case 'G':
	GALAXYCOUNTS=true;;
	arg++;
	break;
      case 'R':
	READILMJLM=true;
	IlmJlmfile=std::string(argv[arg+1]);
	std::cout<<IlmJlmfile<<std::endl;
	arg+=2;
	break;
      case 'C':
	CALCILMJLM=true;
	IlmJlmfile=std::string(argv[arg+1]);
	arg+=2;
	break;
       case 'c':
	CROSS=true;
	Alm2file=std::string(argv[arg+1]);	
	map2file=std::string(argv[arg+2]);
	std::cout<<"calculate cross corr!"<<std::endl;
	arg+=3;
	break;
      case 'N':
	nside=atoi(argv[arg+1]);
	arg+=2;
	break;
       case 'L':
	lmax=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'M':
	MULTIPLE=true;
	Nmaps=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'h':


	std::cout<<BOLDRED<<"Options:"<<RESET<<std::endl;
	std::cout<<std::endl;
	std::cout<<YELLOW<<"-I"<<RESET<<"nput: name of the "<<BLUE<<"'inputfile'"<<RESET<<" in fits"<<std::endl;
	std::cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in ascii format"<<std::endl;
	std::cout<<YELLOW<<"-P"<<RESET<<"eebles: Peebles approximation of Cl is calculated; -R or -C and -m as additional option needed!"<<std::endl;
	std::cout<<YELLOW<<"-m"<<RESET<<"ask: "<<BLUE<<"'inputfile'"<<RESET<<" of mask in fits format"<<std::endl;
	std::cout<<YELLOW<<"-o"<<RESET<<"verdensity:  overdensity input Alm and Map for Peebles"<<std::endl;
	std::cout<<YELLOW<<"-G"<<RESET<<"alaxycounts: if -G the shotnoise of galaxy counts is removed. So input map should be density of counts."<<std::endl;
	std::cout<<YELLOW<<"-R"<<RESET<<"eadIlmJlm: "<<BLUE<<"'inputfile'"<<RESET<<" of Ilm and Jlm in ascii format"<<std::endl;
	std::cout<<YELLOW<<"-C"<<RESET<<"alcIlmJlm: Ilm and Jlm are calculated and written to "<<BLUE<<"'outputfile'"<<RESET<<";"<<std::endl;
	std::cout<<YELLOW<<"-T"<<RESET<<"emperaturmap: "<<BLUE<<"'mapfile'"<<RESET<<" of map in fits format; only for Peebles approximation needed"<<std::endl;
	std::cout<<YELLOW<<"-s"<<RESET<<"hot overdensity: "<<BLUE<<"'mapfile'"<<RESET<<" of map in fits format; with shot noise subtraction"<<std::endl;
	std::cout<<YELLOW<<"-N"<<RESET<<"side:"<<BLUE<<"'nside'"<<RESET<<" of the mask and Ilm and Jlm if Peebles approximation is computed:"<<std::endl;
	std::cout<<YELLOW<<"-L"<<RESET<<"max:"<<BLUE<<"'lmax'"<<RESET<<" of the output Cl if Peebles approximation is computed:"<<std::endl;
	std::cout<<YELLOW<<"-M"<<RESET<<"ultiple: "<<BLUE<<"'nmaps'"<<RESET<<" Number of files; Input file must be in format 'Inputfileroot' so that the whole name is like 'Inputfileroot001.dat'"<<std::endl;
	std::cout<<YELLOW<<"-c"<<RESET<<"rosscorrelation: "<<BLUE<<"'almfile2' 'mapfile2'"<<RESET<<" ;files containing the second alms in fits format and the second map in Healpix map fits format; This option only works for Peebles approx; if Overdensity option chosen, please insert dummy map. "<<std::endl<<std::endl;
	std::cout<<MAGENTA<<"Example:"<<std::endl<<"./Alm2Cl -I input_alm.fits -O output_cl.dat -m mask.fits -C write_Ilm_Jlm_to_this_file.dat -P "<<std::endl<<" computes the Peebles approximation of the power spectrum of the input alms 'input_alm.fits' with mask 'mask.fits' and writes it to output file 'output_cl.dat';"<<RESET<<std::endl<<std::endl;
	START=false;
	arg++;
	break;
      default:
	std::cout<<"bad argument: do Map2Cl -h to get man page. "<<RED<<"Exit now"<<RESET<<std::endl;
	exit(-1);
	break;
      }
    }else{
      std::cout<<"bad argument: do Map2Cl -h to get man page. "<<RED<<"Exit now."<<RESET<<std::endl;
      exit(-1);
    }
  }
  if(!inputfile.empty() && !outputfile.empty()) START=true;

  //constructor for class where analysis is computed
  spherical_shell Shell;
  std::cout<<"Peebles: "<<PEEBLES<<"   Mask:"<<MASK<<" "<<maskfile<<std::endl;
  if(START){
    if(MULTIPLE){
      inputfileroot=inputfile;
      outputfileroot=outputfile;
    }
    
    //loop for multiple files
    for(int i=0; i<Nmaps; i++){

      if(MULTIPLE){
	inputfile=inputfileroot+intToString(i,3)+".dat";
	outputfile=outputfileroot+intToString(i,3)+".fits";
      }

      //for first file calc/read in of: mask, Ilm, Jlm ans first map
     if(i==0){
		Shell.set( inputfile, nside, lmax);
		if(MASK){
		  std::cout<<"read mask now!"<<std::endl;
		  Shell.read_mask(maskfile);
		}
		if(READILMJLM){
		  std::cout<<"read Ilm Jlm now!"<<std::endl;
		  Shell.read_I_lm_J_lm(IlmJlmfile);
		}
		if(CALCILMJLM) Shell.calc_I_lm_J_lm(IlmJlmfile);
    }
     Shell.read_alm(inputfile);     
    if(PEEBLES){
		if(!MASK && (!READILMJLM || !CALCILMJLM)){
		  std::cout<<"Either mask or IlmJlm are not defined! Cannot calculate Peebles approximation, see  Alm2Cl -h to get man page. "<<RED<<"Exit now."<<RESET<<std::endl;
		  exit(-1);
		}
		std::cout<<"Calc PEEbles now!"<<std::endl;
		//std::cout<<mapfile<<"\n"<<std::flush;
		if(CROSS) Shell.calc_PeeblesCross(OVERDENSITY, outputfile, mapfile, GALAXYCOUNTS, Alm2file, map2file);
		else	Shell.calc_PeeblesCl(OVERDENSITY, outputfile, mapfile, GALAXYCOUNTS, SHOTOVERDENSITY);
    }else{
		if(CROSS) Shell.calc_cross(MASK, outputfile,Alm2file, map2file);
		else Shell.calc_cl(MASK, outputfile);
    }
      
    }
  }else std::cout<<"Input or Output file name empty; nothing to do;"<<std::endl;
    
 return 0;
}
