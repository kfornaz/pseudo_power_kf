#include "header.h"
int main(int argc, char *argv[]){


  bool NOISE=false, SUB=false, START=false,CAT=false, ASCII=false, INVERSE=false, MASK=false;
  int arg=1;

  std::string inputfile1, inputfile2, outputfile, maskfile;
  while (arg < argc) {
    if (argv[arg][0] == '-') {
      switch (argv[arg][1]) {
      case 'I':
	inputfile1=std::string(argv[arg+1]);
	inputfile2=std::string(argv[arg+2]);
	arg+=3;
	break;
	case 'S':
	SUB=true;
	arg+=1;
	break;
	case 'i':
	INVERSE=true;
	arg+=1;
	break;
	case 'A':
	ASCII=true;
	arg+=1;
	break;
	case 'C':
	CAT=true;
	arg+=1;
	break;
      case 'O':
 	outputfile=std::string(argv[arg+1]);
	arg+=2;
	break;
      case 'M':
	MASK=true;
	maskfile=std::string(argv[arg+1]);
	arg+=2;
	break;
      case 'h':
	std::cout<<RED<<"Works only for maps with the same resolution Nside!"<<std::endl;
	std::cout<<BOLDRED<<"Options:"<<RESET<<std::endl<<std::endl;
	std::cout<<YELLOW<<"-I"<<RESET<<"nput: name of the "<<BLUE<<"'input map1'"<<RESET<<" name of the "<<BLUE<<"'input map2'"<<RESET<<" in Healpix fits"<<std::endl;
	std::cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'output file'"<<RESET<<" in fits format"<<std::endl;
	std::cout<<YELLOW<<"-S"<<RESET<<"ubtraction: returns map1-map2"<<std::endl;
	std::cout<<YELLOW<<"-C"<<RESET<<"atalogue: returns map1 in ascii format, nothing else is done. Requires dummy names as 'input map2'"<<std::endl;
	std::cout<<YELLOW<<"-I"<<RESET<<"nverse: returns -map1. Requires dummy names as 'input map2"<<std::endl;
	std::cout<<YELLOW<<"-A"<<RESET<<"scii: output is an galaxy catalogue rather then a Healpix fits map."<<std::endl;
	std::cout<<YELLOW<<"-M"<<RESET<<"ask: name of the "<<BLUE<<"'mask file'"<<RESET<<" in Healpix fits format if only certain region should be output"<<std::endl;
	std::cout<<MAGENTA<<"Example:"<<std::endl<<"./AddMAps -I map1.fits map2.fits -O output_map.dat "<<std::endl<<" Adds map2 to map1 and write it into output_map.fits"<<RESET<<std::endl<<std::endl;
	std::cout<<MAGENTA<<"Example:"<<std::endl<<"./AddMAps -I map1.fits dummy -O output_cat.dat -C "<<std::endl<<" Outputs the pixel values of map1 as an ascii catalogue with format RA DEC VALUE"<<RESET<<std::endl<<std::endl;
	std::cout<<MAGENTA<<"Example:"<<std::endl<<"./AddMAps -I map1.fits dummy -O output_cat.dat -C -M mask.fits"<<std::endl<<" Outputs the pixel values defined in the mask (mask=1) of map1 as an ascii catalogue with format RA DEC VALUE"<<RESET<<std::endl<<std::endl;
	std::cout<<GREEN<<"Queries should be directed to wolz.laura@gmail.com"<<RESET<<std::endl;

	START=false;
	arg++;
	break;
 default:
	std::cout<<"bad argument: do AddMaps -h to get man page"<<std::endl;
	exit(-1);
	break;
      }
    }else{
      std::cout<<"bad argument: do AddMaps -h to get man page"<<std::endl;
      exit(-1);
    }
  }
  if(!inputfile1.empty()  && !outputfile.empty()) START=true;
  Healpix_Map<double> map1;
  Healpix_Map<double> map2;
  Healpix_Map<double> mask;

  if(START){
    read_Healpix_map_from_fits(inputfile1,map1, 1, 2);
    std::cout<<!CAT<<"  "<<!INVERSE<<std::endl;
    if(!CAT && !INVERSE){
	std::cout<<"Read second map"<<std::endl;
	 read_Healpix_map_from_fits(inputfile2,map2, 1, 2);
    }else {
	map2.SetNside(map1.Nside(), RING); 
	map2.fill(0.0);
	}
    if(MASK) read_Healpix_map_from_fits(maskfile,mask, 1, 2);
    else{
      mask.SetNside(map1.Nside(), RING);
      mask.fill(1.0);
    }
    if(map1.Nside()!=map2.Nside()){
      std::cout<<"Maps don't have the same resolution! "<<RED<<"Exit now. "<<RESET<<std::endl;
      exit(-1);
    }
if(NOISE) std::cout<<map1.rms()<<"  "<<map2.rms()<<std::endl;

    for(int i=0; i<map1.Npix(); i++){
      if(mask[i]==1){
	if(SUB){
	  map1[i]-=map2[i];
	}else if(INVERSE){
	  map1[i]=-map1[i];
	}else{
	  map1[i]+=map2[i];
	}
      }else map1[i]=0.0;
    }

    if(ASCII){
	pointing point;
      std::ofstream file(&outputfile[0]);
      for(int i=0; i<map1.Npix(); i++){
	point=map1.pix2ang(i);
	  if(mask[i]==1) file<<point.phi<<"  "<<point.theta<<"  "<<map1[i]<<std::endl;
      }
      file.close();
	
    }else  write_Healpix_map_to_fits(outputfile, map1, PLANCK_FLOAT64);
  }else std::cout<<"Input or Output file name empty; nothing to do;"<<std::endl;



return 0;

}
