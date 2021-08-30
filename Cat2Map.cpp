#include "header.h"

#define _GLIBCXX_USE_CXX11_ABI 0

void set_count_map(std::string inputfile, std::string outputfile, int nside, bool OVERDENSITY);

//Options:
//-Input: name of the input file in ascii
//-Output: name of the output file in fits format
//-Density
//-Nside : Nside of the Healpix maps, default=512
int main(int argc, char *argv[]){


  bool OVERDENSITY=false,  START=false;
  std::string inputfile, outputfile;
  std::string inputfileroot, outputfileroot;

  int arg=1;
  int nside=512;
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
      case 'N':
 	nside=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'D':
	OVERDENSITY=true;
	arg++;
	break;
      case 'h':
	std::cout<<BOLDRED<<"Options:"<<RESET<<std::endl<<std::endl;
	std::cout<<YELLOW<<"-I"<<RESET<<"nput: name of the "<<BLUE<<"'input file'"<<RESET<<" in ascii with RA DEC"<<std::endl;
	std::cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'output file'"<<RESET<<" in fits format"<<std::endl;
	std::cout<<YELLOW<<"-D"<<RESET<<"ensity for an overdensity map as output"<<std::endl;
	std::cout<<YELLOW<<"-N"<<RESET<<"side: "<<BLUE<<"'Nside'"<<RESET<<" of the Healpix maps, default=512"<<std::endl;
	std::cout<<MAGENTA<<"Example:"<<std::endl<<" './Cat2Map -I input_catalogue.dat -O output_map.fits -N 1024 '"<<std::endl
	<<" makes a galaxy density map with nside 1024 from the catalogue 'input_catalogue.dat' and writes it to 'output_map.fits';"<<RESET<<std::endl<<std::endl;

	std::cout<<GREEN<<"Queries should be directed to wolz.laura@gmail.com"<<RESET<<std::endl;
	START=false;
	arg++;
	break;
      default:
	std::cout<<"bad argument: do Map2Cat -h to get man page"<<std::endl;
	exit(-1);
	break;
      }
    }else{
      std::cout<<"bad argument: do Map2Cat -h to get man page"<<std::endl;
      exit(-1);
    }
  }
  if(!inputfile.empty() && !outputfile.empty()) START=true;
  
  if(START){  
      
	set_count_map(inputfile, outputfile, nside, OVERDENSITY);
      
    
  }else std::cout<<"Input or Output file name empty; nothing to do;"<<std::endl;
    
  return 0;

}


void set_count_map(std::string inputfile, std::string outputfile, int nside, bool OVERDENSITY){
   pointing point;
   int pixel, count=0;
    Healpix_Map<double> map(log2(nside), RING);
    double tempphi, temptheta, tempcomov, tempmass,temp;
    arr<double> counts_for_map(map.Npix());
    for(int j=0; j<map.Npix(); j++)    counts_for_map[j]=0.0;
    std::ifstream infile;
    char tempstream[1000];
   infile.open(&inputfile[0]);
    if(infile.is_open()){
      infile.getline(tempstream, 1000);
      while(!infile.eof()){
	  infile>>tempphi>>temptheta;
	  point.phi=(tempphi);
	  point.theta=(temptheta);
	
	  pixel=map.ang2pix(point);
	  counts_for_map[pixel]+=1.0;
	
	  count++;
      }
      infile.close();
      std::cout<<"Total galaxy number: "<<count<<std::endl;
    }else{
      std::cout<<"Error with input file! Exit now."<<std::endl;
      exit(-1);
    }
    map.Set(counts_for_map,RING );
    if(OVERDENSITY){
      double average=map.average();
      map.Scale(1.0/average);
      map.Add(-1.0);
    }
    fitshandle opt;
    opt.open(outputfile);
    opt.goto_hdu(0);
    write_Healpix_map_to_fits(opt , map, PLANCK_FLOAT64 );
    opt.close();
}





