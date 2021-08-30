#include "header.h"

spherical_shell::spherical_shell(){}
spherical_shell::~spherical_shell(){}

const double EPS=1.0e-10; 
const double EPSREL=1.0e-9;
const size_t n=10000;
const double LIGHT=3.0e5;


void spherical_shell::set(int l, std::string inputfile){

  read_Healpix_map_from_fits(inputfile, this->hp_map, 1, 2);
  this->scheme=hp_map.Scheme();
  this->nside=hp_map.Nside();
  this->order=round(log2(this->nside));
  this->hp_mask.Set(this->order, this->scheme);
  if(l==-100) this->lmax=2*this->nside+1;
  else this->lmax=l+1;
  this->mmax=this->lmax;
  this->alm.Set(this->lmax, this->mmax);
  this->gal_number=0.0;

  std::cout<<CYAN<<"ORDER: "<<this->order<<"  NSIDE: "<<this->nside<<"  LMAX: "<<this->lmax<<RESET<<std::endl; 

}
void spherical_shell::set(int ns, int l){

  this->scheme=RING;
  this->nside=ns;
  this->order=log2(this->nside);
  this->hp_mask.Set(this->order, this->scheme);
  if(l==-100) this->lmax=2*this->nside;
  else this->lmax=l;
  this->mmax=this->lmax;
  this->alm.Set(this->lmax, this->mmax);
  this->gal_number=0.0;
  this->hp_map.Set(this->order, this->scheme);

  std::cout<<CYAN<<"ORDER: "<<this->order<<"  NSIDE: "<<this->nside<<"  LMAX: "<<this->lmax<<RESET<<std::endl; 

}

void spherical_shell::set(std::string inputfile, int ns){

  this->scheme=RING;
  this->nside=ns;
  this->order=log2((int)this->nside);
  this->hp_map.Set( this->order, this->scheme);
  this->hp_mask.Set(this->order, this->scheme);
  get_almsize(inputfile, this->lmax, this->mmax,2);
  this->alm.Set(this->lmax, this->mmax);
  read_Alm_from_fits(inputfile,this->alm,  this->lmax, this->mmax,2); 
  this->gal_number=0.0;

  std::cout<<CYAN<<"ORDER: "<<this->order<<"  NSIDE: "<<this->nside<<"  LMAX: "<<this->lmax<<RESET<<std::endl; 

}

void spherical_shell::set(std::string inputfile, int ns, int l){

  this->scheme=RING;
  this->nside=ns;
  this->order=log2((int)this->nside);
  this->hp_map.Set( this->order, this->scheme);
  this->hp_mask.Set(this->order, this->scheme);
  get_almsize(inputfile, this->lmax, this->mmax,2);
  this->alm.Set(this->lmax, this->mmax);
  read_Alm_from_fits(inputfile,this->alm,  this->lmax, this->mmax,2); 
  this->gal_number=0.0;
  this->lmax=l;
  this->mmax=this->lmax;
  std::cout<<CYAN<<"ORDER: "<<this->order<<"  NSIDE: "<<this->nside<<"  LMAX: "<<this->lmax<<RESET<<std::endl; 

}





void spherical_shell::set(std::string inputfile){

  get_almsize(inputfile, this->lmax, this->mmax,2);
  this->alm.Set(this->lmax, this->mmax);
  read_Alm_from_fits(inputfile,this->alm,  this->lmax, this->mmax,2);
  this->mmax=this->lmax;
  this->gal_number=0.0;
  std::cout<<CYAN<<"  LMAX: "<<this->lmax<<RESET<<std::endl; 

}

void spherical_shell::read_alm(std::string inputfile){
  read_Alm_from_fits(inputfile,this->alm,  this->lmax, this->mmax,2);

}
void spherical_shell::read_map(std::string inputfile){
  read_Healpix_map_from_fits(inputfile, this->hp_map, 1, 2);
}
void spherical_shell::write_map(std::string outputfile){
  write_Healpix_map_to_fits(outputfile, this->hp_map, PLANCK_FLOAT64);
}

void spherical_shell::write_alm(std::string outputfile){
  write_Alm_to_fits(outputfile, this->alm, this->lmax,this->mmax, PLANCK_FLOAT64);
}

void spherical_shell::read_mask(std::string maskfile){

 this->pixelpos=new int[this->hp_map.Npix()];
 this->mask_area=0;
 Healpix_Map<double> buffer;
 read_Healpix_map_from_fits<float64>(maskfile,buffer, 1, 2 );
 this->hp_mask.Import(buffer,false);
 for(int i=0; i<this->hp_mask.Npix(); i++){
    if(this->hp_mask[i]!=0.0){
      this->pixelpos[this->mask_area]=i;
      this->mask_area++;
    }
  }
 std::cout<<"READ MASK"<<std::endl;
}




void spherical_shell::calc_alm(bool MASK){

  pointing point_here;
  double ylm;
  arr<double>weight(2*hp_map.Nside());
  for(int i=0; i<2*hp_map.Nside(); i++) weight[i]=1.0;
  if(MASK){
    Healpix_Map<double> buffer(this->order, RING);
    buffer.fill(0.0);
    for (int i=0; i<this->hp_map.Npix(); i++)  buffer[i]=this->hp_mask[i]*this->hp_map[i];
    map2alm_iter(buffer, this->alm , 10, weight);
  }else{
    
    map2alm_iter(this->hp_map, this->alm , 10 ,weight);
    /*for (int l=0;l<lmax;l++){
      for (int m=0;m<=l;m++){
     
    
    	for (int k=0; k<this->hp_map.Npix(); k++){
	  point_here=this->hp_map.pix2ang(k);
	  ylm= gsl_sf_legendre_sphPlm(double(l),abs(double(m)),cos(point_here.theta));
     	  this->alm(l,m).re+=(ylm*cos(double(m)*point_here.phi)*this->hp_map[k]);
	  this->alm(l,m).im-=(ylm*sin(double(m)*point_here.phi)*this->hp_map[k]);
    
	}
	this->alm(l,m).re*=(4.0*PI/this->hp_map.Npix());
	this->alm(l,m).im*=(4.0*PI/this->hp_map.Npix());
	std::cout<<l<<" "<<m<<" "<<this->alm(l,m).re<<" "<<this->alm(l,m).im<<std::endl;

	}
	}*/
      }
}





//==========================================================================
//==============================Cl==========================================



void spherical_shell::calc_cl(bool MASK, std::string outputfile){

  arr<double> temp_ps(this->lmax+1);
  double skyfrac;
  double skyfrac_weighted=0.0; //PCM
 
  PowSpec powspec(1, this->lmax);
  if(MASK){
    skyfrac=double(this->mask_area)/double(this->hp_map.Npix());
    for (int i; i< this->mask_area; i++ ){ //PCM
      double mask=this->hp_mask[int(this->pixelpos[i])]; 
      skyfrac_weighted += mask*mask ; 
  }    
  skyfrac_weighted /= double(this->hp_map.Npix()); //PCM  
  }else{
    skyfrac=1.0;
    skyfrac_weighted=1.0; //PCM
  }
  for(int i=0; i<this->hp_map.Npix();i++){
	if(MASK) if(this->hp_mask[i]!=0) this->gal_number+=this->hp_map[i];
	else this->gal_number+=this->hp_map[i];
  }
  double shotnoise=1.0/(this->gal_number/(4.0*pi*skyfrac));
  std::cout<<"Skyfraction: "<<skyfrac_weighted<<std::endl;
  arr<double>winpix(4*this->nside);
  //std::string winpix_directory="/home/pablo/codes/Healpix_3.11/data";
  read_pixwin(winpix_directory,this->nside, winpix);
  

  std::cout<<"Shotnoise: "<<shotnoise<<std::endl;

  std::ofstream file(&outputfile[0]);
  if(file.is_open()){
    file<<"RMS:    Average: "<<std::endl;
    file<<hp_map.rms()<<"   "<<hp_map.average()<<std::endl;
    for(int l=0; l<this->lmax; l++){
      temp_ps[l]=0.0;
      for(int m=0; m<=l; m++){
	temp_ps[l]+=( pow(this->alm(l,m).re, 2.0) + pow(this->alm(l,m).im, 2.0) ); 
	if(m==0) temp_ps[l]/=2.0;   
      }
       
      temp_ps[l]/=((double(l)+0.5)*skyfrac_weighted); //PCM
      temp_ps[l]/=(winpix[l]*winpix[l]);
      //temp_ps[l]/=pow(this->hp_map.average(),2.0) ;
      //temp_ps[l]-=shotnoise;
      file<<l<<std::setw(15)<<temp_ps[l]<<std::endl;
    }
  }else{
    std::cout<<"Cannot open outputfile! "<<RED<<" Exit now."<<RESET<<std::endl;
    exit(-1);
  }
  this->powspec.Set(temp_ps);
}

void spherical_shell::calc_cross(bool MASK, std::string outputfile, std::string Alm2file, std::string map2file){

  arr<double> temp_ps(this->lmax+1);
  double skyfrac;
  double skyfrac_weighted=0.0; //PCM

  Alm<xcomplex <double> > alm2;
  alm2.Set(this->lmax, this->mmax);
  read_Alm_from_fits(Alm2file, alm2,  this->lmax, this->mmax,2);
 
  Healpix_Map<double> map2(9, RING);
  read_Healpix_map_from_fits(map2file,map2, 1,2);
  this->gal_number=0.0;
  int galnum2=0;
  for(int i=0; i<this->hp_map.Npix(); i++){
    if(MASK) if(this->hp_mask[i]!=0) {
      galnum2+=map2[i];
      this->gal_number+=this->hp_map[i];
    }else{galnum2+=map2[i]; this->gal_number+=this->hp_map[i];}
  }

  PowSpec powspec(1, this->lmax);
  if(MASK){
    skyfrac=double(this->mask_area)/double(this->hp_map.Npix());
    for (int i; i< this->mask_area; i++ ){ //PCM
      double mask=this->hp_mask[int(this->pixelpos[i])]; 
      skyfrac_weighted += mask*mask ; 
  }    
  skyfrac_weighted /= double(this->hp_map.Npix()); //PCM  
  }else{
    skyfrac=1.0;
    skyfrac_weighted=1.0; //PCM
  }

  double shotnoise1=1.0/(this->gal_number/(4.0*pi*skyfrac));
  double shotnoise2=1.0/(galnum2/(4.0*pi*skyfrac));
  std::cout<<"Skyfraction: "<<skyfrac_weighted<<std::endl;
  arr<double>winpix(4*this->nside);
  //std::string winpix_directory="/home/lwolz/Libraries/Healpix_donatello/Healpix_2.20a/data/";
  read_pixwin(winpix_directory,this->nside, winpix);
  

  std::cout<<"Shotnoise1: "<<shotnoise1<<std::endl;
  std::cout<<"Shotnoise2: "<<shotnoise2<<std::endl;

  std::ofstream file(&outputfile[0]);
  if(file.is_open()){
    file<<"#RMS1:    Average: "<<std::endl;
    file<<hp_map.rms()<<"   "<<hp_map.average()<<std::endl;
    file<<"#RMS2:    Average: "<<std::endl;
    file<<map2.rms()<<"   "<<map2.average()<<std::endl;
    for(int l=0; l<this->lmax; l++){
      temp_ps[l]=0.0;
      for(int m=0; m<=l; m++){
        temp_ps[l]+= (this->alm(l,m).re * alm2(l,m).re )+( this->alm(l,m).im * alm2(l,m).im ); 
        if(m==0) temp_ps[l]/=2.0;   
      }
       
      temp_ps[l]/=((double(l)+0.5)*skyfrac_weighted); //PCM
      temp_ps[l]/=(winpix[l]*winpix[l]);
      //temp_ps[l]/=pow(this->hp_map.average(),2.0) ;
      //temp_ps[l]-=shotnoise;
      file<<l<<std::setw(15)<<temp_ps[l]<<std::endl;
    }
  }else{
    std::cout<<"Cannot open outputfile! "<<RED<<" Exit now."<<RESET<<std::endl;
    exit(-1);
  }
  this->powspec.Set(temp_ps);
}

float get_Ntot_weighted(const char* nameArq,const char* compararkey)
    {
      std::cout<<"Passei aqui A"<<RED<<std::endl;
      fitsfile *fptr;         
        char card[FLEN_CARD]; 
        //char *compararkey=argv[1];
        char *nr;
        float result;
        

      int status = 0,  nkeys, ii; 

      fits_open_file(&fptr, nameArq, READONLY, &status); // open fits file
      std::cout<<"Passei aqui B"<<RED<<std::endl;
        fits_get_hdrspace(fptr, &nkeys, NULL, &status);// read header
        std::cout<<"Passei aqui C"<<RED<<std::endl;
        
        for (ii = 1; ii <= nkeys; ii++)  { //search for a key
          fits_read_record(fptr, ii, card, &status); 
          
          if(strstr(card,compararkey)!=NULL){
                   nr= strtok (card,"=");//tokenizer for get only number part of key
                   nr= strtok (NULL,"=");  
                   result=atof(nr);        

          }
         
        }
        
        fits_close_file(fptr, &status);
        return result;
     
    }



void spherical_shell::calc_PeeblesCl(bool OVERDENSITY, std::string outputfile, std::string mapfile, bool GALAXYCOUNTS, bool SHOTOVERDENSITY){

  arr<double> Cltemp(this->lmax+1);

  double deltaomega=0.0;
  double complexabs=0.0;
  read_Healpix_map_from_fits(mapfile, this->hp_map, 1,2);
  arr<double>winpix(4*this->nside);
  //std::string winpix_directory="/home/lwolz/Library/Healpix_2.20a/data/";
  read_pixwin(winpix_directory,this->nside, winpix);
  /* std::ofstream winpixfile("/home/lwolz/new/maayane/Winpix_512.dat");
  for(int i =0; i<=this->lmax; i++){
    winpixfile<<i<<"  "<<winpix[i]<<std::endl;
  }
  winpixfile.close(); */
  deltaomega=double(this->mask_area)*4.0*PI/double(this->hp_map.Npix());
  this->gal_number=0.0;
  
  for(int i=0; i<this->hp_map.Npix(); i++){
    if(this->hp_mask[i]!=0) this->gal_number+=this->hp_map[i];
  }
  double  skyfrac=double(this->mask_area)/double(this->hp_map.Npix());
  double shotnoise=1.0/(this->gal_number/(4.0*pi*skyfrac));
  double omega_pix=4.0*pi/this->hp_map.Npix();
  double factor=this->gal_number/deltaomega*omega_pix;


  std::cout<<"galaxy N: "<<this->gal_number<<std::endl<<"mask area: "<<this->mask_area<<std::endl;
  if(GALAXYCOUNTS) std::cout<<"Shotnoise: "<<shotnoise<<std::endl;
  double jlmsum=0.0;
  double testsum=0.0;
  for(int i=0; i<this->lmax; i++) Cltemp[i]=0.0;
  double average=this->gal_number/this->mask_area;
  std::ofstream file(&outputfile[0]);
  if(file.is_open()){
    file<<"RMS:    Average: "<<std::endl;
    file<<hp_map.rms()<<"   "<<average<<std::endl;

    if(OVERDENSITY){
      for(int l=0; l<this->lmax; l++){
        for(int m=0; m<=l; m++){
          Cltemp[l]+=( (pow(this->alm(l,m).re,2.0)+pow(this->alm(l,m).im, 2.0))/ this->J[l][m]);
          if(m==0) Cltemp[l]/=2.0;
        }
        Cltemp[l]/=(double(l)+0.5) ;
        if (SHOTOVERDENSITY){
          std::cout<<"Passei aqui 1"<<RED<<std::endl;
          std::string s = mapfile+"[1]";
          int n = s.length();
          char char_array[n + 1];
          strcpy(char_array, s.c_str());
          std::cout<<"Passei aqui 2"<<RED<<std::endl;
          
          const char *file=char_array;

          std::cout<<"Passei aqui 3"<<RED<<std::endl;

          const char *key="HIERARCH ntot_weighted";

          std::cout<<"Passei aqui 4"<<RED<<std::endl;

          float nr=get_Ntot_weighted(file,key);
          std::cout<<"Passei aqui 5"<<RED<<std::endl;
          //Cltemp[l]-=shotnoise;
          Cltemp[l]-=nr;
          std::cout<<"Passei aqui 6"<<RED<<std::endl;
          
        }
        Cltemp[l]/=(winpix[l]*winpix[l]); //PCM
        file<<l<<std::setw(15)<<Cltemp[l]<<std::endl;
      }
    }else{
      for(int l=0; l<this->lmax; l++){
        jlmsum=0.0;
        for(int m=0; m<=l; m++){
          complexabs=pow(this->alm(l,m).re- factor*this->I[l][m].re, 2.0)+pow(this->alm(l,m).im- factor*this->I[l][m].im, 2.0);
          Cltemp[l]+= (complexabs);// this->J[l][m]);
          jlmsum+=this->J[l][m];
      	  if(m==0) {
            Cltemp[l]/=2.0;
            jlmsum/=2.0;
          }
        }
        Cltemp[l]/=jlmsum;
        //Cltemp[l]/=(double(l)+0.5);
        //	Cltemp[l]/=pow(average,2.0);
        if(GALAXYCOUNTS) Cltemp[l]-=shotnoise;
        Cltemp[l]/=(winpix[l]*winpix[l]);
        file<<l<<std::setw(15)<<Cltemp[l]<<std::endl;
      }
    }
  }else{
    std::cout<<"Cannot open outputfile! "<<RED<<" Exit now."<<RESET<<std::endl;
    exit(-1);
  }
  this->Peebles.Set(Cltemp);

}

void spherical_shell::calc_PeeblesCross(bool OVERDENSITY, std::string outputfile, std::string mapfile, bool GALAXYCOUNTS, std::string Alm2file, std::string map2file){

  arr<double> Cltemp(this->lmax+1);

  double deltaomega=0.0;
  double complexabs=0.0;

  arr<double>winpix(4*this->nside);
//  std::string winpix_directory="/home/lwolz/Library/Healpix_2.20a/data/";
  read_pixwin(winpix_directory,this->nside, winpix);

  Alm<xcomplex <double> > alm2;
  alm2.Set(this->lmax, this->mmax);
  read_Alm_from_fits(Alm2file, alm2,  this->lmax, this->mmax,2);

  deltaomega=double(this->mask_area)*4.0*PI/double(this->hp_map.Npix());

  read_Healpix_map_from_fits(mapfile, this->hp_map, 1,2);
  Healpix_Map<double> map2(9, RING);
  read_Healpix_map_from_fits(map2file,map2, 1,2);
  this->gal_number=0.0;
  int galnum2=0;
  for(int i=0; i<this->hp_map.Npix(); i++){
    galnum2+=map2[i];
    this->gal_number+=this->hp_map[i];
  }
 
  double  skyfrac=double(this->mask_area)/double(this->hp_mask.Npix());
  double shotnoise=1.0/(this->gal_number/(4.0*pi*skyfrac));
  double omega_pix=4.0*pi/this->hp_map.Npix();
  double factor=this->gal_number/deltaomega*omega_pix;
  double factor2=galnum2/deltaomega*omega_pix;

  std::cout<<"galaxy N: "<<this->gal_number<<std::endl<<"mask area: "<<this->mask_area<<std::endl;
  if(GALAXYCOUNTS) std::cout<<"Shotnoise: "<<shotnoise<<std::endl;

  double testsum=0.0;
  for(int i=0; i<this->lmax; i++) Cltemp[i]=0.0;
  double average=this->gal_number/this->mask_area;
  std::ofstream file(&outputfile[0]);
  if(file.is_open()){
    file<<"RMS:    Average: "<<std::endl;
    file<<hp_map.rms()<<"   "<<average<<std::endl;

    if(OVERDENSITY){
      for(int l=0; l<this->lmax; l++){
	for(int m=0; m<=l; m++){
	  xcomplex<double> fac1=this->alm(l,m);
	  xcomplex<double> fac2=alm2(l,m);
	  double val=fac1.re*fac2.re+fac1.im*fac2.im;
	  Cltemp[l]+= (val/ this->J[l][m]);
	  if(m==0) Cltemp[l]/=2.0;
	}
	Cltemp[l]/=(double(l)+0.5) ;
	Cltemp[l]/=(winpix[l]*winpix[l]); //PCM
	file<<l<<std::setw(15)<<Cltemp[l]<<std::endl;
      }
    }else{
      for(int l=0; l<this->lmax; l++){
	for(int m=0; m<=l; m++){
	 
	xcomplex<double> fac1=this->alm(l,m)-factor*this->I[l][m];
	xcomplex<double> fac2=alm2(l,m)-factor2*this->I[l][m];
	double val=fac1.re*fac2.re+fac1.im*fac2.im;

	  Cltemp[l]+= (val/ this->J[l][m]);
      	  if(m==0) Cltemp[l]/=2.0;

	}
	Cltemp[l]/=(double(l)+0.5);
	//	Cltemp[l]/=pow(average,2.0);
	if(GALAXYCOUNTS) Cltemp[l]-=shotnoise;
	Cltemp[l]/=(winpix[l]*winpix[l]);

	file<<l<<std::setw(15)<<Cltemp[l]<<std::endl;
      }
    }
  }else{
    std::cout<<"Cannot open outputfile! "<<RED<<" Exit now."<<RESET<<std::endl;
    exit(-1);
  }
  this->Peebles.Set(Cltemp);

}







//============================================================================
//================   Ilm Jlm ===============================================

void spherical_shell::calc_I_lm_J_lm(std::string IlmJlmfile){

  this->I=new xcomplex<double>*[this->lmax];
  this->J=new double*[this->lmax];
  for(int i=0; i<this->lmax; i++){
    this->I[i]=new xcomplex<double>[this->lmax];
    this->J[i]=new double[this->lmax];
  }
  pointing pointer;
  double ylm=0.0;
  xcomplex<double> Ylm;
  double pixarea=4.0*PI/double(this->hp_map.Npix());
  double deltaomega=double(this->mask_area)*4.0*PI/double(this->hp_map.Npix());
  for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	this->I[l][m].re=0.0;
	this->I[l][m].im=0.0;
	this->J[l][m]=0.0;

      }
  }

  std::cout<<IlmJlmfile<<std::endl;
  std::ofstream file(&IlmJlmfile[0]);
  if(file.is_open()){
for(int l=0; l<this->lmax; l++){
  std::cout<<l<<"   ";
      for(int m=0; m<=l; m++){
	for(int i=0; i<this->mask_area; i++){
	  pointer=this->hp_map.pix2ang(double(this->pixelpos[i]));
	  ylm=gsl_sf_legendre_sphPlm(double(l), abs(double(m)) ,cos(pointer.theta));
	  Ylm.re=ylm*cos(double(m)*pointer.phi);
	  Ylm.im=ylm*sin(double(m)*pointer.phi);     
	  //if(m<0) Ylm=pow(-1.0, double(m))*Ylm.conj();
	  this->I[l][m].re+=(Ylm.re*pixarea);
	  this->I[l][m].im+=(-Ylm.im*pixarea);	
	  this->J[l][m]+=(pow(abs(Ylm), 2.0)*pixarea);  
	}

      }
 }
 
    for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	file<<l<<"  "<<m<<"  "<<I[l][m].re<<"   "<<I[l][m].im<<"   "<<J[l][m]<<std::endl;
      }
    }
    file.close();
  }else{
    std::cout<<"Cannot write to file in Ilm and Jlm! "<<RED<<" Exit now."<<RESET<<std::endl;
    exit(-1);
  }
}


void spherical_shell::calc_I_lm_J_lm_weighted(std::string IlmJlmfile){ //PCM

  this->I=new xcomplex<double>*[this->lmax];
  this->J=new double*[this->lmax];
  for(int i=0; i<this->lmax; i++){
    this->I[i]=new xcomplex<double>[this->lmax];
    this->J[i]=new double[this->lmax];
  }
  pointing pointer;
  double ylm=0.0;
  xcomplex<double> Ylm;
  double pixarea=4.0*PI/double(this->hp_map.Npix());
  double deltaomega=double(this->mask_area)*4.0*PI/double(this->hp_map.Npix());
  for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	this->I[l][m].re=0.0;
	this->I[l][m].im=0.0;
	this->J[l][m]=0.0;

      }
  }

  std::cout<<IlmJlmfile<<std::endl;
  std::ofstream file(&IlmJlmfile[0]);
  if(file.is_open()){
for(int l=0; l<this->lmax; l++){
  std::cout<<l<<"   ";
      for(int m=0; m<=l; m++){
	for(int i=0; i<this->mask_area; i++){
	  pointer=this->hp_map.pix2ang(double(this->pixelpos[i]));
	  double mask=this->hp_mask[int(this->pixelpos[i])]; //PCM
	  ylm=gsl_sf_legendre_sphPlm(double(l), abs(double(m)) ,cos(pointer.theta));
	  Ylm.re=ylm*cos(double(m)*pointer.phi);
	  Ylm.im=ylm*sin(double(m)*pointer.phi);     
	  //if(m<0) Ylm=pow(-1.0, double(m))*Ylm.conj();
	  this->I[l][m].re+=(Ylm.re*pixarea)*mask; // PCM
	  this->I[l][m].im+=(-Ylm.im*pixarea)*mask;  //PCM	
	  this->J[l][m]+=(pow(abs(Ylm), 2.0)*pixarea)*mask*mask;  //PCM  
	}

      }
 }
 
    for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	file<<l<<"  "<<m<<"  "<<I[l][m].re<<"   "<<I[l][m].im<<"   "<<J[l][m]<<std::endl;
      }
    }
    file.close();
  }else{
    std::cout<<"Cannot write to file in Ilm and Jlm! "<<RED<<" Exit now."<<RESET<<std::endl;
    exit(-1);
  }
}

void spherical_shell::calc_I_lm_J_lm(std::string IlmJlmfile, int startl, int stopl){

  this->I=new xcomplex<double>*[this->lmax];
  this->J=new double*[this->lmax];
  for(int i=0; i<this->lmax; i++){
    this->I[i]=new xcomplex<double>[this->lmax];
    this->J[i]=new double[this->lmax];
  }
  pointing pointer;
  double ylm=0.0;
  xcomplex<double> Ylm;
  double pixarea=4.0*PI/double(this->hp_map.Npix());
  double deltaomega=double(this->mask_area)*4.0*PI/double(this->hp_map.Npix());
  for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	this->I[l][m].re=0.0;
	this->I[l][m].im=0.0;
	this->J[l][m]=0.0;

      }
  }

  std::cout<<IlmJlmfile<<std::endl;
  std::ofstream file(&IlmJlmfile[0]);
  if(file.is_open()){
for(int l=startl; l< stopl; l++){
  std::cout<<l<<"   ";
      for(int m=0; m<=l; m++){
	for(int i=0; i<this->mask_area; i++){
	  pointer=this->hp_map.pix2ang(double(this->pixelpos[i]));
	  ylm=gsl_sf_legendre_sphPlm(double(l), abs(double(m)) ,cos(pointer.theta));
	  Ylm.re=ylm*cos(double(m)*pointer.phi);
	  Ylm.im=ylm*sin(double(m)*pointer.phi);     
	  //if(m<0) Ylm=pow(-1.0, double(m))*Ylm.conj();
	  this->I[l][m].re+=(Ylm.re*pixarea);
	  this->I[l][m].im+=(-Ylm.im*pixarea);	
	  this->J[l][m]+=(pow(abs(Ylm), 2.0)*pixarea);  
	}

      }
 }
 
    for(int l=startl; l<stopl ; l++){
      for(int m=0; m<=l; m++){
	file<<l<<"  "<<m<<"  "<<this->I[l][m].re<<"   "<<this->I[l][m].im<<"   "<<this->J[l][m]<<std::endl;
      }
    }
    file.close();
  }else{
    std::cout<<"Cannot write to file in Ilm and Jlm! "<<RED<<" Exit now."<<RESET<<std::endl;
    exit(-1);
  }
}

void spherical_shell::calc_I_lm_J_lm_weighted(std::string IlmJlmfile, int startl, int stopl){ //PCM

  this->I=new xcomplex<double>*[this->lmax];
  this->J=new double*[this->lmax];
  for(int i=0; i<this->lmax; i++){
    this->I[i]=new xcomplex<double>[this->lmax];
    this->J[i]=new double[this->lmax];
  }
  pointing pointer;
  double ylm=0.0;
  xcomplex<double> Ylm;
  double pixarea=4.0*PI/double(this->hp_map.Npix());
  double deltaomega=double(this->mask_area)*4.0*PI/double(this->hp_map.Npix());
  for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	this->I[l][m].re=0.0;
	this->I[l][m].im=0.0;
	this->J[l][m]=0.0;

      }
  }

  std::cout<<IlmJlmfile<<std::endl;
  std::ofstream file(&IlmJlmfile[0]);
  if(file.is_open()){
for(int l=startl; l< stopl; l++){
  std::cout<<l<<"   ";
      for(int m=0; m<=l; m++){
	for(int i=0; i<this->mask_area; i++){
	  pointer=this->hp_map.pix2ang(double(this->pixelpos[i]));
	  double mask=this->hp_mask[int(this->pixelpos[i])]; //PCM
	  ylm=gsl_sf_legendre_sphPlm(double(l), abs(double(m)) ,cos(pointer.theta));
	  Ylm.re=ylm*cos(double(m)*pointer.phi);
	  Ylm.im=ylm*sin(double(m)*pointer.phi);     
	  //if(m<0) Ylm=pow(-1.0, double(m))*Ylm.conj();
	  this->I[l][m].re+=(Ylm.re*pixarea)*mask; // PCM
	  this->I[l][m].im+=(-Ylm.im*pixarea)*mask;  //PCM	
	  this->J[l][m]+=(pow(abs(Ylm), 2.0)*pixarea)*mask*mask;  //PCM  	   
	}

      }
 }
 
    for(int l=startl; l<stopl ; l++){
      for(int m=0; m<=l; m++){
	file<<l<<"  "<<m<<"  "<<this->I[l][m].re<<"   "<<this->I[l][m].im<<"   "<<this->J[l][m]<<std::endl;
      }
    }
    file.close();
  }else{
    std::cout<<"Cannot write to file in Ilm and Jlm! "<<RED<<" Exit now."<<RESET<<std::endl;
    exit(-1);
  }
}




void spherical_shell::read_I_lm_J_lm(std::string IlmJlmfile){ 
  std::ifstream file(IlmJlmfile.c_str());
  int templ, tempm;
  this->I=new xcomplex<double>*[this->lmax];
  this->J=new double*[this->lmax];
  for(int i=0; i<this->lmax; i++){
    this->I[i]=new xcomplex<double>[this->lmax];
    this->J[i]=new double[this->lmax];
  }
  if(file.is_open()){
    for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	file>>templ>>tempm;
        file>>I[templ][tempm].re>>I[templ][tempm].im>>J[templ][tempm];
      }
    }
	std::cout<<"And here 333\n"<<std::flush;
  }else{
    std::cout<<"Cannot read Ilm and Jlm! "<<RED<<" Exit now."<<RESET<<std::endl;
    exit(-1);
  }
}



void spherical_shell::draw_alm(PowSpec Input,  int ll){ 

std::ofstream file;
  file.open("./alm_realisation.dat");
  xcomplex<double> **atemp;
  atemp=new xcomplex<double>*[ll];
  for(int l=0; l<ll; l++) atemp[l]=new xcomplex<double>[2*l+1];
  gsl_rng *random=gsl_rng_alloc(gsl_rng_default);
  int seed=rand();
  gsl_rng_set(random, seed);
  for(int l=0; l<ll;l++){
    atemp[l][l+0].im=0.0;
    atemp[l][l+0].re=(gsl_ran_gaussian(random,sqrt(0.5*Input.tt(l))) );
    file<<l<<" "<<0<<" "<<atemp[l][l+0].re<<std::setw(15)<<atemp[l][l+0].im<<std::setw(15)<<atemp[l][l+0].re<<std::setw(15)<<atemp[l][l+0].im<<std::endl;
    for(int m=1; m<=l; m++){
      atemp[l][l+m].re=(gsl_ran_gaussian(random,sqrt(0.5*Input.tt(l))) );
      atemp[l][l+m].im=(gsl_ran_gaussian(random,sqrt(0.5*Input.tt(l))) );
      atemp[l][-m+l].re=pow(-1.0, m)*atemp[l][l+m].re;
      atemp[l][-m+l].im=-pow(-1.0, m)*atemp[l][l+m].im;  
      if(abs(atemp[l][-m+l])!=abs(atemp[l][l+m])) std::cout<<"alm"<<std::endl;  
      file<<l<<" "<<m<<" "<<atemp[l][l+m].re<<std::setw(15)<<atemp[l][l+m].im<<std::setw(15)<<atemp[l][-m+l].re<<std::setw(15)<<atemp[l][-m+l].im<<std::endl;
    }
  }
  //this->alm.Set(atemp, this->lmax, this->mmax);
  file.close();
}



void spherical_shell::create_map( int ll){

  double ylm=0.0;
  xcomplex<double> **atemp;
  atemp=new xcomplex<double>*[ll];
  for(int l=0; l<ll; l++) atemp[l]=new xcomplex<double>[2*l+1];
  std::ifstream file;
  file.open("/home/lwolz/Documents/TEST/alm_realisation.dat");
  int t;
 for(int l=0; l<ll;l++){
   for(int m=0; m<=l; m++){
     file>>t>>t>>atemp[l][l+m].re>>atemp[l][l+m].im>>atemp[l][-m+l].re>>atemp[l][-m+l].im;
   }
 }
  xcomplex<double> temp,round,  Ylm, Ylm2;
  temp.re=0.0, temp.im=0.0;
  pointing pointer;
  for(int i=0; i<this->hp_map.Npix(); i++){  this->hp_map[i]=0.0;
    pointer=this->hp_map.pix2ang(i);
    temp.re=0.0;
    temp.im=0.0;
    for(int l=0; l<ll;l++){
      //m=0
      round.im=0.0;
      ylm=gsl_sf_legendre_sphPlm(double(l),abs(double(0)),cos(pointer.theta));
      Ylm.re=ylm*cos(double(0)*pointer.phi)/pow(-1.0,double(0)) ;
      Ylm.im=-ylm*sin(double(0)*pointer.phi)/pow(-1.0,double(0)) ;
      temp.re+=(( Ylm.re*atemp[l][0+l].re - Ylm.im*atemp[l][0+l].im));
      temp.im+=((atemp[l][0+l].re*Ylm.im + atemp[l][0+l].im*Ylm.re));
      round.im+=((atemp[l][0+l].re*Ylm.im + atemp[l][0+l].im*Ylm.re));
      
      for(int m=1; m<=l; m++){
	
	  ylm=gsl_sf_legendre_sphPlm(double(l),abs(double(m)),cos(pointer.theta));
	  //negative part -m
	  Ylm.re=ylm*cos(double(m)*pointer.phi)/pow(-1.0,-double(m));
	  Ylm.im=-ylm*sin(double(m)*pointer.phi)/pow(-1.0,-double(m)) ;
	  temp.re+=(( Ylm.re*atemp[l][-m+l].re - Ylm.im*atemp[l][-m+l].im));
	  temp.im+=((atemp[l][-m+l].re*Ylm.im + atemp[l][-m+l].im*Ylm.re));
	  round.im+=((atemp[l][-m+l].re*Ylm.im + atemp[l][-m+l].im*Ylm.re));
	 
	  //positive part +m
	  ylm=gsl_sf_legendre_sphPlm(double(l),abs(double(m)),cos(pointer.theta));

	  Ylm2.re=ylm*cos(double(m)*pointer.phi);
	  Ylm2.im=ylm*sin(double(m)*pointer.phi);
	  temp.re+=(( Ylm2.re*atemp[l][m+l].re - Ylm2.im*atemp[l][m+l].im));
	  temp.im+=((atemp[l][m+l].re*Ylm2.im + atemp[l][m+l].im*Ylm2.re));


	  if(abs(atemp[l][-m+l].re*Ylm.im + atemp[l][-m+l].im*Ylm.re)!=abs(atemp[l][m+l].re*Ylm2.im + atemp[l][m+l].im*Ylm2.re)){
	    std::cout<<"pixel: "<<i<<std::endl;
	    std::cout<<"m: "<<l<<" "<<m<<" "<<((atemp[l][-m+l].re*Ylm.im + atemp[l][-m+l].im*Ylm.re))<<std::endl;
	    std::cout<<"m: "<<l<<" "<<m<<" "<<((atemp[l][m+l].re*Ylm2.im + atemp[l][m+l].im*Ylm2.re))<<std::endl;
	    std::cout<<Ylm<<"  "<<Ylm2<<std::endl;
	    std::cout<<atemp[l][m+l]<<"   "<<atemp[l][-m+l]<<std::endl;
	  }
	  round.im+=((atemp[l][m+l].re*Ylm2.im + atemp[l][m+l].im*Ylm2.re));

      }
      if(round.im>0.001) std::cout<<i<<"  "<<l<<"  "<<round.im<<std::endl;
    }
    if(temp.im>0.001) std::cout<<"Imaginary part : "<<i<<"  "<<temp.im<<"  "<<ylm<<"  "<<pointer<<std::endl;
    this->hp_map[i]=temp.re;
  
  }
}
void spherical_shell::create_map(){

  double ylm=0.0;
  
  xcomplex<double> temp,round,  Ylm, Ylm2;
  temp.re=0.0, temp.im=0.0;
  pointing pointer;
  for(int i=0; i<this->hp_map.Npix(); i++){  this->hp_map[i]=0.0;
    pointer=this->hp_map.pix2ang(i);
    temp.re=0.0;
    temp.im=0.0;
    for(int l=0; l<this->lmax;l++){
      //m=0
      round.im=0.0;
      ylm=gsl_sf_legendre_sphPlm(double(l),abs(double(0)),cos(pointer.theta));
      Ylm.re=ylm*cos(double(0)*pointer.phi)/pow(-1.0,double(0)) ;
      Ylm.im=-ylm*sin(double(0)*pointer.phi)/pow(-1.0,double(0)) ;
      temp.re+=(( Ylm.re*this->alm(l,0+l).re - Ylm.im*this->alm(l,0+l).im));
      temp.im+=((this->alm(l,0+l).re*Ylm.im + this->alm(l,0+l).im*Ylm.re));
      round.im+=((this->alm(l,0+l).re*Ylm.im + this->alm(l,0+l).im*Ylm.re));
      
      for(int m=1; m<=l; m++){
	
	  ylm=gsl_sf_legendre_sphPlm(double(l),abs(double(m)),cos(pointer.theta));
	  //negative part -m
	  Ylm.re=ylm*cos(double(m)*pointer.phi)/pow(-1.0,-double(m));
	  Ylm.im=-ylm*sin(double(m)*pointer.phi)/pow(-1.0,-double(m)) ;
	  temp.re+=(( Ylm.re*this->alm(l,-m+l).re - Ylm.im*this->alm(l,-m+l).im));
	  temp.im+=((this->alm(l,-m+l).re*Ylm.im + this->alm(l,-m+l).im*Ylm.re));
	  round.im+=((this->alm(l,-m+l).re*Ylm.im + this->alm(l,-m+l).im*Ylm.re));
	 
	  //positive part +m
	  ylm=gsl_sf_legendre_sphPlm(double(l),abs(double(m)),cos(pointer.theta));

	  Ylm2.re=ylm*cos(double(m)*pointer.phi);
	  Ylm2.im=ylm*sin(double(m)*pointer.phi);
	  temp.re+=(( Ylm2.re*this->alm(l,m+l).re - Ylm2.im*this->alm(l,m+l).im));
	  temp.im+=((this->alm(l,m+l).re*Ylm2.im + this->alm(l,m+l).im*Ylm2.re));


	  if(abs(this->alm(l,-m+l).re*Ylm.im + this->alm(l,-m+l).im*Ylm.re)!=abs(this->alm(l,m+l).re*Ylm2.im + this->alm(l,m+l).im*Ylm2.re)){
	    std::cout<<"pixel: "<<i<<std::endl;
	    std::cout<<"m: "<<l<<" "<<m<<" "<<((this->alm(l,-m+l).re*Ylm.im + this->alm(l,-m+l).im*Ylm.re))<<std::endl;
	    std::cout<<"m: "<<l<<" "<<m<<" "<<((this->alm(l,m+l).re*Ylm2.im + this->alm(l,m+l).im*Ylm2.re))<<std::endl;
	    std::cout<<Ylm<<"  "<<Ylm2<<std::endl;
	    std::cout<<this->alm(l,m+l)<<"   "<<this->alm(l,-m+l)<<std::endl;
	  }
	  round.im+=((this->alm(l,m+l).re*Ylm2.im + this->alm(l,m+l).im*Ylm2.re));

      }
      if(round.im>0.001) std::cout<<i<<"  "<<l<<"  "<<round.im<<std::endl;
    }
    if(temp.im>0.001) std::cout<<"Imaginary part : "<<i<<"  "<<temp.im<<"  "<<ylm<<"  "<<pointer<<std::endl;
    this->hp_map[i]=temp.re;
  
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//mask

void spherical_shell::calc_mask_lognorm(std::string outputfile){
 
  pointing pointer;
  this->mask_area=0;
  for(int i=0; i<this->hp_mask.Npix(); i++){
    pointer=this->hp_mask.pix2ang(i);
    if(( pointer.phi<=(270.0*PI/180.0) &&  (pointer.phi)>=(110.0*PI/180.0) && pointer.theta>=((20.0)*PI/180.0) &&  pointer.theta<=((95.0)*PI/180.0))){
      this->hp_mask[i]=1.0;
      this->mask_area++;
    }else{
      this->hp_mask[i]=0.0;      
    }
  }
write_Healpix_map_to_fits(outputfile, this->hp_mask, PLANCK_FLOAT64);
}


void spherical_shell::calc_mask_halfsky(std::string outputfile){
 
  pointing pointer;
  this->mask_area=0;
  for(int i=0; i<this->hp_mask.Npix(); i++){
    pointer=this->hp_mask.pix2ang(i);
    if(( pointer.theta<=pi/2)  ){//|| (2.0*PI-pointer.phi)<=this->phi_lim )){
      this->hp_mask[i]=1.0;
      this->mask_area++;
    }else{
      this->hp_mask[i]=0.0;      
    }
  }
write_Healpix_map_to_fits(outputfile, this->hp_mask, PLANCK_FLOAT64);
}

void spherical_shell::calc_mask_fg(std::string outputfile, std::string fgfile, double tlim){
  Healpix_Map <double> fg;
  read_Healpix_map_from_fits(fgfile, fg, 1,2);
  this->mask_area=0;
  this->hp_mask.fill(0.0);
  for(int i=0; i<fg.Npix(); i++){
    if(fg[i]<=tlim) this->hp_mask[i]=1.0;
  }
write_Healpix_map_to_fits(outputfile, this->hp_mask, PLANCK_FLOAT64);
}





void spherical_shell::calc_mask(std::string outputfile, double *theta, double *phi){
 
  pointing pointer;
  this->mask_area=0;
  for(int i=0; i<this->hp_mask.Npix(); i++){
    pointer=this->hp_mask.pix2ang(i);

    if(( pointer.theta>=theta[0] && pointer.theta<=theta[1]) && (pointer.phi>=phi[0] && pointer.phi<=phi[1])  ) {
	
      this->hp_mask[i]=1.0;
      this->mask_area++;
    }else{
      this->hp_mask[i]=0.0;      
    }
  }
  write_Healpix_map_to_fits<float64>(outputfile, this->hp_mask, PLANCK_FLOAT64);

}




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//===============noise=======================


void  noise_map(double *zlim, double fwhm, std::string outputfile, int nside){

  int order=log2(nside);

  Healpix_Map <double> noisemap(order,RING);
  Healpix_Map <double> noisemap2(order,RING);
 double f21cm= 1420.0; //inMHz 
 double f0=f21cm*1e6/(1.0+zlim[1]), f1=f21cm*1e6/(1.0+zlim[2]);
 double deltaf=(f1-f0);
 double nu=f0+deltaf/2;

  double Omega_pix=(4.0*pi)/noisemap.Npix();

  noisemap.fill(0.0);
  noisemap2.fill(0.0);

  double integration=10.0*60.0*300.0;

  double Omega_beam=1.133*pow(fwhm,2.0); //Gaussian beam
 
  std::cout<<"Omega_Gauss= "<<Omega_beam<<std::endl;
  std::cout<<"FWHM= "<<fwhm*180.0/pi<<std::endl;

  int lmax=int(sqrt(12)*noisemap.Nside())-1;
  double frac=0.1;//10 percent of full SKA science requirements
  double AeffoverTsys=2.0e4;
  double rmstemp=(pow(speedOfLight,2.0)) /(pow(nu,2.0)*frac*AeffoverTsys*sqrt(2.0*integration*deltaf)*Omega_beam);
  std::cout<<"frequency: "<<nu<<" A_eff/T_sys: "<<AeffoverTsys<<" t_int: "<<integration<<" Omega_beam: "<<Omega_beam<<std::endl;
  double rmstemp_pix=(pow(speedOfLight,2.0)) /(pow(nu,2.0)*frac*AeffoverTsys*sqrt(2.0*integration*deltaf)*Omega_pix);

  std::cout<<MAGENTA<<"Temperature RMS: "<<rmstemp<<RESET<<std::endl;
 
  std::cout<<MAGENTA<<"Temperature RMS pixel: "<<rmstemp_pix<<RESET<<std::endl;
 arr<double> weight(lmax);
 for(int j=0;j<lmax;j++) weight[j]=1.0;
 Alm< xcomplex <double> > alm(lmax,lmax) ;
 gsl_rng *random=gsl_rng_alloc(gsl_rng_default);
 int seed=time(NULL);
 std::cout<<"Seed: "<<seed<<std::endl;
 gsl_rng_set(random, seed);
 for(int i=0; i<noisemap.Npix(); i++){
 
   noisemap[i]=gsl_ran_gaussian(random, 1.0 );
  
  }
 gsl_rng_free(random);
 alm.SetToZero();
 double rms= noisemap.rms();
 std::cout<<"RMS of Gaussian: "<<rms<<std::endl;
 
 map2alm( noisemap, alm, weight, false) ;  
  smoothWithGauss(alm, fwhm); //fwhm in radian
 PowSpec cl(1, lmax);
 extract_powspec(alm, cl);
 /*std::ofstream file("Powspec_test.dat");
 for(int l=0; l<lmax;l++){
   file<<l<<std::setw(15)<<cl.tt(l)<<std::endl;
 }*/
 alm2map(alm, noisemap2);
 rms= noisemap2.rms();
 
 std::cout<<"RMS after smooth: "<<rms<<std::endl;

 noisemap2.Scale(rmstemp/rms);
 double rmsnew= noisemap2.rms();

std::cout<<"RMS after rescale: "<<rmsnew<<std::endl;

map2alm(noisemap2, alm, weight, false);
smoothWithGauss(alm, -fwhm);
alm2map(alm, noisemap);

rms=noisemap.rms();
std::cout<<"RMS of deconvolved map: "<<rms<<std::endl;


// std::cout<<"RMS after rescale pixel withour smoothing: "<<rmsnew<<std::endl;


 write_Healpix_map_to_fits(outputfile, noisemap2, PLANCK_FLOAT64);

}


void  noise_map(int numb, double fwhm, std::string outputfile, int nside){

  int order=log2(nside);

  Healpix_Map <double> noisemap(order,RING);
  Healpix_Map <double> noisemap2(order,RING);
 //double f21cm= 1420.0; //inMHz 
 double f0=570.0e6, f1=1400.0e6;
 int nmaps=160;
 double deltaf=(f1-f0)/double(nmaps);
 double nu=(f0+double(numb)*deltaf+f0+ double(numb+1)*deltaf)/2.0;


  noisemap.fill(0.0);
  noisemap2.fill(0.0);

  double integration=77.0;

  double Omega_beam=1.133*pow(fwhm,2.0); //Gaussian beam
 
  std::cout<<"Omega_Gauss= "<<Omega_beam<<std::endl;
  std::cout<<"FWHM= "<<fwhm*180.0/pi<<std::endl;

  int lmax=int(sqrt(12)*noisemap.Nside())-1;
  double frac=1.0;//100percent of full SKA science requirements
  double AeffoverTsys=2.0e4;
  double rmstemp=(pow(speedOfLight,2.0)) /(pow(nu,2.0)*frac*AeffoverTsys*sqrt(2.0*integration*deltaf)*Omega_beam);
std::cout<<"Divided by 4   ";
  rmstemp/=4.0; //divide by 4 to decrease noiselevel
   std::cout<<"frequency: "<<nu<<" A_eff/T_sys: "<<AeffoverTsys<<" t_int: "<<integration<<" Omega_beam: "<<Omega_beam<<std::endl;
   std::cout<<MAGENTA<<"Temperature RMS: "<<rmstemp<<RESET<<std::endl;
  arr<double> weight(lmax);
 for(int j=0;j<lmax;j++) weight[j]=1.0;
 Alm< xcomplex <double> > alm(lmax,lmax) ;
 gsl_rng *random=gsl_rng_alloc(gsl_rng_default);
 int seed=time(NULL);
 gsl_rng_set(random, seed);
 for(int i=0; i<noisemap.Npix(); i++){
 
   noisemap[i]=gsl_ran_gaussian(random, 1.0 );
  
  }
 gsl_rng_free(random);
 alm.SetToZero();
 double rms= noisemap.rms();
 
 map2alm( noisemap, alm, weight, false) ;  
  smoothWithGauss(alm, fwhm); //fwhm in radian

 alm2map(alm, noisemap2);
 rms= noisemap2.rms();
 
 noisemap2.Scale(rmstemp/rms);
 double rmsnew= noisemap2.rms();

 std::cout<<rmsnew<<std::setw(25);

 noisemap.Scale(rmstemp/rms);
 rmsnew= noisemap.rms();

 std::cout<<rmsnew<<std::endl;


 write_Healpix_map_to_fits(outputfile, noisemap2, PLANCK_FLOAT64);

//TEST
  int orderdeg=log2(nside/2);
noisemap2.Add(1000.0);
for(int i=0; i<noisemap2.Npix()/2;i++) noisemap2[i]+=1000;
  Healpix_Map <double> noisemapdeg(orderdeg,RING);
  Healpix_Map <double> noisemapup(order,RING);
  std::cout<<"RMS orig: "<<noisemap2.average()<<std::endl;
  noisemapdeg.Import_degrade(noisemap2);
  std::cout<<"RMS degrade: "<<noisemapdeg.average()<<std::endl;
  noisemapup.Import_upgrade(noisemapdeg);
  std::cout<<"RMS upgrade: "<<noisemapdeg.average()<<std::endl;

}





void spherical_shell::calc_mix_matrix(std::string outputfile){
  /*std::cout<<"In calc mixmatrix"<<std::endl;
  std::cout<<"Initializing Wig 3j table..." << std::endl;

  wig_table_init(2*lmax,3);
  wig_temp_init(2*lmax);
  std::cout<<"Initialised Wig 3j table...DONE" << std::endl; */

  gsl_matrix *mixmatrix, *unmixmatrix;
  mixmatrix=gsl_matrix_alloc(  lmax,lmax);
  unmixmatrix=gsl_matrix_alloc(  lmax,lmax);
  double skyfrac=double(this->mask_area)/double(this->hp_map.Npix());
  double skyfrac_weighted=0.0; //PCM
  for (int i; i< this->mask_area; i++ ){
  	  double mask=this->hp_mask[int(this->pixelpos[i])]; 
  	  skyfrac_weighted += mask*mask ; 
  }    
  skyfrac_weighted /= double(this->hp_map.Npix()); //PCM    
  
  //gsl_permutation *p=gsl_permutation_alloc(lmax);
  int *signum;
  signum=new int;
  std::ofstream file(&outputfile[0]);
  double Rll;
  double W[lmax];
  //double wigner=0.0;
  for(int l=0; l<lmax; l++){
    W[l]=0.0;
    for(int m=0; m<=l; m++){  
      W[l]+=pow(abs(this->I[l][m]), 2.0);
      if(m==0) W[l]/=2.0;
    }
    W[l]/=( (double(l)+0.5)) ;
  }
  for(int l=0; l<lmax; l++){
    std::cout<<l <<"\t" <<std::flush;
    for(int ls=0; ls<lmax;ls++){  
      Rll=0.0;
      std::vector<double> array (2 * ((l < ls) ? l : ls ) + 1); //PCM
      uclwig3j::uclwig3j(l,ls,0,array.begin());
      int lss_min = abs(l - ls); 
      int lss = lss_min;
      int lss_max = l + ls + 1; 
      int count = 0; 
      while( lss < lss_max and lss < lmax ){
         double wigner = array [lss - lss_min];
         //printf("%lf, %lf, %lf, %lf\n", W[lss], wigner, skyfrac_weighted, skyfrac );
         Rll+=( (2.0*double(lss)+1.0)*wigner*wigner*W[lss]);
         lss++ ;       
      } //PCM      
      /*for(int lss=0; lss<lmax;lss++){
         wigner=wig3jj(2*l ,2*ls,2*lss ,0, 0, 0);
         Rll+=( (2.0*double(lss)+1.0)*wigner*wigner*W[lss]);        
      }*/
      Rll*=(2.0*double(ls)+1.0)/(4.0*PI)/skyfrac_weighted;//PCM
      file<<l<<std::setw(15)<<ls<<std::setw(15)<<Rll<<std::endl;
      //gsl_matrix_set(mixmatrix, l, ls, Rll);
    }
  }
  std::cout<< "\n" <<std::flush;  
  /* gsl_linalg_LU_decomp(mixmatrix, p, signum);
     gsl_linalg_LU_invert(mixmatrix, p, unmixmatrix);*/
  /*std::ofstream file(&outputfile[0]);
  for(int l=0; l<lmax; l++){
    for(int ls=0; ls<lmax;ls++){
      file<<gsl_matrix_get(mixmatrix, l, ls)<<std::setw(15);
    }
    file<<std::endl;
    }*/
  file.close();

}

void spherical_shell::calc_mix_matrix(std::string outputfile, int lzero){
  
  /* std::cout<<"In calc mixmatrix"<<std::endl;
  std::cout<<"Initializing Wig 3j table..." << std::endl;

  wig_table_init(2*lmax,3);
  wig_temp_init(2*lmax);
  std::cout<<"Initialised Wig 3j table...DONE" << std::endl; */

  gsl_matrix *mixmatrix, *unmixmatrix;
  mixmatrix=gsl_matrix_alloc(  lmax,lmax);
  unmixmatrix=gsl_matrix_alloc(  lmax,lmax);
  double skyfrac=double(this->mask_area)/double(this->hp_map.Npix());
  double skyfrac_weighted=0.0; //PCM
  for (int i; i< this->mask_area; i++ ){
    double mask=this->hp_mask[int(this->pixelpos[i])]; 
    skyfrac_weighted += mask*mask ; 
  } 
  skyfrac_weighted /= double(this->hp_map.Npix()); //PCM     

  //gsl_permutation *p=gsl_permutation_alloc(lmax);
  int *signum;
  signum=new int;
  std::ofstream file(&outputfile[0]);
  double Rll;
  double W[lmax];
  double wigner=0.0;
  for(int l=0; l<lmax; l++){
    W[l]=0.0;
    for(int m=0; m<=l; m++){  
      W[l]+=pow(abs(this->I[l][m]), 2.0);
      if(m==0) W[l]/=2.0;
    }
    W[l]/=( (double(l)+0.5)) ;
  }
  for(int l=lzero; l<=lzero; l++){
    std::cout<<l <<"\t" <<std::flush;
    std::cout<<"What is goung on here!"<<std::flush;
    for(int ls=0; ls<lmax;ls++){  
      Rll=0.0;
      std::vector<double> array (2 * ((l < ls) ? l : ls ) + 1); //PCM
      uclwig3j::uclwig3j(l,ls,0,array.begin());
      int lss_min = abs(l - ls); 
      int lss = lss_min;
      int lss_max = l + ls + 1; 
      int count = 0; 
      while( lss < lss_max and lss < lmax ){
         double wigner = array [lss - lss_min];
         //printf("%lf, %lf, %lf, %lf\n", W[lss], wigner, skyfrac_weighted, skyfrac );
         Rll+=( (2.0*double(lss)+1.0)*wigner*wigner*W[lss]);
         lss++ ;       
      } //PCM  
      /*for(int lss=0; lss<lmax;lss++){
         wigner=wig3jj(2*l ,2*ls,2*lss ,0, 0, 0);
         Rll+=( (2.0*double(lss)+1.0)*wigner*wigner*W[lss]);        
      }*/
      Rll*=(2.0*double(ls)+1.0)/(4.0*PI)/skyfrac_weighted;//PCM
      file<<l<<std::setw(15)<<ls<<std::setw(15)<<Rll<<std::endl;
      //gsl_matrix_set(mixmatrix, l, ls, Rll);
    }
  }
  std::cout<< "\n" <<std::flush;  
  /* gsl_linalg_LU_decomp(mixmatrix, p, signum);
     gsl_linalg_LU_invert(mixmatrix, p, unmixmatrix);*/
  /*std::ofstream file(&outputfile[0]);
  for(int l=0; l<lmax; l++){
    for(int ls=0; ls<lmax;ls++){
      file<<gsl_matrix_get(mixmatrix, l, ls)<<std::setw(15);
    }
    file<<std::endl;
    }*/
  file.close();

}


void spherical_shell::winpix( int whichCl){
 
  arr<double>winpix(4*this->nside);
  //std::string winpix_directory="/home/lwolz/Libraries/Healpix_leo/Healpix_2.20a/data/"; 
  read_pixwin(winpix_directory,this->nside, winpix);
  //for(int i=0; i<4*this->nside; i++) std::cout<<winpix[i]<<"  "<<pow(winpix[i],2.0)<<std::endl;
  for(int l=0; l<this->lmax; l++){
    if(whichCl==0) this->powspec.tt(l)/=(winpix[l]*winpix[l]);
    if(whichCl==1) this->Peebles.tt(l)/=(winpix[l]*winpix[l]);
    if(whichCl==2) this->mix_Cl.tt(l)/=(winpix[l]*winpix[l]);
  } 

}

Healpix_Map<double> spherical_shell::antipixelise(std::string file, Healpix_Map<double> newmap){

  double min=0.0, max=0.0, rand_numb=0.0, rand_density=0.0, meangal_number;
  int rand_pixel=0;
  pointing rand_pointer;
  std::ofstream outfile;
  outfile.open(&file[0]);
  int counter=0;
  planck_rng generator(rand(), rand(), rand(), rand());
  arr<double> data_for_map(hp_map.Npix()); 
  this->hp_map.minmax(min, max);
  if(min<0.0){
    meangal_number=this->gal_number/double(hp_map.Npix());
    for(int i=0; i<hp_map.Npix(); i++){
      this->hp_map[i]=meangal_number*(this->hp_map[i]+1.0);
    }
    this->hp_map.minmax(min, max);
  }

  std::cout<<"Absolute Maximum of map: "<<max<<std::endl;
  std::cout<<"Populate galaxies: "<<this->gal_number<<std::endl;
  do{
    rand_pixel=0.0, rand_density=0.0;
    rand_pixel=int(generator.rand_uni()*double(this->hp_map.Npix()) );
    rand_density=this->hp_map[rand_pixel];
    rand_numb=generator.rand_uni();
    if(rand_density/max>1.0) std::cout<<"ERROR! ANTIPIX not working!"<<std::endl;
    if(rand_density/max > rand_numb){   
      counter++; 
      rand_pointer=this->hp_map.pix2ang(rand_pixel);  
      outfile<<rand_pointer.phi<<std::setw(20)<<(rand_pointer.theta)<<std::endl;
      data_for_map[rand_pixel]+=1.0;
    }

  }while(counter<this->gal_number);
  outfile.close();
  newmap.Set(data_for_map, this->scheme);
  newmap.minmax(min, max);
  std::cout<<"Repix max: "<<max<<std::endl;

  data_for_map.dealloc();
  return newmap;
}









double temperature(double bandwidth, double HImass, double comov, double Omega_pix){

  double A12=2.869e-15;// in s^-1
  double m_Hatom=1.673e-27; //in kg
  double nu_21=speedOfLight/0.21; //in s^-1
  double Megaparsec=3.08568e22; //in m

  //convert comoving distance from Megaparsec to meter
  comov*=Megaparsec;
 
  double factor=(3.0*A12*hPlanck*pow(speedOfLight,2.0))/(32.0*pi*m_Hatom*kBoltzmann);

  double temp=factor*HImass/(pow(comov,2.0)*bandwidth*nu_21*Omega_pix);
  return temp;
}


double hubble(double z){
double hub=0.0; 
 
 double w0=-1.0, wa=0.0, omm=0.3, oml=1.0-omm, hubblenull=100.0*0.70;
 hub=hubblenull*(sqrt((omm*pow((1.0+z),3.0))+
		     (oml*exp(-3.0*wa*(z/(z+1.0)))*pow(1.0+z,3.0*(1.0+wa+w0))) ));
 
return hub;
}


double invhubble(double redshift){

return (1.0/hubble(redshift));
}

struct para {};
/*
double comoving_distance(double redshift){ //[Mpc]

  double (*R)(double)=(&invhubble);
  double temp=0.0;
  gsl_function F;
  para params={};
  F.function=R;
  F.params=&params;
  
  double *result, *abserr;
  abserr=new double;
  result=new double;
  gsl_integration_workspace * work=gsl_integration_workspace_alloc (n);
  gsl_integration_qag (&F ,0.0,redshift, EPS, EPSREL ,n, GSL_INTEG_GAUSS51, work, result,abserr);
  gsl_integration_workspace_free (work);
  temp=*result;
delete result;
delete abserr;
  //temp=(qsimp3)((*R),0.0,redshift,model);
  return (LIGHT*(temp));
  
}
*/





