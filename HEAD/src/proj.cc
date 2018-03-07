
#include"Window.hh"
#include"Timer.hh"

#include "arr.h"
#include "healpix_map.h"
#include "alm.h"
#include "powspec.h"
#include"datatypes.h"
#include "healpix_map_fitsio.h"
#include "alm_healpix_tools.h"
#include "healpix_data_io.h"

#include "powspec_fitsio.h"

#include"string_utils.h"
#include "fitshandle.h"
#include "alm_powspec_tools.h"
#include "lsconstants.h"
#include "announce.h"
#include "paramfile.h"

#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<sstream>
#include<string>
#include<vector>
#include<iterator>
using namespace std;


using GALTYPE = float32;

Timer* timer=nullptr;


template<typename T> struct galaxies{
  arr<T> ra;
  arr<T> dec;
  arr<T> z;
};

template<typename T> class Catalog
{
  
public:
  Catalog(string fn,bool rsd=true){
    fh.open(fn);
    fh.goto_hdu(2);
    //read data arrays
    fh.read_entire_column(2,gal.ra);
    fh.read_entire_column(3,gal.dec);
    fh.read_entire_column(4,gal.z);
    if (rsd) {
      arr<T> dz;
      fh.read_entire_column(5,dz);
      for (size_t i=0;i<dz.size();i++) gal.z[i]+=dz[i];
    }
    fh.close();
  }

  //data
  galaxies<T> gal;
  fitshandle fh;
 
};


template<typename T> class Shell
{

public:

  //constructors
  Shell(const Catalog<T>& c,const Window* w):cat(c),win(w),Nw(0){};

  inline void fillIndex(const uint64& i) {index.push_back(i);}

  void projectMap(int nside){
    map.SetNside(nside,RING); 
    map.fill(0.);
    //loop on gals
#pragma omp parallel for schedule(dynamic,5000)
    //for (auto igal : index) {
    for (size_t i=0;i<index.size();i++) {
      uint32 igal=index[i];
      double theta=degr2rad*(90.-(cat.gal.dec)[igal]);
      double phi=degr2rad*(cat.gal.ra)[igal];
      pointing p(theta,phi);
      p.normalize();
      uint32 ipix=map.ang2pix(p);
      double w=win->weight(cat.gal.z[igal]);
      map[ipix]+=w;
      Nw+=w;
    }
  }
  void writeMap(){
    stringstream os;
    os << "!map_"<< win->zmin() << "-" << win->zmax() <<".fits";
    cout << "writing map " << os.str() <<endl;
    write_Healpix_map_to_fits(os.str(),map,planckType<T>());
  }

  void computeAlm(const int lmax){
    //first rescale the map
    avg=map.average();
    map.Add(-avg);
    map.Scale(1./avg);
    
    arr<double> weight(2*map.Nside());
    weight.fill(1.);

    /*
    arr<double> weight;
    read_weight_ring(string(HEALPIXDATA),map.Nside(),weight);
    for (tsize m=0; m<weight.size(); ++m) weight[m]+=1;
    */

    alm.Set(lmax,lmax);
    map2alm(map,alm,weight);

  }

  //data

  vector<uint64> index;
  const Catalog<T>& cat;
  const Window* win;
  Healpix_Map<double> map;
  double avg;
  Alm<xcomplex<double> > alm;
  double Nw; //weighted # of gals
};


auto main(int argc,char** argv)-> int {

PLANCK_DIAGNOSIS_BEGIN

  announce(argv[0]);

//decode paramfile
 paramfile params(argv[1]);

 const string filein=params.find<string>("filein","XXX");
 const string fileout=params.find<string>("fileout","!cls.fits");

 string mean=params.find<string>("mean","0");

 vector<string> words;
 tokenize(mean,',',words);
 vector<double> zmean;
 for (auto zw : words) zmean.push_back(stringToData<double>(zw));
 cout << "z shells=";
 copy(zmean.begin(),zmean.end(),ostream_iterator<double>(cout,"\t"));
 cout <<endl;

 const double width=params.find<double>("width",0.05);
 const double nsigcut=params.find<double>("n_sigma_cut",3);
 int x_depth =params.find<double>("cross_depth",1); 

 string window_type=params.find<string>("wtype","TopHat");

 vector<Window*> windows;
 for (size_t i=0;i<zmean.size();i++){
   if (window_type=="TopHat") 
     windows.push_back(new UniformWindow(zmean[i]-width,zmean[i]+width));
   else if (window_type=="Gauss")
     windows.push_back(new GaussWindow(zmean[i],width,nsigcut));
   else
     throw string("unknown window="+window_type);
 }

 const int nside=params.find<int>("nside",512); 
 const int lmax=params.find<int>("Lmax",751)-1; 
 const bool rsd=params.find<int>("include_rsd",1)==1;
 const bool remove_SN=params.find<bool>("remove_shotnoise",true);

 timer=new Timer();
 //read catalog
 Catalog<GALTYPE> cat(filein,rsd);
 cout << "Read catalog= " << argv[1] << " with " << cat.gal.z.size()/1e6 << " M entries " << *timer <<endl; 
    
 //create shellss
 vector<Shell<GALTYPE> > shells;
 for (size_t i=0;i<zmean.size();i++) 
   shells.push_back(Shell<GALTYPE>(cat,windows[i]));
    
    
 //singe loop on catalog to fill the shellss
 //#pragma omp parallel for shared(cat,shells)
 for (uint32 i=0;i<cat.gal.z.size();i++){
   double z=cat.gal.z[i];
   for (auto &s : shells){
     if (s.win->in(z)) s.fillIndex(i);
   }
 }
 cout <<"shells index filled" << *timer << endl;
    
 //shells loop
 for (auto& shell : shells){
   cout <<" -shell [" << shell.win->zmin() <<"," <<  shell.win->zmax() << "]: " << shell.index.size()/1e6 << " M galaxies " ;
   shell.projectMap(nside);
   //shell.writeMap();
   double zemin,zemax;
   shell.map.minmax(zemin,zemax);
   cout << "map avg=" << shell.map.average()<< " minmax=[" << zemin << "," << zemax <<"]";
   //attention: ci-dessous change la carte en contraste de densite:
   shell.computeAlm(lmax);
   cout << *timer << endl;
 }

 //combinatorics (from angpow)
 int nWin=shells.size();
 int jMSoff7 = 0;
 if (x_depth < 0 ) { //JEC 11/7/17 x_depth <0 do not compute auto-corr
   x_depth = -x_depth;
   jMSoff7 = 1;
 }     
 //define names
 vector<string> pair_names;
 for(int iMS=0;iMS<nWin;iMS++){
   std::stringstream ss1; ss1 << iMS;
   for(int jMS=iMS+jMSoff7; jMS<=std::min(iMS+x_depth,nWin-1); jMS++){
     std::stringstream ss2; ss2 << jMS;
     pair_names.push_back(ss1.str()+ss2.str());
   }
 }

 //computing+writing Cls
 fitshandle fout;
 fout.create(fileout);

 //header
 std::vector<fitscolumn> cols;
 cols.push_back(fitscolumn("ell", "",1,planckType<double>()));
 for (size_t i=0;i<pair_names.size();i++)
   cols.push_back(fitscolumn("cl"+pair_names[i],"",1,planckType<double>()));
 fout.insert_bintab(cols);

 //columns
 int icol=1;
 arr<int> ell(lmax+1);
 std::iota (&ell[0],&ell[0]+lmax+1,0);
 fout.write_column(icol++,ell);
 
 for(int iMS=0;iMS<nWin;iMS++){
   for(int jMS=iMS+jMSoff7; jMS<=std::min(iMS+x_depth,nWin-1); jMS++){
     cout << iMS <<"x" << jMS <<endl;
     PowSpec powspec;
     iMS==jMS ? extract_powspec(shells[iMS].alm,powspec)
     : extract_crosspowspec(shells[iMS].alm,shells[jMS].alm,powspec);

     arr<double> cl=powspec.tt();

     //shot noise for auto-spectrea
     if (iMS==jMS && remove_SN)
       for (auto &cell:cl) cell-=(4*M_PI)/shells[iMS].Nw;

     fout.write_column(icol++,cl);
     


   }
 }
 
 //keys
 

 //loop on shells
 fout.set_key("removeSN",remove_SN,"has the shot noise been subtracted?");
 for (size_t i=0;i<shells.size();i++){
   double Ntot=shells[i].index.size();
   fout.add_comment("NEW SHELL********************************************* ");
   fout.set_key("Ngal"+dataToString(i),Ntot,"number of galaxies in shell");
   fout.set_key("Nw"+dataToString(i),shells[i].Nw,"weighted number of galaxies in shell");
   const Window* win=shells[i].win;
   fout.set_key("wintype"+dataToString(i),win->type(),"window type");
   double z1=(win->zmin()+win->zmax())/2;
   fout.set_key("zmean"+dataToString(i),z1,"mean redshift window");
   fout.set_key("zmin"+dataToString(i),win->zmin(),"min redshift");
   fout.set_key("zmax"+dataToString(i),win->zmax(),"max redshift");
 }
 cout << "Writing -> " << fileout << endl;
 cout << "partial time=" << timer->partial() << " Total=" << timer->total() << endl;
 delete timer;

PLANCK_DIAGNOSIS_END
}
