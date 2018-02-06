
#include"Window.hh"

#include "healpix_map.h"
#include "alm.h"
#include "powspec.h"
#include"datatypes.h"
#include "healpix_map_fitsio.h"
#include "alm_healpix_tools.h"
#include "healpix_data_io.h"

#include "powspec_fitsio.h"

#include"string_utils.h"
#include "arr.h"
#include "fitshandle.h"
#include "alm_powspec_tools.h"
#include"lsconstants.h"
#include "announce.h"

#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<sstream>

#include"Timer.hh"
#include <thread>

#include<string>
#include<vector>
using namespace std;


using GALTYPE = float64;

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
    readcol(2,gal.ra);
    readcol(3,gal.dec);
    readcol(4,gal.z);
    if (rsd) {
      arr<T> dz;
      readcol(5,dz);
      for (size_t i=0;i<dz.size();i++) gal.z[i]+=dz[i];
    }
    fh.close();
      
  }
  void read_ra(){
    fh.read_entire_column(2,gal.ra);
  }

  void readcol(int icol,arr<T>& arr)
  {
     fh.read_entire_column(icol,arr);
  }

  galaxies<T> gal;
  fitshandle fh;
 
};


template<typename T> class Shell
{

public:

  //constructors
  Shell(const Catalog<T>& c,const Window* w):cat(c),win(w){};

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
      map[ipix]+=win->weight(cat.gal.z[igal]);
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
  Healpix_Map<T> map;
  double avg;
  Alm<xcomplex<double> > alm;

};


auto main(int argc,char** argv)-> int {

PLANCK_DIAGNOSIS_BEGIN

  announce(argv[0]);
 vector<double> zmean={0.15,0.25,0.35,0.45};
 double width=0.05;
 double nsigcut=3;
 int x_depth = 1;

 string window_type="TopHat";
 vector<Window*> windows;
 for (size_t i=0;i<zmean.size();i++){
   if (window_type=="TopHat") 
     windows.push_back(new UniformWindow(zmean[i]-width,zmean[i]+width));
   else if (window_type=="Gauss")
     windows.push_back(new GaussWindow(zmean[i],width,nsigcut));
   else
     throw string("unknown window="+window_type);
 }

 const int nside=512;
 const int lmax=750;
 const char* fileout="!cls.fits";
 bool rsd=true;

 timer=new Timer();
 //read catalog
 Catalog<GALTYPE> cat(argv[1],rsd);
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
 cols.push_back(fitscolumn("ell", "",1,planckType<int>()));
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
     if (iMS==jMS)
       extract_powspec(shells[iMS].alm,powspec);
     else
       extract_crosspowspec(shells[iMS].alm,shells[jMS].alm,powspec);
 
       fout.write_column(icol++,powspec.tt());
     


   }
 }
 
 //keys
 

 //loop on shells
 for (size_t i=0;i<shells.size();i++){
   double Ntot=shells[i].index.size();
   fout.add_comment("NEW SHELL********************************************* ");
   fout.set_key("Ngal"+dataToString(i),Ntot,"number of galaxies in shell");
   fout.set_key("Nbar"+dataToString(i),Ntot/(4*M_PI),"mean number per steradian");
   const Window* win=shells[i].win;
   fout.set_key("wintype"+dataToString(i),win->type(),"window type");
   double z1=(win->zmin()+win->zmax())/2;
   fout.set_key("zmean"+dataToString(i),z1,"mean redshift window");
   fout.set_key("zmin"+dataToString(i),win->zmin(),"min redshift");
   fout.set_key("zmax"+dataToString(i),win->zmax(),"max redshift");
 }

 cout << "partial time=" << timer->partial() << " Total=" << timer->total() << endl;
 delete timer;

PLANCK_DIAGNOSIS_END
}
