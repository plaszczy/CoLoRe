
#include"Window.hh"

#include "healpix_map.h"
#include "alm.h"
#include "powspec.h"
#include"datatypes.h"
#include "healpix_map_fitsio.h"
#include "alm_healpix_tools.h"
#include"healpix_data_io.h"
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
    fitshandle fh;
    fh.open(fn);
    fh.goto_hdu(2);
    //read data arrays
    fh.read_entire_column(2,gal.ra);
    fh.read_entire_column(3,gal.dec);
    fh.read_entire_column(4,gal.z);
    if (rsd) {
      arr<float> dz;
      fh.read_entire_column(5,dz);
      for (size_t i=0;i<dz.size();i++) gal.z[i]+=dz[i];
    }
    cout << "Read catalog= " << fn << " with " << gal.z.size()/1e6 << " M entries " << *timer <<endl; 
    fh.close();
      
  }

  galaxies<T> gal;

 
};


template<typename T> class Shell
{

public:

  //constructors
  Shell(const Catalog<T>& c,const Window* w):cat(c),win(w){};

  inline void fillIndex(const uint64& i) {index.push_back(i);}

  void projectMap(int nside){
    map.SetNside(nside,RING); 
    //loop on gals
    for (auto igal : index) {
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
    os << "map_"<< win->zmin() << "-" << win->zmax() <<".fits";
    write_Healpix_map_to_fits(os.str(),map,planckType<T>());
  }

  void computeAlm(const int lmax){
    //first rescale the map
    double avg=map.average();
    map.Add(-avg);
    map.Scale(1./avg);
    arr<double> weight;
    read_weight_ring(string(HEALPIXDATA),map.Nside(),weight);
    map2alm_iter(map,alm,0,weight);

  }



  vector<uint64> index;
  const Catalog<T>& cat;
  const Window* win;
  Healpix_Map<T> map;
  Alm<xcomplex<double> > alm;

};


auto main(int argc,char** argv)-> int {

  announce(argv[0]);
  vector<double> zmean={0.15,0.25,0.35,0.45};
  double width=0.05;
  const int nside=512;
  const int lmax=750;


  bool rsd=true;

  try{
    

    timer=new Timer();
    //read catalog
    Catalog<GALTYPE> cat(argv[1],rsd);
    
    //create shellss
    vector<Shell<GALTYPE> > shells;
    for (size_t i=0;i<zmean.size();i++) 
      shells.push_back(Shell<GALTYPE>(cat,new UniformWindow(zmean[i]-width,zmean[i]+width)));
    
    
    //singe loop on catalog to fill the shellss
    for (uint64 i=0;i<cat.gal.z.size();i++){
      double z=cat.gal.z[i];
      for (auto &s : shells){
	if (s.win->in(z)) s.fillIndex(i);
      }
    }
    cout <<"shells index filled " << *timer << endl;
    
    //shells loop
     for (auto shell : shells){
       shell.projectMap(nside);
       //s.writeMap();
       shell.computeAlm(lmax);
       cout <<"shells [" << shell.win->zmin() <<"," <<  shell.win->zmax() << "]: " << shell.index.size()/1e6 << " M galaxies " << *timer << endl;
     }

     //combine alms for auto/cross spectra
     //for (size_t i=0;i<shells.size()
	    


  }
  catch (PlanckError& e){   
    cerr << e.what()<<endl;
    return 1;
  }
  return 0;

}
