
#include"Window.hh"

#include "healpix_map.h"
#include "alm.h"
#include "powspec.h"
#include"datatypes.h"
#include "healpix_map_fitsio.h"
#include "alm_healpix_tools.h"
#include "arr.h"
#include "fitshandle.h"
#include "powspec_fitsio.h"
#include"lsconstants.h"
#include "announce.h"

#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include<string>
#include<vector>
using namespace std;


using GALTYPE = float32;

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
    cout << "OK: read catalog= " << fn << " with " << gal.z.size() << " entries" << endl;
    fh.close();
      
  }

  galaxies<T> gal;

 
};


template<typename T> class Shell
{

public:

  //constructors
  Shell(const Catalog<T>& c,const Window* w):cat(c),win(w){};

  inline void fillIndex(const size_t& i) {index.push_back(i);}

  void projectMap(int nside){
    map.SetNside(nside,RING); 
    //loop on gals
    for (auto i : index) {
      size_t igal=index[i];
      double theta=degr2rad*(90.-(cat.gal.dec)[igal]);
      double phi=degr2rad*(cat.gal.ra)[igal];
      //cout << "gal=" << igal << ": " << theta << "\t" << phi <<endl; 
      pointing p(theta,phi);
      p.normalize();
      long ipix=map.ang2pix(p);
      map[ipix]+=win->weight(cat.gal.z[igal]);
    }
  }
  void writeMap(){
    stringstream os;
    os << "mapgal_"<< win->zmin() << "-" << win->zmax() <<".fits";
    write_Healpix_map_to_fits(os.str(),map,planckType<T>());
  }

  vector<size_t> index;
  const Catalog<T>& cat;
  const Window* win;
  Healpix_Map<T> map;


};


auto main(int argc,char** argv)-> int {

  announce(argv[0]);
  vector<double> zmean={0.15,0.25,0.35,0.45};
  double width=0.05;
  const int nside=512;

  bool rsd=true;

  try{
    
    //read catalog
    Catalog<GALTYPE> cat(argv[1],rsd);
    
    //create shells
    vector<Shell<GALTYPE> > shell;
    for (size_t i=0;i<zmean.size();i++) 
      shell.push_back(Shell<GALTYPE> (cat,new UniformWindow(zmean[i]-width,zmean[i]+width)));


    //fill shell index looping once on the catalog
    for (size_t i=0;i<cat.gal.z.size();i++){
      double z=cat.gal.z[i];
      for (auto &s : shell)
	if (s.win->in(z)) s.fillIndex(i);
    }
    
    //shell loop
    for (auto s : shell){
      cout <<" shell [" << s.win->zmin() <<"," <<  s.win->zmax() << "]: " << s.index.size() << " gals" << endl;
      s.projectMap(nside);
      s.writeMap();
    }



  }
  catch (PlanckError& e){   
    cerr << e.what()<<endl;
    return 1;
  }
  return 0;

}
