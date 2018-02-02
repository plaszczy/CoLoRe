//--------------------------------------------------------------------------
//
// Description:
// 	class Window :
// standard z selection window
//
//
// Author List:
//	Stephane Plaszczynski (plaszczy@lal.in2p3.fr)
//
// History (add to end):
//	creation:   vendredi 2 février 2018, 10:45:19 (UTC+0100) 
//
//------------------------------------------------------------------------

#ifndef Window_hh
#define Window_hh

#include<cmath>

class Window
{

public:

  //constructors
  Window(double zmin,double zmax):_zmin(zmin),_zmax(zmax){};

  inline double zmin() const { return _zmin;}
  inline double zmax() const { return _zmax;}
  inline bool in(const double& z) const {return z>=_zmin && z<=_zmax;}

  virtual double weight(const double& z)=0;


private:
  double _zmin,_zmax;
  
};

class UniformWindow : public Window 
{
public:
  UniformWindow(double zmin,double zmax):Window(zmin,zmax){}
  inline double weight(const double& z){return 1;}
};

class GaussWindow : public Window 
{
public:
  GaussWindow(double mean,double sigma,double nsigcut=3):Window(mean-nsigcut*sigma,mean+nsigcut*sigma),_zmean(mean),_sigma(sigma){}
  inline double weight(const double& z){return std::exp(std::pow((z-_zmean)/_sigma,2)/2);}
private:
  double _zmean,_sigma;
};







#endif

