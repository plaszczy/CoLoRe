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
#include<string>
#include<iostream>

class Window
{

public:

  //constructors
  Window(double zmin,double zmax):_zmin(zmin),_zmax(zmax){};

  inline double zmin() const { return _zmin;}
  inline double zmax() const { return _zmax;}
  inline bool in(const double& z) const {return z>=_zmin && z<=_zmax;}

  virtual double weight(const double& z) const =0;
  virtual std::string type() const=0;

protected:
  double _zmin,_zmax;
  
};

class UniformWindow : public Window 
{
public:
  UniformWindow(double zmin,double zmax):Window(zmin,zmax){}
  inline double weight(const double& z) const {return 1.;}  
  std::string type() const {return "TopHat";}

};

class GaussWindow : public Window 
{
public:
  GaussWindow(double mean,double sigma,double nsigcut=3):Window(mean-nsigcut*sigma,mean+nsigcut*sigma),_zmean(mean),_sigma(sigma)
  {
    double w2=sqrt(M_PI)/2*_sigma*(std::erf((_zmax-mean)/_sigma)+std::erf((mean-_zmin)/_sigma));
    _autonorm=sqrt(w2)/(_zmax-_zmin);
    std::cout << "gaussian norm=" << _autonorm << std::endl;
  }
  inline double weight(const double& z) const 
  {
    return _autonorm*std::exp(-std::pow((z-_zmean)/_sigma,2)/2);
  }
  
std::string type() const {return "Gauss";}

private:
  double _zmean,_sigma,_autonorm;
};







#endif

