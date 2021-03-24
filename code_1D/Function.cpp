#include "Function.h"
#include "DataFile.h"
#include "termcolor.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>

Function::Function()
{
}

Function::Function(DataFile* DF, Mesh* mesh):
  _DF(DF), _mesh(mesh), _xmin(mesh->getxMin()), _xmax(mesh->getxMax()), _g(_DF->getGravityAcceleration()), _nCells(mesh->getNumberOfCells()), _cellCenters(mesh->getCellCenters())
{
}

void Function::Initialize(DataFile* DF, Mesh* mesh)
{
  // Initialisation
  _DF = DF;
  _mesh = mesh;
  _xmin = mesh->getxMin();
  _xmax = mesh->getxMax();
  _g = DF->getGravityAcceleration();
  _i = 0;
  _nCells = mesh->getNumberOfCells();
  _cellCenters = mesh->getCellCenters();
  this->Initialize();
}

void Function::Initialize()
{
  // Logs de début
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Building topography and initial condition..." << std::endl;

  // Resize la condition initiale, la topographie et le terme source
  _Sol0.resize(_nCells, 2);
  _topography.resize(_nCells, 2);
  _source.resize(_nCells, 2);

  // Initialise la topographie
  _topography.col(0) = _cellCenters;
  if (_DF->getTopographyType() == "FlatBottom")
    {
      _topography.col(1).setZero();
    }
  else if (_DF->getTopographyType() == "LinearUp")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _topography(i,1) = 0.05 * (_cellCenters(i) - _xmin);
        }
    }
  else if (_DF->getTopographyType() == "LinearDown")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _topography(i,1) = 0.05 * (_xmax - _cellCenters(i));
        }
    }
  else if (_DF->getTopographyType() == "SineLinearUp")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _topography(i,1) = 0.05 * (_cellCenters(i) - _xmin) + 0.05*sin(20*M_PI*_cellCenters(i)/(_xmax - _xmin));
        }
    }
  else if (_DF->getTopographyType() == "SineLinearDown")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _topography(i,1) = 0.05 * (_xmax - _cellCenters(i)) + 0.05*sin(20*M_PI*_cellCenters(i)/(_xmax - _xmin));
        }
    }
  else if (_DF->getTopographyType() == "EllipticBump")
    {
      _topography.col(1).setZero();
      double bumpCenter(_xmin + 0.75 * (_xmax - _xmin));
      double bumpHeight(0.8);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double bumpFunction(bumpHeight - pow(_cellCenters(i) - bumpCenter,2));
          if (bumpFunction > 0)
            {
              _topography(i,1) = bumpFunction;
            }
        }
    }
  else if (_DF->getTopographyType() == "File")
    {
      const std::string topoFile(_DF->getTopographyFile());
      std::ifstream topoStream(topoFile);
      std::string line, properLine;
      double dummy;
      int i(0);
      if (!topoStream.is_open())
        {
          std::cout << termcolor::red << "ERROR::TOPOGRAPHY : Unable to open the topography file : " << topoFile << std::endl;
          std::cout << termcolor::reset << "====================================================================================================" << std::endl;
          exit(-1);
        }
      else
        {
          std::cout << "Building the topography from file : " << topoFile << std::endl;
        }
      while(getline(topoStream, line))
        {
          properLine = regex_replace(line, std::regex(",") , std::string(" "));
          std::stringstream ss(properLine);
          ss >> dummy >> _topography(i,1);
          ++i;
        }
    }
  else
    {
      std::cout << termcolor::red << "ERROR::TOPOGRAPHY : Case not implemented" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
  std::cout << termcolor::green << "SUCCESS::TOPOGRAPHY : Topography was successfully built." << std::endl;
  std::cout << termcolor::reset;

  // Initialise la condition initiale
  if (_DF->getScenario() == "ConstantWaterHeight")
    {
      _Sol0.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = 3.;
        }
    }
  else if (_DF->getScenario() == "RestingLake")
    {
      _Sol0.col(1).setZero();
      double H(3.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = std::max(H - _topography(i,1), 0.);
        }
    }
  else if (_DF->getScenario() == "DamBreak")
    {
      _Sol0.col(1).setZero();
      double Hg(1.), Hd(2.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (_cellCenters(i) < 0.5*(_xmax + _xmin))
            {
              _Sol0(i,0) = std::max(Hg - _topography(i,1), 0.);
            }
          else
            {
              _Sol0(i,0) = std::max(Hd - _topography(i,1), 0.);
            }
        }
    }
  else if (_DF->getScenario() == "SineWave")
    {
      _Sol0.col(1).setZero();
      double H(2.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = std::max(H + 0.2 * cos(M_PI * _cellCenters(i)) - _topography(i,1), 0.);
        }
    }
  else if (_DF->getScenario() == "SinePerturbation")
    {
      _Sol0.col(1).setZero();
      double H(2.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (-1 < _cellCenters(i) && _cellCenters(i) < 1)
            {
              _Sol0(i,0) = std::max(H + 0.2 * cos(M_PI * _cellCenters(i)) - _topography(i,1), 0.);
            }
          else
            {
              _Sol0(i,0) = std::max(1.8 - _topography(i,1),0.);
            }
        }
    }
    /*else if (_DF->getScenario() == "LaSalie")
    {
      _Sol0.col(1).setZero();
      double dx(_DF->getSpaceStep());
      for (int i(0) ; i < _nCells ; ++i)
      {
        _Sol0(i,0) =
      }

    }*/
  else
    {
      std::cout << termcolor::red << "ERROR::SCENARIO : Case not implemented" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
  std::cout << termcolor::green << "SUCCESS::SCENARIO : Initial Conditions was successfully built." << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;

  _i = 0; // On initialise i au cas ou
}




void Function::buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  // Construit le terme source en fonction de la topographie.
  if (_DF->getTopographyType() == "FlatBottom")
    {
      _source.setZero();
    }
  else if (_DF->getTopographyType() == "LinearUp")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _source(i,1) = - _g * Sol(i,0) * 0.05;
        }
    }
  else if (_DF->getTopographyType() == "LinearDown")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _source(i,1) = _g * Sol(i,0) * 0.05;
        }
    }
  else if (_DF->getTopographyType() == "SineLinearUp")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _source(i,1) = - _g * Sol(i,0) * (0.05 + M_PI/(_xmax - _xmin)*cos(20*M_PI*_cellCenters(i)/(_xmax - _xmin)));
        }
    }
  else if (_DF->getTopographyType() == "SineLinearDown")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _source(i,1) = - _g * Sol(i,0) * (-0.05 + M_PI/(_xmax - _xmin)*cos(20*M_PI*_cellCenters(i)/(_xmax - _xmin)));
        }
    }
  else if (_DF->getTopographyType() == "EllipticBump")
    {
      _source.setZero();
      double bumpCenter(_xmin + 0.75 * (_xmax - _xmin));
      double bumpHeight(0.8);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double bumpFunction(bumpHeight - pow(_cellCenters(i) - bumpCenter,2));
          if (bumpFunction > 0)
            {
              _source(i,1) = _g * Sol(i,0) * 2. * (_cellCenters(i) - bumpCenter);
            }
        }
    }

  // Pour un fichier de topographie, la dérivée est approchée par une formule de
  // différence finie centrée d'ordre deux à l'intérieur, et par une formule
  // décentrée d'ordre deux sur les bords
  else if (_DF->getTopographyType() == "File")
    {
      double dx(_mesh->getSpaceStep());
      _source.col(0).setZero();
      _source(0,1) = - _g * Sol(0,0) * (-_topography(2) + 4.*_topography(1) - 3.*_topography(0))/(2.*dx);
      for (int i(1) ; i < _nCells - 1 ; ++i)
        {
          _source(i,1) = - _g * Sol(i,0) * (_topography(i+1) - _topography(i-1))/(2. * dx);
        }
      _source(_nCells - 1, 1) = - _g * Sol(_nCells - 1,0) * (3.*_topography(_nCells - 1) - 4.*_topography(_nCells - 2) + _topography(_nCells - 3))/(2.*dx);
    }
  else
    {
      std::cout << termcolor::red << "ERROR::SOURCETERM : Case not implemented." << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
}




Eigen::Vector2d Function::dirichletFunction(double x, double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& donnees, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  Eigen::Vector2d h(0.,0.);
  double a, b, c, dx, dt, u1, u2, h1, h2, x1, u_xe;
  int i_max;
  dx = _DF->getDx();
  dt = _DF->getTimeStep();
  x1 = _DF->getXmin() + dx/2; // x1 est le milieu de la première maille
  i_max = donnees.size() - _i;

  if (_DF->getScenario() == "LaSalie") // Cas d'étude pratique
  {
    //std::ifstream HauteurFile("experimental_data/water_height_35.txt", ios::in);
    double temps1(donnees(_i,0)), temps2(donnees(_i+1,0)), hauteur1(donnees(_i,1)), hauteur2(donnees(_i+1,1));
    while ((not(temps1 < t <= temps2)) & (_i < i_max)) // On cherche a trouver les temps connus des capteurs tq temps1 < t <= temps2
    {
      _i++;
      temps1 = donnees(_i,0);
      temps2 = donnees(_i+1,0);
    }
    if (_i == i_max)
    {
      std::cout << "Aucun pas de temps ne correspond --> pas trop eleve ?" << std::endl;
      return h;
    }
    hauteur1 = donnees(_i,1);
    hauteur2 = donnees(_i+1,1);
    h(0) = (t - temps1)*(hauteur2 - hauteur1)/(temps2 - temps1) - hauteur1;

    h1 = Sol(0,0);
    h2 = Sol(1,0);
    u1 = Sol(0,1);
    u2 = Sol(1,1);
    a = pow(1 + dt/dx * (u2 - u1), 2);
    b = 2*dt*(u1 - x1/dx * (u2 - u1)) * (1 + dt/dx * (u2 - u1)) - dt*dt*_g*(h2 - h1)/dx;
    c = pow(dt*u1 - dt/dx * x1 * (u2 - u1), 2) - dt*dt * _g * (h1 - x1/dx * (h2 - h1));

    double racine;
    racine = FindRacine(a, b, c);
    if (racine < 0)
    {
      std::cout << "Pb dans les conditions aux bords de Neumann" << std::endl;
      h(1) = 0; // C'est n'importe quoi et on le sait
    }
    else
    {
      u_xe = u1 + (racine - x1)*(u2 - u1)/dx;
    }

    double beta_moins_xe_tn, beta_moins_0_tnplus1;
    double source_terme_xe; // IL FAUT RECUPERER LE TERME SOURCE EN XE PAR INTERPOLATION
    source_terme_xe = FindSourceX(racine);
    beta_moins_xe_tn = u_xe - 2*sqrt(_g*h(0));
    beta_moins_0_tnplus1 = beta_moins_xe_tn - _g*dt*source_terme_xe;
    h(1) = h(0)*(beta_moins_0_tnplus1 + 2*sqrt(_g*h(0)));
  }
  return h;
}




Eigen::Vector2d Function::neumannFunction(double x, double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& donnees, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  Eigen::Vector2d h(0.,0.);
  double a, b, c, dx, dt, u1, u2, h1, h2, x1, u_xe;
  int i_max;
  dx = _DF->getDx();
  dt = _DF->getTimeStep();
  x1 = _DF->getXmin() + dx/2; // x1 est le milieu de la première maille
  i_max = donnees.size() - _i;

  if (_DF->getScenario() == "LaSalie") // Cas d'étude pratique
  {
    //std::ifstream HauteurFile("experimental_data/water_height_35.txt", ios::in);
    double temps1(donnees(_i,0)), temps2(donnees(_i+1,0)), hauteur1(donnees(_i,1)), hauteur2(donnees(_i+1,1));
    while ((not(temps1 < t <= temps2)) & (_i < i_max)) // On cherche a trouver les temps connus des capteurs tq temps1 < t <= temps2
    {
      _i++;
      temps1 = donnees(_i,0);
      temps2 = donnees(_i+1,0);
    }
    if (_i == i_max)
    {
      std::cout << "Aucun pas de temps ne correspond --> pas trop eleve ?" << std::endl;
      return h;
    }
    hauteur1 = donnees(_i,1);
    hauteur2 = donnees(_i+1,1);
    h(0) = (t - temps1)*(hauteur2 - hauteur1)/(temps2 - temps1) - hauteur1;

    h1 = Sol(0,0);
    h2 = Sol(1,0);
    u1 = Sol(0,1);
    u2 = Sol(1,1);
    a = pow(1 + dt/dx * (u2 - u1), 2);
    b = 2*dt*(u1 - x1/dx * (u2 - u1)) * (1 + dt/dx * (u2 - u1)) - dt*dt*_g*(h2 - h1)/dx;
    c = pow(dt*u1 - dt/dx * x1 * (u2 - u1), 2) - dt*dt * _g * (h1 - x1/dx * (h2 - h1));

    double racine;
    racine = FindRacine(a, b, c);
    if (racine < 0)
    {
      std::cout << "Pb dans les conditions aux bords de Neumann" << std::endl;
      h(1) = 0; // C'est n'importe quoi et on le sait
    }
    else
    {
      u_xe = u1 + (racine - x1)*(u2 - u1)/dx;
    }

    double beta_moins_xe_tn, beta_moins_0_tnplus1;
    double source_terme_xe; // IL FAUT RECUPERER LE TERME SOURCE EN XE PAR INTERPOLATION
    source_terme_xe = FindSourceX(racine);
    beta_moins_xe_tn = u_xe - 2*sqrt(_g*h(0));
    beta_moins_0_tnplus1 = beta_moins_xe_tn - _g*dt*source_terme_xe;
    h(1) = h(0)*(beta_moins_0_tnplus1 + 2*sqrt(_g*h(0)));
  }
  return h;
}


double Function::FindRacine(double a, double b, double c) // Calcule les racines de ax^2 + bx + c
{
  double discriminant;
  discriminant = b*b - 4*a*c;
  if (discriminant < 0)
  {
    std::cout << "Pas de racine au polynome -> pb au niveau des CdB" << std::endl;
    return -1;
  }
  else if (discriminant == 0)
  {
    double r(-b/(2*a));

    if (r > 0) {return r;}
    else
    {
      std::cout << "Probleme au niveau des racines : on trouve x_e negatif" << std::endl;
      return -1;
    }
  }
  else
  {
    double r1, r2;
    r1 = (-b - sqrt(discriminant))/(2*a);
    r2 = (-b + sqrt(discriminant))/(2*a);
    if (r2*r1 < 0) {return r2;}
    else
    {
      std::cout << "Probleme au niveau des racines : on trouve x_e negatif" << std::endl;
      return -1;
    }
  }
}

double Function::FindSourceX(double x) // Donne le terme source en x par interpolation
{
  int i(0);
  double source1, source2, source;
  double x1(_source(i,0)); //mesure en x
  double x2;
  while ((x1 < x)||(i < _source.col(0).size()))
  {
    i++;
    x1 = _source(i,0);
  }
  x2 = _source(i+1,0);
  source1 = _source(i,1);
  source2 = _source(i+1,1);
  source = (x - x1)*(source2 - source1)/(x2 - x1) - source1;
  return source;
}
