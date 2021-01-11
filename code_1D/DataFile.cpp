#include "DataFile.h"
#include "termcolor.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <regex>

DataFile::DataFile()
{
}

DataFile::DataFile(const std::string& fileName):
  _fileName(fileName), _scenario("none")
{
}

void DataFile::Initialize(const std::string& fileName)
{
  _fileName = fileName;
  _scenario = "none";
}

std::string DataFile::cleanLine(std::string &line)
{
  std::string res = line;

  // Remove everything after a possible #
  res = regex_replace(res, std::regex("#.*$"), std::string(""));
  // Replace tabulation(s) by space(s)
  res = regex_replace(res, std::regex("\t"), std::string(" "), std::regex_constants::match_any);
  // Replace multiple spaces by 1 space
  res = regex_replace(res, std::regex("\\s+"), std::string(" "), std::regex_constants::match_any);
  // Remove any leading spaces
  res = regex_replace(res, std::regex("^ *"), std::string(""));

  return res;
}

// Lit le fichier de paramètres ligne par ligne et affecte la valeur
// adéquate à chaque paramètre. Des vérifications sont faites et la valeur
// de certain paramètres peut changer selon celles d'autres paramètres. Dans ce
// cas là, un warning est affiché sur le terminal.
void DataFile::readDataFile()
{
  // Open the data file
  std::ifstream data_file(_fileName.data());
  if (!data_file.is_open())
    {
      std::cout << termcolor::red << "ERROR::DATAFILE : Unable to open file " << _fileName << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
  else
    {
      std::cout << "====================================================================================================" << std::endl;
      std::cout << "Reading data file " << _fileName << std::endl;
    }
  // Pour stocker chaque ligne
  std::string line;
  // Run through the data_file to find the parameters
  while (getline(data_file, line))
    {
      // Clean line
      std::string proper_line(cleanLine(line));
      if (proper_line.find("TimeScheme") != std::string::npos)
        {
          data_file >> _timeScheme;
        }
      if (proper_line.find("NumericalFlux") != std::string::npos)
        {
          data_file >> _numericalFlux;
        }
      if (proper_line.find("ResultsDir") != std::string::npos)
        {
          data_file >> _resultsDir;
        }
      if (proper_line.find("xmin") != std::string::npos)
        {
          data_file >> _xmin;
        }
      if (proper_line.find("xmax") != std::string::npos)
        {
          data_file >> _xmax;
        }
      if (proper_line.find("dx") != std::string::npos)
        {
          data_file >> _dx;
        }
      if (proper_line.find("InitialTime") != std::string::npos)
        {
          data_file >> _initialTime;
        }      
      if (proper_line.find("FinalTime") != std::string::npos)
        {
          data_file >> _finalTime;
        }
      if (proper_line.find("TimeStep") != std::string::npos)
        {
          data_file >> _timeStep;
        }
      if (proper_line.find("CFL") != std::string::npos)
        {
          data_file >> _CFL;
        }
      if (proper_line.find("GravityAcceleration") != std::string::npos)
        {
          data_file >> _g;
        }
      if (proper_line.find("SaveFrequency") != std::string::npos)
        {
          data_file >> _saveFrequency;
        }
      if (proper_line.find("Scenario") != std::string::npos)
        {
          data_file >> _scenario;
        }
      if (proper_line.find("IsTopography") != std::string::npos)
        {
          data_file >> _isTopography;
        }
      if (proper_line.find("TopographyType") != std::string::npos)
        {
          data_file >> _topographyType;
        }
      if (proper_line.find("TopographyFile") != std::string::npos)
        {
          data_file >> _topographyFile;
        }
    }

  // Création et nettoyage du dossier de résultats
  std::cout << "Creating the results directory..." << std::endl;
  system(("mkdir -p ./" +_resultsDir).c_str());
  system(("rm -f ./" +_resultsDir + "/solution*").c_str());
  system(("cp -r ./" + _fileName + " ./" + _resultsDir + "/parameters.txt").c_str());

  // Logs
  std::cout << termcolor::green << "SUCCESS::DATAFILE : Results directory created successfully !" << std::endl;
  std::cout << termcolor::reset;

  // Pour le scénario LaSalie, impose la topographie
  if (_scenario == "LaSalie")
    {
      _isTopography = true;
      _topographyType = "File";
      _topographyFile = "topography_la_salie.csv";
    }

  // Si pas de topo --> impose un fond plat
  if (_isTopography == false)
    {
      _topographyType = "FlatBottom";
    }

  // Ajustement des paramètres du maillage
  // Si fichier de topographie, fait correspondre le maillage aux données
  if (_topographyType == "File")
    {
      std::ifstream topoStream(_topographyFile.data());
      if (!topoStream.is_open())
        {
          std::cout << termcolor::red << "ERROR::DATAFILE : Unable to open the topography file " << _topographyFile << std::endl;
          std::cout << termcolor::reset;
          exit(-1);
        }
      else
        {
          std::cout << "Reading topography file " << _topographyFile << std::endl;
        }
      std::cout << termcolor::magenta << "WARNING::DATAFILE : Adjusting the mesh parameters to fit with the topography." << std::endl;
      std::cout << termcolor::reset;
      // Réinitialise les paramètres du maillage
      _xmin = 0.;
      _xmax = 0.;
      _Nx = 1;
      double dummy;
      std::string properLine;
      // Lit la première ligne pour récupérer _xmin
      getline(topoStream, line);
      properLine = regex_replace(line, std::regex(","), std::string(" "));
      std::stringstream ss(properLine);
      ss >> _xmin >> dummy;
      // Lit le reste du fichier pour récupérer le nombre de cellules, _xmax et dx
      while (getline(topoStream, line))
        {
          properLine = regex_replace(line, std::regex(","), std::string(" "));
          std::stringstream ss(properLine);
          ss >> _xmax >> dummy;
          ++_Nx;
        }
      _dx = (_xmax - _xmin)/_Nx;
      // Affiche les nouvelles valeurs
      std::cout << "New values : " << std::endl;
      std::cout << "  |xmin           = " << _xmin << std::endl;
      std::cout << "  |xmax           = " << _xmax << std::endl;
      std::cout << "  |dx             = " << _dx << std::endl;
    }
  // Sinon, ajuste uniquement dx pour coller avec le domaine
  else
    {
      // Ajustement du pas d'espace pour correspondre au domaine
      std::cout << termcolor::magenta << "WARNING::DATAFILE : Adjusting dx to fit within the spatial domain." << std::endl;
      _Nx = int(ceil((_xmax - _xmin)/_dx));
      _dx = (_xmax - _xmin)/_Nx;
      // Affiche la nouvelle valeur
      std::cout << termcolor::reset << "New value : dx = " << _dx << std::endl;
    }
  // Logs de succès
  std::cout << termcolor::green << "SUCCESS::DATAFILE : File read successfully" << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}


// Affiche les paramètres sur le terminal
void DataFile::printData() const
{
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Printing parameters of " << _fileName << std::endl;
  std::cout << "Mesh              = Generated" << std::endl;
  std::cout << "  |xmin           = " << _xmin << std::endl;
  std::cout << "  |xmax           = " << _xmax << std::endl;
  std::cout << "  |dx             = " << _dx << std::endl;
  std::cout << "Time Scheme       = " << _timeScheme << std::endl;
  std::cout << "Initial time      = " << _initialTime << std::endl;
  std::cout << "Final time        = " << _finalTime << std::endl;
  std::cout << "Time step         = " << _timeStep << std::endl;
  std::cout << "Gravity           = " << _g << std::endl;
  std::cout << "Numerical Flux    = " << _numericalFlux << std::endl;
  std::cout << "Results directory = " << _resultsDir << std::endl;
  std::cout << "Save Frequency    = " << _saveFrequency << std::endl;
  std::cout << "Scenario          = " << _scenario << std::endl;
  std::cout << "Topography        = " << _topographyType << std::endl;
  if (_topographyType == "File")
    {
      std::cout << "Topography file   = " << _topographyFile << std::endl;
    }
  std::cout << "====================================================================================================" << std::endl << std::endl;
}
