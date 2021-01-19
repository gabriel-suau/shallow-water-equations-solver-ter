#ifndef DATA_FILE_H
#define DATA_FILE_H

#include <iostream>
#include <string>

class DataFile
{
private:
  // Nom du fichier de paramètres
  std::string _fileName;

  // Paramètres généraux
  std::string _scenario;
  std::string _resultsDir;

  // Mesh parameters
  std::string _meshFile;

  // Numerical Flux
  std::string _numericalFlux;

  // Time parameters
  std::string _timeScheme;
  double _initialTime;
  double _finalTime;
  double _timeStep;
  double _CFL;

  // Gravity Acceleration
  double _g;

  // Fréquence de sauvegarde de la solution
  int _saveFrequency;

  // Topography
  bool _isTopography;
  std::string _topographyType;
  std::string _topographyFile;

public:
  // Constructeurs
  DataFile();
  DataFile(const std::string& fileName);

  // Destructeur
  ~DataFile() = default;

  // Initialise l'objet
  void Initialize(const std::string& fileName);

  // Lit le fichier
  void readDataFile();

  // Nettoyer une ligne du fichier
  std::string cleanLine(std::string &line);

  // Getters
  const std::string& getFileName() const {return _fileName;};
  const std::string& getScenario() const {return _scenario;};
  const std::string& getResultsDirectory() const {return _resultsDir;};
  const std::string& getMeshFile() const {return _meshFile;};
  const std::string& getNumericalFlux() const {return _numericalFlux;};
  const std::string& getTimeScheme() const {return _timeScheme;};
  double getInitialTime() const {return _initialTime;};
  double getFinalTime() const {return _finalTime;};
  double getTimeStep() const {return _timeStep;};
  double getCFL() const {return _CFL;};
  double getGravityAcceleration() const {return _g;};
  int getSaveFrequency() const {return _saveFrequency;};
  bool isTopography() const {return _isTopography;};
  const std::string& getTopographyType() const {return _topographyType;};
  const std::string& getTopographyFile() const {return _topographyFile;};

  // Affichage des paramètres
  void printData() const;
};

#endif // DATA_FILE_H
