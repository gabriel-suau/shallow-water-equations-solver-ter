/*!
 * @file DataFile.h
 *
 * Handles the reading of the data file.
 *
 * @authors Gabriel Suau, Remi Pegouret, Lucas Trautmann
 *
 * @version 0.1.0
 *
 * @copyright © 2021 Gabriel Suau
 * @copyright © 2021 Remi Pegouret
 * @copyright © 2021 Lucas Trautmann
 * 
 * @copyright This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * @copyright This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * @copyright You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef DATA_FILE_H
#define DATA_FILE_H

#include "Eigen/Eigen/Dense"

#include <iostream>
#include <string>
#include <vector>

class DataFile
{
private:
  std::string _fileName;

  std::string _scenario;
  std::string _resultsDir;

  std::string _meshFile;

  std::string _numericalFlux;

  // Time parameters
  std::string _timeScheme;
  double _initialTime;
  double _finalTime;
  double _timeStep;
  double _CFL;

  double _g;

  int _saveFrequency;

  // Topography
  bool _isTopography;
  std::string _topographyType;
  std::string _topographyFile;

  // Boundary conditions
  int _nBoundaries;
  Eigen::VectorXi _boundaryConditionReference;
  std::vector<std::string> _boundaryConditionType;
  
public:
  DataFile();
  DataFile(const std::string& fileName);

  ~DataFile() = default;

  void Initialize(const std::string& fileName);

  void readDataFile();

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
  const Eigen::VectorXi& getBoundaryConditionReference() const {return _boundaryConditionReference;};
  const std::vector<std::string>& getBoundaryConditionType() const {return _boundaryConditionType;};

  // Print the parameters
  void printData() const;
};

#endif // DATA_FILE_H
