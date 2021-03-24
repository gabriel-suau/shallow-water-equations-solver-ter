/*!
 * @file Mesh.h
 *
 * Defines a Mesh class (and other cool classes).
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

#ifndef MESH_H
#define MESH_H

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"
#include "DataFile.h"
#include <fstream>
#include <vector>

//--------------------------------------------------//
//---------------------Vertices---------------------//
//--------------------------------------------------//
class Vertex
{
private:
  // Indice du sommet
  int _index;
  // Coordonnées du sommet
  Eigen::Vector2d _coordinates;
  
public:
  // Constructeurs
  Vertex();
  Vertex(double x, double y, int index);

  // Destructeur
  ~Vertex() = default;
  
  // Getters
  int getIndex() const {return _index;};
  const Eigen::Vector2d& getCoordinates() const {return _coordinates;};

  // Printer (for debugging purposes)
  void print() const;
};


//-----------------------------------------------//
//---------------------Edges---------------------//
//-----------------------------------------------//
class Edge
{
private:
  // Indice de l'arête
  int _index;
  // Références des cellules adjacentes à l'arête
  int _c1, _c2;
  // Longueur de l'arête
  double _edgeLength;
  // Type de condition aux limites : Dirichlet, Neumann ou None (si arête intérieure)
  std::string _boundaryCondition;
  
  // Références des sommets de l'arête
  Eigen::Vector2i _verticesIndex;
  // Normale unitaire à l'arête
  Eigen::Vector2d _edgeNormal;
  // Coordonnées du centre de l'arête
  Eigen::Vector2d _edgeCenter;

public:
  // Constructeurs
  Edge();
  Edge(int vertex1, int vertex2, int index, std::string boundaryCondition);
  
  // Destructeur
  ~Edge() = default;

  // Getters
  int getIndex() const {return _index;};
  int getC1() const {return _c1;};
  int getC2() const {return _c2;};
  double getLength() const {return _edgeLength;};
  const Eigen::Vector2i& getVerticesIndex() const {return _verticesIndex;};
  const Eigen::Vector2d& getNormal() const {return _edgeNormal;};
  const Eigen::Vector2d& getCenter() const {return _edgeCenter;};
  const std::string& getBoundaryCondition() const {return _boundaryCondition;};

  // Add neighbour triangle
  void addNeighbourCell(int c)
  {
    if (_c1 == -1)
      {
        _c1 = c;
      }
    else
      {
        _c2 = c;
      }
  }
  
  // Printer (for debugging purposes)
  void print() const;
};


//------------------------------------------------------------//
//---------------------Generic Cell class---------------------//
//------------------------------------------------------------//
class Cell
{
private:
  // Indice de la cellule
  int _index;
    // Nombre de sommets
  int _nbVertices;
  // Aire de la cellule
  double _area;
  // Indices des sommets
  Eigen::VectorXi _verticesIndex;
  // Coordonnées du centre de gravité
  Eigen::Vector2d _center;

public:
  // Constructeurs
  Cell();
  Cell(const Eigen::VectorXi& verticesIndex, int index);

  // Destructeur
  ~Cell() = default;

  // Getters
  int getIndex() const {return _index;};
  int getNumberOfVertices() const {return _nbVertices;};
  double getArea() const {return _area;};
  const Eigen::VectorXi& getVerticesIndex() const {return _verticesIndex;};
  const Eigen::Vector2d& getCenter() const {return _center;};

  // Printer (for debugging purposes=)
  void print() const;
};

//----------------------------------------------//
//---------------------Mesh---------------------//
//----------------------------------------------//
class Mesh
{
private:
  // Pointeur vers le fichier de paramètres
  DataFile* _DF;
  // Nom du fichier de maillage
  std::string _meshFile;

  // Sommets
  int _numberOfVertices;
  std::vector<Vertex> _vertices;

  // Cells
  int _numberOfCells;
  int _numberOfVerticesPerCell;
  std::string _cellType;
  std::vector<Cell> _cells;
  Eigen::Matrix<double, Eigen::Dynamic, 2> _cellsCenter;
  Eigen::VectorXd _cellsArea;
  Eigen::VectorXd _cellsPerimeter;

  // Arêtes
  int _numberOfEdges; 
  std::vector<Edge> _edges;
  Eigen::Matrix<double, Eigen::Dynamic, 2> _edgesCenter;
  Eigen::Matrix<double, Eigen::Dynamic, 2> _edgesNormal;
  Eigen::VectorXd _edgesLength;

  // Conditions aux limites
  Eigen::VectorXi _boundaryConditionReference;
  std::vector<std::string> _boundaryConditionType;
  
public:
  // Constructeurs
  Mesh();
  Mesh(DataFile* DF);

  // Destructeur
  ~Mesh() = default;

  // Initialisation
  void Initialize(DataFile* DF);
  void Initialize();

  // Getters
  
  // Mesh File
  const std::string& getMeshFile() const {return _meshFile;};

  // Vertices
  int getNumberOfVertices() const {return _numberOfVertices;};
  const std::vector<Vertex>& getVertices() const {return _vertices;};

  // Cells
  int getNumberOfCells() const {return _numberOfCells;};
  int getNumberOfVerticesPerCell() const {return _numberOfVerticesPerCell;};
  const std::string& getCellType() const {return _cellType;};
  const std::vector<Cell>& getCells() const {return _cells;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getCellsCenter() const {return _cellsCenter;};
  const Eigen::VectorXd& getCellsArea() const {return _cellsArea;};
  const Eigen::VectorXd& getCellsPerimeter() const {return _cellsPerimeter;};

  // Edges
  int getNumberOfEdges() const {return _numberOfEdges;};
  const std::vector<Edge>& getEdges() const {return _edges;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getEdgesCenter() const {return _edgesCenter;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getEdgesNormal() const {return _edgesNormal;};
  const Eigen::VectorXd& getEdgesLength() const {return _edgesLength;};

  // Useful methods
  void buildCellsCenterAndAreaAndPerimeter();
  void buildEdgesNormalAndLengthAndCenter();
  
  // Printer (for information purposes)
  void printParameters() const;

protected:
  // Add an Edge (must not be public for obvious reasons)
  void addEdge(const Edge& edge, int nt, std::vector<int>& headMinv, std::vector<int>& nextEdge, int& nbEdge);
};

#endif // MESH_H
