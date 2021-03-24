/*!
 * @file Mesh.cpp
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

#include "Mesh.h"
#include "DataFile.h"
#include "termcolor.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <fstream>
#include <vector>

//--------------------------------------------------//
//---------------------Vertices---------------------//
//--------------------------------------------------//
Vertex::Vertex():
  _index(-1), _coordinates(0.,0.)
{
}

Vertex::Vertex(double x, double y, int index):
  _index(index), _coordinates(x,y)
{
}

void Vertex::print() const
{
  std::cout << "(x,y,index) = (" << _coordinates(0)  << "," << _coordinates(1) << "," << _index << ")" << std::endl;
}


//-----------------------------------------------//
//---------------------Edges---------------------//
//-----------------------------------------------//
Edge::Edge():
  _index(-1), _verticesIndex(-1,-1)
{
}

Edge::Edge(int vertex1, int vertex2, int index, std::string boundaryCondition):
  _index(index), _c1(-1), _c2(-1), _boundaryCondition(boundaryCondition)
{
  if (vertex1 > vertex2)
    {
      _verticesIndex(0) = vertex2;
      _verticesIndex(1) = vertex1;
    }
  else
    {
      _verticesIndex(0) = vertex1;
      _verticesIndex(1) = vertex2;
    }
}

void Edge::print() const
{
  std::cout << "(vertex1,vertex2,index) = (" << _verticesIndex(0)  << "," << _verticesIndex(1) << "," << _index << ")" << std::endl;
}

//------------------------------------------------------------//
//---------------------Generic Cell class---------------------//
//------------------------------------------------------------//
Cell::Cell()
{
}

Cell::Cell(const Eigen::VectorXi& verticesIndex, int index):
  _index(index), _nbVertices(verticesIndex.size()), _verticesIndex(verticesIndex)
{  
}

void Cell::print() const
{
  std::cout << "Cell (index,vertices) = (" << _index;
  for (int i(0) ; i < _nbVertices ; ++i)
    {
      std::cout << "," << _verticesIndex(i);
    }
  std::cout << ")" << std::endl;
}

//----------------------------------------------//
//---------------------Mesh---------------------//
//----------------------------------------------//
Mesh::Mesh()
{
}

Mesh::Mesh(DataFile* DF):
  _DF(DF), _meshFile(_DF->getMeshFile()), _boundaryConditionReference(_DF->getBoundaryConditionReference()), _boundaryConditionType(_DF->getBoundaryConditionType())
{
}

void Mesh::Initialize(DataFile* DF)
{
  _DF = DF;
  _meshFile = _DF->getMeshFile();
  _boundaryConditionReference = _DF->getBoundaryConditionReference();
  _boundaryConditionType = _DF->getBoundaryConditionType();
  this->Initialize();
}

// Ajoute une arête
void Mesh::addEdge(const Edge& edge, int nc, std::vector<int>& headMinv, std::vector<int>& nextEdge, int& nbEdge)
{
  int vertex1 = edge.getVerticesIndex()(0);
  int vertex2 = edge.getVerticesIndex()(1);
  int ref = edge.getIndex();
  std::string BCType = edge.getBoundaryCondition();

  bool exist = false;
  // we look at the list of edges leaving from n1
  // if we find the same edge than n1->n2 we add the edge
  for (int e(headMinv[vertex1]) ; e != -1 ; e = nextEdge[e])
    {
      if (_edges[e].getVerticesIndex()(1) == vertex2)
        {
          if (nc >= 0)
            {
              _edges[e].addNeighbourCell(nc);
            }
          exist = true;
        }
    }

  // if the edge has not been found, we create it
  if (!exist)
    {
      // we initialize the edge
      _edges[nbEdge] = Edge(vertex1, vertex2, ref, BCType);
      if (nc >= 0)
        {
          _edges[nbEdge].addNeighbourCell(nc);
        }
      // we update the arrays next_edge and head_minv
      nextEdge[nbEdge] = headMinv[vertex1];
      headMinv[vertex1] = nbEdge;
      nbEdge++;
    }
}

// Calcule les centres, les aires et les périmètres des cellules
void Mesh::buildCellsCenterAndAreaAndPerimeter()
{
  _cellsCenter.resize(_numberOfCells,2);
  _cellsArea.resize(_numberOfCells);
  _cellsPerimeter.resize(_numberOfCells);
  // Boucle sur les cellules
  for (int i(0) ; i < _numberOfCells ; ++i)
    {
      // Récupère les coordonnées des sommets de la cellule
      const Eigen::VectorXi& verticesIndex(_cells[i].getVerticesIndex());
      int nbVertices(verticesIndex.size());
      
      // Calcul du centre
      for (int j(0) ; j < nbVertices ; ++j)
        {
          double x(_vertices[verticesIndex(j)].getCoordinates()[0]);
          double y(_vertices[verticesIndex(j)].getCoordinates()[1]);
          _cellsCenter(i,0) += x;
          _cellsCenter(i,1) += y;
        }
      _cellsCenter.row(i) /= nbVertices;

      // Calul du périmètre et de l'aire
      
      // // Pour tout polygone convexe
      // for (int j(0) ; j < nbVertices - 1  ; ++j)
      //   {
      //     double x1(_vertices[verticesIndex(j)].getCoordinates()[0]);
      //     double y1(_vertices[verticesIndex(j)].getCoordinates()[1]);
      //     double x2(_vertices[verticesIndex(j+1)].getCoordinates()[0]);
      //     double y2(_vertices[verticesIndex(j+1)].getCoordinates()[1]);
      //     _cellsPerimeter(i) += sqrt(pow(x2-x1,2) + pow(y2-y1,2));
      //     _cellsArea(i) += (x2+x1)*(y2-y1);
      //     // Verification de l'aiere
      //     //std::cout << i << " " << j << " " << _cellsArea(i) << std::endl;
          
          
      //   }
      // double x1(_vertices[nbVertices-1].getCoordinates()[0]);
      // double y1(_vertices[nbVertices-1].getCoordinates()[1]);
      // double x2(_vertices[0].getCoordinates()[0]);
      // double y2(_vertices[0].getCoordinates()[1]);
      // _cellsPerimeter(i) += sqrt(pow(x2-x1,2) + pow(y2-y1,2));
      // _cellsArea(i) += (x2+x1)*(y2-y1);
      // _cellsArea(i) = abs(0.5 * _cellsArea(i));
      
      // Pour des triangles //
      // Calcul de l'aire
      double x1(_vertices[verticesIndex(0)].getCoordinates()[0]);
      double y1(_vertices[verticesIndex(0)].getCoordinates()[1]);
      double x2(_vertices[verticesIndex(1)].getCoordinates()[0]);
      double y2(_vertices[verticesIndex(1)].getCoordinates()[1]);
      double x3(_vertices[verticesIndex(2)].getCoordinates()[0]);
      double y3(_vertices[verticesIndex(2)].getCoordinates()[1]);
       double l12(sqrt(pow(x1-x2,2) + pow(y1-y2,2)));
       double l13(sqrt(pow(x1-x3,2) + pow(y1-y3,2)));
       double l23(sqrt(pow(x2-x3,2) + pow(y2-y3,2)));
       double p(0.5 * (l12 + l13 + l23));
      _cellsArea(i) = sqrt(p*(p-l12)*(p-l13)*(p-l23));
      _cellsArea(i) = 0.5 * abs((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
      _cellsPerimeter(i) = 2.*p;
    }
}

// Calcule le centre, la longueur et la normale de chaque arête
void Mesh::buildEdgesNormalAndLengthAndCenter()
{
  _edgesCenter.resize(_numberOfEdges,2);
  _edgesLength.resize(_numberOfEdges);
  _edgesNormal.resize(_numberOfEdges,2);

  for (int i(0) ; i < _numberOfEdges ; ++i)
    {
      // Calcul de la longueur
      int vertex1(_edges[i].getVerticesIndex()(0));
      int vertex2(_edges[i].getVerticesIndex()(1));
      double x1(_vertices[vertex1].getCoordinates()(0));
      double y1(_vertices[vertex1].getCoordinates()(1));
      double x2(_vertices[vertex2].getCoordinates()(0));
      double y2(_vertices[vertex2].getCoordinates()(1));
      _edgesLength(i) = sqrt(pow(x1-x2,2) + pow(y1-y2,2));
      
      // Calcul du centre
      _edgesCenter(i,0) = 0.5 * (x1 + x2);
      _edgesCenter(i,1) = 0.5 * (y1 + y2);
      
      // Calcul du vecteur (centre de la cellule c1 to centre de l'arête)
      int c1(_edges[i].getC1());
      Eigen::Vector2d diff(_edgesCenter.row(i) - _cellsCenter.row(c1));
      // Calcul de la normale dans un sens arbitraire
      _edgesNormal(i,0) = y1 - y2;
      _edgesNormal(i,1) = x2 - x1;
      // Produit scalaire entre la normale et le vecteur diff
      double scalar(_edgesNormal.row(i).dot(diff));
      // Forcer la normale à sortir de T1
      if (scalar < 0.)
        {
          _edgesNormal.row(i) = -_edgesNormal.row(i);
        }
      // Normalisation
      _edgesNormal.row(i) /= _edgesLength(i);
    }
}

// Build the mesh from the mesh file
void Mesh::Initialize()
{
  std::cout << "====================================================================================================" << std::endl;

  std::ifstream meshStream(_meshFile, std::ios::in);

  // Vérifie que le fichier est bien ouvert.
  if (!meshStream.is_open())
    {
      std::cout << termcolor::red << "ERROR::MESH : Unable to open the mesh file : " << _meshFile << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
  else
    {
      std::cout << "Generating a 2D mesh from file : " << _meshFile << std::endl;
    }

  std::string line;
  int dimension(3), nBoundaryEdges(0);
  std::vector<Edge> boundaryEdges;

  // Parcours les lignes du fichier de maillage
  while (getline(meshStream, line))
    {
      // Dimension du maillage
      if (line.find("Dimension") != std::string::npos)
        {
          meshStream >> dimension;
        }
      // Création des sommets du maillage
      else if (line.find("Vertices") != std::string::npos)
        {
          meshStream >> _numberOfVertices;
          _vertices.resize(_numberOfVertices);
          for (int i(0) ; i < _numberOfVertices ; ++i)
            {
              double x(0.), y(0.), z(0.);
              int index(0);
              meshStream >> x >> y >> z >> index;
              _vertices[i] = Vertex(x, y, index);
            }
        }
      // Création des arêtes extérieures du maillage
      else if (line.find("Edges") != std::string::npos)
        {
          meshStream >> nBoundaryEdges;
          boundaryEdges.resize(nBoundaryEdges);
          for (int i(0) ; i < nBoundaryEdges ; ++i)
            {
              int vertex1(0), vertex2(0);
              int index(0);
              meshStream >> vertex1 >> vertex2 >> index;
              --vertex1; --vertex2;
              std::string BCType("none");
              for (int i(0) ; i < _boundaryConditionReference.size() ; ++i)
                {
                  if (index == _boundaryConditionReference[i])
                    {
                      BCType = _boundaryConditionType[i];
                    }
                }
              if (BCType == "none")
                {
                  std::cout << termcolor::red << "ERROR::MESH : Problem with boundary conditions in your mesh (reference or types are wrong)" << std::endl;
                  std::cout << termcolor::reset;
                  exit(-1);
                }
              boundaryEdges[i] = Edge(vertex1, vertex2, index, BCType);
            }
        }
      // Création des cellules du maillage
      // Triangles
      else if (line.find("Triangles") != std::string::npos)
        {
          _numberOfVerticesPerCell = 3;
          _cellType = "Triangles";
          meshStream >> _numberOfCells;
          _cells.resize(_numberOfCells);
          for (int i(0) ; i < _numberOfCells ; ++i)
            {
              Eigen::Vector3i verticesIndex;
              int vertex1, vertex2, vertex3;
              int index;
              meshStream >> vertex1 >> vertex2 >> vertex3 >> index;
              --vertex1; --vertex2; --vertex3;
              verticesIndex << vertex1, vertex2, vertex3;
              _cells[i] = Cell(verticesIndex, index);
            }
        }
      // Quadrilatères
      else if (line.find("Quadrilaterals") != std::string::npos)
        {
          _numberOfVerticesPerCell = 4;
          _cellType = "Quadrilaterals";
          meshStream >> _numberOfCells;
          _cells.resize(_numberOfCells);
          for (int i(0) ; i < _numberOfCells ; ++i)
            {
              Eigen::Vector4i verticesIndex;
              int vertex1, vertex2, vertex3, vertex4;
              int index;
              meshStream >> vertex1 >> vertex2 >> vertex3 >> vertex4 >> index;
              --vertex1; --vertex2; --vertex3, --vertex4;
              verticesIndex << vertex1, vertex2, vertex3, vertex4;
              _cells[i] = Cell(verticesIndex, index);
            }
        }
    }

  // Création du vecteur _edges.
  if (_cellType == "Triangles")
    {
      _numberOfEdges = (3*_numberOfCells + nBoundaryEdges)/2; 
    }
  else if (_cellType == "Quadrilaterals")
    {
      _numberOfEdges = (4*_numberOfCells + nBoundaryEdges)/2;
    }
  else
    {
      std::cout << termcolor::red << "ERROR::MESH : Cell type not supported !" << std::endl;
      std::cout << termcolor::reset << "Supported types : Triangles, Quadrilaterals" << std::endl;
      std::cout << "====================================================================================================" << std::endl << std::endl;
    }
  _edges.resize(_numberOfEdges);

  std::vector<int> headMinv(_numberOfVertices, -1);
  std::vector<int> nextEdge(_numberOfEdges, -1);

  // Ajout des arêtes extérieures
  int nbEdges(0);
  for (int i(0) ; i < nBoundaryEdges ; ++i)
    {
      addEdge(boundaryEdges[i], -1, headMinv, nextEdge, nbEdges);
    }

  // Ajout des arêtes intérieures
  for (int i(0); i < _numberOfCells; i++)
    {
      const Eigen::VectorXi& nv(_cells[i].getVerticesIndex());
      for (int j = 0; j < _numberOfVerticesPerCell; j++)
        {
          Edge edge(nv(j), nv((j+1)%_numberOfVerticesPerCell), 0, "none");
          addEdge(edge, i, headMinv, nextEdge, nbEdges);
        }
    }

  buildCellsCenterAndAreaAndPerimeter();
  buildEdgesNormalAndLengthAndCenter();
  
  std::cout << termcolor::green << "SUCCESS::MESH : Mesh generated succesfully !" << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}

// Printer (for information purposes)
void Mesh::printParameters() const
{
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Printing parameters of the mesh..." << std::endl;
  std::cout << "Mesh file           = " << _meshFile << std::endl;
  std::cout << "Number of vertices  = " << _numberOfVertices << std::endl;
  std::cout << "Number of edges     = " << _numberOfEdges << std::endl;
  std::cout << "Number of Cells     = " << _numberOfCells << std::endl;
  std::cout << "Cells type          = " << _cellType << std::endl;
  std::cout << "====================================================================================================" << std::endl << std::endl;
}
