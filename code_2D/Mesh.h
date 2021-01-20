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

  // Printer
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
  // Références des triangles adjacents à l'arête
  int _t1, _t2;
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
  int getT1() const {return _t1;};
  int getT2() const {return _t2;};
  double getLength() const {return _edgeLength;};
  const Eigen::Vector2i& getVerticesIndex() const {return _verticesIndex;};
  const Eigen::Vector2d& getNormal() const {return _edgeNormal;};
  const Eigen::Vector2d& getCenter() const {return _edgeCenter;};

  // Printer
  void print() const;
};

//---------------------------------------------------//
//---------------------Triangles---------------------//
//---------------------------------------------------//
class Triangle
{
private:
  // Indice du triangle
  double _index;
  // Aire du triangle
  double _area;
  // Indices des sommets
  Eigen::Vector3i _verticesIndex;
  // Coordonnées du centre de gravité
  Eigen::Vector2d _center;

public:
  // Constructeurs
  Triangle();
  Triangle(int vertex1, int vertex2, int vertex3, int index);

  // Destructeur
  ~Triangle() = default;

  //Getters
  double getIndex() const {return _index;};
  double getArea() const {return _area;};
  const Eigen::Vector3i& getVerticesReference() const {return _verticesIndex;};
  const Eigen::Vector2d& getCenter() const {return _center;};

  // Printer
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
  // Triangles
  int _numberOfTriangles;
  std::vector<Triangle> _triangles;
  Eigen::Matrix<double, Eigen::Dynamic, 2> _trianglesCenter;
  Eigen::VectorXd _trianglesArea;
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
  const std::string& getMeshFile() const {return _meshFile;};
  int getNumberOfVertices() const {return _numberOfVertices;};
  const std::vector<Vertex>& getVertices() const {return _vertices;};
  int getNumberOfTriangles() const {return _numberOfTriangles;};
  const std::vector<Triangle>& getTriangles() const {return _triangles;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getTrianglesCenter() const {return _trianglesCenter;};
  const Eigen::VectorXd& getTrianglesArea() const {return _trianglesArea;};
  int getNumberOfEdges() const {return _numberOfEdges;};
  const std::vector<Edge>& getEdges() const {return _edges;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getEdgesCenter() const {return _edgesCenter;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getEdgesNormal() const {return _edgesNormal;};
  const Eigen::VectorXd& getEdgesLength() const {return _edgesLength;};

  // Printer
  void printParameters() const;
};

#endif // MESH_H
