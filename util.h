#ifndef _UTIL_H_
#define _UTIL_H_
#endif

#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/readPLY.h>
#include "Grid.cpp"
#include <igl/unproject_onto_mesh.h>

using namespace std;

void draw_curve(igl::opengl::glfw::Viewer& viewer, const MatrixXd& V);

void draw_points(igl::opengl::glfw::Viewer& viewer, const MatrixXd& V);

MatrixXd createCircle(float radius, int numberPoints);




RowVector3d get_MousePositionCoord(igl::opengl::glfw::Viewer& viewer, MatrixXd& V, MatrixXi& F);
int get_ClosestVertex(MatrixXd& V, float x, float y);


void createRectangleMouse(MatrixXd& Vertices, MatrixXi& Faces, float size);