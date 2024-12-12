#include "Util.h"
#include <cmath>
#include <igl/unproject_onto_mesh.h>

// Function to draw a curve using edges
void draw_curve(igl::opengl::glfw::Viewer& viewer, const MatrixXd& V) {
    for (unsigned i = 0; i < V.rows() - 1; ++i)
        viewer.data().add_edges(
            V.row(i),
            V.row(i + 1),
            Eigen::RowVector3d(1, 1, 0));
    viewer.data().add_edges(
        V.row(V.rows() - 1),
        V.row(0),
        Eigen::RowVector3d(1, 1, 0));
}

// Function to draw points from the list of points V
void draw_points(igl::opengl::glfw::Viewer& viewer, const MatrixXd& V) {
    viewer.data(0).add_points(V, Eigen::RowVector3d(0, 0.8, 0));
}

// Function to create a circle with a given radius and number of points
MatrixXd createCircle(float radius, int numberPoints) {
    float da = 2 * M_PI / numberPoints;
    int i = 0;
    MatrixXd Vc(numberPoints, 3);
    for (double angle = 0; i < numberPoints; angle += da) {
        float x = radius * cos(angle);
        float y = radius * sin(angle);
        Vc(i, 0) = x;
        Vc(i, 1) = y;
        Vc(i, 2) = 0;
        i++;
    }
    return Vc;
}


RowVector3d get_MousePositionCoord(igl::opengl::glfw::Viewer& viewer, MatrixXd& V, MatrixXi& F)
{
	int fid;
	Eigen::Vector3f bc;
	// Cast a ray in the view direction starting from the mouse position
	double x = viewer.current_mouse_x;
	double y = viewer.core().viewport(3) - viewer.current_mouse_y;
	if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
		viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
	{
		// 3d position of hit
		const RowVector3d m3 = V.row(F(fid, 0)) * bc(0) + V.row(F(fid, 1)) * bc(1) + V.row(F(fid, 2)) * bc(2);
		return m3;
	}
	
	RowVector3d non;
	non << -1, -1, -1;
	return non;

}

int get_ClosestVertex(MatrixXd& V, float x, float y){
	int closest = 0;
	for (int i = 0; i < V.rows(); i++){
		if (std::sqrt(std::pow(V.row(i)(0) - x, 2) + std::pow(V.row(i)(1) - y, 2)) < 0.3){
			return(i);
		}
	}
	return -1;
}

void createRectangleMouse(MatrixXd& Vertices, MatrixXi& Faces, float size)
{
	Vertices = MatrixXd(5, 3);
	Faces = MatrixXi(4, 3);

	Vertices << -size, -size, 0.0,
		size, -size, 0.0,
		size, size, 0.0,
		-size, size, 0.0,
		0, 0, 0;

		Faces << 0,1,4,
			     1,2,4,
			     2,3,4,
			     3,0,4;
}

void Deform(const std::vector<std::vector<float>> &weights, MatrixXd &vertices, const MatrixXd &cage) {
    for (int i = 0; i < vertices.rows(); i++) {
        RowVectorXd point(3);
        point << 0, 0, 0;
        for (int j = 0; j < weights[i].size(); j++) {
            point(0) += weights[i][j] * cage(j, 0);
            point(1) += weights[i][j] * cage(j, 1);
        }
        vertices.row(i) = point;
    }
}

MatrixXd createSphere(float radius, int numLatitude, int numLongitude) {
    int totalPoints = (numLatitude + 1) * (numLongitude + 1); // Include poles and rings
    MatrixXd Vs(totalPoints, 3);

    int index = 0;
    for (int lat = 0; lat <= numLatitude; ++lat) {
        // Latitude angle (from 0 to pi)
        double theta = M_PI * lat / numLatitude;
        double sinTheta = sin(theta);
        double cosTheta = cos(theta);

        for (int lon = 0; lon <= numLongitude; ++lon) {
            // Longitude angle (from 0 to 2*pi)
            double phi = 2 * M_PI * lon / numLongitude;
            double sinPhi = sin(phi);
            double cosPhi = cos(phi);

            // Calculate 3D position
            double x = radius * sinTheta * cosPhi;
            double y = radius * sinTheta * sinPhi;
            double z = radius * cosTheta;

            // Store in matrix
            Vs(index, 0) = x;
            Vs(index, 1) = y;
            Vs(index, 2) = z;
            ++index;
        }
    }

    return Vs;
}
