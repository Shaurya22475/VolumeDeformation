#include "Util.h"

#include <ctime>
#include <chrono>


MatrixXd cageVertices;
MatrixXi cageFaces;

MatrixXd meshVertices;
MatrixXi meshFaces;

MatrixXd mouseView;
MatrixXi mousePoints;

std::vector<std::vector<float>> weights;

bool click = false;
int clickedVertex;



bool onMousemove(igl::opengl::glfw::Viewer &viewer, int mouse_x, int mouse_y) {
    if (click) {
        // Update clicked vertex position
        Vector3d Pos = get_MousePositionCoord(viewer, mouseView, mousePoints);
        cageVertices.row(clickedVertex) = Pos;

        // Clear and redraw cage vertices
        viewer.data(0).clear();
        viewer.data(0).add_points(cageVertices, Eigen::RowVector3d(1, 0, 0));
        draw_curve(viewer, cageVertices);

        // Deform the mesh based on cage movement
        Deform(weights, meshVertices, cageVertices);

        // Prepare colors based on the clicked vertex's influence
        Eigen::MatrixXd colors(meshVertices.rows(), 3);
        for (int i = 0; i < meshVertices.rows(); ++i) {
            float influence = weights[i][clickedVertex];
            colors.row(i) = Eigen::RowVector3d(std::pow((influence),0.35), std::pow((influence),0.35), std::pow((influence),0.35)); // Example: red to green gradient
            // colors.row(i) = Eigen::RowVector3d(0, 0, std::min(1.0,pow(influence,0.35))); // Example: red to green gradient
        }

        // Visualize mesh vertices with gradient colors
        viewer.data(0).add_points(meshVertices, colors);

        // Optionally redraw the mesh edges (curve)
        draw_curve(viewer, meshVertices);
    }
    return true;
}


bool onMouseDown(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {
    Vector3d Pos;
    Pos = get_MousePositionCoord(viewer, mouseView, mousePoints);
    clickedVertex = get_ClosestVertex(cageVertices, Pos(0), Pos(1));
    if (clickedVertex == -1)
        return false;
    click = true;
    return false;
}

bool onMouseup(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {
    if (click) {
        Deform(weights, meshVertices, cageVertices);
        viewer.data(0).add_points(meshVertices, Eigen::RowVector3d(1, 0, 0));
        draw_curve(viewer, meshVertices);
    }

    click = false;
    return false;
}

int main(int argc, char *argv[])
{
	createRectangleMouse(mouseView, mousePoints, 20.0);
	cageVertices = createCircle(5, 10);
	
	meshVertices = createCircle(3, 30);
	
	
	Grid G(7);
	G.setupCage(cageVertices);
	G.setupMesh(meshVertices);
	
    G.Flood_Fill(0, 0); 
    for (int x = 0; x < G.getsize(); x++) { 
        for (int y = 0; y < G.getsize(); y++) {
            if (G.getGridValue(x,y) == UNVISITED) { 
                G.setGridValue(x,y, INTERIOR);  
                G.addInternalPoint({x, y});
            }
        }
    }

	for (int i = 0; i < cageVertices.rows();i++){
		std::cout<<"processing vertex number : "<< i <<'\n';
		G.LaplacianSmooth(0.000001, i);
	}
	
	G.computeWeights();

	
	weights = G.get_weights();

	igl::opengl::glfw::Viewer viewer;
	// viewer.core().background_color = Eigen::RowVector4f(0.0f, 0.0f, 0.0f, 1.0f);


	draw_points(viewer, cageVertices); 
	draw_curve(viewer, cageVertices);
	draw_points(viewer, meshVertices);
	draw_curve(viewer, meshVertices);
	viewer.core(0).align_camera_center(cageVertices, cageFaces);


	viewer.callback_mouse_move = &onMousemove;
	viewer.callback_mouse_down = &onMouseDown;
	viewer.callback_mouse_up = &onMouseup;

	viewer.launch(); 	
}
