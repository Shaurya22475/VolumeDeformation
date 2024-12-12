	#include "Util.h"
	#include <ctime>
	#include <chrono>

	MatrixXd cageVertices;
	MatrixXi cageFaces;

	MatrixXd meshVertices;
	MatrixXi meshFaces;

	MatrixXd Vmouse;
	MatrixXi Fmouse;

	bool parallel = false;

	std::vector<std::vector<float>> weights;

	bool click = false;
	int clickedVertex;

	void updateMesh(const std::vector<std::vector<float>> &weights, MatrixXd &vertices, const MatrixXd &cage){
		for (int i = 0; i < vertices.rows(); i++){
			RowVectorXd point(3);
			point << 0,0,0;
			for (int j = 0; j < weights[i].size(); j++){
				point(0) += weights[i][j]*cage(j,0);
				point(1) += weights[i][j]*cage(j,1);
			}
			vertices.row(i) = point;
		}
	}

	bool mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y){
		if (click)
			{
				Vector3d Pos;
				Pos = get_MousePositionCoord(viewer, Vmouse, Fmouse);
				cageVertices.row(clickedVertex) = Pos;
				viewer.data(0).clear();
				viewer.data().clear();
				viewer.data(0).add_points(cageVertices, Eigen::RowVector3d(1, 0, 0));
				draw_curve(viewer, cageVertices);
				// Drawing mesh
				updateMesh(weights, meshVertices, cageVertices);
				viewer.data(0).add_points(meshVertices, Eigen::RowVector3d(1, 0, 0));
				draw_curve(viewer, meshVertices);
			}
			return true;
	};

	bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier){
			Vector3d Pos;
			Pos = get_MousePositionCoord(viewer, Vmouse, Fmouse);
			clickedVertex = get_ClosestVertex(cageVertices ,Pos(0), Pos(1));
			if(clickedVertex ==-1)
				return false;
			click = true;
			return false;
		};

	bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier){
		if(click){
			updateMesh(weights, meshVertices, cageVertices);
			viewer.data(0).add_points(meshVertices, Eigen::RowVector3d(1, 0, 0));
			draw_curve(viewer, meshVertices);
		}

		click = false;
		return false;
	};


int main(int argc, char *argv[])
{
	createRectangleMouse(Vmouse, Fmouse, 20.0);
	cageVertices = createCircle(5, 10);
	meshVertices = createCircle(3, 30);
	
	
	Grid G(7);
	G.setupCage(cageVertices);
	G.setupMesh(meshVertices);
	
    G.Flood_Fill(0, 0); 
    for (int x = 0; x < G.getsize(); x++) { 
        for (int y = 0; y < G.getsize(); y++) {
            if (G.getGridValue(x,y) == UNTYPED) { 
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

	draw_points(viewer, cageVertices); 
	draw_curve(viewer, cageVertices);
	draw_points(viewer, meshVertices);
	draw_curve(viewer, meshVertices);
	viewer.core(0).align_camera_center(cageVertices, cageFaces);


	viewer.callback_mouse_move = &mouse_move;
	viewer.callback_mouse_down = &mouse_down;
	viewer.callback_mouse_up = &mouse_up;

	viewer.launch(); 	
}
