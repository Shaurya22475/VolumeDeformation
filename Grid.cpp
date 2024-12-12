#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <Eigen/Core>
#include <stack>

using namespace Eigen;

#define UNTYPED 0
#define BOUNDARY 1
#define INTERIOR 2
#define EXTERIOR 3

const float xMin = -10;
const float yMin = -10;
const float xMax = 10;
const float yMax = 10;


class Grid
{
  private:
    int size;
    float step; 
    Eigen::MatrixXi G; // Grid size*size. G(x,y) represents a Region: UNTYPED, BOUNDARY, INTERIOR, EXTERIOR
    Eigen::MatrixXd cageV; // Cage vertices
    Eigen::MatrixXd meshV; // Mesh vertices

    std::vector<Eigen::MatrixXd> Harmonics;
    std::vector<std::vector<float>> weights;

    std::vector<std::pair<int, int>> internalPoints;

  public:
    Grid(const int &s) {
        size = std::pow(2, s);
        step = (xMax - xMin) / (size - 1);
        G = Eigen::MatrixXi::Zero(size, size);
    }

    int getsize(){
        return size;
    }

    // Getter for G(x, y)
    int getGridValue(int x, int y) const {
        if (x < 0 || x >= size || y < 0 || y >= size) {
            throw std::out_of_range("Grid index out of range");
        }
        return G(x, y);
    }

    // Setter for G(x, y)
    void setGridValue(int x, int y, int value) {
        if (x < 0 || x >= size || y < 0 || y >= size) {
            throw std::out_of_range("Grid index out of range");
        }
        G(x, y) = value;
    }

    void addInternalPoint(const std::pair<int, int>& point) {
        internalPoints.push_back(point);
    }

    void setupCage(const MatrixXd &cage) {
        cageV = cage;

        Harmonics.resize(cage.rows(), Eigen::MatrixXd::Zero(size, size));

        for (int i = 0; i < cage.rows(); i++) {
            int x0, y0, x1, y1;
            x0 = std::round((cage.row(i)(0) - xMin) / step);
            y0 = std::round((cage.row(i)(1) - yMin) / step);
            if (i < cage.rows() - 1) {
                x1 = std::round((cage.row(i + 1)(0) - xMin) / step);
                y1 = std::round((cage.row(i + 1)(1) - yMin) / step);
            } else {
                x1 = std::round((cage.row(0)(0) - xMin) / step);
                y1 = std::round((cage.row(0)(1) - yMin) / step);
            }
            RasterizeLine(x0, y0, x1, y1, i);
        }
    }

    void setupMesh(const MatrixXd &mesh) {
        meshV = mesh;
    }

    void computeWeights() {
        const int numVertices = meshV.rows();
        const int numHarmonics = Harmonics.size();

        weights.reserve(numVertices);


        auto convertToGrid = [this](float coord, float minCoord) {
            return std::round((coord - minCoord) / step);
        };

        for (int i = 0; i < numVertices; ++i) {
            std::vector<float> temp;
            temp.reserve(numHarmonics);

            const float x = meshV.row(i)(0);
            const float y = meshV.row(i)(1);
            const int gridX = convertToGrid(x, xMin);
            const int gridY = convertToGrid(y, yMin);

            for (int j = 0; j < numHarmonics; ++j) {
                temp.emplace_back(Harmonics[j](gridX, gridY));
            }

            weights.emplace_back(std::move(temp));
        }
    }


    void RasterizeLine(int x0, int y0, int x1, int y1, int cageIdx) {

        Harmonics[cageIdx](x0, y0) = 1;

        int dx = abs(x1 - x0);
        int sx = x0 < x1 ? 1 : -1;
        int dy = abs(y1 - y0);
        int sy = y0 < y1 ? 1 : -1;
        int err = (dx > dy ? dx : -dy) / 2;
        int e2;
        int pointno = 0;


        std::vector<std::pair<int, int>> pts;

        while (1) {
            G(x0, y0) = BOUNDARY;
            if (x0 == x1 && y0 == y1) break;
            e2 = err;
            if (e2 > -dx) {
                err -= dy;
                x0 += sx;
            }
            if (e2 < dy) {
                err += dx;
                y0 += sy;
            }

            pts.push_back(std::pair<int,int>(x0, y0));

        }
        float t = 1.0/pts.size();
        float sigma = pts.size() / 4.0; // Standard deviation

        // Gaussian Interpolation.
        for (int j = 0; j < pts.size(); j++) {
            float distance = j * t;
            float weight = exp(-distance * distance / (2 * sigma * sigma));

            Harmonics[cageIdx](pts[j].first, pts[j].second) += weight;

            if (cageIdx < Harmonics.size() - 1) {
                Harmonics[cageIdx + 1](pts[j].first, pts[j].second) += (1 - weight);
            } else {
                Harmonics[0](pts[j].first, pts[j].second) += (1 - weight);
            }
        }
    }


    

    void Fill_Grid_Regions() {
        Flood_Fill(0, 0);
        for (int x = 0; x < size; x++) {
            for (int y = 0; y < size; y++) {
                if (G(x, y) == UNTYPED) {
                    G(x, y) = INTERIOR;
                    internalPoints.push_back({x,y});
                }
            }
        }
    }

    void updateMesh(const MatrixXd &cage){
		for (int i = 0; i < meshV.rows(); i++){
			RowVectorXd point(3);
			point << 0,0,0;
			for (int j = 0; j < cage.rows(); j++){
				point(0) += weights[i][j]*cage(j,0);
				point(1) += weights[i][j]*cage(j,1);
			}
			meshV.row(i) = point;
		}
	}

    void LaplacianSmooth(float tolerance = 0.00001, int harmonicIndex = 0) {
        float maxDifference = 1.0f; 
        Eigen::MatrixXf updatedHarmonic = Harmonics[harmonicIndex].cast<float>(); 

        while (maxDifference > tolerance) {
            maxDifference = 0.0f;

            for (const auto& point : internalPoints) {
                int gridX = point.first;
                int gridY = point.second;

                float neighborAverage = 0.0f;
                int count = 0;

                if (gridX + 1 < Harmonics[harmonicIndex].rows()) {
                    neighborAverage += Harmonics[harmonicIndex](gridX + 1, gridY);
                    count++;
                }
                if (gridX - 1 >= 0) {
                    neighborAverage += Harmonics[harmonicIndex](gridX - 1, gridY);
                    count++;
                }
                if (gridY + 1 < Harmonics[harmonicIndex].cols()) {
                    neighborAverage += Harmonics[harmonicIndex](gridX, gridY + 1);
                    count++;
                }
                if (gridY - 1 >= 0) {
                    neighborAverage += Harmonics[harmonicIndex](gridX, gridY - 1);
                    count++;
                }

                if (count > 0) {
                    neighborAverage /= count;
                }

                
                float difference = std::abs(neighborAverage - Harmonics[harmonicIndex](gridX, gridY));
                maxDifference = std::max(maxDifference, difference);

            
                updatedHarmonic(gridX, gridY) = neighborAverage;
            }

        
            Harmonics[harmonicIndex] = updatedHarmonic.cast<double>();
        }
    }

    std::vector<std::vector<float>> get_weights(){
		return weights;
	}

    void Flood_Fill(int startX, int startY) {
        std::stack<std::pair<int, int>> stack;
        stack.push({startX, startY});
        
        while (!stack.empty()) {
            int x = stack.top().first;
            int y = stack.top().second;
            stack.pop();

    
            if (x < 0 || x >= size || y < 0 || y >= size || G(x, y) != UNTYPED) {
                continue;
            }
            G(x, y) = EXTERIOR;

    
            stack.push({x + 1, y});
            stack.push({x - 1, y});
            stack.push({x, y + 1});
            stack.push({x, y - 1});
        }
    }

    

    void Print_Grid() {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                std::cout << G(i, j) << "  ";
            }
            std::cout << "\n";
        }
    }
};
