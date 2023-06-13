#include "ShallowWaterSimulation.h"

#include "DataStructure/Matrix3.h"
#include "Utils/FastNoiseLit.h"
#include "Utils/Utils.h"

namespace ShallowWater {

void visualize(std::vector<GridCell> &grid, int N, int resized) {

    Matrix3<float> values(N, N);
    for (int i = 0; i < N*N; i++)
        values[i] = grid[i].height;

    if (resized > 0) {
        values = values.resize(resized, resized, 1);
    }

    std::cout << values.displayAsPlot() << std::endl;
    std::cout << values.min() << "   " << values.max() << "\n\n\n";
    std::cout << std::endl;
}

void visualize_diff(std::vector<GridCell> initialgrid, std::vector<GridCell> grid, int N, int resized) {
    std::vector<GridCell> diff(N*N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            diff[i * N + j] = initialgrid[i * N + j] - grid[i * N + j];
        }
    }
    ShallowWater::visualize(diff, N, resized);
}
GridCell operator-(const GridCell& cell1, const GridCell& cell2) {
    GridCell res;
    res.height = cell1.height - cell2.height;
    res.water_height = cell1.water_height - cell2.water_height;
    return res;
}

Pipe::Pipe(GridCell *cell1, GridCell *cell2)
    : cell1(cell1), cell2(cell2), flow_rate(0), sediments(0) {}


void Pipe::update_flow_rate() {
    // Calculate the height difference between the two cells
    float height_difference = cell1->total_height() - cell2->total_height();

    // Only allow flow from a higher cell to a lower cell
    if (height_difference > 0) {
        flow_rate = height_difference;
    } else {
        flow_rate = 0;
    }
}

GridCell *Pipe::downstream() {
    if (flow_rate > 0) {
        return cell2;
    } else {
        return cell1;
    }
}

GridCell::GridCell() : GridCell(0) {}

GridCell::GridCell(float height)
    : height(height), water_height(0) {}

void GridCell::add_pipe(Pipe *pipe) {
    pipes.push_back(pipe);
}

void GridCell::update_water_height(float dt) {
    for (Pipe* pipe : pipes) {
        if (pipe->cell1 == this) {
            water_height -= pipe->flow_rate * dt;
        } else {
            water_height += pipe->flow_rate * dt;
        }
    }
    water_height = std::max(0.f, water_height);
}

void GridCell::evaporate(float dt) {
    float qtt = 1.0;
    water_height -= qtt * dt;
    water_height = std::max(0.f, water_height);
}

void GridCell::deposit(float dt) {
    for (Pipe* pipe : pipes) {
        if (pipe->downstream() == this) {
            float qtt = pipe->sediments * dt;
            height += qtt;
            pipe->sediments -= qtt;
        }
    }
}

void GridCell::abrasion(float dt) {
    float K = 1.f;// 0.1;

    float totalAbrasion = 0.f;
    for (Pipe* pipe : pipes) {
        if (pipe->downstream() == this) {
            float qtt = pipe->flow_rate * K * dt;
            totalAbrasion += qtt;
        }
    }
    totalAbrasion = std::min(totalAbrasion, float(height));
    height -= totalAbrasion;
    for (Pipe* pipe : pipes) {
        pipe->sediments += totalAbrasion / float(pipes.size());
    }
}

float GridCell::total_height() {
    return height + water_height;
}

void ShallowWaterSimulation::run(float dt, int num_steps) {
    for (int t = 0; t < num_steps; t++) {
        // Update the flow rates
        for (GridCell& cell : grid) {
            for (Pipe* pipe : cell.pipes) {
                pipe->update_flow_rate();
            }
        }

        // Update the water heights
        for (GridCell& cell : grid) {
            cell.update_water_height(dt);
            cell.abrasion(dt);
//            cell.deposit(dt);
            cell.evaporate(0.1);
        }
//        visualize(grid, N);
    }
}

void ShallowWaterSimulation::init(int N)
{
    this->N = N;
    // Create a grid of cells with random heights
    FastNoiseLite noise;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float height = (1.f + noise.GetNoise((float)i, (float)j)) * .5f * 100.f;//rand() % 100;  // Random height between 0 and 99
            this->grid.push_back(GridCell(height));
        }
    }

    auto initialgrid = this->grid;

    // Connect each cell to its neighbors with pipes
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            GridCell& cell = this->grid[i * N + j];
            if (i > 0) {
                // Connect to the cell above
                Pipe* pipe = new Pipe(&this->grid[i * N + j], &this->grid[(i - 1) * N + j]);
                cell.add_pipe(pipe);
                this->grid[(i - 1) * N + j].add_pipe(pipe);
                pipe->dir = Pipe::V;
            }
            if (j > 0) {
                // Connect to the cell to the left
                Pipe* pipe = new Pipe(&this->grid[i * N + j], &this->grid[i * N + (j-1)]);
                cell.add_pipe(pipe);
                this->grid[i * N + (j - 1)].add_pipe(pipe);
                pipe->dir = Pipe::H;
            }
        }
    }
}

}
