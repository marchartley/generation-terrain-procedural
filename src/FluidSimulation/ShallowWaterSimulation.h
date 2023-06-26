#ifndef SHALLOWWATERSIMULATION_H
#define SHALLOWWATERSIMULATION_H

#include <vector>
#include <string>
#include <iostream>

namespace ShallowWater {

class Pipe;
class GridCell;
class ShallowWaterSimulation;

// Define the grid size and parameters
//const int N = 3; // grid size
const float dt = 0.01; // time step
const float g = 9.8; // gravity constant
const float k = 0.01; // evaporation rate
const float e = 0.1; // erosion rate
const float c = 0.5; // sediment capacity
const float d = 0.01; // diffusion rate
const float s = 0.01; // smoothing factor
class GridCell;

class Pipe {
public:
    GridCell* cell1;  // The first cell connected by the pipe
    GridCell* cell2;  // The second cell connected by the pipe
    float flow_rate; // Flow rate of water in the pipe
    float sediments;
    enum DIRECTION {H, V};
    DIRECTION dir;

    Pipe(GridCell* cell1, GridCell* cell2);

    // Update the flow rate based on the difference in water height
    void update_flow_rate();

    GridCell* downstream();
};

class GridCell {
public:
    float height;            // Height of the terrain in the cell
    float water_height;      // Height of the water in the cell
    std::vector<Pipe*> pipes; // Pipes connecting this cell to its neighbors

    GridCell();
    GridCell(float height);

    // Add a pipe connecting this cell to a neighbor
    void add_pipe(Pipe* pipe);

    // Update the water height based on the flow rates of the pipes
    void update_water_height(float dt);

    void evaporate(float dt);

    void deposit(float dt);

    void abrasion(float dt);

    // Get the total height of the cell, including the water
    float total_height();

    friend GridCell operator-(const GridCell& cell1, const GridCell& cell2);
};


class ShallowWaterSimulation {
public:
    std::vector<GridCell> grid; // The grid of cells representing the terrain

    void run(float dt, int num_steps);
    void init(int N);

    int N;
};

void visualize(std::vector<GridCell> &grid, int N, int resized = -1);
void visualize_diff(std::vector<GridCell> initialgrid, std::vector<GridCell> grid, int N, int resized = -1);

}

#endif // SHALLOWWATERSIMULATION_H
