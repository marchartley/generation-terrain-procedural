layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

void main() {
    ivec3 grid_pos = ivec3(gl_WorkGroupID);
    uint index = get_grid_index(grid_pos);
    ivec3 wall = ivec3(200, 2, 200);
    cell[index].type = (grid_pos.x < wall.x && grid_pos.y < wall.y && grid_pos.z < wall.z ? SOLID : AIR);
    cell[index].vel = vec3(0);
}
