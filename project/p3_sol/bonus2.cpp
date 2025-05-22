#include <iostream>
#include <vector>
#include <string>
#include "sokoban.hpp"

int main()
{
    std::vector<std::string> grid;
    read_map(grid);
    std::cout << solve(grid) << std::endl;
    return 0;
}