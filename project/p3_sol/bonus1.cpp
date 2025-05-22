#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem> // C++17

using namespace std;
namespace fs = std::filesystem;

void writeMapStateToFile(const vector<vector<char>>& grid, ofstream& outFile, int step = -1, char move = '\0') {
    if (step == -1) {
        outFile << "Initial state:\n";
    } else {
        string direction;
        switch (move) {
            case 'U': direction = "Up"; break;
            case 'D': direction = "Down"; break;
            case 'L': direction = "Left"; break;
            case 'R': direction = "Right"; break;
        }
        outFile << "Step " << step << " (Move " << move << " - " << direction << "):\n";
    }
    
    for (const auto& row : grid) {
        for (char cell : row) {
            outFile << cell;
        }
        outFile << '\n';
    }
    outFile << '\n';
}

void simulateGameSteps(const vector<string>& inputMap, const string& moves, const string& outputFile) {
    vector<vector<char>> grid;
    grid.reserve(inputMap.size());
    for (const auto& mapRow : inputMap) {
        vector<char> currentRow(mapRow.begin(), mapRow.end());
        grid.push_back(move(currentRow));
    }

    size_t playerRow = 0, playerCol = 0;
    for (size_t r = 0; r < grid.size(); ++r) {
        for (size_t c = 0; c < grid[0].size(); ++c) {
            if (grid[r][c] == 'S') {
                playerRow = r;
                playerCol = c;
                break;
            }
        }
    }
    fs::create_directories(fs::path(outputFile).parent_path());

    ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        cerr << "无法打开文件: " << outputFile << endl;
        return;
    }

    writeMapStateToFile(grid, outFile);

    for (size_t i = 0; i < moves.length(); ++i) {
        char move = moves[i];
        size_t newRow = playerRow, newCol = playerCol;
        
        switch (move) {
            case 'U': newRow--; break;
            case 'D': newRow++; break;
            case 'L': newCol--; break;
            case 'R': newCol++; break;
        }
        if (grid[newRow][newCol] == 'B' || grid[newRow][newCol] == 'R') {
            size_t boxNewRow = newRow + (newRow - playerRow);
            size_t boxNewCol = newCol + (newCol - playerCol);
            if (grid[boxNewRow][boxNewCol] == 'T') {
                grid[boxNewRow][boxNewCol] = 'B';
            } else {
                grid[boxNewRow][boxNewCol] = 'B';
            }
        }

        // 更新玩家位置
        char originalCell = grid[playerRow][playerCol];
        grid[playerRow][playerCol] = (originalCell == 'R' ? 'T' : '.');
        grid[newRow][newCol] = 'S';
        
        playerRow = newRow;
        playerCol = newCol;

        // 写入当前状态
        writeMapStateToFile(grid, outFile, i + 1, move);
    }

    outFile.close();
}

int main() {
    vector<string> inputMap = {
        "########",
        "#.RT..T#",
        "#.B.BT.#",
        "#T.##B.#",
        "#...BS.#",
        "#..#..B#",
        "#TT....#",
        "########"
    };

    string moves = "LRDDRUUUULDDDRDLLLRUULLULUURLDDDDUURURR";
    string outputFile = "detail/big_10_detail.out";

    simulateGameSteps(inputMap, moves, outputFile);
    return 0;
}