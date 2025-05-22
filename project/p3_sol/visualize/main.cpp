#include <SFML/Graphics.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

using namespace std;

const int TILE_SIZE = 32; 

void drawState(sf::RenderWindow& window, const vector<string>& state) {
    window.clear();
    for (int i = 0; i < state.size(); ++i) {
        for (int j = 0; j < state[i].size(); ++j) {
            sf::RectangleShape tile(sf::Vector2f(TILE_SIZE, TILE_SIZE));
            if (state[i][j] == '#') {
                tile.setFillColor(sf::Color::Black); // 墙壁用黑色表示
            } else if (state[i][j] == 'S') {
                tile.setFillColor(sf::Color::Green); // 玩家用绿色表示
            } else if (state[i][j] == 'B') {
                tile.setFillColor(sf::Color::Blue); // 箱子用蓝色表示
            } else if (state[i][j] == 'T') {
                tile.setFillColor(sf::Color::Yellow); // 目标位置用黄色表示
            } else if (state[i][j] == 'R') {
                tile.setFillColor(sf::Color::Cyan); // 箱子在目标位置用青色表示
            } else {
                tile.setFillColor(sf::Color::White); // 空地用白色表示
            }
            tile.setPosition(j * TILE_SIZE, i * TILE_SIZE);
            window.draw(tile);
        }
    }
    window.display();
    this_thread::sleep_for(chrono::milliseconds(1000)); 
}

vector<vector<string>> readStepsFromFile(const string& filename) {
    ifstream file(filename);
    vector<vector<string>> steps;
    string line;
    vector<string> currentStep;
    bool readingMap = false;

    while (getline(file, line)) {
        if (line.empty()) {
            if (readingMap && !currentStep.empty()) {
                steps.push_back(currentStep);
                currentStep.clear();
            }
            readingMap = false;
        } else if (line[0] == '#') {
            readingMap = true;
            currentStep.push_back(line);
        }
    }
    if (!currentStep.empty()) {
        steps.push_back(currentStep);
    }

    return steps;
}

int main() {
    vector<vector<string>> steps = readStepsFromFile("big_10_detail.out");

    if (steps.empty()) {
        cerr << "Failed to read steps from file." << endl;
        return 1;
    }

    sf::RenderWindow window(sf::VideoMode(steps[0][0].size() * TILE_SIZE, steps[0].size() * TILE_SIZE), "Sokoban Animation");

    for (const auto& state : steps) {
        drawState(window, state);
    }

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
    }

    return 0;
}