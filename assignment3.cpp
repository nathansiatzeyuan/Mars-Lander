#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

int main() {
    // declare variables
    int x_position;
    double m, k, x, v, t_max, dt, t, a, i;
    vector<double> t_list, x_list, v_list;

    // mass, spring constant, initial position and velocity
    m = 1;
    k = 1;
    x = 0;
    v = 1;
    i = 0;

    // simulation time and timestep
    t_max = 100;
    dt = 0.1;

    // Euler integration
    for (t = 0; t <= t_max; t = t + dt) {

        // append current state to trajectories
        t_list.push_back(t);

        if (i < 2) {
            a = -k * x / m;
            x = x + dt * v;
            v = v + dt * a;
            x_list.push_back(x);
            v_list.push_back(v);
            std::cout << x << endl;
        }
        else {
            x_position = (x_list.size()) - 1;
            //std::cout << x_position << endl;
            x = x_list[x_position];
            a = -k * x / m;
            x = (2 * x) - (x_list[x_position - 1]) + (pow(dt, 2)) * a;
            v = (x - (x_list[x_position])) / dt;
            x_list.push_back(x);
            v_list.push_back(v);
            std::cout << x << endl;
        };

        i += 1;
    }

    // Define the directory path where you want to save the file
    string directoryPath = "C:/Users/Dell/OneDrive/Documents/Cambridge University/Engineering Course/First Year/Mars Lander/lander/";

    // Define the complete file path by concatenating the directory and filename
    string filePath = directoryPath + "trajectories.txt";

    // Write the trajectories to file
    ofstream fout(filePath);

    if (fout) { // file opened successfully
        for (int i = 0; i < t_list.size(); i = i + 1) {
            fout << t_list[i] << ' ' << x_list[i] << ' ' << v_list[i] << endl;
        }
    }
    else { // file did not open successfully
        std::cout << "Could not open trajectory file for writing" << endl;
    }
}

