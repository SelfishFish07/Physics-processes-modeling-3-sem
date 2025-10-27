#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "json.hpp"

using json = nlohmann::json;

struct Config {
    double dt;
    double x0;
    double v0;
    double w;
    int steps;
    std::string output_path;
};

Config load_config(const std::string &path) {
    std::ifstream f(path);
    
    json j;
    f >> j;

    Config c;
    c.dt = j.value("dt", 0.001);
    c.x0 = j.value("x0", 10.0);
    c.v0 = j.value("v0", 0.0);
    c.w = j.value("w", 10.0);
    c.steps = j.value("steps", 10000);
    c.output_path = j.value("output_path", "data.txt");

    return c;
}

int main(int argc, char *argv[]) {
    Config cfg = load_config(argv[1]);

    std::vector<double> t = {0.0};
    std::vector<double> x = {cfg.x0};
    std::vector<double> v = {cfg.v0};

    double dt = cfg.dt;
    double w = cfg.w;

    for (int i = 0; i < cfg.steps; ++i) {
        double k1_x = v[i];
        double k1_v = -w * w * std::sin(x[i]);

        double k2_x = v[i] + dt / 2 * k1_v;
        double k2_v = -w * w *std::sin (x[i] + dt / 2 * k1_x);

        double k3_x = v[i] + dt / 2 * k2_v;
        double k3_v = -w * w * std::sin(x[i] + dt / 2 * k2_x);

        double k4_x = v[i] + dt * k3_v;
        double k4_v = -w * w * std::sin(x[i] + dt * k3_x);

        x.push_back(x[i] + dt / 6 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x));
        v.push_back(v[i] + dt / 6 * (k1_v + 2 * k2_v + 2 * k3_v + k4_v));
        t.push_back((i + 1) * dt);
    }

    std::ofstream file(cfg.output_path);

    for (size_t i = 0; i < t.size(); ++i)
        file << t[i] << " " << x[i] << " " << v[i] << "\n";

    return 0;
}
