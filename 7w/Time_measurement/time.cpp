#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <chrono>
#include "json.hpp"

using json = nlohmann::json;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

struct Config {
    double dt;
    double x0;
    double v0;
    double w;
    double k;
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
    c.k = j.value("k", 0.1);
    c.steps = j.value("steps", 10000);
    c.output_path = j.value("output_path", "data.txt");
    return c;
}

std::pair<double,double> f(double x, double v, const Config& cfg){
    double dt = cfg.dt;
    double w = cfg.w;
    double k = cfg.k;
    std::pair<double,double> res;
    res.first = v;
    res.second = -(w*w*x+2*k*v);
    return res;
}

int main(int argc, char *argv[]) {
    Config cfg = load_config(argv[1]);

    double dt = cfg.dt;
    double w = cfg.w;
    double k = cfg.k;

    std::vector<double> t = {0.0};
    std::vector<double> x = {cfg.x0};
    std::vector<double> v = {cfg.v0};

    auto t1 = high_resolution_clock::now();

    for (int i = 0; i < cfg.steps; ++i) {
        double k1_x = f(x[i],v[i],cfg).first;
        double k1_v = f(x[i],v[i],cfg).second;

        double k2_x = f(x[i] + dt / 2 * k1_x,v[i] + dt / 2 * k1_v,cfg).first;
        double k2_v = f(x[i] + dt / 2 * k1_x,v[i] + dt / 2 * k1_v,cfg).second;

        double k3_x = f(x[i] + dt / 2 * k2_x,v[i] + dt / 2 * k2_v,cfg).first;
        double k3_v = f(x[i] + dt / 2 * k2_x,v[i] + dt / 2 * k2_v,cfg).second;

        double k4_x = f(x[i] + dt * k3_x,v[i] + dt * k3_v,cfg).first;
        double k4_v = f(x[i] + dt * k3_x,v[i] + dt * k3_v,cfg).second;

        x.push_back(x[i] + dt / 6 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x));
        v.push_back(v[i] + dt / 6 * (k1_v + 2 * k2_v + 2 * k3_v + k4_v));
        t.push_back((i + 1) * dt);
    }

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;


    std::ofstream file(cfg.output_path);


    file << ms_double.count() <<"\n";

    return 0;
}
