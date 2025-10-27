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
    std::string output_path_rk4;
    std::string output_path_h;
    std::string output_path_e;
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
    c.output_path_rk4 = j.value("output_path_rk4", "data.txt");
    c.output_path_h = j.value("output_path_h", "data.txt");
    c.output_path_e = j.value("output_path_e", "data.txt");

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

void rk4(const Config& cfg){
    double dt = cfg.dt;
    double w = cfg.w;
    double k = cfg.k;

    std::vector<double> t = {0.0};
    std::vector<double> x = {cfg.x0};
    std::vector<double> v = {cfg.v0};

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


    std::ofstream file(cfg.output_path_rk4);

    for (size_t i = 0; i < t.size(); ++i)
        file << t[i] << " " << x[i] << " " << v[i] << "\n";

}

void heun(const Config& cfg){
    double dt = cfg.dt;
    double w = cfg.w;
    double k = cfg.k;

    std::vector<double> t = {0.0};
    std::vector<double> x = {cfg.x0};
    std::vector<double> v = {cfg.v0};

    for (int i = 0; i<=cfg.steps; ++i){
        x.push_back(x[i]+dt*f(x[i],v[i],cfg).first);
        v.push_back(v[i]+dt*f(x[i],v[i],cfg).second);
        x[i+1]=x[i]+dt*f((x[i]+x[i+1])/2,(v[i]+v[i+1])/2,cfg).first;
        v[i+1]=v[i]+dt*f((x[i]+x[i+1])/2,(v[i]+v[i+1])/2,cfg).second;
        t.push_back((i+1)*dt);
    }

    std::ofstream file(cfg.output_path_h);

    for (size_t i = 0; i < t.size(); ++i)
        file << t[i] << " " << x[i] << " " << v[i] << "\n";
}

void euler(const Config& cfg){
    double dt = cfg.dt;
    double w = cfg.w;
    double k = cfg.k;

    std::vector<double> t = {0.0};
    std::vector<double> x = {cfg.x0};
    std::vector<double> v = {cfg.v0};

    for (int i = 0; i<=cfg.steps; ++i){
        x.push_back(x[i]+dt*f(x[i],v[i],cfg).first);
        v.push_back(v[i]+dt*f(x[i],v[i],cfg).second);
        t.push_back((i+1)*dt);
    }

    std::ofstream file_e(cfg.output_path_e);

    for (size_t i = 0; i < t.size(); ++i)
        file_e << t[i] << " " << x[i] << " " << v[i] << "\n";
}

int main(int argc, char *argv[]) {
    Config cfg = load_config(argv[1]);
    rk4(cfg);
    heun(cfg);
    euler(cfg);

    return 0;
}
