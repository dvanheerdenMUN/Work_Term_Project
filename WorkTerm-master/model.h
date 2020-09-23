#ifndef MODEL_H
#define MODEL_H

#include <Eigen/Core>
#include <vector>

class Model
{
public:
    Model();
    bool isBoundary;

    void AddSector(double ux, double uy, double uz,
                   double vx, double vy, double vz,
                   double sx, double sy, double txy);

    //void Evaluate();    // calculate the result

    struct Sector
    {
        Eigen::Vector2d u, v; // edges
        Eigen::Vector2d u_normalized, v_normalized; // edges
        Eigen::Vector2d u_p, v_p; // edges
        Eigen::Matrix2d stress;
        Eigen::Vector2d t0, t1;
        double angle0, angle1;
        double angle_span;
    };

    std::vector<Sector> fan;

    struct Result
    {
        Eigen::Vector2d traction[2];    // traction of the "left" and "right" side of the fracture plane
        double angle_fwd, angle_bwd;
        double phi[2];
        double theta[2];
        double t0_normal, t0_tangential;
        double t1_normal, t1_tangential;
        double trac_normal, trac_tangential;
        Eigen::Vector2d tn, tn_perp;
    };

    static const unsigned number_of_ponts = 200;
    Result results[number_of_ponts];


    Eigen::Vector2d dir;
    unsigned evaluation_time;
    double max_angle;
    unsigned fwd_sector, bwd_sector;
    unsigned idx_max_sector;
    double max_normal_traction;
    double approx_radius;

    //friend struct getTraj;

    void getMax(double min_w, double max_w);

};

//struct getTraj
//{
//        Model M;
//        double operator()(double& angle_fwd);

//};

#endif // MODEL_H
