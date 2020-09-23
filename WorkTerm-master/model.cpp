#include "model.h"
#include "boost/math/tools/minima.hpp"
#include <iostream>
#include <cmath>
#include <cfloat>
#include <Eigen/Geometry>

using namespace std;

Model::Model()
{

}

void Model::AddSector(double ux, double uy, double uz,
               double vx, double vy, double vz,
               double sx, double sy, double txy)
{
    Sector s;
    s.u << ux, uy;
    s.v << vx, vy;
    s.stress << sx, txy, txy, sy;
    s.u_normalized = s.u.normalized();
    s.v_normalized = s.v.normalized();

    s.u_p << -s.u_normalized.y(), s.u_normalized.x();
    s.v_p << -s.v_normalized.y(), s.v_normalized.x();

    fan.push_back(s);
}

struct funcdouble
{
    Model M;
    double operator()(double const& angle_fwd)
    {
//double getTraj::operator()(double& angle_fwd)
//{
        auto get_angle = [](Eigen::Vector2d u, Eigen::Vector2d v)
        {
            double dot = u.dot(v)/(u.norm()*v.norm());
            if(dot > 1) dot = 1.0;
            else if(dot < -1.0) dot = -1.0;
            return acos(dot);
        };

        // precompute t0, t1 for each sector
        double end_angle =  0;
        for(Model::Sector &s : M.fan)
        {
            s.angle0 = end_angle;
            s.angle_span = get_angle(s.u,s.v);
            end_angle += s.angle_span;
            s.angle1 = end_angle;

            s.t0 = s.stress * s.u_p;
            s.t1 = s.stress * s.v_p;
        }

        // evaluate the function for different angles and record the maximum
        M.max_normal_traction = -DBL_MAX;
        M.max_angle = end_angle;

        M.dir = Eigen::Vector2d::Zero();

        // discretize (CCW)
            Model::Result ssr;
            ssr.tn = Eigen::Vector2d::Zero();
            ssr.traction[0] = ssr.traction[1] = Eigen::Vector2d::Zero();

            //angle_fwd = (double)end_angle/M.number_of_ponts;
            ssr.angle_fwd = angle_fwd;

            double angle_bwd = angle_fwd+end_angle/2;
            if (angle_bwd >= end_angle) angle_bwd -= end_angle;
            ssr.angle_bwd = angle_bwd;

            // integrate traction
            int sector = (M.isBoundary || angle_fwd < angle_bwd) ? 0 : 1;

            for (std::size_t f=0; f < M.fan.size(); f++)
            {
                Model::Sector &fp = M.fan[f];

                if (angle_fwd >= fp.angle0 && angle_fwd < fp.angle1)
                {
                    ssr.phi[0] = angle_fwd - fp.angle0;
                    ssr.theta[0] = fp.angle1 - angle_fwd;


                    double ratio = ssr.phi[0]/fp.angle_span;
                    ssr.tn = fp.u_normalized*(1-ratio) + fp.v_normalized*ratio;
                    ssr.tn.normalize();

                    ssr.tn_perp = fp.u_p*(1-ratio) + fp.v_p*ratio;
                    ssr.tn_perp.normalize();
                    Eigen::Vector2d tmult = fp.stress * ssr.tn_perp;

                    ssr.traction[sector] += tmult - fp.t0;
                    sector = 1-sector;
                    ssr.traction[sector] += fp.t1 - tmult;
                }
                else if (!M.isBoundary && angle_bwd >= fp.angle0 && angle_bwd < fp.angle1)
                {
                    ssr.phi[1] = angle_bwd - fp.angle0;
                    ssr.theta[1] = fp.angle1 - angle_bwd;

                    double ratio = ssr.phi[1]/fp.angle_span;
    //                Eigen::Vector2d tn = fp.u_normalized*(1-ratio) + fp.v_normalized*ratio;

                    Eigen::Vector2d tn_perp = fp.u_p*(1-ratio) + fp.v_p*ratio;
                    tn_perp.normalize();

                    Eigen::Vector2d tmult = fp.stress * tn_perp;

                    ssr.traction[sector] += tmult - fp.t0;
                    sector = 1-sector;
                    ssr.traction[sector] += fp.t1 - tmult;
                }
                else
                {
                    ssr.traction[sector] += fp.t1 - fp.t0;
                }
            }   // nFans


            ssr.t0_tangential = ssr.traction[0].dot(ssr.tn);
            ssr.t1_tangential = ssr.traction[1].dot(ssr.tn);
            ssr.t0_normal = ssr.tn_perp.dot(ssr.traction[0]);
            ssr.t1_normal = ssr.tn_perp.dot(ssr.traction[1]);

            ssr.trac_normal = ssr.t0_normal-ssr.t1_normal;
            ssr.trac_tangential = ssr.t0_tangential-ssr.t1_tangential;

            if(!M.isBoundary)
            {
                ssr.trac_normal/=2;
                ssr.trac_tangential/=2;
            }

            if(M.max_normal_traction < ssr.trac_normal) {
                M.max_normal_traction = ssr.trac_normal;
                M.dir = ssr.tn;
            }

            return -ssr.trac_normal;
        }
};

void Model::getMax(double min_w, double max_w)
{
    using boost::math::tools::brent_find_minima;

        boost::uintmax_t max_iter = 30;
        int bits = std::numeric_limits<double>::digits;

        std::pair<double, double> r = brent_find_minima(funcdouble(), min_w, max_w, bits, max_iter);

        std::cout << "x at minimum = " << r.first << ", f(" << r.first << ") = " << r.second << std::endl;
}



