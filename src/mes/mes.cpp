#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <limits>
#include <cmath>

#include "mes.h"
#include "../solver/cholesky_ldl.h"

namespace Fem {
    // --------------Triangles----------------
    constexpr double bc_pc_xi_tr[6] = {-0.57735f, 0.57735f, 0.57735f, -0.57735f, -1.0f, -1.0f};
    constexpr double bc_pc_eta_tr[6] = {-1.0f, -1.0f, -0.57735f, 0.57735f, 0.57735, -0.57735};
    constexpr double bc_pc_w_tr[6] = {1, 1, 1, 1, 1, 1};

    constexpr double pc_xi_tr[3] = {-2.0f/3.0f, -2.0f/3.0f, 1.0f/3.0f};
    constexpr double pc_eta_tr[3] = {-2.0f/3.0f, 1.0f/3.0f, -2.0f/3.0f};
    constexpr double pc_w_tr[3] = {2.0f/3.0f, 2.0f/3.0f, 2.0f/3.0f};

    double N1_tr(double xi, double eta) {
        return -(xi+eta)*0.5f;
    }
    double N2_tr(double xi, double eta) {
        return 0.5f*(1+xi);
    }
    double N3_tr(double xi, double eta) {
        return 0.5f*(1+eta);
    }

    double Triangle::dN1dx_tr(std::vector<Node> &nodes) {
        double x1, x2, x3, y1, y2, y3;

        x1 = nodes[this->node_ids[0]].x;
        y1 = nodes[this->node_ids[0]].y;
        x2 = nodes[this->node_ids[1]].x;
        y2 = nodes[this->node_ids[1]].y;
        x3 = nodes[this->node_ids[2]].x;
        y3 = nodes[this->node_ids[2]].y;

        const double half_area = ((x2-x1)*(y3-y1))-((x3-x1)*(y2-y1));

        return ((y2-y3)/half_area);
    }

    double Triangle::dN2dx_tr(std::vector<Node> &nodes) {
        double x1, x2, x3, y1, y2, y3;

        x1 = nodes[this->node_ids[0]].x;
        y1 = nodes[this->node_ids[0]].y;
        x2 = nodes[this->node_ids[1]].x;
        y2 = nodes[this->node_ids[1]].y;
        x3 = nodes[this->node_ids[2]].x;
        y3 = nodes[this->node_ids[2]].y;

        const double half_area = ((x2-x1)*(y3-y1))-((x3-x1)*(y2-y1));

        return ((y3-y1)/half_area);
    }

    double Triangle::dN3dx_tr(std::vector<Node> &nodes) {
        double x1, x2, x3, y1, y2, y3;

        x1 = nodes[this->node_ids[0]].x;
        y1 = nodes[this->node_ids[0]].y;
        x2 = nodes[this->node_ids[1]].x;
        y2 = nodes[this->node_ids[1]].y;
        x3 = nodes[this->node_ids[2]].x;
        y3 = nodes[this->node_ids[2]].y;

        return ((y1-y2)/(((x2-x1)*(y3-y1))-((x3-x1)*(y2-y1))));
    }

    double Triangle::dN1dy_tr(std::vector<Node> &nodes) {
        double x1, x2, x3, y1, y2, y3;

        x1 = nodes[this->node_ids[0]].x;
        y1 = nodes[this->node_ids[0]].y;
        x2 = nodes[this->node_ids[1]].x;
        y2 = nodes[this->node_ids[1]].y;
        x3 = nodes[this->node_ids[2]].x;
        y3 = nodes[this->node_ids[2]].y;

        return ((x3-x2)/(((x2-x1)*(y3-y1))-((x3-x1)*(y2-y1))));
    }

    double Triangle::dN2dy_tr(std::vector<Node> &nodes) {
        double x1, x2, x3, y1, y2, y3;

        x1 = nodes[this->node_ids[0]].x;
        y1 = nodes[this->node_ids[0]].y;
        x2 = nodes[this->node_ids[1]].x;
        y2 = nodes[this->node_ids[1]].y;
        x3 = nodes[this->node_ids[2]].x;
        y3 = nodes[this->node_ids[2]].y;

        return ((x1-x3)/(((x2-x1)*(y3-y1))-((x3-x1)*(y2-y1))));
    }

    double Triangle::dN3dy_tr(std::vector<Node> &nodes) {
        double x1, x2, x3, y1, y2, y3;

        x1 = nodes[this->node_ids[0]].x;
        y1 = nodes[this->node_ids[0]].y;
        x2 = nodes[this->node_ids[1]].x;
        y2 = nodes[this->node_ids[1]].y;
        x3 = nodes[this->node_ids[2]].x;
        y3 = nodes[this->node_ids[2]].y;

        return ((x2-x1)/(((x2-x1)*(y3-y1))-((x3-x1)*(y2-y1))));
    }

    //---------------Quads-----------------

    std::vector<double> pc_xi_q = {-1.0/std::sqrt(3), 1.0/std::sqrt(3), 1.0/std::sqrt(3), -1.0/std::sqrt(3)};
    std::vector<double> pc_eta_q = {-1.0/std::sqrt(3), -1.0/std::sqrt(3), 1.0/std::sqrt(3), 1.0/std::sqrt(3)};
    std::vector<double> pc_weights_q = {1,1,1,1};

    std::vector<double> bc_xi_q = {-1.0/std::sqrt(3), 1.0/std::sqrt(3), 1,1, 1.0/std::sqrt(3), -1.0/std::sqrt(3), -1, -1};
    std::vector<double> bc_eta_q = {-1,-1, -1.0/std::sqrt(3), 1.0/std::sqrt(3), 1, 1, 1.0/std::sqrt(3), -1.0/std::sqrt(3)};
    std::vector<double> bc_weights_q = {1,1,1,1,1,1,1,1};

    double dN1dxi_q(double eta) { return -0.25 * (1 - eta); }
    double dN2dxi_q(double eta) { return 0.25 * (1 - eta); }
    double dN3dxi_q(double eta) { return 0.25 * (1 + eta); }
    double dN4dxi_q(double eta) { return -0.25 * (1 + eta); }

    double dN1deta_q(double xi) { return -0.25 * (1 - xi); }
    double dN2deta_q(double xi) { return -0.25 * (1 + xi); }
    double dN3deta_q(double xi) { return 0.25 * (1 + xi); }
    double dN4deta_q(double xi) { return 0.25 * (1 - xi); }

    double N1_q(double xi, double eta) { return 0.25 * ((1 - xi) * (1 - eta)); }
    double N2_q(double xi, double eta) { return 0.25 * ((1 + xi) * (1 - eta)); }
    double N3_q(double xi, double eta) { return 0.25 * ((1 + xi) * (1 + eta)); }
    double N4_q(double xi, double eta) { return 0.25 * ((1 - xi) * (1 + eta)); }

    Matrix Quad::jacobian_mat(Quad &element, std::vector<Node> &nodes, double pc_xi, double pc_eta)
    {
        Matrix jacobian = Matrix::Zero(2,2);

        jacobian(0,0) = dN1dxi_q(pc_eta)*nodes[element.node_ids[0]].x+Fem::dN2dxi_q(pc_eta)*nodes[element.node_ids[1]].x+Fem::dN3dxi_q(pc_eta)*nodes[element.node_ids[2]].x+Fem::dN4dxi_q(pc_eta)*nodes[element.node_ids[3]].x;
        jacobian(0,1) = dN1dxi_q(pc_eta)*nodes[element.node_ids[0]].y+Fem::dN2dxi_q(pc_eta)*nodes[element.node_ids[1]].y+Fem::dN3dxi_q(pc_eta)*nodes[element.node_ids[2]].y+Fem::dN4dxi_q(pc_eta)*nodes[element.node_ids[3]].y;
        jacobian(1,0) = dN1deta_q(pc_xi)*nodes[element.node_ids[0]].x+Fem::dN2deta_q(pc_xi)*nodes[element.node_ids[1]].x+Fem::dN3deta_q(pc_xi)*nodes[element.node_ids[2]].x+Fem::dN4deta_q(pc_xi)*nodes[element.node_ids[3]].x;
        jacobian(1,1) = dN1deta_q(pc_xi)*nodes[element.node_ids[0]].y+Fem::dN2deta_q(pc_xi)*nodes[element.node_ids[1]].y+Fem::dN3deta_q(pc_xi)*nodes[element.node_ids[2]].y+Fem::dN4deta_q(pc_xi)*nodes[element.node_ids[3]].y;

        return jacobian;
    }

    double Quad::det_jacobian(Fem::Matrix jacobian_mat)
    {
        return jacobian_mat(0,0)*jacobian_mat(1,1) - jacobian_mat(1,0)*jacobian_mat(0,1);
    }

    Matrix Quad::inv_jacobian_mat(Fem::Matrix jacobian_mat)
    {
        Fem::Matrix inv_jacobian = Matrix::Zero(2,2);
        double det_J = this->det_jacobian(jacobian_mat);
        inv_jacobian(0,0) = jacobian_mat(1,1)/det_J;
        inv_jacobian(0,1) = -jacobian_mat(0,1)/det_J;
        inv_jacobian(1,0) = -jacobian_mat(1,0)/det_J;
        inv_jacobian(1,1) = jacobian_mat(0,0)/det_J;

        return inv_jacobian;
    }

    std::vector<Node> load_nodes(const std::string& file_name)
    {
        std::fstream file;
        std::vector<Node> result;

        file.open(file_name);
        if(!file.good()){
            std::cout << "Couldn't load text file\n";
            return result;
        }

        std::string line;
        bool node_selection = false;

        while (std::getline(file, line)) {
            if (line == "*Nodes") {
                node_selection = true;
                continue;
            }

            if (line == "*Triangles" || line == "*BC" || line == "*Quads") {
                node_selection = false;
            }

            if (node_selection) {
                int id;
                double x, y;
                char comma;
                std::istringstream iss(line);
                if (iss >> id >> comma >> x >> comma >> y) {
                    result.emplace_back(x, y);
                }
            }
        }

        file.close();
        return result;
    }

    std::vector<Triangle> load_triangles(const std::string& file_name)
    {
        std::fstream file;
        std::vector<Triangle> result;

        file.open(file_name);
        if(!file.good()){
            std::cout << "Couldn't load text file\n";
            return result;
        }

        std::string line;
        bool element_selection = false;

        while (std::getline(file, line)) {
            if (line == "*Triangles") {
                element_selection = true;
                continue;
            }

            if (line == "*Nodes" || line == "*BC" || line == "*Quads") {
                element_selection = false;
            }

            if (element_selection) {
                int id;
                int node_id1, node_id2, node_id3;
                char comma;
                std::istringstream iss(line);
                if (iss >> id >> comma >> node_id1 >> comma >> node_id2 >> comma >> node_id3) {
                    result.emplace_back(node_id1, node_id2, node_id3);
                }
            }
        }

        file.close();
        return result;
    }

    std::vector<Quad> load_quad_elements(std::string file_name)
    {
        std::fstream file;
        std::vector<Fem::Quad> result;

        file.open(file_name);
        if(!file.good()){
            std::cout << "Couldn't load text file\n";
            return result;
        }

        std::string line;
        bool element_selection = false;

        while (std::getline(file, line)) {
            if (line == "*Quads") {
                element_selection = true;
                continue;
            }

            if (line == "*Nodes" || line == "*BC" || line == "*Triangles") {
                element_selection = false;
            }

            if (element_selection) {
                int id;
                int node_id1, node_id2, node_id3, node_id4;
                char comma;
                std::istringstream iss(line);
                if (iss >> id >> comma >> node_id1 >> comma >> node_id2 >> comma >> node_id3 >> comma >> node_id4) {
                    result.emplace_back(id, node_id1, node_id2, node_id3, node_id4);
                }
            }
        }

        file.close();

        return result;
    }

    GlobalData load_configuration(const std::string& file_name)
    {
        std::fstream file;
        GlobalData data;

        file.open(file_name);
        if(!file.good()){
            std::cout << "Couldn't load text file\n";
            return data;
        }

        std::string line;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string key;
            double value;

            if (iss >> key >> value) {
                if (key == "SimulationTime") data.total_time = value;
                else if (key == "SimulationStepTime") data.time_step = value;
                else if (key == "Conductivity") data.conductivity = value;
                else if (key == "InitialTemp") data.init_temperature = value;
                else if (key == "Density") data.density = value;
                else if (key == "SpecificHeat") data.specific_heat = value;
                else if (key == "Nodes_number") data.node_number = value;
                else if (key == "Triangle_number") data.trian_number = value;
                else if (key == "Quad_number") data.quad_number= value;
            }
        }

        file.close();
        return data;
    }

    std::vector<Node> load_bc(const std::string& file_name, std::vector<Node>& nodes)
    {
        std::fstream file;

        file.open(file_name);
        if(!file.good()){
            std::cout << "Couldn't load text file\n";
            return nodes;
        }

        std::string line;
        bool bc_selection = false;

        while (std::getline(file, line)) {
            if (line == "*BC") {
                bc_selection = true;
                continue;
            }

            if (line == "*Nodes" || line == "*Triangles") {
                bc_selection = false;
            }

            if (bc_selection) {
                int node_id;
                double flux, alfa, t_ext, dir_temp;
                char comma;
                std::istringstream iss(line);
                if (iss >> node_id >> comma >> dir_temp >> comma >> flux >> comma >> alfa >> comma >> t_ext) {
                    BC_node temp_bc(node_id, flux, alfa, t_ext, dir_temp);
                    nodes[node_id].bc = temp_bc;
                }
                else {
                    std::cout << "Failed to parse line: " << line << "\n";
                }
            }
        }

        file.close();
        return nodes;
    }

    void print_config(GlobalData configuration)
    {
        std::cout<<"-----------------Configuration-----------------\n";
        std::cout<<"Total time = "<<configuration.total_time<<"\n";
        std::cout<<"Time step = "<<configuration.time_step<<"\n";
        std::cout<<"Conductivity = "<<configuration.conductivity<<"\n";
        std::cout<<"Initial temperature = "<<configuration.init_temperature<<"\n";
        std::cout<<"Density = "<<configuration.density<<"\n";
        std::cout<<"Specific heat = "<<configuration.specific_heat<<"\n";
        std::cout<<"Number of nodes = "<<configuration.node_number<<"\n";
        std::cout<<"Number of triangles = "<<configuration.trian_number<<"\n";
        std::cout<<"Number of quads = "<<configuration.quad_number<<"\n";
        std::cout<<"-----------------Configuration-----------------\n";
    }


    double dist(Node a, Node b){
        return sqrt(pow( a.x - b.x,2)+pow(a.y - b.y,2));
    }

    Matrix calc_local_H_tr(Triangle &local_el, std::vector<Node> &nodes, const double conductivity)
    {
        Matrix dNdx = Matrix::Zero(3,1);
        Matrix dNdy = Matrix::Zero(3,1);

        dNdx(0,0) = local_el.dN1dx_tr(nodes);
        dNdx(1,0) = local_el.dN2dx_tr(nodes);
        dNdx(2,0) = local_el.dN3dx_tr(nodes);
        dNdy(0,0) = local_el.dN1dy_tr(nodes);
        dNdy(1,0) = local_el.dN2dy_tr(nodes);
        dNdy(2,0) = local_el.dN3dy_tr(nodes);

        double x1 = nodes[local_el.node_ids[0]].x;
        double y1 = nodes[local_el.node_ids[0]].y;
        double x2 = nodes[local_el.node_ids[1]].x;
        double y2 = nodes[local_el.node_ids[1]].y;
        double x3 = nodes[local_el.node_ids[2]].x;
        double y3 = nodes[local_el.node_ids[2]].y;

        double D = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
        double Area = std::abs(D) / 2.0f;

        Matrix H_local = ((dNdx*dNdx.transpose()) + (dNdy*dNdy.transpose())) * (conductivity * Area);

        return H_local;
    }

    Matrix calc_local_Hbc_tr(Triangle &local_el, std::vector<Node> &nodes)
    {
        Matrix Hbc = Matrix::Zero(3, 3);

        for (int i = 0; i < 3; i++) {
            int id1 = local_el.node_ids[i];
            int id2 = local_el.node_ids[(i + 1) % 3];

            if (! (nodes[id1].bc.exist && nodes[id2].bc.exist)) continue;

            const double alfa = nodes[id1].bc.alfa;
            const double detJ = 0.5f*dist(nodes[id1], nodes[id2]);

            for (int pc=0; pc<2; ++pc) {
                const int pc_aktualny = i*2+pc;
                Matrix Ns = Matrix::Zero(3,1);

                Ns(0,0) = N1_tr(bc_pc_xi_tr[pc_aktualny], bc_pc_eta_tr[pc_aktualny]);
                Ns(1,0) = N2_tr(bc_pc_xi_tr[pc_aktualny], bc_pc_eta_tr[pc_aktualny]);
                Ns(2,0) = N3_tr(bc_pc_xi_tr[pc_aktualny], bc_pc_eta_tr[pc_aktualny]);

                Hbc += (Ns*Ns.transpose()) * (bc_pc_w_tr[pc_aktualny]*detJ*alfa);
            }
        }

        return Hbc;
    }

    Matrix calc_p_vec_tr(const Triangle &local_el, const std::vector<Node> &nodes)
    {
        Matrix p_vec = Matrix::Zero(3, 1);

        for (int i = 0; i < 3; i++) {
            const int id1 = local_el.node_ids[i];
            const int id2 = local_el.node_ids[(i + 1) % 3];

            if (!(nodes[id1].bc.exist && nodes[id2].bc.exist)) continue;

            const double detJ = 0.5f*dist(nodes[id1], nodes[id2]);

            const double flux = nodes[id1].bc.flux;

            const double alfa = nodes[id1].bc.alfa;
            const double t_ext = nodes[id2].bc.t_ext;

            Matrix pvec_edge = Matrix::Zero(3,1);

            for (int pc=0; pc<2; ++pc) {
                int pc_aktualny = i*2+pc;
                Matrix Ns = Matrix::Zero(3,1);

                Ns(0,0) = N1_tr(bc_pc_xi_tr[pc_aktualny], bc_pc_eta_tr[pc_aktualny]);
                Ns(1,0) = N2_tr(bc_pc_xi_tr[pc_aktualny], bc_pc_eta_tr[pc_aktualny]);
                Ns(2,0) = N3_tr(bc_pc_xi_tr[pc_aktualny], bc_pc_eta_tr[pc_aktualny]);

                for (int j=0; j<3; ++j) {
                    pvec_edge(j,0) += Ns(j,0)*bc_pc_w_tr[pc_aktualny];
                }
            }

            for (int k =0; k<3; ++k) {
                p_vec(k,0) += pvec_edge(k,0)*alfa*detJ*t_ext;
            }

            for (int k=0; k<3; ++k) {
                p_vec(k,0) += pvec_edge(k,0)*detJ*flux;
            }
        }

        return p_vec;
    }

    Matrix calc_c_tr(const Triangle &local_el, const std::vector<Node> &nodes, const double density, const double specific_heat, const int c_lump=0)
    {
        Matrix c_matrix = Matrix::Zero(3, 3);

        const double x1 = nodes[local_el.node_ids[0]].x;
        const double y1 = nodes[local_el.node_ids[0]].y;
        const double x2 = nodes[local_el.node_ids[1]].x;
        const double y2 = nodes[local_el.node_ids[1]].y;
        const double x3 = nodes[local_el.node_ids[2]].x;
        const double y3 = nodes[local_el.node_ids[2]].y;

        double detJ = 0.25f *((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
        detJ = std::abs(detJ);

        for (int pc =0; pc<3; ++pc) {
            Matrix Ns = Matrix::Zero(3,1);
            Ns(0,0) = N1_tr(pc_xi_tr[pc], pc_eta_tr[pc]);
            Ns(1,0) = N2_tr(pc_xi_tr[pc], pc_eta_tr[pc]);
            Ns(2,0) = N3_tr(pc_xi_tr[pc], pc_eta_tr[pc]);

            c_matrix += (Ns*Ns.transpose()) * (pc_w_tr[pc]*density*specific_heat*detJ);
        }

        if (c_lump == 1) {
            Matrix c_lumped = Matrix::Zero(3, 3);
            for(int i=0; i<3; ++i) {
                c_lumped(i,i) = c_matrix.row(i).sum();
            }
            c_matrix = c_lumped;
        }

        return c_matrix;
    }

    void aggregate_tr(Matrix &Global, const Triangle& element, Matrix &Local){
        for(int i=0; i<3;i++){
            for(int j=0; j<3;j++){
                const int glob_i = element.node_ids[i];
                const int glob_j = element.node_ids[j];
                Global(glob_i,glob_j) += Local(i,j);
            }
        }
    }

    void aggregate_p_vec_tr(Matrix &P_vec, const Triangle& element, Matrix &Local){
        for(int i=0; i<3; i++){
            const int glob_i = element.node_ids[i];
            P_vec(glob_i,0) += Local(i,0);
        }
    }

    Matrix calc_local_H_q(Fem::Quad &element, std::vector<Fem::Node> &nodes, double conductivity)
    {
        Fem::Matrix H = Matrix::Zero(4,4);

        for(int i=0; i<4; i++){
            Fem::Matrix jacob = element.jacobian_mat(element, nodes, Fem::pc_xi_q[i], Fem::pc_eta_q[i]);
            double det_J = element.det_jacobian(jacob);

            Fem::Matrix inv_jacob = element.inv_jacobian_mat(jacob);

            Fem::Matrix dNdx = Matrix::Zero(4,1);
            Fem::Matrix dNdy = Matrix::Zero(4,1);

            dNdx(0,0) = inv_jacob(0,0)*Fem::dN1dxi_q(Fem::pc_eta_q[i]) + inv_jacob(0,1)*Fem::dN1deta_q(Fem::pc_xi_q[i]);
            dNdx(1,0) = inv_jacob(0,0)*Fem::dN2dxi_q(Fem::pc_eta_q[i]) + inv_jacob(0,1)*Fem::dN2deta_q(Fem::pc_xi_q[i]);
            dNdx(2,0) = inv_jacob(0,0)*Fem::dN3dxi_q(Fem::pc_eta_q[i]) + inv_jacob(0,1)*Fem::dN3deta_q(Fem::pc_xi_q[i]);
            dNdx(3,0) = inv_jacob(0,0)*Fem::dN4dxi_q(Fem::pc_eta_q[i]) + inv_jacob(0,1)*Fem::dN4deta_q(Fem::pc_xi_q[i]);

            dNdy(0,0) = inv_jacob(1,0)*Fem::dN1dxi_q(Fem::pc_eta_q[i]) + inv_jacob(1,1)*Fem::dN1deta_q(Fem::pc_xi_q[i]);
            dNdy(1,0) = inv_jacob(1,0)*Fem::dN2dxi_q(Fem::pc_eta_q[i]) + inv_jacob(1,1)*Fem::dN2deta_q(Fem::pc_xi_q[i]);
            dNdy(2,0) = inv_jacob(1,0)*Fem::dN3dxi_q(Fem::pc_eta_q[i]) + inv_jacob(1,1)*Fem::dN3deta_q(Fem::pc_xi_q[i]);
            dNdy(3,0) = inv_jacob(1,0)*Fem::dN4dxi_q(Fem::pc_eta_q[i]) + inv_jacob(1,1)*Fem::dN4deta_q(Fem::pc_xi_q[i]);

            H += ((dNdx*dNdx.transpose()) + (dNdy*dNdy.transpose())) * (conductivity*det_J*Fem::pc_weights_q[i]);
        }

        return H;
    }

    Matrix calc_local_Hbc_q(Fem::Quad &element, std::vector<Fem::Node> &nodes)
    {
        Fem::Matrix H_bc = Matrix::Zero(4,4);

        for(int edge=0; edge<4; edge++){
            Fem::Matrix Hbc_edge = Matrix::Zero(4,4);
            int id1 = element.node_ids[edge];
            int id2 = element.node_ids[(edge+1)%4];

            if(!(nodes[id1].bc.exist && nodes[id2].bc.exist)){
                continue;
            }

            double alfa = nodes[id1].bc.alfa;

            double det_J = 0.5*dist(nodes[id1], nodes[id2]);

            for(int pc=0; pc<2; pc++){
                int pc_on_egde = edge*2+pc;
                Fem::Matrix Ns = Matrix::Zero(4,1);

                Ns(0,0) = Fem::N1_q(Fem::bc_xi_q[pc_on_egde], Fem::bc_eta_q[pc_on_egde]);
                Ns(1,0) = Fem::N2_q(Fem::bc_xi_q[pc_on_egde], Fem::bc_eta_q[pc_on_egde]);
                Ns(2,0) = Fem::N3_q(Fem::bc_xi_q[pc_on_egde], Fem::bc_eta_q[pc_on_egde]);
                Ns(3,0) = Fem::N4_q(Fem::bc_xi_q[pc_on_egde], Fem::bc_eta_q[pc_on_egde]);

                Hbc_edge += (Ns*Ns.transpose()) * Fem::bc_weights_q[pc_on_egde];
            }

            H_bc += Hbc_edge * (det_J*alfa);
        }

        return H_bc;
    }

    Matrix calc_P_q(Fem::Quad &element, std::vector<Fem::Node> &nodes)
    {
        Fem::Matrix p_vec = Matrix::Zero(4,1);

        for(int edge=0; edge<4; edge++){
            Fem::Matrix P_edge = Matrix::Zero(4,1);
            int id1 = element.node_ids[edge];
            int id2 = element.node_ids[(edge+1)%4];

            if(!(nodes[id1].bc.exist && nodes[id2].bc.exist)){
                continue;
            }

            double temperature = nodes[id1].bc.t_ext;
            double alfa = nodes[id1].bc.alfa;

            double flux = nodes[id1].bc.flux;

            double det_J = 0.5*dist(nodes[id1], nodes[id2]);

            for(int pc=0; pc<2; pc++){
                int pc_on_egde = edge*2+pc;
                Fem::Matrix Ns = Matrix::Zero(4,1);

                Ns(0,0) = Fem::N1_q(Fem::bc_xi_q[pc_on_egde], Fem::bc_eta_q[pc_on_egde]);
                Ns(1,0) = Fem::N2_q(Fem::bc_xi_q[pc_on_egde], Fem::bc_eta_q[pc_on_egde]);
                Ns(2,0) = Fem::N3_q(Fem::bc_xi_q[pc_on_egde], Fem::bc_eta_q[pc_on_egde]);
                Ns(3,0) = Fem::N4_q(Fem::bc_xi_q[pc_on_egde], Fem::bc_eta_q[pc_on_egde]);

                for(int i=0; i<4;i++){
                    P_edge(i,0) += Ns(i,0)*Fem::bc_weights_q[pc_on_egde]*temperature;
                }
            }

            for(int i=0; i<4;i++){
                p_vec(i,0) += P_edge(i,0)*det_J*alfa;
            }

            for(int i=0; i<4;i++){
                p_vec(i,0) += P_edge(i,0)*det_J*flux;
            }
        }

        return p_vec;
    }

    Matrix calc_local_C_q(Fem::Quad &element, std::vector<Fem::Node> &nodes, double rho, double c)
    {
        Fem::Matrix C = Matrix::Zero(4,4);

        for(int i=0; i<4; i++){
            Fem::Matrix jacob = element.jacobian_mat(element, nodes, Fem::pc_xi_q[i], Fem::pc_eta_q[i]);
            double det_J = element.det_jacobian(jacob);

            Fem::Matrix base_functions_values = Matrix::Zero(4,1);

            base_functions_values(0,0) = Fem::N1_q(Fem::pc_xi_q[i], Fem::pc_eta_q[i]);
            base_functions_values(1,0) = Fem::N2_q(Fem::pc_xi_q[i], Fem::pc_eta_q[i]);
            base_functions_values(2,0) = Fem::N3_q(Fem::pc_xi_q[i], Fem::pc_eta_q[i]);
            base_functions_values(3,0) = Fem::N4_q(Fem::pc_xi_q[i], Fem::pc_eta_q[i]);

            C += (base_functions_values*base_functions_values.transpose()) * (rho*c*det_J*Fem::pc_weights_q[i]);
        }

        return C;
    }

    void aggregate_q(Fem::Matrix &Global, Fem::Quad element, Fem::Matrix &Local)
    {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int glob_i = element.node_ids[i];
                int glob_j = element.node_ids[j];

                Global(glob_i,glob_j) += Local(i,j);
            }
        }
    }

    void aggregate_p_vec_q(Fem::Matrix &P_vec, Fem::Quad element, Fem::Matrix &Local)
    {
        for(int i=0; i<4; i++){
            int glob_i = element.node_ids[i];
            P_vec(glob_i,0) += Local(i,0);
        }
    }

    void write_to_vtu_file(int step, const std::vector<Fem::Node> &nodes,
    const std::vector<double> &temp, const std::vector<Fem::Triangle> &elements,
    const std::vector<Fem::Quad> &quads)
    {
        std::ostringstream fname;
        fname << "sol_" << step << ".vtu";
        const std::string path = "Data/" + fname.str();

        std::ofstream out(path);
        if (!out.is_open()) {
            std::perror(("Nie można otworzyć pliku " + path).c_str());
            return;
        }

        size_t total_cells = elements.size() + quads.size();

        out << R"(<?xml version="1.0"?>)"
            << "\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            "  <UnstructuredGrid>\n"
            "    <Piece NumberOfPoints=\"" << nodes.size()
            << "\" NumberOfCells=\"" << total_cells << "\">\n";

        out << "      <Points>\n"
            "        <DataArray type=\"double32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        out << std::setprecision(6);
        for (const auto& n : nodes) {
            out << "          " << n.x << " " << n.y << " 0\n";
        }
        out << "        </DataArray>\n"
            "      </Points>\n";

        out << "      <Cells>\n"
            "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

        for (const auto& e : elements) {
            out << "          " << e.node_ids[0] << " " << e.node_ids[1] << " " << e.node_ids[2] << "\n";
        }
        for (const auto& q : quads) {
            out << "          " << q.node_ids[0] << " " << q.node_ids[1] << " " << q.node_ids[2] << " " << q.node_ids[3] << "\n";
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        int offset = 0;
        for (size_t i = 0; i < elements.size(); ++i) {
            offset += 3;
            out << "          " << offset << "\n";
        }
        for (size_t i = 0; i < quads.size(); ++i) {
            offset += 4;
            out << "          " << offset << "\n";
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (size_t i = 0; i < elements.size(); ++i) {
            out << "          5\n";
        }
        for (size_t i = 0; i < quads.size(); ++i) {
            out << "          9\n";
        }
        out << "        </DataArray>\n"
            "      </Cells>\n";

        out << "      <PointData Scalars=\"Temperature\">\n"
            "        <DataArray type=\"double32\" Name=\"Temperature\" format=\"ascii\">\n";
        for (double t : temp) {
            out << "          " << t << "\n";
        }
        out << "        </DataArray>\n"
            "      </PointData>\n";

        out << "    </Piece>\n"
            "  </UnstructuredGrid>\n"
            "</VTKFile>\n";

        out.close();
    }

    void Solution::apply_dirichlet_symmetric(Matrix &H, std::vector<double> &P) {
        for (int d = 0; d < conf.node_number; ++d) {
            if (nodes[d].bc.exist && nodes[d].bc.dir_temp != 0.0) {
                double T_presc = nodes[d].bc.dir_temp;

                for (int i = 0; i < conf.node_number; ++i) {
                    if (i != d) {
                        P[i] -= H(i,d) * T_presc;
                        H(i,d) = 0.0;
                        H(d,i) = 0.0;
                    }
                }

                for (int j = 0; j < conf.node_number; ++j) {
                    if (j != d) {
                        H(d,j) = 0.0;
                    }
                }

                P[d] = +H(d,d) * T_presc;
            }
        }
    }

    double Solution::stability_dt_max() {
        double dt_max = std::numeric_limits<double>::max();
        for (int j = 0; j < conf.node_number; j++) {
            // Dirichlet and disconnected nodes are held/overwritten, not integrated.
            if (nodes[j].bc.exist && nodes[j].bc.dir_temp != 0.0) {
                continue;
            }

            if (Global_C(j,j) == 0.0) {
                continue;
            }

            double row_abs = 0.0;
            for (int k = 0; k < conf.node_number; k++) {
                row_abs += std::abs(Global_H(j,k));
            }
            if (row_abs == 0.0) continue;

            dt_max = std::min(dt_max, 2.0 * Global_C(j,j) / row_abs);
        }
        return dt_max;
    }

    void Solution::solve(const bool write_vtu, bool print_conf) {

        std::cout << "\n==================================================\n";
        std::cout << "Solving non-stationary (Transient) heat equation...\n";
        std::cout << "==================================================\n";

        std::vector<double> t0(conf.node_number);
        std::vector<double> t1(conf.node_number);

        for(int i=0; i<conf.node_number;i++){
            if (nodes[i].bc.exist && nodes[i].bc.dir_temp != 0.0) {
                t0[i] = nodes[i].bc.dir_temp;
            } else {
                t0[i] = conf.init_temperature;
            }
        }

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        if (this->solver_type == "implicit_euler") {
            Fem::Matrix Global_static = Matrix::Zero(conf.node_number, conf.node_number);
            for(int i=0; i<conf.node_number; i++){
                for(int j=0; j<conf.node_number; j++){
                    Global_static(i,j) = Global_H(i,j) + (Global_C(i,j) / conf.time_step);
                }
            }

            std::cout<<"Beginning time integration (Implicit Euler)...\n";
            for(int i=conf.time_step; i<=conf.total_time; i+= conf.time_step){
                Fem::Matrix Global_active = Global_static;

                for(int j=0; j<conf.node_number; j++){
                    double rhs = 0;
                    for(int k = 0; k<conf.node_number; k++){
                        rhs += (Global_C(j,k) / conf.time_step) * t0[k];
                    }
                    rhs += Global_P(j,0);
                    t1[j] = rhs;
                }

                apply_dirichlet_symmetric(Global_active, t1);

                t1 = cholesky_ldl(Global_active, t1);

                std::cout << "\nMIN: " << *std::min_element(t1.begin(), t1.end()) << " MAX: " << *std::max_element(t1.begin(), t1.end()) << std::endl;

                if(write_vtu){
                    write_to_vtu_file((int)i, nodes, t1, triangles, quads);
                }

                showProgress(i, conf.total_time);

                t0=t1;
            }
            std::cout << "\nMIN: " << *std::min_element(t1.begin(), t1.end()) << " MAX: " << *std::max_element(t1.begin(), t1.end()) << std::endl;
        }

        else if (this->solver_type == "explicit_euler") {
            const double dt_max = stability_dt_max();
            std::cout << "Max stable dt = " << dt_max << " (using dt = " << conf.time_step << ")\n";
            if (conf.time_step > dt_max) {
                std::cout << "[WARNING] dt exceeds stability limit. Simulation will be UNSTABLE!\n";
            }

            std::cout << "Beginning time integration (Explicit Euler)\n";
            int steps = static_cast<int>(conf.total_time / conf.time_step);
            for(int i = 0; i <= steps; i++) {
                for (int j = 0; j < conf.node_number; j++) {

                    if (nodes[j].bc.exist && nodes[j].bc.dir_temp != 0.0) {
                        t1[j] = nodes[j].bc.dir_temp;
                        continue;
                    }

                    if (Global_C(j,j) == 0.0) {
                        t1[j] = t0[j];
                        continue;
                    }

                    double H_t0 = 0.0;
                    for (int k = 0; k < conf.node_number; k++) {
                        H_t0 += Global_H(j,k) * t0[k];
                    }

                    double nawias = H_t0 - Global_P(j,0);
                    double wspolczynnik = conf.time_step / Global_C(j,j);

                    t1[j] = t0[j] - wspolczynnik * nawias;
                }

                if(write_vtu){
                    write_to_vtu_file((int)i, nodes, t1, triangles, quads);
                }

                showProgress(i, conf.total_time);
                t0 = t1;
            }

            std::cout << "\nMIN: " << *std::min_element(t1.begin(), t1.end())<< " MAX: " << *std::max_element(t1.begin(), t1.end()) << std::endl;
        }
        else if (this->solver_type == "crank-nicolson") {
            Fem::Matrix Global_static = Matrix::Zero(conf.node_number, conf.node_number);
            for(int i=0; i<conf.node_number; i++){
                for(int j=0; j<conf.node_number; j++){
                    Global_static(i,j) = 0.5*Global_H(i,j) + (Global_C(i,j) / conf.time_step);
                }
            }

            std::cout<<"Beginning time integration (Crank-Nicolson)...\n";
            for(int i=conf.time_step; i<=conf.total_time; i+= conf.time_step){
                Fem::Matrix Global_active = Global_static;

                for(int j=0; j<conf.node_number; j++){
                    double rhs = 0;
                    for(int k = 0; k<conf.node_number; k++){
                        rhs += ((Global_C(j,k) / conf.time_step) -0.5*Global_H(j,k))* t0[k];
                    }
                    rhs += Global_P(j,0);
                    t1[j] = rhs;
                }
                apply_dirichlet_symmetric(Global_active, t1);

                t1 = cholesky_ldl(Global_active, t1);

                std::cout << "\nMIN: " << *std::min_element(t1.begin(), t1.end()) << " MAX: " << *std::max_element(t1.begin(), t1.end()) << std::endl;

                if(write_vtu){
                    write_to_vtu_file((int)i, nodes, t1, triangles, quads);
                }

                showProgress(i, conf.total_time);

                t0=t1;
            }
            std::cout << "\nMIN: " << *std::min_element(t1.begin(), t1.end()) << " MAX: " << *std::max_element(t1.begin(), t1.end()) << std::endl;
        }
        else {
            std::cout<<"Unknown solver type\n";
        }

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        std::cout << "Solving time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
        std::cout << "==================================================\n\n";

    }

    void Solution::solve_stationary(const bool write_vtu) {
        std::cout << "\n==================================================\n";
        std::cout << "Solving stationary (steady-state) heat equation...\n";
        std::cout << "==================================================\n\n";

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        Fem::Matrix A = this->Global_H;
        std::vector<double> B(conf.node_number);
        for (int i = 0; i < conf.node_number; ++i) {
            B[i] = this->Global_P(i,0);
        }

        apply_dirichlet_symmetric(A, B);

        std::vector<double> T = cholesky_ldl(A, B);

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        std::cout << "MIN Temp: " << *std::min_element(T.begin(), T.end())
                  << " | MAX Temp: " << *std::max_element(T.begin(), T.end()) << std::endl;
        std::cout << "Stationary solving time = "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

        if (write_vtu) {
            write_to_vtu_file(0, nodes, T, triangles, quads);
        }
        std::cout << "==================================================\n\n";
    }

}

void fem_solve(const std::string& solver_type) {
    Fem::Solution solution("Data/fem_data.txt", solver_type);

    if (solver_type == "stationary") {
        solution.solve_stationary(true);
    } else {
        solution.solve(true, true);
    }
}