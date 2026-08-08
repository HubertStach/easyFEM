#include "geometry.h"

#include <algorithm>

#include "raylib.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "rlgl.h"
#include "imgui.h"
#include "../progress_bar.h"

static Color color_for_temp(float val, float min_val, float max_val) {
    if (std::fabs(max_val - min_val) < 0.0001f) return GREEN;
    float t = std::clamp((val - min_val) / (max_val - min_val), 0.0f, 1.0f);
    unsigned char r, g, b;
    if (t < 0.5f) {
        float lt = t * 2.0f;
        r = (unsigned char)(255 * lt); g = r; b = 255;
    } else {
        float lt = (t - 0.5f) * 2.0f;
        r = 255; g = (unsigned char)(255 * (1.0f - lt)); b = g;
    }
    return Color{r, g, b, 255};
}

geo::Node::Node()
{
    this->x = 0.0;
    this->y = 0.0;
    this->bc.is_bc = true;
}

geo::Node::Node(float x, float y)
{
    this->x = x;
    this->y = y;
    this->bc.is_bc = true;
}

geo::Node::Node(const float x, const float y, const bool is_bc) {
    this->x = x;
    this->y = y;
    this->bc.is_bc = is_bc;

    if (!is_bc) {
        this->color = GREEN;
    }
}

geo::Edge::Edge()
{
    this->node_ids[0] = -1;
    this->node_ids[1] = -1;
}

geo::Edge::Edge(int n1, int n2)
{
    this->node_ids[0] = n1;
    this->node_ids[1] = n2;
}

geo::Triangle::Triangle()
{
    this->node_ids[0] = -1;
    this->node_ids[1] = -1;
    this->node_ids[2] = -1;
}

geo::Triangle::Triangle(int n1, int n2, int n3)
{
    this->node_ids[0] = n1;
    this->node_ids[1] = n2;
    this->node_ids[2] = n3;
}

geo::Quad::Quad() = default;

geo::Quad::Quad(int n1, int n2, int n3, int n4) {
    this->node_ids[0] = n1;
    this->node_ids[1] = n2;
    this->node_ids[2] = n3;
    this->node_ids[3] = n4;
}

geo::Mesh::Mesh() {

}


void geo::Mesh::add_point(float x, float y, bool is_bc)
{
    for(Node n : this->nodes){
        if(n.x == x && n.y == y){
            return;
        }
    }

    Node temp_node(x, y);
    temp_node.bc.is_bc = is_bc;
    if (!is_bc) {
        temp_node.color = GREEN;
    }

    this->nodes.push_back(temp_node);
    this->initial_bc_nodes.push_back(temp_node);

    size_t s = this->nodes.size();

    if (s == 2) {
        geo::Edge e(0, 1);
        e.bc_edge.is_bc = true;
        this->edges.push_back(e);
    }
    else if (s == 3) {
        geo::Edge e1(1, 2); // B->C
        e1.bc_edge.is_bc = true;
        this->edges.push_back(e1);

        geo::Edge e2(2, 0);
        e2.bc_edge.is_bc = true;
        this->edges.push_back(e2);
    }
    else if (s > 3) {
        this->edges.pop_back();

        geo::Edge e1(static_cast<int>(s) - 2, static_cast<int>(s) - 1);
        e1.bc_edge.is_bc = true;
        this->edges.push_back(e1);

        geo::Edge e2(static_cast<int>(s) - 1, 0);
        e2.bc_edge.is_bc = true;
        this->edges.push_back(e2);
    }
}

void geo::Mesh::pop_point()
{
    if(this->mesh_created){
        this->triangles.clear();
        this->nodes.clear();
        this->edges.clear();
        this->initial_bc_nodes.clear();
        this->quads.clear();
        mesh_created = false;
        return;
    }

    if(this->nodes.empty()) return;

    size_t s = this->nodes.size();

    this->nodes.pop_back();
    this->initial_bc_nodes.pop_back();

    if (s == 2) {
        this->edges.clear();
    }
    else if (s == 3) {
        this->edges.pop_back();
        this->edges.pop_back();
    }
    else if (s > 3) {
        this->edges.pop_back();
        this->edges.pop_back();

        geo::Edge closing(static_cast<int>(s) - 2, 0);
        closing.bc_edge.is_bc = true;
        this->edges.push_back(closing);
    }
}

void geo::Mesh::create_edges() {
    this->edges.clear();

    std::vector<std::pair<int, int>> temp_edges;

    temp_edges.reserve(this->triangles.size() * 3 + this->quads.size() * 4);

    auto add_edge = [&](int n1, int n2) {
        temp_edges.emplace_back(std::min(n1, n2), std::max(n1, n2));
    };

    //zbieramy krawędzie trójkątów
    for (const auto& tr : this->triangles) {
        add_edge(tr.node_ids[0], tr.node_ids[1]);
        add_edge(tr.node_ids[1], tr.node_ids[2]);
        add_edge(tr.node_ids[2], tr.node_ids[0]);
    }

    //zbieramy krawędzie czworokątów
    for (const auto& quad : this->quads) {
        add_edge(quad.node_ids[0], quad.node_ids[1]);
        add_edge(quad.node_ids[1], quad.node_ids[2]);
        add_edge(quad.node_ids[2], quad.node_ids[3]);
        add_edge(quad.node_ids[3], quad.node_ids[0]);
    }

    //sortujemy krawędzie
    std::sort(temp_edges.begin(), temp_edges.end());

    //usuwamy duplikaty
    temp_edges.erase(std::unique(temp_edges.begin(), temp_edges.end()), temp_edges.end());

    this->edges.reserve(temp_edges.size());
    for (const auto& ep : temp_edges) {
        this->edges.emplace_back(ep.first, ep.second);
    }
}


//------------------- wyświetlanie -------------------

void geo::Mesh::draw_nodes(const float size) const {
    for(const Node &n:this->nodes){
        Color color=n.color;
        Vector2 pos;
        pos.x = n.x;
        pos.y = n.y;

        if (n.bc.initialised) {
            color = ORANGE;
        }

        DrawCircleV(pos, size, color);
    }
}

void geo::Mesh::draw_edges()
{
    for (const auto& edge : this->edges) {
        if (edge.node_ids[0] < 0 || edge.node_ids[0] >= (int)nodes.size() || edge.node_ids[1] < 0 || edge.node_ids[1] >= (int)nodes.size())
        {
            continue;
        }

        const geo::Node& n1 = this->nodes[edge.node_ids[0]];
        const geo::Node& n2 = this->nodes[edge.node_ids[1]];

        Vector2 start_pos = { n1.x, n1.y };
        Vector2 end_pos   = { n2.x, n2.y };

        Color col = BLUE;
        if (edge.bc_edge.initialised) {
            col = RED;
        }

        DrawLineV(start_pos, end_pos, col);
    }
}


void geo::Mesh::draw_tr(Color c)
{
    for (const geo::Triangle& tr : this->triangles) {

        const geo::Node& node1 = this->nodes[tr.node_ids[0]];
        const geo::Node& node2 = this->nodes[tr.node_ids[1]];
        const geo::Node& node3 = this->nodes[tr.node_ids[2]];

        const Vector2 pos1 = { node1.x, node1.y };
        const Vector2 pos2 = { node2.x, node2.y };
        const Vector2 pos3 = { node3.x, node3.y };

        DrawLineV(pos1, pos2, c);
        DrawLineV(pos2, pos3, c);
        DrawLineV(pos3, pos1, c);
    }
}

void geo::Mesh::draw_q(Color c) {
    for (const geo::Quad& q : this->quads) {

        const geo::Node& node1 = this->nodes[q.node_ids[0]];
        const geo::Node& node2 = this->nodes[q.node_ids[1]];
        const geo::Node& node3 = this->nodes[q.node_ids[2]];
        const geo::Node& node4 = this->nodes[q.node_ids[3]];

        const Vector2 pos1 = { node1.x, node1.y };
        const Vector2 pos2 = { node2.x, node2.y };
        const Vector2 pos3 = { node3.x, node3.y };
        const Vector2 pos4 = { node4.x, node4.y };

        DrawLineV(pos1, pos2, c);
        DrawLineV(pos2, pos3, c);
        DrawLineV(pos3, pos4, c);
        DrawLineV(pos4, pos1, c);
    }
}

void geo::Mesh::draw_q_grad(std::vector<double> &temp, float max, float min) const {
    rlBegin(RL_TRIANGLES);

    for (const geo::Quad& q : this->quads) {
        const geo::Node& node1 = this->nodes[q.node_ids[0]];
        const geo::Node& node2 = this->nodes[q.node_ids[1]];
        const geo::Node& node3 = this->nodes[q.node_ids[2]];
        const geo::Node& node4 = this->nodes[q.node_ids[3]];

        float t1 = temp[q.node_ids[0]];
        float t2 = temp[q.node_ids[1]];
        float t3 = temp[q.node_ids[2]];
        float t4 = temp[q.node_ids[3]];

        float cx = (node1.x + node2.x + node3.x + node4.x) / 4.0f;
        float cy = (node1.y + node2.y + node3.y + node4.y) / 4.0f;
        float tc = (t1 + t2 + t3 + t4) / 4.0f;

        Color col1 = color_for_temp(t1, min, max);
        Color col2 = color_for_temp(t2, min, max);
        Color col3 = color_for_temp(t3, min, max);
        Color col4 = color_for_temp(t4, min, max);
        Color colC = color_for_temp(tc, min, max);

        //rysujemy trójkąty 2 krotnie aby uniknąć back-cullingu -> trójkąty się nie wyświetlają bo och wektor normalny
        //jest w drugą stronę od kamery/ekranu
        auto draw_triangle_both_ways = [&](const geo::Node& nA, Color cA, const geo::Node& nB, Color cB) {
            // Rysowanie z jedną rotacją
            rlColor4ub(cA.r, cA.g, cA.b, cA.a); rlVertex2f(nA.x, nA.y);
            rlColor4ub(cB.r, cB.g, cB.b, cB.a); rlVertex2f(nB.x, nB.y);
            rlColor4ub(colC.r, colC.g, colC.b, colC.a); rlVertex2f(cx, cy);

            // Rysowanie z odwrotną rotacją (gwarantuje widoczność 2D)
            rlColor4ub(cA.r, cA.g, cA.b, cA.a); rlVertex2f(nA.x, nA.y);
            rlColor4ub(colC.r, colC.g, colC.b, colC.a); rlVertex2f(cx, cy);
            rlColor4ub(cB.r, cB.g, cB.b, cB.a); rlVertex2f(nB.x, nB.y);
        };

        draw_triangle_both_ways(node1, col1, node2, col2);
        draw_triangle_both_ways(node2, col2, node3, col3);
        draw_triangle_both_ways(node3, col3, node4, col4);
        draw_triangle_both_ways(node4, col4, node1, col1);
    }

    // Kończymy rysowanie
    rlEnd();
}

void geo::Mesh::draw_tr_grad(std::vector<double> &temp, float max, float min) const {
    rlBegin(RL_TRIANGLES);

    for (const geo::Triangle& tr : this->triangles) {
        const geo::Node& node1 = this->nodes[tr.node_ids[0]];
        const geo::Node& node2 = this->nodes[tr.node_ids[1]];
        const geo::Node& node3 = this->nodes[tr.node_ids[2]];

        float t1 = temp[tr.node_ids[0]];
        float t2 = temp[tr.node_ids[1]];
        float t3 = temp[tr.node_ids[2]];

        Color col1 = color_for_temp(t1, min, max);
        Color col2 = color_for_temp(t2, min, max);
        Color col3 = color_for_temp(t3, min, max);

        // --- RYSOWANIE TRÓJKĄTA (Zgodnie z kierunkiem wskazówek zegara / przeciw) ---
        // Należy zdefiniować kolor ZANIM zdefiniuje się wierzchołek

        // Trójkąt pierwszy

        rlColor4ub(col1.r, col1.g, col1.b, col1.a);
        rlVertex2f(node1.x, node1.y);

        rlColor4ub(col2.r, col2.g, col2.b, col2.a);
        rlVertex2f(node2.x, node2.y);

        rlColor4ub(col3.r, col3.g, col3.b, col3.a);
        rlVertex2f(node3.x, node3.y);

        // Trójkąt drugi (Odwrotna kolejność wierzchołków, odpowiada to rysowaniu
        // dwustronnemu z oryginalnego kodu - DrawTriangle(pos1, pos3, pos2))


        rlColor4ub(col1.r, col1.g, col1.b, col1.a);
        rlVertex2f(node1.x, node1.y);

        rlColor4ub(col3.r, col3.g, col3.b, col3.a);
        rlVertex2f(node3.x, node3.y);

        rlColor4ub(col2.r, col2.g, col2.b, col2.a);
        rlVertex2f(node2.x, node2.y);

    }

    // Kończymy rysowanie
    rlEnd();
}


//---------------Ładowanie siatki------------------

void geo::Mesh::load_nodes(const std::string &filepath) {
    std::vector<Node> loaded_nodes;
    std::fstream file;

    file.open(filepath);
    if(!file.good()){
        std::cout << "Couldn't load text file\n";
        return ;
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
            float x, y;
            char comma;
            std::istringstream iss(line);
            if (iss >> id >> comma >> x >> comma >> y) {
                loaded_nodes.emplace_back(x*100, y*100, false);
            }
        }
    }
    file.close();
    this->nodes = loaded_nodes;
}

void geo::Mesh::load_tr_elements(const std::string &filepath) {
    std::vector<Triangle> tr_elements;
    std::fstream file;

    file.open(filepath);
    if(!file.good()){
        std::cout << "Couldn't load text file\n";
        return ;
    }

    std::string line;
    bool tr_selection = false;

    while (std::getline(file, line)) {
        if (line == "*Triangles") {
            tr_selection = true;
            continue;
        }

        if (line == "*Nodes" || line == "*BC" || line == "*Quads") {
            tr_selection = false;
        }

        if (tr_selection) {
            int id;
            int node_id1, node_id2, node_id3;
            char comma;
            std::istringstream iss(line);
            if (iss >> id >> comma >> node_id1 >> comma >> node_id2 >> comma >> node_id3) {
                tr_elements.emplace_back(node_id1, node_id2, node_id3);
            }
        }
    }
    file.close();
    this->triangles = tr_elements;
}

void geo::Mesh::load_q_elements(const std::string &filepath) {
    std::vector<Quad> q_elements;
    std::fstream file;

    file.open(filepath);
    if(!file.good()){
        std::cout << "Couldn't load text file\n";
        return ;
    }

    std::string line;
    bool tr_selection = false;

    while (std::getline(file, line)) {
        if (line == "*Quads") {
            tr_selection = true;
            continue;
        }

        if (line == "*Nodes" || line == "*BC" || line=="*Triangles") {
            tr_selection = false;
        }

        if (tr_selection) {
            int id;
            int node_id1, node_id2, node_id3, node_id4;
            char comma;
            std::istringstream iss(line);
            if (iss >> id >> comma >> node_id1 >> comma >> node_id2 >> comma >> node_id3 >> comma >>node_id4) {
                q_elements.emplace_back(node_id1, node_id2, node_id3, node_id4);
            }
        }
    }
    file.close();
    this->quads = q_elements;
}

void geo::Mesh::load_bcs(const std::string &filepath) {
    std::fstream file;

    file.open(filepath);
    if(!file.good()){
        std::cout << "Couldn't load text file\n";
        return ;
    }

    std::string line;
    bool bc_selection = false;

    while (std::getline(file, line)) {
        if (line == "*BC") {
            bc_selection = true;
            continue;
        }

        if (line == "*Nodes" || line == "*Triangles" || line == "*Quads") {
            bc_selection = false;
        }

        if (bc_selection) {
            int node_id;
            float flux, alfa, t_ext, dir_temp;
            char comma;
            std::istringstream iss(line);
            if (iss >> node_id >> comma >> dir_temp >> comma >> flux >> comma >> alfa >> comma >> t_ext) {
                this->nodes[node_id].bc.flux = flux;
                this->nodes[node_id].bc.alfa = alfa;
                this->nodes[node_id].bc.t_ext = t_ext;
                this->nodes[node_id].bc.dir_temp = dir_temp;
                this->nodes[node_id].bc.is_bc = true;
                this->nodes[node_id].bc.initialised = true;

            }
            else {
                std::cout << "Failed to parse line: " << line << "\n"; // Debug output
            }
        }
    }

    file.close();
}

//
void geo::Mesh::reconstruct_edges() {
    this->edges.clear();

    std::map<std::pair<int, int>, int> edge_count;

    for (const geo::Triangle& tr : this->triangles) {
        for (int i = 0; i < 3; ++i) {
            int n1 = tr.node_ids[i];
            int n2 = tr.node_ids[(i + 1) % 3];

            int min_n = std::min(n1, n2);
            int max_n = std::max(n1, n2);

            edge_count[{min_n, max_n}]++;
        }
    }

    for (const geo::Quad& q : this->quads) {
        for (int i = 0; i < 4; ++i) {
            int n1 = q.node_ids[i];
            int n2 = q.node_ids[(i + 1) % 4];

            int min_n = std::min(n1, n2);
            int max_n = std::max(n1, n2);

            edge_count[{min_n, max_n}]++;
        }
    }

    for (const auto& pair : edge_count) {
        if (pair.second == 1) {
            int n1 = pair.first.first;
            int n2 = pair.first.second;

            geo::Edge new_edge(n1, n2);

            if (this->nodes[n1].bc.is_bc && this->nodes[n2].bc.is_bc) {
                new_edge.bc_edge.is_bc = true;

                if (this->nodes[n1].bc.initialised && this->nodes[n2].bc.initialised) {
                    new_edge.bc_edge.initialised = true;
                    new_edge.bc_edge.flux = this->nodes[n1].bc.flux;
                    new_edge.bc_edge.alfa = this->nodes[n1].bc.alfa;
                    new_edge.bc_edge.t_ext = this->nodes[n1].bc.t_ext;
                }
            }

            this->edges.push_back(new_edge);
        }
    }

    this->initial_edges = this->edges;
}


void geo::Mesh::load_mesh_from_txt(const std::string& filepath) {
    this->load_nodes(filepath);
    this->load_tr_elements(filepath);
    this->load_q_elements(filepath);
    this->load_bcs(filepath);
    this->reconstruct_edges();
    this->mesh_created = true;
}

std::vector<int> geo::Mesh::get_continuous_edges(int clicked_edge_idx) {
    std::vector<int> result;

    if (clicked_edge_idx < 0 || clicked_edge_idx >= static_cast<int>(this->edges.size())) {
        return result;
    }

    result.push_back(clicked_edge_idx);

    //zaznaczamy tam które już odwiedziliśmy
    std::vector<bool> visited(this->edges.size(), false);
    visited[clicked_edge_idx] = true;

    geo::Edge start_edge = this->edges[clicked_edge_idx];
    int n1_id = start_edge.node_ids[0];
    int n2_id = start_edge.node_ids[1];

    // znormalizowany wektor kieurnkowy
    float dx0 = this->nodes[n2_id].x - this->nodes[n1_id].x;
    float dy0 = this->nodes[n2_id].y - this->nodes[n1_id].y;
    float len0 = std::sqrt(dx0 * dx0 + dy0 * dy0);

    if (len0 < 1e-6f) return result;

    float ux = dx0 / len0;
    float uy = dy0 / len0;

    auto walk = [&](int start_node_id) {
        int current_node = start_node_id;

        while (true) {
            int next_edge_idx = -1;
            int next_node = -1;

            for (size_t i = 0; i < this->edges.size(); ++i) {
                if (visited[i]) continue;

                //sąsiad w prawo
                if (this->edges[i].node_ids[0] == current_node) {
                    next_edge_idx = i;
                    next_node = this->edges[i].node_ids[1];
                    break;
                }
                //sąsiad w lewo
                else if (this->edges[i].node_ids[1] == current_node) {
                    next_edge_idx = i;
                    next_node = this->edges[i].node_ids[0];
                    break;
                }
            }

            if (next_edge_idx == -1) break;

            float dx = this->nodes[next_node].x - this->nodes[current_node].x;
            float dy = this->nodes[next_node].y - this->nodes[current_node].y;
            float len = std::sqrt(dx * dx + dy * dy);

            if (len < 1e-6f) break;

            float vx = dx / len;
            float vy = dy / len;

            float cross = ux * vy - uy * vx;

            constexpr float tolerance = 0.01f;

            if (std::abs(cross) < tolerance) {
                result.push_back(next_edge_idx);
                visited[next_edge_idx] = true;

                current_node = next_node;
            } else {
                break;
            }
        }
    };

    walk(n2_id);
    walk(n1_id);

    return result;
}



//-------------------triangulacja siatki---------------

float geo::len(geo::Node A, geo::Node B){
    return sqrt(pow(B.x-A.x,2)+pow(B.y-A.y,2));
}

float geo::tr_size(geo::Triangle &tr, std::vector<geo::Node> nodes){
    const geo::Node A = nodes[tr.node_ids[0]];
    const geo::Node B = nodes[tr.node_ids[1]];
    const geo::Node C = nodes[tr.node_ids[2]];

    const float a = len(A, B);
    const float b = len(B, C);
    const float c = len(C, A);

    return 0.25f*(std::sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)));
}

void geo::Mesh::interpolate_bc_points(float spacing)
{
    size_t parent_edge_count = this->edges.size();
    if(parent_edge_count == 0) return;

    std::vector<geo::Node> new_nodes;
    std::vector<geo::Edge> new_edges;

    int global_node_index = 0;

    for (size_t i = 0; i < parent_edge_count; i++) {
        const geo::Edge& parent_edge = this->edges[i];

        geo::Node& A = this->nodes[parent_edge.node_ids[0]];
        geo::Node& B = this->nodes[parent_edge.node_ids[1]];

        float len = sqrt(pow(B.x-A.x,2) + pow(B.y-A.y,2));

        int N = std::round(len/spacing);
        if (N < 1) N = 1;

        auto x_u = [&](float u) { return A.x + (B.x - A.x) * u; };
        auto y_u = [&](float u) { return A.y + (B.y - A.y) * u; };

        float x_diff = (B.x-A.x)*(1.0f/static_cast<float>(N));
        float y_diff = (B.y-A.y)*(1.0f/static_cast<float>(N));
        float temp_bc_max_len = std::sqrt(std::pow(x_diff,2) + std::pow(y_diff,2));

        if (this->max_bc_len <= temp_bc_max_len) {
            this->max_bc_len = temp_bc_max_len;
        }

        for (int j = 0; j < N; j++) {
            float u = static_cast<float>(j) / static_cast<float>(N);

            // 1. Tworzenie węzła
            geo::Node temp_node(x_u(u), y_u(u));
            temp_node.bc.is_bc = true;

            // Przypisanie BC do węzła
            if (parent_edge.bc_edge.initialised) {
                temp_node.bc = parent_edge.bc_edge;
                temp_node.bc.initialised = true;
            }
            new_nodes.push_back(temp_node);

            geo::Edge temp_edge(global_node_index, global_node_index + 1);
            temp_edge.bc_edge.is_bc = true;

            if (parent_edge.bc_edge.initialised) {
                temp_edge.bc_edge = parent_edge.bc_edge;
                temp_edge.bc_edge.initialised = true;
            }

            new_edges.push_back(temp_edge);

            global_node_index++;
        }
    }

    // Podmiana węzłów
    this->nodes = new_nodes;

    if (!new_edges.empty()) {
        new_edges.back().node_ids[1] = 0;
    }

    // Podmiana krawędzi
    this->edges = new_edges;
}

std::vector<geo::Node> geo::Mesh::super_triangle()
{
    float max_x, max_y, min_x, min_y;

    max_x = nodes[0].x;
    max_y = nodes[0].y;
    min_x = nodes[0].x;
    min_y = nodes[0].y;

    for(auto&it:nodes){
        if(it.x > max_x){
            max_x = it.x;
        }

        if(it.x < min_x){
            min_x = it.x;
        }

        if(it.y > max_y){
            max_y = it.y;
        }

        if(it.y < min_y){
            min_y = it.y;
        }
    }

    float w = max_x - min_x;
    float h = max_y - min_y;

    Node p1(min_x - w*0.1f, min_y-h);
    Node p2(min_x+w*1.7f, min_y+h*0.5f);
    Node p3(min_x-w*0.1f, min_y+h*2.0f);

    std::vector<geo::Node> result = {p1, p2, p3};
    return result;
}

bool geo::Mesh::inside_circumcircle(geo::Triangle tr, int node_id)
{
    geo::Node A = this->nodes[tr.node_ids[0]];
    geo::Node B = this->nodes[tr.node_ids[1]];
    geo::Node C = this->nodes[tr.node_ids[2]];
    geo::Node D = this->nodes[node_id];

    double adx = A.x - D.x;
    double ady = A.y - D.y;
    double bdx = B.x - D.x;
    double bdy = B.y - D.y;
    double cdx = C.x - D.x;
    double cdy = C.y - D.y;

    double bcdet = bdx * cdy - bdy * cdx;
    double cadet = cdx * ady - cdy * adx;
    double abdet = adx * bdy - ady * bdx;

    double alift = adx * adx + ady * ady;
    double blift = bdx * bdx + bdy * bdy;
    double clift = cdx * cdx + cdy * cdy;

    double det = alift * bcdet + blift * cadet + clift * abdet;

    double area_ABC = (B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y);

    if (area_ABC == 0.0) {
        return false;
    }

    return (det * area_ABC > 0);
}

void geo::Mesh::circumcenter(const geo::Triangle &tr, geo::Node &center) {
    geo::Node A = this->nodes[tr.node_ids[0]];
    geo::Node B = this->nodes[tr.node_ids[1]];
    geo::Node C = this->nodes[tr.node_ids[2]];

    const double D = 2*(A.x*(B.y-C.y)+B.x*(C.y-A.y)+C.x*(A.y-B.y));

    if (std::abs(D) < 1e-10f) {
        return;
    }

    const double X = (1/D)*((pow(A.x,2)+pow(A.y,2))*(B.y-C.y) + (pow(B.x,2)+pow(B.y,2))*(C.y-A.y)+
        (pow(C.x,2)+pow(C.y,2))*(A.y-B.y));

    const double Y = (1/D)*((pow(A.x,2)+pow(A.y,2))*(C.x-B.x) + (pow(B.x,2)+pow(B.y,2))*(A.x-C.x)+
        (pow(C.x,2)+pow(C.y,2))*(B.x-A.x));

    geo::Node temp(X, Y, false);

    center = temp;
}

bool geo::Mesh::same_triangle(geo::Triangle A, geo::Triangle B )
{

    int matches = 0;

    for (const int node_id : A.node_ids) {
        for (const int j : B.node_ids) {
            if ((node_id == j) && (node_id == j)) {
                matches++;
                break;
            }
        }
    }

    return matches == 3;
}

bool geo::Mesh::is_boundary_edge(std::vector<geo::Triangle> &triangles, geo::Edge edge)
{
    int count = 0;

    for (const auto& triangle : triangles) {
        for (int i=0; i<3; i++) {
            if ((triangle.node_ids[i] == edge.node_ids[0] && triangle.node_ids[(i+1)%3] == edge.node_ids[1]) ||
                (triangle.node_ids[(i+1)%3] == edge.node_ids[0] && triangle.node_ids[i] == edge.node_ids[1]))
            {
                count++;
                if (count > 1) return false;
            }
        }
    }
    return count == 1;
}

//algorytm bowyera-watsona
void geo::Mesh::triangulate()
{
    size_t original_nodes_count = this->nodes.size();
    if (original_nodes_count < 3) return;

    std::vector<geo::Node> original_nodes = this->nodes;

    std::vector<geo::Node> st = this->super_triangle();
    this->nodes.insert(this->nodes.begin(), st.begin(), st.end());

    std::vector<geo::Triangle> triangulation;
    triangulation.emplace_back(0, 1, 2);

    for (size_t i = 0; i < original_nodes_count; ++i) {
        int node_id = (int)(i + 3);

        // collect bad triangle indices (O(n))
        std::vector<size_t> bad_idx;
        for (size_t t = 0; t < triangulation.size(); ++t) {
            if (inside_circumcircle(triangulation[t], node_id))
                bad_idx.push_back(t);
        }

        // boundary edges: those belonging to exactly one bad triangle (O(k log k))
        std::map<std::pair<int,int>, int> edge_count;
        for (size_t t : bad_idx) {
            const int* n = triangulation[t].node_ids;

            for (int j = 0; j < 3; ++j)
            {
                edge_count[{std::min(n[j], n[(j+1)%3]), std::max(n[j], n[(j+1)%3])}]++;
            }
        }

        std::vector<bool> is_bad(triangulation.size(), false);
        for (size_t t : bad_idx) is_bad[t] = true;
        size_t w = 0;
        for (size_t t = 0; t < triangulation.size(); ++t)
            if (!is_bad[t]) triangulation[w++] = triangulation[t];
        triangulation.resize(w);

        for (auto& [e, cnt] : edge_count){
            if (cnt == 1){
                triangulation.emplace_back(e.first, e.second, node_id);
            }
        }
    }

    this->nodes = original_nodes;

    this->triangles.clear();
    for (const auto& tr : triangulation) {
        if (tr.node_ids[0] >= 3 && tr.node_ids[1] >= 3 && tr.node_ids[2] >= 3){
            this->triangles.emplace_back(tr.node_ids[0]-3, tr.node_ids[1]-3, tr.node_ids[2]-3);
        }
    }
}

bool geo::Mesh::point_in_mesh(float x, float y) {
    bool inside = false;
    size_t count = this->initial_bc_nodes.size();

    for (size_t i = 0, j = count - 1; i < count; j = i++) {
        float xi = this->initial_bc_nodes[i].x;
        float yi = this->initial_bc_nodes[i].y;
        float xj = this->initial_bc_nodes[j].x;
        float yj = this->initial_bc_nodes[j].y;

        bool intersect = ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);

        if (intersect) {
            inside = !inside;
        }
    }

    return inside;
}


// źródło algorytmu
// Owen, Steven. (2000). A Survey of Unstructured Mesh Generation Technology. 7th International Meshing Roundtable. 3.
// rodział 2.2.1 Point insertion
void geo::Mesh::create_nodes(float alpha, float beta) {
    // Krok 1 (z artykułu): Obliczenie funkcji dystrybucji punktów (dp_o) dla brzegu.
    std::vector<float> dp(this->nodes.size(), 0.0f);
    std::vector<int> edge_count(this->nodes.size(), 0);

    for (const auto& edge : this->edges) {
        int n1 = edge.node_ids[0];
        int n2 = edge.node_ids[1];
        float length = geo::len(this->nodes[n1], this->nodes[n2]);

        dp[n1] += length;
        dp[n2] += length;
        edge_count[n1]++;
        edge_count[n2]++;
    }

    // Wyciągnięcie średniej dla każdego węzła początkowego
    for (size_t i = 0; i < this->nodes.size(); ++i) {
        if (edge_count[i] > 0) {
            dp[i] /= static_cast<float>(edge_count[i]);
        } else {
            dp[i] = 1.0f; // Zabezpieczenie (fallback), jeśli węzeł "wisi" w powietrzu
        }
    }

    // Krok 2: Generacja początkowej triangulacji z punktów brzegowych
    std::cout << "Initial triangulation...\n";
    this->triangulate();

    // Główna pętla wstawiania punktów (Sweep)
    bool points_added = true;
    int current_iter = 0;
    int max_iter = 200;

    std::cout << "Inserting interior points...\n";
    while (points_added && current_iter < max_iter) {
        points_added = false;
        current_iter++;

        // Bufory na punkty zaakceptowane w TYM KONKRETNYM przejściu (Krok 3)
        std::vector<geo::Node> new_nodes_this_sweep;
        std::vector<float> new_dp_this_sweep;

        // Krok 4: Pętla po wszystkich wygenerowanych trójkątach
        for (const geo::Triangle& tr : this->triangles) {
            int id1 = tr.node_ids[0];
            int id2 = tr.node_ids[1];
            int id3 = tr.node_ids[2];

            geo::Node n1 = this->nodes[id1];
            geo::Node n2 = this->nodes[id2];
            geo::Node n3 = this->nodes[id3];

            // Krok 4a: Definicja potencjalnego punktu Q (Środek ciężkości)
            float cx = (n1.x + n2.x + n3.x) / 3.0f;
            float cy = (n1.y + n2.y + n3.y) / 3.0f;

            // Upewniamy się, że punkt znajduje się fizycznie wewnątrz naszej figury
            if (!this->point_in_mesh(cx, cy)) {
                continue;
            }

            geo::Node Q(cx, cy, false);

            // Krok 4b: Interpolacja gęstości dp_Q z węzłów trójkąta (średnia arytmetyczna z 3 punktów)
            float dp_Q = (dp[id1] + dp[id2] + dp[id3]) / 3.0f;

            // Krok 4c: TEST ALFA (Odległość od Q do wierzchołków jego własnego trójkąta)
            float d1 = geo::len(Q, n1);
            float d2 = geo::len(Q, n2);
            float d3 = geo::len(Q, n3);

            // Jeśli punkt jest zbyt blisko KTÓREGOKOLWIEK z istniejących wierzchołków trójkąta -> odrzuć
            if (d1 < alpha * dp_Q || d2 < alpha * dp_Q || d3 < alpha * dp_Q) {
                continue;
            }

            // Krok 4c: TEST BETA (Strefa wykluczenia dla innych punktów dodawanych w TEJ SAMEJ iteracji)
            bool rejected_by_beta = false;
            for (size_t j = 0; j < new_nodes_this_sweep.size(); ++j) {
                float dist_to_Pj = geo::len(Q, new_nodes_this_sweep[j]);

                // Jeśli nowy punkt leży w promieniu ochronnym innego nowo-zaproponowanego punktu -> odrzuć
                if (dist_to_Pj < beta * dp_Q) {
                    rejected_by_beta = true;
                    break;
                }
            }

            if (rejected_by_beta) {
                continue;
            }

            // Krok 4d: Akceptacja punktu
            new_nodes_this_sweep.push_back(Q);
            new_dp_this_sweep.push_back(dp_Q); // Zachowujemy jego "skalę" dla przyszłych generacji
            points_added = true;
        }

        // Krok 5 i 6: Dodanie nowych węzłów do globalnej listy i retriangulacja Delaunaya
        if (points_added) {
            for (size_t i = 0; i < new_nodes_this_sweep.size(); ++i) {
                this->nodes.push_back(new_nodes_this_sweep[i]);
                dp.push_back(new_dp_this_sweep[i]);
            }

            this->triangulate();
        }

        showProgress(current_iter, max_iter);
    }
    std::cout << "\nDone. Total nodes: " << this->nodes.size() << "\n";
}

void geo::Mesh::remove_outside_triangles() {
    std::vector<Triangle> kept;

    for (const auto& tr : triangles) {
        double cx = (nodes[tr.node_ids[0]].x + nodes[tr.node_ids[1]].x + nodes[tr.node_ids[2]].x) / 3.0;
        double cy = (nodes[tr.node_ids[0]].y + nodes[tr.node_ids[1]].y + nodes[tr.node_ids[2]].y) / 3.0;

        if (this->point_in_mesh(cx, cy)) {
            kept.push_back(tr);
        }
    }

    this->triangles = kept;
}

std::vector<geo::Edge> geo::Mesh::find_missing_bondaries() {
    struct SortedEdge {
        int a, b;
        SortedEdge(int n1, int n2) {
            a = std::min(n1, n2);
            b = std::max(n1, n2);
        }
        // Operator potrzebny dla std::set
        bool operator<(const SortedEdge& other) const {
            if (a != other.a) return a < other.a;
            return b < other.b;
        }
    };

    std::set<SortedEdge> existing_edges;
    std::vector<Edge> omited_edges;

    for (const auto& tr : this->triangles) {
        existing_edges.insert(SortedEdge(tr.node_ids[0], tr.node_ids[1]));
        existing_edges.insert(SortedEdge(tr.node_ids[1], tr.node_ids[2]));
        existing_edges.insert(SortedEdge(tr.node_ids[2], tr.node_ids[0]));
    }

    for (const auto& boundary_edge : this->edges) {
        SortedEdge search_edge(boundary_edge.node_ids[0], boundary_edge.node_ids[1]);

        if (existing_edges.find(search_edge) == existing_edges.end()) {
            omited_edges.push_back(boundary_edge);
        }
    }

    return omited_edges;
}

void geo::Mesh::recreate_mising_triangles() {
    std::vector<Edge> missing = this->find_missing_bondaries();
    std::vector<Triangle> recreated_tr;

    for (auto& missing_edge : missing) {
        int n1 = missing_edge.node_ids[0];
        int n2 = missing_edge.node_ids[1];

        float mid_x = (nodes[n1].x + nodes[n2].x) / 2.0f;
        float mid_y = (nodes[n1].y + nodes[n2].y) / 2.0f;

        int best_n3 = -1;
        float min_dist_sq = 1e15f;

        for (int i = 0; i < (int)nodes.size(); i++) {
            if (i == n1 || i == n2) continue;

            float d_sq = distSq(mid_x, mid_y, nodes[i].x, nodes[i].y);

            if (d_sq < min_dist_sq) {
                float ctx = (nodes[n1].x + nodes[n2].x + nodes[i].x) / 3.0f;
                float cty = (nodes[n1].y + nodes[n2].y + nodes[i].y) / 3.0f;

                if (this->point_in_mesh(ctx, cty)) {
                    min_dist_sq = d_sq;
                    best_n3 = i;
                }
            }
        }

        if (best_n3 != -1) {
            recreated_tr.emplace_back(n1, n2, best_n3);
        }
    }
    for (auto& tr : recreated_tr) {
        this->triangles.push_back(tr);
    }
}

void geo::Mesh::create_mesh(float spacing, float alfa, float beta)
{
    std::cout << "==================================================\n";
    std::cout << "Generating mesh...\n";
    std::cout << "==================================================\n\n";

    if(this->mesh_created){
        this->nodes.clear();
        this->triangles.clear();
        this->edges.clear();
        this->nodes = this->initial_bc_nodes;
        this->edges = this->initial_edges;
        this->mesh_created = false;
    } else {
        this->initial_edges = this->edges;
    }

    std::cout << "Interpolating boundary points...\n";
    this->interpolate_bc_points(spacing);

    this->create_nodes(alfa, beta);

    std::cout << "Removing outside triangles...\n";
    this->remove_outside_triangles();

    std::cout << "Repairing missing boundaries...\n";
    this->recreate_mising_triangles();

    this->mesh_created = true;
    std::cout << "Mesh complete. Triangles: " << this->triangles.size()
              << ", Nodes: " << this->nodes.size() << "\n";
    std::cout << "==================================================\n";
}

//dobieranie warunków brzegowych
int geo::get_node_clicked(const std::vector<geo::Node>& nodes, float x_pos, float y_pos) {
    constexpr float clickRadius = 0.1f;

    for (int i = 0; i < static_cast<int>(nodes.size()); i++) {
        if (std::abs(nodes[i].x - x_pos) < clickRadius && std::abs(nodes[i].y - y_pos) < clickRadius) {
            return i;
        }
    }

    return -1;
}

float geo::distSq(float x1, float y1, float x2, float y2) {
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

int geo::get_edge_clicked(const std::vector<geo::Edge>& edges, const std::vector<geo::Node>& nodes, float x_pos, float y_pos) {
    constexpr float clickThreshold = 0.1f;
    constexpr float clickThresholdSq = clickThreshold * clickThreshold;

    for (int i = 0; i < static_cast<int>(edges.size()); i++) {
        const auto& n1 = nodes[edges[i].node_ids[0]];
        const auto& n2 = nodes[edges[i].node_ids[1]];

        const float x1 = n1.x; float y1 = n1.y;
        const float x2 = n2.x; float y2 = n2.y;

        float dx = x2 - x1;
        float dy = y2 - y1;

        if (dx == 0 && dy == 0) {
            if (distSq(x_pos, y_pos, x1, y1) <= clickThresholdSq) return i;
            continue;
        }

        float t = ((x_pos - x1) * dx + (y_pos - y1) * dy) / (dx * dx + dy * dy);

        t = std::max(0.0f, std::min(1.0f, t));

        float closestX = x1 + t * dx;
        float closestY = y1 + t * dy;

        float dist = distSq(x_pos, y_pos, closestX, closestY);

        if (dist <= clickThresholdSq) {
            return i;
        }
    }

    return -1;
}

