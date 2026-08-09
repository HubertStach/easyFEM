#pragma once

#include "raylib.h"
#include "imgui.h"
#include <vector>
#include <string>
#include <future>

#include "mes/mes.h"
#include "mesh/geometry.h"
#include "visual/visualisation.h"

struct MainWindow {
    const int screenWidth = 1280;
    const int screenHeight = 720;
    float toolbarWidth = 300.0f;

    Camera2D camera = { 0 };
    float camera_base_zoom = 80.0f;
    float mouseSensitivity = 8.0f;

    bool showXY = true;
    bool show_grid_bool = true;
    bool creatingMesh = false;
    float spacing = 1.0f;
    float mesh_alfa = 0.7f;
    float mesh_beta = 0.5f;
    bool mesh_created = false;
    bool problem_solved = false;

    // Solver & Viz
    geo::Mesh mesh;
    Visualisation vis;
    Fem::GlobalData configuration;
    bool loading_visual = false;
    std::future<void> fem_future;
    bool solving = false;
    bool auto_play = false;
    bool display_nodes = true;
    bool display_triangles = true;
    bool display_quads = true;
    int current_solver = 1;
    std::string solver_type_str = "implicit_euler";

    // BCs
    int bc_edge_clicked = -1;
    std::vector<int> selected_edges;
    float bc_flux = 0.0f;
    float bc_alfa = 0.0f;
    float bc_text = 0.0f;
    float bc_dir_temp = 0.0f;

    MainWindow();
    void Run();

private:
    void ApplyModernDarkStyle();
    void DrawTabMesh();
    void DrawTabBCs();
    void DrawTabSolver();
    void DrawTabVisualisation();
};