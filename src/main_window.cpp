#include "../../_deps/imgui-src/imgui.h"
#include "../../_deps/rlimgui-src/rlImGui.h"
#include "../../_deps/rlimgui-src/rlImGuiColors.h"

#include "raylib.h"
#include "../include/rlgl.h"
#include "../include/raymath.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <thread>

#include "app_ui/axis.h"
#include "app_ui/grid.h"

#include "mes/mes.h"
#include "mesh/geometry.h"
#include "saving_data.h"
#include "visual/visualisation.h"
#include "main_window.h"

#include <iomanip>

MainWindow::MainWindow() {
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    InitWindow(screenWidth, screenHeight, "easyFEM - Pipeline UI");
    SetTargetFPS(60);
    rlImGuiSetup(true);

    ApplyModernDarkStyle();

    camera.offset = { screenWidth / 2.0f + 0.5f * toolbarWidth, screenHeight / 2.0f };
    camera.target = { 0.0f, 0.0f };
    camera.zoom = camera_base_zoom;

    clean_vtu_files();
}

void MainWindow::ApplyModernDarkStyle() {
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 6.0f;
    style.ChildRounding = 4.0f;
    style.FrameRounding = 4.0f;
    style.TabRounding = 4.0f;
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0.08f, 0.08f, 0.09f, 1.00f);
    style.Colors[ImGuiCol_ChildBg] = ImVec4(0.11f, 0.11f, 0.13f, 1.00f);
}

void MainWindow::Run() {
    while (!WindowShouldClose()) {
        int currentWidth = GetRenderWidth();
        int currentHeight = GetRenderHeight();

        Vector2 worldPos = GetScreenToWorld2D(GetMousePosition(), camera);
        worldPos.x = round(worldPos.x / 0.1f) * 0.1f;
        worldPos.y = round(worldPos.y / 0.1f) * 0.1f;
        Vector2 screenPos = GetWorldToScreen2D(worldPos, camera);

        BeginDrawing();
        ClearBackground({20, 20, 25, 255}); // Tło viewportu

        rlImGuiBegin();

        // ROOT WINDOW
        ImGui::SetNextWindowPos({ 0,0 }, ImGuiCond_Always);
        ImGui::SetNextWindowSize({ (float)currentWidth, (float)currentHeight }, ImGuiCond_Always);
        ImGui::Begin("Root", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoBackground);

        // LEWY PANEL STEROWANIA
        ImGui::BeginChild("Toolbar", { toolbarWidth, 0 }, true);
        ImGui::TextColored({ 0.3f, 0.6f, 1.0f, 1.0f }, "SIMULATION STEPS");
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::BeginTabBar("Steps")) {
            if (ImGui::BeginTabItem("1. Mesh")) { DrawTabMesh(); ImGui::EndTabItem(); }
            if (ImGui::BeginTabItem("2. BCs")) { DrawTabBCs(); ImGui::EndTabItem(); }
            if (ImGui::BeginTabItem("3. Solver")) { DrawTabSolver(); ImGui::EndTabItem(); }
            if (ImGui::BeginTabItem("4. Results")) { DrawTabVisualisation(); ImGui::EndTabItem(); }
            ImGui::EndTabBar();
        }

        ImGui::SetCursorPosY(ImGui::GetWindowHeight() - 110);
        ImGui::Separator();
        ImGui::Checkbox("Show XY Axis", &showXY);
        ImGui::Checkbox("Show Grid", &show_grid_bool);
        if (ImGui::Button("Reset Camera View", ImVec2(-1, 0))) { camera.zoom = camera_base_zoom; camera.target = { 0,0 }; }

        ImGui::EndChild();

        ImGui::SameLine();

        // GŁÓWNY PANEL VIEWPORTU
        ImGui::BeginChild("MainPanel", { 0,0 }, false, ImGuiWindowFlags_NoBackground);
        ImVec2 panelMin = ImGui::GetCursorScreenPos();
        ImVec2 panelSize = ImGui::GetContentRegionAvail();
        bool panelHover = ImGui::IsWindowHovered();
        ImGui::EndChild();

        ImGui::End(); // End Root
        rlImGuiEnd();

        // --- OBSŁUGA KAMERY (FIX: Działa zawsze, niezależnie od trybu) ---
        camera.offset = { panelMin.x + panelSize.x / 2.0f, panelMin.y + panelSize.y / 2.0f };
        if (panelHover) {
            // Przesuwanie widoku (LPM) - teraz działa ZAWSZE
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                Vector2 delta = GetMouseDelta();
                camera.target.x -= delta.x / camera.zoom;
                camera.target.y -= delta.y / camera.zoom;
            }
            // Zoom (Kółko)
            float wheel = GetMouseWheelMove();
            if (wheel != 0.0f) {
                camera.zoom += wheel * mouseSensitivity;
                if (camera.zoom < 5.0f) camera.zoom = 5.0f;
            }
        }

        // --- RYSOWANIE SYMULACJI ---
        rlEnableScissorTest();
        rlScissor((int)panelMin.x, currentHeight - (int)(panelMin.y + panelSize.y), (int)panelSize.x, (int)panelSize.y);

        BeginMode2D(camera);
            if (show_grid_bool) show_adaptive_grid_dim(camera);

            mesh.draw_edges();

            // Zaznaczone krawędzie (MAGENTA)
            if (bc_edge_clicked != -1 && !selected_edges.empty()) {
                for (int id : selected_edges) {
                    int n1 = mesh.edges[id].node_ids[0], n2 = mesh.edges[id].node_ids[1];
                    DrawLineEx({ (float)mesh.nodes[n1].x, (float)mesh.nodes[n1].y },
                               { (float)mesh.nodes[n2].x, (float)mesh.nodes[n2].y },
                               3.0f * (1.0f / camera.zoom), MAGENTA);
                }
            }

            if (loading_visual && vis.solved) {
                std::vector<double> current_temps;
                for (size_t i = 0; i < mesh.nodes.size(); i++)
                    current_temps.push_back(i < (size_t)vis.temps_matrix.rows() ? vis.temps_matrix(i, vis.current_step) : 0.0);
                mesh.draw_tr_grad(current_temps, vis.max_temp, vis.min_temp);
                mesh.draw_q_grad(current_temps, vis.max_temp, vis.min_temp);
            }

            if (display_triangles || !loading_visual) mesh.draw_tr(WHITE);
            if (display_quads || !loading_visual) mesh.draw_q(WHITE);
            if (display_nodes || !loading_visual) mesh.draw_nodes(3 * (1 / camera.zoom));
        EndMode2D();

        // --- LEGENDA ---
        if (loading_visual && vis.solved) {
            float legW = 25, legH = std::min(400.0f, panelSize.y - 100.0f);
            float legX = panelMin.x + panelSize.x - 110, legY = panelMin.y + (panelSize.y - legH) / 2.0f;
            DrawRectangleGradientV((int)legX, (int)legY, (int)legW, (int)legH/2, RED, WHITE);
            DrawRectangleGradientV((int)legX, (int)(legY + legH/2), (int)legW, (int)legH/2, WHITE, BLUE);
            DrawRectangleLines((int)legX, (int)legY, (int)legW, (int)legH, RAYWHITE);
            for (int i = 0; i < 6; i++) {
                float t = 1.0f - (i / 5.0f);
                float val = vis.min_temp + t * (vis.max_temp - vis.min_temp);
                DrawTextEx(GetFontDefault(), TextFormat("%.1f", val), { legX + 35, legY + i*(legH/5) - 7 }, 16, 1, RAYWHITE);
            }
        }

        // --- INPUT: KLIKANIE KRAWĘDZI I TWORZENIE SIATKI ---
        // Wybieranie krawędzi (Tylko gdy nie tworzymy siatki)
        if (panelHover && IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && !creatingMesh) {
            int hit = geo::get_edge_clicked(mesh.edges, mesh.nodes, worldPos.x, worldPos.y);
            if (hit != -1) {
                bc_edge_clicked = hit;
                if (IsKeyDown(KEY_LEFT_SHIFT)) selected_edges = mesh.get_continuous_edges(hit);
                else { selected_edges.clear(); selected_edges.push_back(hit); }
                if (mesh.edges[hit].bc_edge.initialised) {
                    bc_flux = mesh.edges[hit].bc_edge.flux; bc_alfa = mesh.edges[hit].bc_edge.alfa; bc_text = mesh.edges[hit].bc_edge.t_ext;
                }
            }
        }

        // Tryb tworzenia siatki
        if (creatingMesh) {
            DrawCircleV(screenPos, 4, RED);
            DrawTextEx(GetFontDefault(), TextFormat("[%0.1f, %0.1f]", worldPos.x, worldPos.y), { screenPos.x + 10, screenPos.y - 20 }, 18, 1, RAYWHITE);
            if (IsKeyPressed(KEY_E)) {
                mesh.add_point(worldPos.x, worldPos.y);
            }

            if (IsKeyPressed(KEY_Q)) {
                if(mesh_created){
                    if(problem_solved){
                        vis.solved = false;
                        loading_visual = false;
                        problem_solved = false;
                    }
                    mesh.triangles.clear();
                    mesh.edges.clear();
                    selected_edges.clear();
                    bc_edge_clicked = -1;
                }

                mesh.pop_point();
                mesh_created = false;
            }

        }

        rlDisableScissorTest();
        if (showXY) DrawCoordinateAxes(camera, panelMin, panelSize);

        EndDrawing();
    }
    rlImGuiShutdown();
    CloseWindow();
}

// ================= WIDGETY ZAKŁADEK =================

void MainWindow::DrawTabMesh() {
    ImGui::Spacing();
    if (ImGui::Button("Load .txt File", ImVec2(-1, 0))) {
        mesh.load_mesh_from_txt(); mesh_created = true;
    }
    if (ImGui::Button("Load .inp File", ImVec2(-1, 0))) {
        load_inp_mesh(mesh); mesh.create_edges(); mesh_created = true;
    }
    ImGui::Separator();
    ImGui::Checkbox("Enable Point Placement", &creatingMesh);
    ImGui::TextDisabled("E: Place | Q: Remove last");
    ImGui::Text("Boundary node spacing");
    ImGui::SliderFloat("##", &spacing, 0.1f, 5.0f);

    ImGui::Text("Mesh alfa parameter");
    ImGui::SliderFloat("alfa", &mesh_alfa, 0.1f, 2.0f);

    ImGui::Text("Mesh beta parameter");
    ImGui::SliderFloat("beta", &mesh_beta, 0.1f, 5.0f);

    if (ImGui::Button("GENERATE MESH", ImVec2(-1, 35))) {
        if (mesh.nodes.size() >= 3){
            mesh_created = true;
            mesh.create_mesh(spacing, mesh_alfa, mesh_beta);
        }
    }
    if (ImGui::Button("Clear All", ImVec2(-1, 0))) { mesh.nodes.clear(); mesh.edges.clear(); mesh_created = false; bc_edge_clicked = -1; }
}

void MainWindow::DrawTabBCs() {
    ImGui::Spacing();
    if (!mesh_created) { ImGui::TextWrapped("Generate mesh first."); return; }
    if (bc_edge_clicked != -1) {
        ImGui::TextColored({1,0,1,1}, "Selected: %zu edges", selected_edges.size());
        ImGui::InputFloat("Temperature", &bc_dir_temp);
        ImGui::InputFloat("Flux", &bc_flux);
        ImGui::InputFloat("Alpha", &bc_alfa);
        ImGui::InputFloat("T_ext", &bc_text);
        if (ImGui::Button("Save BC Values", ImVec2(-1, 30))) {
            for (int id : selected_edges) {
                int n1 = mesh.edges[id].node_ids[0], n2 = mesh.edges[id].node_ids[1];
                mesh.nodes[n1].bc.initialised = mesh.nodes[n2].bc.initialised = mesh.edges[id].bc_edge.initialised = true;
                mesh.nodes[n1].bc.flux = mesh.nodes[n2].bc.flux = mesh.edges[id].bc_edge.flux = bc_flux;
                mesh.nodes[n1].bc.alfa = mesh.nodes[n2].bc.alfa = mesh.edges[id].bc_edge.alfa = bc_alfa;
                mesh.nodes[n1].bc.t_ext = mesh.nodes[n2].bc.t_ext = mesh.edges[id].bc_edge.t_ext = bc_text;
                mesh.nodes[n1].bc.dir_temp = mesh.nodes[n2].bc.dir_temp = mesh.edges[id].bc_edge.dir_temp = bc_dir_temp;
            }
        }
    } else ImGui::TextDisabled("No edge selected. Click one!");
}

void MainWindow::DrawTabSolver() {
    ImGui::Spacing();

    static int solver_mode = 1; // 0 = Stationary (Ustalony), 1 = Transient (Nieustalony)

    ImGui::TextColored({ 0.3f, 0.6f, 1.0f, 1.0f }, "ANALYSIS TYPE");
    ImGui::Separator();
    ImGui::RadioButton("Stationary", &solver_mode, 0);
    ImGui::SameLine();
    ImGui::RadioButton("Transient", &solver_mode, 1);
    ImGui::Separator();

    if (solver_mode == 1) {
        ImGui::Text("Time Parameters:");
        ImGui::InputDouble("Sim Time", &configuration.total_time);
        ImGui::InputDouble("Step dT", &configuration.time_step);
        ImGui::Separator();

        ImGui::Text("Time Integration Scheme:");
        ImGui::RadioButton("Explicit", &current_solver, 0);
        ImGui::RadioButton("Implicit", &current_solver, 1);
        ImGui::RadioButton("Crank-Nicolson", &current_solver, 2);

        if (current_solver == 0) {
            solver_type_str = "explicit_euler";

            if (ImGui::Button("Check stability") && mesh_created) {
                // Build the FEM matrices so the limit reflects element shape AND
                // convection BCs, not just mesh geometry.
                save_fem_data(mesh, configuration);
                Fem::Solution probe("Data/fem_data.txt", "explicit_euler");
                const double dt_max = probe.stability_dt_max();

                std::cout << "Max stable dt = " << std::setprecision(5) << dt_max
                          << " (using dt = " << configuration.time_step << ")\n";
                if (configuration.time_step > dt_max) {
                    std::cout << "[WARNING] dt exceeds stability limit. Simulation will be UNSTABLE!\n";
                }
            }
        }

        else if (current_solver == 1) {
            solver_type_str = "implicit_euler";
        }
        else {
            solver_type_str = "crank-nicolson";
        }

    }
    else {
        solver_type_str = "stationary";
        ImGui::TextDisabled("Time parameters and schemes are disabled in stationary mode.");
    }

    ImGui::Separator();
    ImGui::Text("Initial/Reference temperature");
    ImGui::InputDouble("C", &configuration.init_temperature);

    if (ImGui::CollapsingHeader("Material data")) {
        ImGui::Text("Density");
        ImGui::InputDouble("kg/m^3", &configuration.density);

        ImGui::Text("Conductivity");
        ImGui::InputDouble("W/m*K", &configuration.conductivity);

        ImGui::Text("Specific heat");
        ImGui::InputDouble("J/kg*K", &configuration.specific_heat);
    }

    if (ImGui::Button("RUN SIMULATION", ImVec2(-1, 40))) {
        if (mesh_created) {
            save_fem_data(mesh, configuration);
            std::thread fem_thread(fem_solve, solver_type_str);
            fem_thread.join();
            vis.init_visualisation(mesh);
            problem_solved = true;
            loading_visual = true;
        }
    }
}

void MainWindow::DrawTabVisualisation() {
    ImGui::Spacing();
    if (!problem_solved) { ImGui::TextDisabled("Run solver first."); return; }
    ImGui::Checkbox("Show Temperature Heatmap", &loading_visual);
    ImGui::InputFloat("Range Min", &vis.min_temp);
    ImGui::InputFloat("Range Max", &vis.max_temp);
    ImGui::SliderInt("Current Step", &vis.current_step, 0, (int)vis.time_ids.size() - 1);
    ImGui::Checkbox("Animation On", &auto_play);
    if (auto_play) vis.current_step = (vis.current_step + 1) % vis.time_ids.size();
    ImGui::Separator();
    ImGui::Checkbox("Show Nodes", &display_nodes);
    ImGui::Checkbox("Show Triangles", &display_triangles);
}