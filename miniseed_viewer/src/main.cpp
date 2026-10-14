// miniseed_viewer: quick-look viewer for miniSEED files.
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

#include <GLFW/glfw3.h>

#include "app.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

#include "../third_party/fonts/roboto_medium.h"

static std::vector<std::string> g_dropped;

static void on_drop(GLFWwindow*, int count, const char** paths) {
    for (int i = 0; i < count; i++) g_dropped.push_back(paths[i]);
}

// Test hook: if MSV_SCRIPT names a text file, its lines are played back as mouse
// input, one step per frame, so the program can be driven without a person.
//   frames N | move X Y | down | up | rdown | rup | wheel DY | click X Y | shot FILE.ppm | quit
struct Script {
    std::vector<std::string> lines;
    size_t pos = 0;
    int wait = 0;
    bool active = false;
};

static void save_ppm(GLFWwindow* window, const std::string& path) {
    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    std::vector<unsigned char> px((size_t)w * h * 3), flipped((size_t)w * h * 3);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, px.data());
    for (int y = 0; y < h; y++) std::memcpy(&flipped[(size_t)y * w * 3], &px[(size_t)(h - 1 - y) * w * 3], (size_t)w * 3);
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return;
    std::fprintf(f, "P6\n%d %d\n255\n", w, h);
    std::fwrite(flipped.data(), 1, flipped.size(), f);
    std::fclose(f);
}

// Called once per frame, after drawing and before the buffers are swapped.
static void run_script(Script& sc, ImGuiIO& io, GLFWwindow* window) {
    if (!sc.active) return;
    if (sc.wait > 0) {
        sc.wait--;
        return;
    }
    while (sc.pos < sc.lines.size()) {
        char cmd[32] = "", arg[512] = "";
        float a = 0, b = 0;
        const std::string line = sc.lines[sc.pos++];
        if (std::sscanf(line.c_str(), "%31s %511s", cmd, arg) < 1) continue;
        std::sscanf(line.c_str(), "%*s %f %f", &a, &b);
        std::string c = cmd;
        if (c == "frames") { sc.wait = (int)a; return; }
        else if (c == "shot") save_ppm(window, arg);
        else if (c == "move") io.AddMousePosEvent(a, b);
        else if (c == "down") { io.AddMouseButtonEvent(0, true); return; }
        else if (c == "up") { io.AddMouseButtonEvent(0, false); return; }
        else if (c == "rdown") { io.AddMouseButtonEvent(1, true); return; }
        else if (c == "rup") { io.AddMouseButtonEvent(1, false); return; }
        else if (c == "wheel") { io.AddMouseWheelEvent(0, a); return; }
        else if (c == "click") {  // move, then press and release on separate frames
            io.AddMousePosEvent(a, b);
            sc.lines.insert(sc.lines.begin() + sc.pos, {"frames 1", "down", "up"});
            return;
        }
        else if (c == "quit") { glfwSetWindowShouldClose(window, 1); sc.active = false; return; }
    }
    sc.active = false;
}

static void on_error(int code, const char* msg) { std::fprintf(stderr, "GLFW error %d: %s\n", code, msg); }

int main(int argc, char** argv) {
    std::vector<std::string> files;
    for (int i = 1; i < argc; i++) {
        if (!std::strcmp(argv[i], "-h") || !std::strcmp(argv[i], "--help")) {
            std::printf(
                "usage: miniseed_viewer file1.mseed [file2.mseed ...]\n"
                "Files can also be dropped onto the window.\n");
            return 0;
        }
        files.push_back(argv[i]);
    }

    glfwSetErrorCallback(on_error);
    if (!glfwInit()) return 1;
#if defined(__APPLE__)
    const char* glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#else
    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
#endif
    GLFWwindow* window = glfwCreateWindow(1400, 850, "miniseed_viewer", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    glfwSetDropCallback(window, on_drop);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.IniFilename = nullptr;  // do not leave window-layout files lying around
    float xs = 1, ys = 1;
    glfwGetWindowContentScale(window, &xs, &ys);
    io.Fonts->AddFontFromMemoryCompressedBase85TTF(roboto_medium_compressed_data_base85, 15.0f * xs);
    ImGui::StyleColorsDark();
    ImGui::GetStyle().ScaleAllSizes(xs);
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    {
        Script script;
        if (const char* sp = std::getenv("MSV_SCRIPT")) {
            std::FILE* f = std::fopen(sp, "r");
            char buf[600];
            while (f && std::fgets(buf, sizeof buf, f)) {
                std::string l = buf;
                while (!l.empty() && (l.back() == '\n' || l.back() == '\r')) l.pop_back();
                if (!l.empty()) script.lines.push_back(l);
            }
            if (f) std::fclose(f);
            script.active = true;
        }
        int busy_frames = 0;
        msv::App app(files);
        while (!glfwWindowShouldClose(window)) {
            // Redraw quickly while something is happening, rest otherwise (saves battery).
            if (script.active)
                glfwPollEvents();
            else
                glfwWaitEventsTimeout(busy_frames > 0 ? 0.01 : 0.2);
            if (!g_dropped.empty()) {
                app.add_files(g_dropped);
                g_dropped.clear();
            }
            io.FontGlobalScale = app.cfg.ui_scale;

            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::SetNextWindowSize(io.DisplaySize);
            ImGui::Begin("##root", nullptr,
                         ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings |
                             ImGuiWindowFlags_NoBringToFrontOnFocus);
            app.frame();
            ImGui::End();
            ImGui::Render();

            int w, h;
            glfwGetFramebufferSize(window, &w, &h);
            glViewport(0, 0, w, h);
            glClearColor(0.08f, 0.08f, 0.09f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            run_script(script, io, window);
            glfwSwapBuffers(window);

            bool active = io.MouseDelta.x != 0 || io.MouseDelta.y != 0 || io.MouseWheel != 0 || io.MouseWheelH != 0 ||
                          io.InputQueueCharacters.Size > 0 || ImGui::IsAnyMouseDown() || app.busy();
            for (int k = ImGuiKey_NamedKey_BEGIN; k < ImGuiKey_NamedKey_END && !active; k++)
                if (ImGui::IsKeyDown((ImGuiKey)k)) active = true;
            busy_frames = active ? 6 : std::max(0, busy_frames - 1);
        }
        app.save_all();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
