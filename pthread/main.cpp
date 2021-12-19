#include <graphic/graphic.hpp>
#include <imgui_impl_sdl.h>
#include <cstring>
#include <chrono>
#include <hdist/hdist.hpp>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>

template<typename ...Args>
void UNUSED(Args &&... args [[maybe_unused]]) {}

using namespace std;

// pthread_mutex_t mutex_p;
pthread_cond_t mutex_threshold;
#define THREADS_NUM 4

struct Arguments{
    hdist::State * state;
    hdist::Grid * grid;
    bool * stable_flag;
    int pid;
    int local_bodies;
    int bodies;
};

ImColor temp_to_color(double temp) {
    auto value = static_cast<uint8_t>(temp / 100.0 * 255.0);
    return {value, 0, 255 - value};
}

void *local_process(void *arg_ptr){
    auto arguments = static_cast<Arguments *>(arg_ptr);
    // pthread_mutex_lock(&mutex_p);
    bool local_stable;
    local_stable = hdist::calculate(*(arguments->state), *(arguments->grid),  
                                        arguments->pid+1, arguments->local_bodies,
                                        arguments-> bodies);
    // pthread_mutex_unlock(&mutex_p);   
    *(arguments->stable_flag) &= local_stable;                            
    delete arguments;
    return nullptr;
}


int main(int argc, char **argv) {
    UNUSED(argc, argv);
    bool first = true;
    bool finished = false;
    static hdist::State current_state, last_state;
    static std::chrono::high_resolution_clock::time_point begin, end;
    static const char* algo_list[2] = { "jacobi", "sor" };
    graphic::GraphicContext context{"Assignment 4"};
    auto grid = hdist::Grid{
            static_cast<size_t>(current_state.room_size),
            current_state.border_temp,
            current_state.source_temp,
            static_cast<size_t>(current_state.source_x),
            static_cast<size_t>(current_state.source_y)};
    context.run([&](graphic::GraphicContext *context [[maybe_unused]], SDL_Window *) {
        auto io = ImGui::GetIO();
        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f));
        ImGui::SetNextWindowSize(io.DisplaySize);
        ImGui::Begin("Assignment 4", nullptr,
                     ImGuiWindowFlags_NoMove
                     | ImGuiWindowFlags_NoCollapse
                     | ImGuiWindowFlags_NoTitleBar
                     | ImGuiWindowFlags_NoResize);
        ImDrawList *draw_list = ImGui::GetWindowDrawList();
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate,
                    ImGui::GetIO().Framerate);
        ImGui::DragInt("Room Size", &current_state.room_size, 10, 200, 1600, "%d");
        ImGui::DragFloat("Block Size", &current_state.block_size, 0.01, 0.1, 10, "%f");
        ImGui::DragFloat("Source Temp", &current_state.source_temp, 0.1, 0, 100, "%f");
        ImGui::DragFloat("Border Temp", &current_state.border_temp, 0.1, 0, 100, "%f");
        ImGui::DragInt("Source X", &current_state.source_x, 1, 1, current_state.room_size - 2, "%d");
        ImGui::DragInt("Source Y", &current_state.source_y, 1, 1, current_state.room_size - 2, "%d");
        ImGui::DragFloat("Tolerance", &current_state.tolerance, 0.01, 0.01, 1, "%f");
        ImGui::ListBox("Algorithm", reinterpret_cast<int *>(&current_state.algo), algo_list, 2);

        if (current_state.algo == hdist::Algorithm::Sor) {
            ImGui::DragFloat("Sor Constant", &current_state.sor_constant, 0.01, 0.0, 20.0, "%f");
        }

        if (current_state.room_size != last_state.room_size) {
            grid = hdist::Grid{
                    static_cast<size_t>(current_state.room_size),
                    current_state.border_temp,
                    current_state.source_temp,
                    static_cast<size_t>(current_state.source_x),
                    static_cast<size_t>(current_state.source_y)};
            first = true;
        }

        if (current_state != last_state) {
            last_state = current_state;
            finished = false;
        }

        if (first) {
            first = false;
            finished = false;
            begin = std::chrono::high_resolution_clock::now();
        }

        if (!finished) {
            // finished = hdist::calculate(current_state, grid);
            bool stable_flag =  true;
            pthread_attr_t attr;
            pthread_attr_init(&attr);
            pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

            pthread_t threads[THREADS_NUM];
            int local_size = ceil(double(current_state.room_size)/THREADS_NUM);

            for(int m=0; m<THREADS_NUM; m++){
                    pthread_create(&threads[m], nullptr, local_process, new Arguments{
                        .state = &current_state,
                        .grid = &grid,
                        .stable_flag = &stable_flag,
                        .pid = m,
                        .local_bodies = local_size,
                        .bodies = current_state.room_size
                    });
            }
            for(auto & n : threads){
                pthread_join(n, nullptr);
            }
            pthread_attr_destroy(&attr);
            finished = stable_flag;

            if (finished) end = std::chrono::high_resolution_clock::now();
        } else {
            ImGui::Text("stabilized in %ld ns", std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
        }

        const ImVec2 p = ImGui::GetCursorScreenPos();
        float x = p.x + current_state.block_size, y = p.y + current_state.block_size;
        for (size_t i = 0; i < current_state.room_size; ++i) {
            for (size_t j = 0; j < current_state.room_size; ++j) {
                auto temp = grid[{i, j}];
                auto color = temp_to_color(temp);
                draw_list->AddRectFilled(ImVec2(x, y), ImVec2(x + current_state.block_size, y + current_state.block_size), color);
                y += current_state.block_size;
            }
            x += current_state.block_size;
            y = p.y + current_state.block_size;
        }
        ImGui::End();
    });
}
