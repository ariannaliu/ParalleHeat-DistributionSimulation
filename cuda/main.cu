#include <graphic/graphic.hpp>
#include <imgui_impl_sdl.h>
#include <cstring>
#include <chrono>
#include <hdist/hdist.hpp>
#include <iostream>
#include <random>
#include <utility>
#include <time.h>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>


using namespace std;

template<typename ...Args>
void UNUSED(Args &&... args [[maybe_unused]]) {}



__device__
int getBlockId() {
  return blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
}

__device__
int getLocalThreadId() {
  return (threadIdx.z * (blockDim.x * blockDim.y)) + (threadIdx.y * blockDim.x) + threadIdx.x;
}

__device__
int getThreadId() {
  int blockId = getBlockId();
  int localThreadId = getLocalThreadId();
  return blockId * (blockDim.x * blockDim.y * blockDim.z) + localThreadId;
}




// __shared__ bool local_stable = true;
// __device__ bool *local_stable;

__host__ __device__
int get_index_in_array(size_t i, size_t j, int length){
    return i * length + j;
}


__device__
int get_buffer(int current_buffer){
    if (current_buffer == 0){
        return 0;
    }
    return 1;
}

__device__
int get_buffer_alt(int current_buffer){
    if (current_buffer == 0){
        return 1;
    }
    return 0;
}

__device__
int switch_buffer(int current_buffer){
    if (current_buffer == 0){
        return 1;
    }
    return 0;
}

__device__
void update_single(int i, int j, 
                   int room_size, float block_size,
                int source_x, int source_y,
                float source_temp, float border_temp,
                float tolerance, float sor_constant,
                int algo,
                double *data0, double *data1, int length, int current_buffer,
                bool & stable, double & temp){
    if (i == 0 || j == 0 || i == room_size - 1 || j == room_size - 1) {
        temp = border_temp;
    } else if (i == source_x && j == source_y) {
        temp = source_temp;
    } else {
        double sum;
        if(current_buffer == 0){
            sum = data0[get_index_in_array(i + 1, j, length)] + data0[get_index_in_array(i - 1, j, length)]
                        + data0[get_index_in_array(i, j+1, length)] + data0[get_index_in_array(i, j-1, length)];
        }else{
            sum = data1[get_index_in_array(i + 1, j, length)] + data1[get_index_in_array(i - 1, j, length)]
                        + data1[get_index_in_array(i, j+1, length)] + data1[get_index_in_array(i, j-1, length)];
        }
        if(algo == 0){
            temp =  0.25 * sum;
        }else{
            if(current_buffer == 0){
                temp = data0[get_index_in_array(i, j, length)] + (1.0 / sor_constant)*(sum - 4.0 * data0[get_index_in_array(i, j, length)]);
            }else{
                temp = data1[get_index_in_array(i, j, length)] + (1.0 / sor_constant)*(sum - 4.0 * data1[get_index_in_array(i, j, length)]);
            }
        }
    }
    if(current_buffer == 0){
        stable = std::fabs(data0[get_index_in_array(i, j, length)] - temp) < tolerance;   
    }else{
        stable = std::fabs(data1[get_index_in_array(i, j, length)] - temp) < tolerance;   
    }
}

__global__
void calculate(int room_size, float block_size,
                int source_x, int source_y,
                float source_temp, float border_temp,
                float tolerance, float sor_constant,
                int algo,
                double *data0, double *data1, int length, int current_buffer, int local_size, bool *devicestable){
    bool stabilized = true;
    *devicestable = true;
    
    int thread_id = getThreadId();
    printf("Thread_id is %d \n", thread_id);
    if(algo == 0){
        for (size_t i = (thread_id)*local_size; i < (thread_id+1)*local_size; ++i) {
            for (size_t j = 0; j < room_size; ++j) {
                if(get_index_in_array(i, j, length)<room_size*room_size){

                    bool t_stable;
                    double t_temp;
                    update_single(i,j,room_size, block_size, 
                                source_x, source_y, source_temp, border_temp,
                                tolerance, sor_constant, algo, data0, data1, length, current_buffer,
                                t_stable, t_temp);
                    stabilized &= t_stable;
                    if(current_buffer == 0){
                        data1[get_index_in_array(i, j, length)] = t_temp;
                    }else{
                        data0[get_index_in_array(i, j, length)] = t_temp;
                    }
                }
            }
        }
        current_buffer = switch_buffer(current_buffer);
    }else{
        for (auto k : {0, 1}) {
            for (size_t i = (thread_id)*local_size; i < min(room_size, (thread_id+1)*local_size); ++i) {
                for (size_t j = 0; j < room_size; j++) {
                    if(get_index_in_array(i, j, length)<room_size*room_size){
                        bool t_stable;
                        double t_temp;
                        if (k == ((i + j) & 1)) {
                            update_single(i,j,room_size, block_size, 
                                                        source_x, source_y, source_temp, border_temp,
                                                        tolerance, sor_constant, algo, data0, data1, length, current_buffer,
                                                        t_stable, t_temp);
                            stabilized &= t_stable;

                            if(current_buffer == 0){
                                data1[get_index_in_array(i, j, length)] = t_temp;
                            }else{
                                data0[get_index_in_array(i, j, length)] = t_temp;
                            }
                        } else {
                            // grid[{alt, i, j}] = grid[{i, j}];
                            if(current_buffer == 0){
                                data1[get_index_in_array(i, j, length)] = data0[get_index_in_array(i, j, length)];
                            }else{
                                data0[get_index_in_array(i, j, length)] = data1[get_index_in_array(i, j, length)];
                            }
                        }
                    }
                }
            }
            current_buffer = switch_buffer(current_buffer);
        }
    }
    *devicestable &= stabilized;
}





















ImColor temp_to_color(double temp) {
    auto value = static_cast<uint8_t>(temp / 100.0 * 255.0);
    return {value, 0, 255 - value};
}

void init_Grid(size_t size, double border_temp,double source_temp,
                size_t x, size_t y, double *data0, double *data1, 
                int length){
    for (size_t i = 0; i < length; ++i){
        for (size_t j = 0; j < length; ++j){
            int index = get_index_in_array(i, j, length);
            if (i == 0 || j == 0 || i == length - 1 || j == length - 1) {
                data0[index] = border_temp;
            }else if (i == x && j == y){
                data0[index] = source_temp;
            }else{
                data0[index] = 0;
            }
            
        }
    }
}

int main(int argc, char **argv) {
    // UNUSED(argc, argv);
    int thread_num;
    if (argc < 2){
        // if user did not provide the size of array
        // the defualt value is set to be 100
        thread_num = 4;
    }else{
        thread_num = atoi(argv[1]);
    }
    bool first = true;
    bool finished = false;
    // static hdist::State current_state, last_state;
    int room_size_C = 300, room_size_S = 300;
    float block_size_C = 2, block_size_S = 2;
    int source_x_C = room_size_C / 2, source_x_S = room_size_S / 2;
    int source_y_C = room_size_C / 2, source_y_S = room_size_S / 2;
    float source_temp_C = 100, source_temp_S = 100;
    float border_temp_C = 36, border_temp_S = 36;
    float tolerance_C = 0.02, tolerance_S = 0.02;
    float sor_constant_C = 4.0, sor_constant_S = 4.0;
    int algo_C = 0, algo_S = 0;


    static std::chrono::high_resolution_clock::time_point begin, end;
    static const char* algo_list[2] = { "jacobi", "sor" };
    graphic::GraphicContext context{"Assignment 4"};

    double * Hostdata0 = new double[room_size_C*room_size_C];
    double * Hostdata1 = new double[room_size_C*room_size_C];
    bool *Hostdevicestable = new bool;

    double *data0;
    cudaMalloc(&data0, sizeof(double) * room_size_C*room_size_C);
    double *data1;
    cudaMalloc(&data1, sizeof(double) * room_size_C*room_size_C);
    bool *devicestable;
    cudaMalloc(&devicestable, sizeof(bool));

    int length = room_size_C;
    int current_buffer = 0;
    init_Grid(
            static_cast<size_t>(room_size_C),
            border_temp_C, source_temp_C,
            static_cast<size_t>(source_x_C), static_cast<size_t>(source_y_C),
            Hostdata0, Hostdata1, length);
    
    cudaMemcpy(data0, Hostdata0, sizeof(double)*room_size_C*room_size_C, cudaMemcpyHostToDevice);
    cudaMemcpy(data1, Hostdata1, sizeof(double)*room_size_C*room_size_C, cudaMemcpyHostToDevice);
    
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
        ImGui::DragInt("Room Size", &room_size_C, 10, 200, 1600, "%d");
        ImGui::DragFloat("Block Size", &block_size_C, 0.01, 0.1, 10, "%f");
        ImGui::DragFloat("Source Temp", &source_temp_C, 0.1, 0, 100, "%f");
        ImGui::DragFloat("Border Temp", &border_temp_C, 0.1, 0, 100, "%f");
        ImGui::DragInt("Source X", &source_x_C, 1, 1, room_size_C - 2, "%d");
        ImGui::DragInt("Source Y", &source_y_C, 1, 1, room_size_C - 2, "%d");
        ImGui::DragFloat("Tolerance", &tolerance_C, 0.01, 0.01, 1, "%f");
        ImGui::ListBox("Algorithm", reinterpret_cast<int *>(&algo_C), algo_list, 2);

        if (algo_C == 1) {
            ImGui::DragFloat("Sor Constant", &sor_constant_C, 0.01, 0.0, 20.0, "%f");
        }

        if (room_size_C != room_size_S) {
            delete[] Hostdata0;
            delete[] Hostdata1;
            Hostdata0 = new double[room_size_C*room_size_C];
            Hostdata1 = new double[room_size_C*room_size_C];
            first = true;
            length = room_size_C;
            current_buffer = 0;
            init_Grid(
                    static_cast<size_t>(room_size_C),
                    border_temp_C,
                    source_temp_C,
                    static_cast<size_t>(source_x_C),
                    static_cast<size_t>(source_y_C), Hostdata0, Hostdata1, length);
            cudaMemcpy(data0, Hostdata0, sizeof(double)*(length*length), cudaMemcpyHostToDevice);
            cudaMemcpy(data1, Hostdata1, sizeof(double)*(length*length), cudaMemcpyHostToDevice);
        }

        // if (current_state != last_state) {
        //     last_state = current_state;
        //     finished = false;
        // }
        if (room_size_C != room_size_S){
            room_size_S = room_size_C;
            finished = false;
        }
        if (block_size_C != block_size_S){
            block_size_S = block_size_C;
            finished = false;
        }
        if (source_x_C != source_x_S){
            source_x_S = source_x_C;
            finished = false;
        }
        if (source_y_C != source_y_S){
            source_y_S = source_y_C;
            finished = false;
        }
        if (source_temp_C != source_temp_S){
            source_temp_S = source_temp_C;
            finished = false;
        }
        if (border_temp_C != border_temp_S){
            border_temp_S = border_temp_C;
            finished = false;
        }
        if (tolerance_C != tolerance_S){
            tolerance_S = tolerance_C;
            finished = false;
        }
        if (sor_constant_C != sor_constant_S){
            sor_constant_S = sor_constant_C;
            finished = false;
        }
        if(algo_C != algo_S){
            algo_S = algo_C;
            finished = false;
        }

        if (first) {
            first = false;
            finished = false;
            begin = std::chrono::high_resolution_clock::now();
        }

        int local_size = ceil(double(room_size_C)/thread_num);
        bool Hoststable = true;
        if (!finished) {
            // std::cout<<"Length of data1 "<<sizeof(data1)<<std::endl;
            // std::cout<<"Here1"<<std::endl;
            
            calculate<<<thread_num,1>>>(room_size_C, block_size_C,
                    source_x_C, source_y_C,
                    source_temp_C, border_temp_C,
                    tolerance_C, sor_constant_C,
                    algo_C,
                    data0, data1, length, current_buffer, local_size, devicestable);
            cudaDeviceSynchronize();
            cudaMemcpy(Hostdata0, data0, sizeof(double)*room_size_C*room_size_C, cudaMemcpyDeviceToHost);
            cudaMemcpy(Hostdata1, data1, sizeof(double)*room_size_C*room_size_C, cudaMemcpyDeviceToHost);
            cudaMemcpy(Hostdevicestable, devicestable, sizeof(bool), cudaMemcpyDeviceToHost);
            
            finished = *Hostdevicestable;
            if (finished) end = std::chrono::high_resolution_clock::now();
        } else {
            ImGui::Text("stabilized in %ld ns", std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
        }

        const ImVec2 p = ImGui::GetCursorScreenPos();
        float x = p.x + block_size_C, y = p.y + block_size_C;
        // std::cout<<"Here5"<<std::endl;
        for (size_t i = 0; i < room_size_C; ++i) {
            for (size_t j = 0; j < room_size_C; ++j) {
                double temp;
                if(current_buffer == 0){
                    temp = Hostdata0[get_index_in_array(i, j, length)];
                }else{
                    temp = Hostdata1[get_index_in_array(i, j, length)];
                }
                auto color = temp_to_color(temp);
                draw_list->AddRectFilled(ImVec2(x, y), ImVec2(x + block_size_C, y + block_size_C), color);
                y += block_size_C;
            }
            x += block_size_C;
            y = p.y + block_size_C;
        }
        ImGui::End();
    });
}
