#include <graphic/graphic.hpp>
#include <imgui_impl_sdl.h>
#include <cstring>
#include <chrono>
#include <hdist/hdist.hpp>
#include <mpi.h>
#include <iostream>
#include <sys/time.h>
#include <cmath>

MPI_Status status;
MPI_Request request;
int total_size;
template<typename ...Args>
void UNUSED(Args &&... args [[maybe_unused]]) {}

ImColor temp_to_color(double temp) {
    auto value = static_cast<uint8_t>(temp / 100.0 * 255.0);
    return {value, 0, 255 - value};
}

int main(int argc, char **argv) {
    // UNUSED(argc, argv);
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_size);
    int slave_size = total_size - 1;

    if (rank == 0){
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
            
            int room_size = current_state.room_size;
            int local_size = ceil(double(room_size)/slave_size);

            int slave_flage;
            if(!finished){
                slave_flage = 1;
            }else slave_flage = 0;
            for(int i=0; i<slave_size; i++){
                MPI_Send(&slave_flage, 1, MPI_INT, (i+1), 0, MPI_COMM_WORLD);
            }
            //MPI_Barrier(MPI_COMM_WORLD);


            if (!finished) {
                // implement here!!
                // parameter for size
                int para_Size[2];
                para_Size[0] = room_size;
                para_Size[1] = local_size;

                // parameter for current_state which will be passed from master
                int para_State_int[4];
                para_State_int[0] = current_state.room_size;
                para_State_int[1] = current_state.source_x;
                para_State_int[2] = current_state.source_y;
                if(current_state.algo == hdist::Algorithm::Jacobi){
                    para_State_int[3] = 0;
                }else para_State_int[3] = 1;

                float para_State_float[5];
                para_State_float[0] = current_state.block_size;
                para_State_float[1] = current_state.source_temp;
                para_State_float[2] = current_state.border_temp;
                para_State_float[3] = current_state.tolerance;
                para_State_float[4] = current_state.sor_constant;

                // parameter for Grid structure
                int m_length = (grid.data0).size();
                double g_data0[m_length];
                double g_data1[m_length];
                for (int m=0; m<m_length; m++){
                    g_data0[m] = grid.data0[m];
                    g_data1[m] = grid.data1[m];
                }
                int para_Grid[2];
                para_Grid[0] = grid.current_buffer;
                para_Grid[1] = grid.length;

                
                // pass information
                

                for(int i=0; i<slave_size; i++){
                    MPI_Send(&(para_Size[0]), 2, MPI_INT, (i+1), (i+1)*10, MPI_COMM_WORLD);
                    MPI_Send(&(para_State_int[0]), 4, MPI_INT, (i+1), (i+1)*10+1, MPI_COMM_WORLD);
                    MPI_Send(&(para_State_float[0]), 5, MPI_FLOAT, (i+1), (i+1)*10+2, MPI_COMM_WORLD);
                    MPI_Send(&(g_data0[0]), m_length, MPI_DOUBLE, (i+1), (i+1)*10+3, MPI_COMM_WORLD);
                    MPI_Send(&(g_data1[0]), m_length, MPI_DOUBLE, (i+1), (i+1)*10+4, MPI_COMM_WORLD);
                    MPI_Send(&(para_Grid[0]), 2, MPI_INT, (i+1), (i+1)*10+5, MPI_COMM_WORLD);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                

                // MPI_Isend(&s_finished_int, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &request);
                // MPI_Isend(&(s_grid.data0[0]), 1, MPILocalVector, 0, 2, MPI_COMM_WORLD, &request);
                // MPI_Isend(&(s_grid.data1[0]), 1, MPILocalVector, 0, 3, MPI_COMM_WORLD, &request);
                
                int slave_rank;
                bool slave_finish = true;
                int slave_finish_int;
                for (int i=0; i<slave_size; i++){
                    MPI_Recv(&slave_finish_int, 1, MPI_INT,MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,&status);
                    if(slave_finish_int == 0){
                        slave_finish &= false;
                    }

                  
                    MPI_Recv(&(grid.current_buffer), 1, MPI_INT,MPI_ANY_SOURCE, 2, MPI_COMM_WORLD,&status);

                    MPI_Recv(&(g_data0[0]), m_length, MPI_DOUBLE,MPI_ANY_SOURCE, 3, MPI_COMM_WORLD,&status);
                    slave_rank = status.MPI_SOURCE;
                    for (int j=0;j<local_size*room_size; j++){
                        if(((slave_rank-1)*local_size*room_size + j) < room_size*room_size){
                            grid.data0[(slave_rank-1)*local_size*room_size + j] = g_data0[(slave_rank-1)*local_size*room_size + j];
                        }
                    }

                    MPI_Recv(&(g_data1[0]), m_length, MPI_DOUBLE,MPI_ANY_SOURCE, 4, MPI_COMM_WORLD,&status);
                    slave_rank = status.MPI_SOURCE;
                    for (int j=0;j<local_size*room_size; j++){
                        if(((slave_rank-1)*local_size*room_size + j) < room_size*room_size){
                            grid.data1[(slave_rank-1)*local_size*room_size + j] = g_data1[(slave_rank-1)*local_size*room_size + j];
                        }
                    }
                }

                // finished = hdist::calculate(current_state, grid);
                finished = slave_finish;
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
    else{
        // slave process
        while(true){
            
            // judge whether finished or not
            int s_flage = 1;
            MPI_Recv(&s_flage, 1, MPI_INT,0,0,MPI_COMM_WORLD, &status);
            if(s_flage == 0) break;

            // Receive Size Information
            int s_para_Size[2];
            MPI_Recv(&(s_para_Size[0]), 2, MPI_INT,0,rank*10,MPI_COMM_WORLD, &status);
            int s_room_size = s_para_Size[0];
            int s_local_size = s_para_Size[1];

            // Receive Information

            hdist::State s_current_state;
            
            int s_para_State_int[4];
            float s_para_State_float[5];
            int s_para_Grid[2];
            MPI_Recv(&(s_para_State_int[0]), 4, MPI_INT,0,rank*10+1,MPI_COMM_WORLD, &status);
            MPI_Recv(&(s_para_State_float[0]), 5, MPI_FLOAT,0,rank*10+2,MPI_COMM_WORLD, &status);
            // update local current_state
            s_current_state.room_size = s_para_State_int[0];
            s_current_state.source_x = s_para_State_int[1];
            s_current_state.source_y = s_para_State_int[2];
            if(s_para_State_int[3] == 0){
                s_current_state.algo = hdist::Algorithm::Jacobi;
            }else s_current_state.algo = hdist::Algorithm::Sor;

            s_current_state.block_size = s_para_State_float[0];
            s_current_state.source_temp = s_para_State_float[1];
            s_current_state.border_temp = s_para_State_float[2];
            s_current_state.tolerance = s_para_State_float[3];
            s_current_state.sor_constant = s_para_State_float[4];

            // create slave grid;
            auto s_grid = hdist::Grid{
                static_cast<size_t>(s_current_state.room_size),
                s_current_state.border_temp,
                s_current_state.source_temp,
                static_cast<size_t>(s_current_state.source_x),
                static_cast<size_t>(s_current_state.source_y)};

            int s_length = (s_grid.data0).size();
            double s_data0[s_length];
            double s_data1[s_length];
            MPI_Recv(&(s_data0[0]), s_length, MPI_DOUBLE,0,rank*10+3,MPI_COMM_WORLD, &status);
            MPI_Recv(&(s_data1[0]), s_length, MPI_DOUBLE,0,rank*10+4,MPI_COMM_WORLD, &status);
            MPI_Recv(&(s_para_Grid[0]), 2, MPI_INT,0,rank*10+5,MPI_COMM_WORLD, &status);


            // update local grid
            s_grid.current_buffer = s_para_Grid[0];
            s_grid.length = s_para_Grid[1];
            for (int m=0; m<s_length; m++){
                s_grid.data0[m] = s_data0[m];
                s_grid.data1[m] = s_data1[m];
            }
            
            bool s_finished = false;
            s_finished = hdist::calculate(s_current_state, s_grid, s_local_size, s_room_size, rank);   // rank starts from 1

            // send back s_finished and grid parameters
            // finished: flase = 0; true = 1
            int s_finished_int = 0;
            if (s_finished){
                s_finished_int = 1;
            }

            for (int n=0; n<s_length; n++){
                s_data0[n] = s_grid.data0[n];
                s_data1[n] = s_grid.data1[n];
            }

            MPI_Isend(&s_finished_int, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &request);
            MPI_Isend(&(s_grid.current_buffer), 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &request);
            MPI_Isend(&(s_data0[0]), s_length, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &request);
            MPI_Isend(&(s_data1[0]), s_length, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &request);
            // MPI_Isend(&(s_grid.length), 1, MPI_INT, 0, 5, MPI_COMM_WORLD, &request);
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
    return 0;
}
