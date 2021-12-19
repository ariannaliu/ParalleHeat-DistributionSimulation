#pragma once

#include <vector>
#include <cmath>

// namespace hdist {

    // enum class Algorithm : int {
    //     Jacobi = 0,
    //     Sor = 1
    // };

    // struct State {
    //     int room_size = 300;
    //     float block_size = 2;
    //     int source_x = room_size / 2;
    //     int source_y = room_size / 2;
    //     float source_temp = 100;
    //     float border_temp = 36;
    //     float tolerance = 0.02;
    //     float sor_constant = 4.0;
    //     Algorithm algo = hdist::Algorithm::Jacobi;

    //     bool operator==(const State &that) const = default;
    // };

    // struct Alt {
    // };

    // constexpr static inline Alt alt{};

    // struct Grid {
        // std::vector<double> data0, data1;
        // size_t current_buffer = 0;
        // size_t length;

        // explicit Grid(size_t size,
        //               double border_temp,
        //               double source_temp,
        //               size_t x,
        //               size_t y)
        //         : data0(size * size), data1(size * size), length(size) {
        //     for (size_t i = 0; i < length; ++i) {
        //         for (size_t j = 0; j < length; ++j) {
        //             if (i == 0 || j == 0 || i == length - 1 || j == length - 1) {
        //                 this->operator[]({i, j}) = border_temp;
        //             } else if (i == x && j == y) {
        //                 this->operator[]({i, j}) = source_temp;
        //             } else {
        //                 this->operator[]({i, j}) = 0;
        //             }
        //         }
        //     }
        // }

        // std::vector<double> &get_current_buffer() {
        //     if (current_buffer == 0) return data0;
        //     return data1;
        // }

        // double &operator[](std::pair<size_t, size_t> index) {
        //     return get_current_buffer()[index.first * length + index.second];
        // }

        // double &operator[](std::tuple<Alt, size_t, size_t> index) {
        //     return current_buffer == 1 ? data0[std::get<1>(index) * length + std::get<2>(index)] : data1[
        //             std::get<1>(index) * length + std::get<2>(index)];
        // }

        // void switch_buffer() {
        //     current_buffer = !current_buffer;
        // }
    // };

    // struct UpdateResult {
    //     bool stable;
    //     double temp;
    // };

    // UpdateResult update_single(size_t i, size_t j, Grid &grid, const State &state) {
    //     UpdateResult result{};
    //     if (i == 0 || j == 0 || i == state.room_size - 1 || j == state.room_size - 1) {
    //         result.temp = state.border_temp;
    //     } else if (i == state.source_x && j == state.source_y) {
    //         result.temp = state.source_temp;
    //     } else {
    //         auto sum = (grid[{i + 1, j}] + grid[{i - 1, j}] + grid[{i, j + 1}] + grid[{i, j - 1}]);
    //         switch (state.algo) {
    //             case Algorithm::Jacobi:
    //                 result.temp = 0.25 * sum;
    //                 break;
    //             case Algorithm::Sor:
    //                 result.temp = grid[{i, j}] + (1.0 / state.sor_constant) * (sum - 4.0 * grid[{i, j}]);
    //                 break;
    //         }
    //     }
    //     result.stable = std::fabs(grid[{i, j}] - result.temp) < state.tolerance;
    //     return result;
    // }

    // bool calculate(const State &state, Grid &grid) {
    //     bool stabilized = true;

    //     switch (state.algo) {
    //         case Algorithm::Jacobi:
    //             for (size_t i = 0; i < state.room_size; ++i) {
    //                 for (size_t j = 0; j < state.room_size; ++j) {
    //                     auto result = update_single(i, j, grid, state);
    //                     stabilized &= result.stable;
    //                     grid[{alt, i, j}] = result.temp;
    //                 }
    //             }
    //             grid.switch_buffer();
    //             break;
    //         case Algorithm::Sor:
    //             for (auto k : {0, 1}) {
    //                 for (size_t i = 0; i < state.room_size; i++) {
    //                     for (size_t j = 0; j < state.room_size; j++) {
    //                         if (k == ((i + j) & 1)) {
    //                             auto result = update_single(i, j, grid, state);
    //                             stabilized &= result.stable;
    //                             grid[{alt, i, j}] = result.temp;
    //                         } else {
    //                             grid[{alt, i, j}] = grid[{i, j}];
    //                         }
    //                     }
    //                 }
    //                 grid.switch_buffer();
    //             }
    //     }
    //     return stabilized;
    // };


// } // namespace hdist

// __device__ __managed__ bool local_stable;

// int get_index_in_array(size_t i, size_t j, int length){
//     return i * length + j;
// }

// int get_buffer(int current_buffer){
//     if (current_buffer == 0){
//         return 0;
//     }
//     return 1;
// }

// int get_buffer_alt(int current_buffer){
//     if (current_buffer == 0){
//         return 1;
//     }
//     return 0;
// }

// int switch_buffer(int current_buffer){
//     if (current_buffer == 0){
//         return 1;
//     }
//     return 0;
// }

// void update_single(int i, int j, 
//                    int room_size, float block_size,
//                 int source_x, int source_y,
//                 float source_temp, float border_temp,
//                 float tolerance, float sor_constant,
//                 int algo,
//                 double data0[], double data1[], int length, int current_buffer,
//                 bool & stable, double & temp){
//     if (i == 0 || j == 0 || i == room_size - 1 || j == room_size - 1) {
//         temp = border_temp;
//     } else if (i == source_x && j == source_y) {
//         temp = source_temp;
//     } else {
//         if(current_buffer == 0){
//             auto sum = data0[get_index_in_array(i + 1, j, length)] + data0[get_index_in_array(i - 1, j, length)]
//                         + data0[get_index_in_array(i, j+1, length)] + data0[get_index_in_array(i, j-1, length)];
//         }else{
//             auto sum = data1[get_index_in_array(i + 1, j, length)] + data1[get_index_in_array(i - 1, j, length)]
//                         + data1[get_index_in_array(i, j+1, length)] + data1[get_index_in_array(i, j-1, length)];
//         }
//         if(algo == 0){
//             temp =  0.25 * sum;
//         }else{
//             if(current_buffer == 0){
//                 temp = data0[get_index_in_array(i, j, length)] + (1.0 / sor_constant)*(sum - 4.0 * data0[get_index_in_array(i, j, length)]);
//             }else{
//                 temp = data1[get_index_in_array(i, j, length)] + (1.0 / sor_constant)*(sum - 4.0 * data1[get_index_in_array(i, j, length)]);
//             }
//         }
//     }
//     if(current_buffer == 0){
//         stable = std::fabs(data0[get_index_in_array(i, j, length)] - temp) < state.tolerance;   
//     }else{
//         stable = std::fabs(data1[get_index_in_array(i, j, length)] - temp) < state.tolerance;   
//     }
// }

// void calculate(int room_size, float block_size,
//                 int source_x, int source_y,
//                 float source_temp, float border_temp,
//                 float tolerance, float sor_constant,
//                 int algo,
//                 double data0[], double data1[], int length, int current_buffer, int local_size){
//     bool stabilized = true;
//     int thread_id = getThreadId();
    
//     if(algo == 0){
//         for (size_t i = (thread_id)*local_size; i < min(room_size, (thread_id+1)*local_size); ++i) {
//             for (size_t j = 0; j < room_size; ++j) {
//                 bool t_stable;
//                 double t_temp;
//                 update_single(i,j,room_size, block_size, 
//                             source_x, source_y, source_temp, border_temp,
//                             tolerance, sor_constant, algo, data0, data1, length, current_buffer,
//                             t_stable, t_temp);
//                 stabilized &= t_stable;
//                 if(current_buffer == 0){
//                     data1[get_index_in_array(i, j, length)] = t_temp;
//                 }else{
//                     data0[get_index_in_array(i, j, length)] = t_temp;
//                 }
//             }
//         }
//         current_buffer = switch_buffer(current_buffer);
//     }else{
//         for (auto k : {0, 1}) {
//             for (size_t i = (thread_id)*local_size; i < min(room_size, (thread_id+1)*local_size); ++i) {
//                 for (size_t j = 0; j < room_size; j++) {
//                     bool t_stable;
//                     double t_temp;
//                     if (k == ((i + j) & 1)) {
//                         auto result = update_single(i,j,room_size, block_size, 
//                                                     source_x, source_y, source_temp, border_temp,
//                                                     tolerance, sor_constant, algo, data0, data1, length, current_buffer,
//                                                     t_stable, t_temp);
//                         stabilized &= t_stable;

//                         if(current_buffer == 0){
//                             data1[get_index_in_array(i, j, length)] = t_temp;
//                         }else{
//                             data0[get_index_in_array(i, j, length)] = t_temp;
//                         }
//                     } else {
//                         // grid[{alt, i, j}] = grid[{i, j}];
//                         if(current_buffer == 0){
//                             data1[get_index_in_array(i, j, length)] = data0[get_index_in_array(i, j, length)];
//                         }else{
//                             data0[get_index_in_array(i, j, length)] = data1[get_index_in_array(i, j, length)];
//                         }
//                     }
//                 }
//             }
//             current_buffer = switch_buffer(current_buffer);
//         }
//     }
//     local_stable &= stabilized;
// }