#pragma once

#include <stdint.h>
#include <algorithm>
//#include <string>
#include <iostream>
//#include <array>

using Line = uint8_t;
using MyPair = std::pair<size_t, size_t>;
using Bip = std::vector<Line>;


inline bool cmd_option_exists(int argc, char*  argv[], const std::string& option){
    return std::find(argv, argv + argc, option) != argv + argc;
}

inline char* get_cmd_option(int argc, char*  argv[], const std::string & option, int paramNr = 1){
    char ** itr = std::find(argv, argv + argc, option) + paramNr-1;
    if (itr != argv + argc && ++itr != argv + argc){
        return *itr;
    }
    return 0;
}


template<typename T>
constexpr T ipow(T base, T exp){
    T result = 1;
    for (;;){
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (exp == 0)
            break;
        base *= base;
    }
    return result;
}



template<typename T>
inline void print_vec(std::vector<T> const& vec){
    std::cout << "[ " << std::endl << "\t";
    for (int a : vec)
        std::cout << a << " ";
    std::cout << std::endl << "]" << std::endl << std::endl;
}







