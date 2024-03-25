#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <iterator>
#include <limits>
//#include <filesystem>

#include <gmp.h>
#include <gmpxx.h>
#include <boost/dynamic_bitset.hpp>

#include "num_of_prt_arr.h"
#include "bipermutations.h"
#include "utils.h"

int main(int argc, char* argv[]){     
    Bip bip;
    std::string bipPath = "";
    std::string outPath = "";
    int numberOfThreads = 1;
    int debug = 1;

    //std::cout << std::filesystem::current_path() << std::endl;

    if(argc < 2 || cmd_option_exists(argc, argv, "-h") || cmd_option_exists(argc, argv, "--help")){
        std::cout << "usage: ./prt_arr_counting <bip-path> [-o <out-path>] [-p <part>] [-n <num of parts>] [-d <debug level>] [--stats]" << std::endl;
        return 0;
    }

    std::vector<size_t> inBip;
    bipPath = argv[1];
    std::ifstream bipFile(bipPath);
    if(!bipFile.is_open()){
        std::cout << "ERROR: unable to open file at path: " << bipPath << std::endl;
        return -1;
    }
    std::istream_iterator<size_t> itr_start(bipFile), itr_end;
    std::copy(itr_start, itr_end, std::back_inserter(inBip));
    bipFile.close();
    size_t lineLimit = std::numeric_limits<Line>::max()-1;
    if(std::any_of(inBip.begin(), inBip.end(), [lineLimit](auto a){return a > lineLimit;})){
        std::cout << "ERROR: line labels over " << lineLimit << " are not supported. change definition of 'Line' in utils.h to raise it" << std::endl;
        return -1;
    }


    bip.resize(inBip.size());
    std::copy(inBip.begin(), inBip.end(), bip.begin());

    if(cmd_option_exists(argc, argv, "-o"))
        outPath = std::string(get_cmd_option(argc, argv, "-o"));

    if(cmd_option_exists(argc, argv, "-t"))
        numberOfThreads = atoi(get_cmd_option(argc, argv, "-t"));
    
    if(cmd_option_exists(argc, argv, "-d")) 
        debug = atoi(get_cmd_option(argc, argv, "-d"));

    if(cmd_option_exists(argc, argv, "-m")) 
        MBYTES_LIMIT = atoi(get_cmd_option(argc, argv, "-m"));

    if(debug >= 1){
        std::cout 
            << "inpath:\t\t" << bipPath << std::endl
            << "outpath:\t" << outPath << std::endl
            << "threads:\t" << numberOfThreads << std::endl;

        std::cout << "bipermutation:\t[ ";
        for (Line a : bip) std::cout << (int)a << " ";
        std::cout << "]" << std::endl;
    }

    SAVED_VALUES_DISTRIBUTION.resize(bip.size() / 2);
    ACCESSED_VALUES_DISTRIBUTION.resize(bip.size() / 2);

    size_t maxLevel =  bip.size() / 2+1;

    
    for(uint i = 5; i < knownValues.size(); ++i){
        std::vector<Line> v(2*i);
        std::iota(v.begin(), v.end(), 0);
        std::copy_n(v.begin(), i, v.begin()+i);
        size_t vHash = BipHasher{}(v) % (NUMBER_OF_MUTEXES/maxLevel);
        PREVIOUS_VALUES[(maxLevel-1)*NUMBER_OF_MUTEXES/maxLevel + vHash][v] = knownValues[i];
    }
    

    auto begin = std::chrono::system_clock::now();
    const std::clock_t cpuBegin = std::clock();
    mpz_class N = num_of_prt_arr_with_checks(bip, maxLevel, numberOfThreads, debug);
    auto end = std::chrono::system_clock::now();
    const std::clock_t cpuEnd = std::clock();

   std::chrono::duration<float> realTime = end - begin; //in float seconds
   float cpuTime = float(cpuEnd - cpuBegin)/CLOCKS_PER_SEC;

    if(cmd_option_exists(argc, argv, "--stats")){
        std::cout << "N\t\tsaved\t\taccessed: \n";
        for (size_t i = 5; i < SAVED_VALUES_DISTRIBUTION.size(); ++i)
            std::cout << i << "\t\t" << SAVED_VALUES_DISTRIBUTION[i] << "\t\t" << ACCESSED_VALUES_DISTRIBUTION[i] << "\n";
        std::cout << std::endl;
    }

    if(debug >= 2) 
        std::cout << "values saved:\t" << NUMBER_SAVED << "\n";

    if(debug >= 1){
        std::cout << "CPU time:\t" << cpuTime << "s" << std::endl;
        std::cout << "Real time:\t" << realTime.count() << "s" << std::endl;
        std::cout << "count:\t" << N << "\n";
    }

    if(!outPath.empty()){
        if(debug >= 1)
            std::cout << "save results at " << outPath << std::endl;
            
        std::ofstream outFile(outPath);
        outFile << "{ 'bip' : '";
        for(size_t i = 0; i < bip.size()-1; ++i)
            outFile << (int) bip[i] << " ";
        outFile << (int) bip[bip.size()-1] << "', ";
        outFile << " 'count' : "        << N
                << " , 'real_time' : "  << realTime.count()
                << " , 'threads' : "    << numberOfThreads
                << " , 'cores' : "      << std::thread::hardware_concurrency()
                << " , 'cpu_time' : "   << cpuTime
                << " , 'max_memory' : " << ((MAX_MBYTES_STORED)? MAX_MBYTES_STORED : (50 + 14*numberOfThreads + TOTAL_BYTES_STORED/1048576))
                << " }" << std::endl;
        outFile.close();
    }
    
    return 0;
}
