#pragma once

#include <stdint.h>
#include <iostream>
#include <vector>
//#include <algorithm>
//#include <stdexcept>
#include <unordered_map>
//#include <cmath>
#include <numeric>
#include <cassert>
#include <bitset>
#include <shared_mutex>
#include <thread>
#include <future>
#include <atomic>
#include <chrono>

#include <gmp.h>
#include <gmpxx.h>
#include <boost/dynamic_bitset.hpp>

#include "prt_ord.h"
#include "utils.h"

#define NUMBER_OF_MUTEXES 100003

struct BipHasher { //from www.cse.yorku.ca/~oz/hash.html
	size_t operator()(Bip const& bip) const {
		size_t hash = 5381;
		for(auto c : bip) hash = ((hash << 5) + hash) + c;
		return hash;
	}
};

extern std::shared_mutex SHARED_MUTEXES[NUMBER_OF_MUTEXES];
extern std::unordered_map<Bip, mpz_class, BipHasher> PREVIOUS_VALUES[NUMBER_OF_MUTEXES];
extern std::vector<size_t> SAVED_VALUES_DISTRIBUTION;
extern std::vector<size_t> ACCESSED_VALUES_DISTRIBUTION;
extern std::atomic_size_t NUMBER_SAVED;
extern std::atomic_size_t TOTAL_BYTES_STORED;
extern size_t MBYTES_LIMIT;
extern size_t MAX_MBYTES_STORED;


mpz_class 
num_of_prt_arr(
	Bip bip,
	size_t const level,
	size_t const maxLevel,
	size_t const numberOfThreads,
	int const debug = 0
	);

mpz_class 
num_of_prt_arr_with_checks(
	Bip bip,
	size_t const maxLevel,
	size_t const numberOfThreads = 1,
	int const debug = 0
	);

uint8_t 
num_of_prt_arr_4_lines(Bip const& bip);


void 
reduce_bip(Bip & bip);

Bip 
get_canonical(Bip const& bip, std::vector<MyPair>const& lineIndeces);

std::vector<MyPair> 
get_all_indeces(Bip const& bip);

std::vector<boost::dynamic_bitset<>> 
get_adj_matrix(Bip const& bip);

size_t 
select_line(std::vector<MyPair>const& allLineIndeces, std::vector<boost::dynamic_bitset<>>const& adjMatrix);

PartialOrder 
get_partial_order(Line splitLine, MyPair splitIdx, Bip const& bip, std::vector<boost::dynamic_bitset<>>const& adjMatrix);

std::pair<std::vector<Line>, std::vector<Line>> 
get_half_bips(MyPair splitIdx, Bip const& bip, std::vector<Line>const& crossingLines);

mpz_class 
sum_recursion(
	std::vector<Line> leftHalf, std::vector<Line> rightHalf, 
	//std::vector<MyPair>& fixedlLeftLineIndeces, std::vector<MyPair>& fixedRightLineIndeces,
	size_t const leftSize, 
	size_t const rightSize,
	size_t const cutsize,
	std::vector<std::vector<Line>>const& topSortings,
	size_t const level,
	size_t const maxLevel,
	size_t const numberOfThreads,
	size_t const thread = 0, 
	int const debug = 0
	);


//bool intersect(size_t a, size_t b, Bip const& bip);
//bool intersect(MyPair const& aIdx, MyPair const& bIdx);
//std::tuple<int, int> get_indeces(int a, Bip const& bip);
bool is_bip(Bip const& bip);





//extern "C" const char * num_prt_arr_char_p(char * b, size_t s, size_t part, size_t parts, int debug);


/*

using fn_type = mpz_class(Bip const&, std::vector<MyPair>const&, size_t, size_t, bool, bool);
static constexpr size_t MAX_NUMBER_OF_LINES = 30;

template<std::size_t... I>
constexpr auto make_helper(std::index_sequence<I...>) {
    return std::array<fn_type *, sizeof...(I)> { calc_num_of_prt_arr<I>... };
}
template <std::size_t MAX>
constexpr auto make_lut() {
    return make_helper(std::make_index_sequence<MAX>{});
}

static constexpr auto CALC_NUM_OF_PRT_ARR_LUT = make_lut<MAX_NUMBER_OF_LINES+1>();
*/

/*
struct VectorHasher {
#define djb2
#ifdef djb2
//from www.cse.yorku.ca/~oz/hash.html
	size_t operator()(Bip const& bip) const {
		size_t hash = 5381;
		for(auto c : bip)
			hash = ((hash << 5) + hash) + c;
		
		return hash;
	}
#endif
#ifdef alt
//form https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector/72073933#72073933
	std::size_t operator()(std::vector<Line> const& vec) const {
		std::size_t seed = vec.size();
		for (auto x : vec) {
			x = ((x << 16) ^ x) * 0x45d9f3b;
			x = ((x >> 16) ^ x) * 0x45d9f3b;
			x = (x >> 16) ^ x;
			seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
#endif
#ifdef manf
	static size_t power(size_t base, size_t exponent) {
		size_t tmp = 1;
		while (exponent) {
			if (exponent % 2) tmp *= base;
			base *= base;
			exponent /= 2;
		}
		return tmp;
	}

	std::size_t operator()(std::vector<size_t> const& vec) const {
		std::size_t seed = 0;
		size_t n = vec.size() / 2;
		size_t visited = -1;
		for (int i = 0; i < 2 * n; i++) {
			size_t x = vec[i];
			if (x > visited) {
				visited = x;
			}
			else { // already seen this value
				seed += power(x, i);
			}
		}
		return seed;
	}
#endif
#undef djb2
#undef alt
#undef manf
};
*/
