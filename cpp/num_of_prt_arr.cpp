#include "num_of_prt_arr.h"

std::shared_mutex SHARED_MUTEXES[NUMBER_OF_MUTEXES];
std::unordered_map<Bip, mpz_class, BipHasher> PREVIOUS_VALUES[NUMBER_OF_MUTEXES];
size_t BYTES_STORED[NUMBER_OF_MUTEXES] = {};

std::vector<size_t> SAVED_VALUES_DISTRIBUTION;
std::vector<size_t> ACCESSED_VALUES_DISTRIBUTION;
std::atomic_size_t NUMBER_SAVED = 0;
size_t NUMBER_RECOMPUTED = 0;
std::atomic_size_t SORTING_INDEX = 0;

std::atomic_size_t TOTAL_BYTES_STORED = 0;
size_t MBYTES_LIMIT = 0;
size_t MAX_MBYTES_STORED = 0;



void free_memory(size_t const numberOfThreads, size_t const maxLevel){
	for (size_t i = 0; i < NUMBER_OF_MUTEXES; ++i){
		SHARED_MUTEXES[i].lock();
	}
	MAX_MBYTES_STORED = std::max(MAX_MBYTES_STORED, 50 + 14*numberOfThreads + TOTAL_BYTES_STORED/1048576); 
	for (size_t l = maxLevel - 1; l > 0 && TOTAL_BYTES_STORED/1048576 > MBYTES_LIMIT/2 - 23 - 7*numberOfThreads; --l){
		for (size_t i = 0; i < NUMBER_OF_MUTEXES/maxLevel; ++i){
			size_t const index = l*NUMBER_OF_MUTEXES/maxLevel + i;
			NUMBER_SAVED -= PREVIOUS_VALUES[index].size();
			PREVIOUS_VALUES[index].clear();
			TOTAL_BYTES_STORED -= BYTES_STORED[index];
			BYTES_STORED[index] = 0;
		}
	}
	for (size_t i = 0; i < NUMBER_OF_MUTEXES; ++i){
		SHARED_MUTEXES[i].unlock();
	}
}


mpz_class num_of_prt_arr(
	Bip bip,
	size_t const level,
	size_t const maxLevel,
	size_t const numberOfThreads, 
	int const debug
	) {
	size_t numberOfLines = bip.size() / 2;

	if (numberOfLines >= 9){ 
		reduce_bip(bip);
		numberOfLines = bip.size() / 2;
		assert(is_bip(bip));
	}

	if (numberOfLines <= 2) return 1;
	if (numberOfLines == 3)	return (bip[0] == bip[3] && bip[1] == bip[4]) ? 2 : 1;
	if (numberOfLines == 4) return num_of_prt_arr_4_lines(bip);

	std::vector<MyPair> lineIndeces = get_all_indeces(bip);

	auto canonicalBip = get_canonical(bip, lineIndeces);
	size_t bipHash = BipHasher{}(canonicalBip) % (NUMBER_OF_MUTEXES/maxLevel);
	bool already_computing = false;
	{	
		const std::shared_lock<std::shared_mutex> lock(SHARED_MUTEXES[level*NUMBER_OF_MUTEXES/maxLevel + bipHash]);
		if (PREVIOUS_VALUES[level*NUMBER_OF_MUTEXES/maxLevel + bipHash].find(canonicalBip) != PREVIOUS_VALUES[level*NUMBER_OF_MUTEXES/maxLevel + bipHash].end()) {
			already_computing = true;
			ACCESSED_VALUES_DISTRIBUTION[numberOfLines]++;
			mpz_class r = PREVIOUS_VALUES[level*NUMBER_OF_MUTEXES/maxLevel + bipHash][canonicalBip];
			if(r != 0)
				return r;
		}
	}
	if(already_computing){
		for(size_t i = 0; i < 100; ++i){
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			{
				const std::shared_lock<std::shared_mutex> lock(SHARED_MUTEXES[level*NUMBER_OF_MUTEXES/maxLevel + bipHash]);
				mpz_class r = PREVIOUS_VALUES[level*NUMBER_OF_MUTEXES/maxLevel + bipHash][canonicalBip];
				if(r != 0)
					return r;
			}
		}
	}

	{
		const std::unique_lock<std::shared_mutex> lock(SHARED_MUTEXES[level*NUMBER_OF_MUTEXES/maxLevel + bipHash]);
		PREVIOUS_VALUES[level*NUMBER_OF_MUTEXES/maxLevel + bipHash][canonicalBip] = 0;
	}
	
	auto const adjMatrix = get_adj_matrix(bip);

	size_t const splitLine = select_line(lineIndeces, adjMatrix);
	if(debug >= 2 && level == maxLevel-1) std::cout << "Split at line " << splitLine << std::endl << std::endl;

	PartialOrder partialOrder = get_partial_order(splitLine, lineIndeces[splitLine], bip, adjMatrix);
	std::vector<std::vector<Line>> const topSortings = partialOrder.all_linear_extensions();

	auto [leftHalf, rightHalf] = get_half_bips(lineIndeces[splitLine], bip, partialOrder.lineLabels);

	size_t const cutSize = topSortings[0].size();
	size_t const leftSize = leftHalf.size() + cutSize;
	size_t const rightSize = rightHalf.size() + cutSize;

	leftHalf.resize(leftSize);
	rightHalf.resize(rightSize);

	if (level == maxLevel-1 && numberOfThreads){
		std::vector<std::future<mpz_class>> results(numberOfThreads);
		for(size_t i = 0; i < numberOfThreads; ++i){
			results[i] = std::async(sum_recursion,
									leftHalf, rightHalf,  
									leftSize, rightSize, 
									cutSize,
									topSortings,
									level,
									maxLevel,
									numberOfThreads,
									i,
									debug
									);
		}
		mpz_class total = 0;
		for(size_t i = 0; i < numberOfThreads; ++i){
			total += results[i].get();
		}
		return total;
	}else{
		mpz_class total = sum_recursion(
			leftHalf, rightHalf, 
			leftSize, rightSize, 
			cutSize,
			topSortings,
			level,
			maxLevel,
			numberOfThreads
			);

		if(level != maxLevel-1){
			const std::unique_lock<std::shared_mutex> lock(SHARED_MUTEXES[level*NUMBER_OF_MUTEXES/maxLevel + bipHash]);
			if (PREVIOUS_VALUES[level*NUMBER_OF_MUTEXES/maxLevel + bipHash][canonicalBip] != 0)
				NUMBER_RECOMPUTED++;
			else{
				NUMBER_SAVED++;
				SAVED_VALUES_DISTRIBUTION[numberOfLines]++;
				size_t const index = level*NUMBER_OF_MUTEXES/maxLevel + bipHash;
				PREVIOUS_VALUES[index][canonicalBip] = total;
				size_t const storage = canonicalBip.size()*sizeof(Line) + mpz_sizeinbase(total.get_mpz_t(),2)/8 + 100;
				TOTAL_BYTES_STORED += storage;
				BYTES_STORED[index] += storage;
			}
		}

		if(MBYTES_LIMIT && 50 + 14*numberOfThreads + TOTAL_BYTES_STORED/1048576 > MBYTES_LIMIT){
			std::async(free_memory, numberOfThreads, maxLevel);
		}
		return total;
	}
}


mpz_class sum_recursion(
	std::vector<Line> leftHalf, std::vector<Line> rightHalf, 
	size_t const leftSize, 
	size_t const rightSize,
	size_t const cutSize,
	std::vector<std::vector<Line>>const& topSortings,
	size_t const level,
	size_t const maxLevel,
	size_t const numberOfThreads,
	size_t const thread, 
	int const debug
	){

	mpz_class total = 0;

	size_t const printNumber = topSortings.size() / 10 + 1;

	for (size_t i = 0; i < topSortings.size(); i++) {
		if (level == maxLevel-1 && ((i%numberOfThreads) != thread)) continue;
		if(debug >= 2){
			if((SORTING_INDEX++%printNumber) == 0 || debug >= 3) {			
				std::cout << "Sorting \t" << SORTING_INDEX << "/" << topSortings.size()
					<< "\nValues Saved: \t" << NUMBER_SAVED
					<< "\nRecomputed: \t" << NUMBER_RECOMPUTED
					<< "\nMbytes: \t" << 50 + 14*numberOfThreads + TOTAL_BYTES_STORED/1048576;
				if (numberOfThreads == 1)
					std::cout << "\nCurrent Total: \t" << total;
				std::cout << std::endl << std::endl;
			}
		}
		
		

		std::vector<Line>const& sorting = topSortings[i];

		std::copy(sorting.rbegin(), sorting.rend(), leftHalf.end() - cutSize);	
		std::copy(sorting.begin(), sorting.end(), rightHalf.end() - cutSize);

		mpz_class leftN;
		mpz_class rightN;

		if 		(leftSize <= 2) leftN = 1;
		else if (leftSize == 3) leftN = (leftHalf[0] == leftHalf[3] && leftHalf[1] == leftHalf[4]) ? 2 : 1;
		else if (leftSize == 4) leftN = num_of_prt_arr_4_lines(leftHalf);
		else 					leftN = num_of_prt_arr(leftHalf, level-1, maxLevel, numberOfThreads);

		if 		(rightSize <= 2) rightN = 1;
		else if (rightSize == 3) rightN = (rightHalf[0] == rightHalf[3] && rightHalf[1] == rightHalf[4]) ? 2 : 1;
		else if (rightSize == 4) rightN = num_of_prt_arr_4_lines(rightHalf);
		else 					 rightN = num_of_prt_arr(rightHalf, level-1, maxLevel, numberOfThreads);

		total += leftN * rightN;
	}

	return total;
}


mpz_class num_of_prt_arr_with_checks(Bip bip, size_t const maxLevel, size_t const numberOfThreads, int const debug) {
	for (Line a : bip) {
		int number = std::count_if(bip.begin(), bip.end(), [a](Line b){return a == b;});
		if (number != 2) {
			std::cout << "label " << a << " appears only once" << std::endl;
			throw std::invalid_argument("label appears only once");
		}
	}

	for (size_t i = 0; i < bip.size() / 2; i++) {
		bool valueAppeared = std::any_of(bip.begin(), bip.end(), [i](Line a){return a==i;});
		if (!valueAppeared) {
			std::cout << "label " << i << "does not appear" << std::endl;
			throw std::invalid_argument("a value does not appear");
		}
	}

	return num_of_prt_arr(bip, maxLevel-1, maxLevel, numberOfThreads, debug);
}


uint8_t num_of_prt_arr_4_lines(Bip const& bip) {
	std::bitset<4*4> A = 0;
	std::bitset<4*4> seenOnce = 0; //only 4 needed
	std::bitset<4*4> mask = 0b1111;

	for(size_t a : bip){
		if(seenOnce[a]){ 
			seenOnce[a] = 0;
			A ^= seenOnce << 4*a;
		}
		else{
			A &= (mask << 4*a).flip();
			A |= seenOnce << 4*a;
			seenOnce[a] = 1;
		} 
	}

	uint8_t d[4];

	d[0] = (A & (mask << 4*0)).count();
	d[1] = (A & (mask << 4*1)).count();
	d[2] = (A & (mask << 4*2)).count();
	d[3] = (A & (mask << 4*3)).count();
	
	uint8_t twiceNumOfCr = d[0] + d[1] + d[2] + d[3];

	if(twiceNumOfCr <= 2*2) return 1;
	if(twiceNumOfCr == 3*2){
		if(d[0] == 0 || d[1] == 0 || d[2] == 0 || d[3] == 0) return 2;
		else return 1;
	}
	if(twiceNumOfCr == 4*2){
		if(d[0] == 1 || d[1] == 1 || d[2] == 1 || d[3] == 1) return 2;
		else return 1;
	}
	if(twiceNumOfCr == 5*2) return 3;
	else return 8; //if(twiceNumOfCr == 6*2) return 8;
}


void reduce_bip(Bip & bip){
	size_t numberOfLines = bip.size() / 2;

	auto M = get_adj_matrix(bip);
	std::vector<bool> to_keep(numberOfLines, false);

	for(size_t i = 0; i < numberOfLines; ++i){
		for(size_t j = i+1; j < numberOfLines; ++j){
			if(!M[i][j]) continue;
			for(size_t k = j+1; k < numberOfLines; ++k){
				if(M[i][k] && M[j][k]){
					to_keep[i] = true;
					to_keep[j] = true;
					to_keep[k] = true;
				}
			}
		}
	}

	size_t oldNumberOfLines = numberOfLines;
	numberOfLines = 0;
	for(size_t i = 0; i < oldNumberOfLines; ++i){
		numberOfLines += to_keep[i];
	}

	size_t sentinel = 0;
	for(size_t i = 0; i < 2*numberOfLines;){
		if(to_keep[bip[sentinel]])
			bip[i++] = bip[sentinel++];
		else
			sentinel++;
	}

	bip.resize(2*numberOfLines);

	std::vector<Line> labels(oldNumberOfLines, -1);
	Line currentLabel = 0;

	for(auto& a : bip){
		if(labels[a] == (Line)-1){ 
			labels[a] = currentLabel;
			currentLabel++;
		}
		a = labels[a];
	}
}


std::vector<MyPair> get_all_indeces(Bip const& bip) {
	size_t const numberOfLines = bip.size() / 2;
	std::vector<MyPair> allIndeces(numberOfLines, {-1, -1});

	for (size_t i = 0; i < bip.size(); i++) {
		assert(bip[i] < allIndeces.size());
		if (allIndeces[bip[i]].first == -1) 
			allIndeces[bip[i]].first = i;
		
		else 
			allIndeces[bip[i]].second = i;
	}
	return allIndeces;
}


std::vector<boost::dynamic_bitset<>> get_adj_matrix(Bip const& bip) {
	std::vector<boost::dynamic_bitset<>> adjMatrix(bip.size() / 2);
	boost::dynamic_bitset<> seenOnce(bip.size() / 2, 0);

	for(size_t a : bip){
		if(seenOnce[a]){ 
			seenOnce[a] = 0;
            adjMatrix[a] ^= seenOnce;
		}
		else{
			adjMatrix[a] = seenOnce;	
			seenOnce[a] = 1;
		} 
	}
	return adjMatrix;
}


size_t select_line(std::vector<MyPair>const& allLineIndeces, std::vector<boost::dynamic_bitset<>> const& adjMatrix) {
	size_t const numberOfLines = allLineIndeces.size();
	size_t const size = numberOfLines*2;
	size_t minScore = size;
	size_t minScore2 = size;
	Line minLabel = 0;
	//bool foundALine = false;

	for (Line a = 0; a < numberOfLines; a++) {
		//size_t degree = std::reduce(adjMatrix[a].begin(), adjMatrix[a].end(), (size_t)0);
		size_t degree = adjMatrix[a].count();

		auto const&[x1, x2] = allLineIndeces[a];

		size_t da = x2 - x1;
		size_t score = std::max(da, size - da) + degree;
		size_t score2 = std::min(da, size - da) + degree;

		if (score < minScore || (score == minScore && score2 < minScore2)) {
			minScore = score;
			minScore2 = score2;
			minLabel = a;
		}
	}

	return minLabel;
}


PartialOrder get_partial_order(Line splitLine, MyPair splitIdx, Bip const& bip, std::vector<boost::dynamic_bitset<>>const& adjMatrix) {
	//size_t degree = std::reduce(adjMatrix[splitLine].begin(), adjMatrix[splitLine].end(), (size_t)0);
	size_t degree = adjMatrix[splitLine].count();

	PartialOrder partialOrder(degree);
	std::vector<Line> crossingLines;
	crossingLines.reserve(bip.size() / 2);

	for (size_t i = splitIdx.first + 1; i < splitIdx.second; ++i) {
		if (adjMatrix[splitLine][bip[i]]) {
			for (size_t j = 0; j < crossingLines.size(); ++j) {
				if (!adjMatrix[bip[i]][crossingLines[j]]) {
					partialOrder.add_edge(crossingLines.size(), j); //j needs to come above the new element, i.e. AFTER so we need new > j
				}
			}
			crossingLines.push_back(bip[i]);
		}
	}

	partialOrder.lineLabels = crossingLines;

	return partialOrder;
}


std::pair<std::vector<Line>, std::vector<Line>> get_half_bips(MyPair splitIdx, Bip const& bip, std::vector<Line>const& crossingLines) {
	//size_t size = bip.size();
	auto kItr1 = bip.begin() + splitIdx.first;
	auto kItr2 = bip.begin() + splitIdx.second;
	std::vector<Line> rightHalf(kItr1 + 1, kItr2);
	std::vector<Line> leftHalf(kItr2 + 1, bip.end());
	leftHalf.reserve(leftHalf.size() + std::distance(bip.begin(), kItr1));
	leftHalf.insert(leftHalf.end(), bip.begin(), kItr1);

	//renormalize: crossing lines should have labels 0,...,h in the order they appear in crossingLines 
	//other lines get other labels in arbitrary order

	std::unordered_map<size_t, size_t> crossLabels; //should be vector
	std::unordered_map<size_t, size_t> newLabels;
	size_t currentLabel = 0;
	for (size_t const& a : crossingLines) {
		crossLabels[a] = currentLabel;
		currentLabel++;
	}

	newLabels = crossLabels;
	currentLabel = crossingLines.size();
	for (Line& a : leftHalf) {
		if (newLabels.find(a) == newLabels.end()) {
			newLabels[a] = currentLabel;
			a = currentLabel;
			currentLabel++;
		}
		else {
			a = newLabels[a];
		}
	}

	newLabels = crossLabels;
	currentLabel = crossingLines.size();
	for (Line& a : rightHalf) {
		if (newLabels.find(a) == newLabels.end()) {
			newLabels[a] = currentLabel;
			a = currentLabel;
			currentLabel++;
		}
		else {
			a = newLabels[a];
		}
	}

	return { leftHalf, rightHalf };
}


Bip get_canonical(Bip const& bip, std::vector<MyPair>const& lineIndeces){
	size_t const size = bip.size();
	size_t const numberOfLines = size/2;
	

	//get promising indeces
	size_t minScore = size;
	/**/
	std::vector<MyPair> promisingIndeces(numberOfLines);
	size_t promisingIndecesSize = 0;
	promisingIndeces.resize(numberOfLines);
	for (auto const& [x, y] : lineIndeces) {
		size_t dx = y - x;
		size_t score = std::min(dx, size - dx);

		bool notWorse = (score <= minScore);
		bool better = (score < minScore);

		promisingIndecesSize = better ? 0 : promisingIndecesSize;
		minScore = better ? score : minScore;

		promisingIndeces[promisingIndecesSize] = (score == dx) ? std::make_pair(x, y) : std::make_pair(y, x); // dieser schreibzugriff ist egal! solang promisingIndecesSize nicht erhöht wird zählt das sowieso nicht
		promisingIndecesSize = notWorse ? promisingIndecesSize + 1 : promisingIndecesSize;
	}
	promisingIndeces.resize(promisingIndecesSize);

	/*
	std::vector<MyPair> promisingIndeces;
	promisingIndeces.reserve(numberOfLines);
	for (auto &[x, y] : lineIndeces) {
		size_t score = std::min(y - x, size - (y - x));
		if (score < minScore) {
			promisingIndeces.clear();
			promisingIndeces.push_back({x,y});
			minScore = score;
		}
		else if (score == minScore) {
			promisingIndeces.push_back({x,y});// order needs to change here depending on y-x vs size - (y-x)
		}
	}
*/
	//get minimal bip
	Bip bestBip = bip;

	//forwards
	std::vector<Line> newLabels(numberOfLines, -1);
	Bip currentBip;
	currentBip.reserve(size);
	//bool isSmaller;
	Line label;

	for (auto const&[x, _] : promisingIndeces) {
		if (x == 0) continue;
		bool isSmaller = false;
		label = 0;

		for (size_t i = 0; i < size; i++) {
			Line a = bip[(i + x) % size];
			if (newLabels[a] == (Line)-1) {
				newLabels[a] = label;
				label++;
			}

			if (isSmaller || newLabels[a] == bestBip[i]) {
				currentBip.push_back(newLabels[a]);
			}
			else if (newLabels[a] < bestBip[i]) {
				currentBip.push_back(newLabels[a]);
				isSmaller = true;
			}
			else {
				break;
			}
		}

		if (isSmaller)
			bestBip = currentBip;

		std::fill(newLabels.begin(), newLabels.end(), -1);
		currentBip.resize(0);
	}

	//backwards
	for (auto const&[_, y] : promisingIndeces) {
		bool isSmaller = false;
		label = 0;

		for (size_t i = 0; i < size; i++) {
			Line a = bip.rbegin()[(i + size - y - 1) % size];
			if (newLabels[a] == (Line)-1) {
				newLabels[a] = label;
				label++;
			}

			if (isSmaller || newLabels[a] == bestBip[i]) {
				currentBip.push_back(newLabels[a]);
			}
			else if (newLabels[a] < bestBip[i]) {
				currentBip.push_back(newLabels[a]);
				isSmaller = true;
			}
			else {
				break;
			}
		}

		if (isSmaller)
			bestBip = currentBip;

		std::fill(newLabels.begin(), newLabels.end(), -1);
		currentBip.resize(0);
	}

	return bestBip;
}


/*
bool intersect(MyPair const& aIdx, MyPair const& bIdx) {
	auto& [a, b] = aIdx;
	auto& [u, v] = bIdx;

	return (a < u) != (a < v) != (b < u) != (b < v);
}*/

/*
bool intersect(size_t a, size_t b, Bip const& bip) {
	if (a == b) {
		return false;
	}
	size_t state = -1;

	for (size_t c : bip) {
		if (c == a) {
			if (state == a)
				return false;
			else
				state = a;
		}
		else if (c == b) {
			if (state == b)
				return false;
			else
				state = b;
		}
	}
	return true;
}*/


/*
std::tuple<int, int> get_indeces(int a, Bip const& bip) {
	int idx1 = -1;
	int idx2 = -1;

	for (size_t i = 0; i < bip.size(); i++) {
		if (bip[i] == a) {
			if (idx1 == -1) {
				idx1 = i;
			}
			else {
				idx2 = i;
				break;
			}
		}
	}

	return std::make_tuple(idx1, idx2);
}
*/


bool is_bip(Bip const& bip) {
	for (size_t a : bip) {
		int number = std::count_if(bip.begin(), bip.end(), [a](Line b){return a == b;});
		if (number != 2) {
			return false;
		}
	}

	return true;
}


/*
const char * num_prt_arr_char_p(char * b, size_t s, size_t part, size_t parts, int debug = 0){
	Bip bip(b, b + s);
	if(debug >= 1){
		std::cout << "compute number of prt arr for \n";
		std::cout << "[ ";
        for (size_t a : bip) std::cout << a << " ";
        std::cout << "]" << std::endl;
	}


	SAVED_VALUES_DISTRIBUTION.resize(bip.size() / 2);
    ACCESSED_VALUES_DISTRIBUTION.resize(bip.size() / 2);
	mpz_class N = num_of_prt_arr_with_checks(bip, part, parts, true, debug);

	if(debug >= 1){
		std::cout << "result: " << N << std::endl;
	}

	std::string N_string;
	N_string = N.get_str();
    return N_string.c_str();
}
*/
