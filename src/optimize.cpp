#include "headers.h"

/* ===== Function Description:
	Constructor of the 'Optimizer' class.
*/
Optimizer::Optimizer(int precision, float cost_single, bool is_same) : _r(precision), _is_same(is_same), _cost_single(cost_single) {
	_heights		= vector<int>(_r, 0);
	_n_carry		= vector<int>(_r, 0);
	_n_counter		= vector<int>(_r, 0);
	_n_split_from	= vector<int>(_r, 0);
	_n_split_to		= vector<int>(_r, 0);
	_counter_sizes	= vector<vector<int>>(_r);
	_bit_table		= vector<vector<Bit>>(_r);
}

/* ===== Function Description:
	Find the peaks and update the '_max_height' variable.
*/
void Optimizer::updatePeaks(vector<int>& peaks) {		
	peaks.clear();
	_max_height = *max_element(_heights.begin(), _heights.end());
	if (_max_height == 0) return;

	for (int i = _r - 1; i >= 0; i--) {
		if (_heights[i] == _max_height) {
			peaks.push_back(i);
		}
	}
}

/* ===== Function Description:	
	Main synthesis process.
*/
pair<float, int> Optimizer::optimize(bool to_print_info) {
	float total_cost = 0;
	int n_adder = 0;
	for (Gate* gate : _gate_list) {
		if (gate->getTypeStr() == "cp") {
			total_cost += COST_TOFFOLI;
		}
	}
	if (_is_same) total_cost += _cost_single * _r;

	vector<int> peaks, remaining;
	for (int ith_iter = 0; ; ++ith_iter) {
		if (to_print_info) 
			printInfo("iteration " + to_string(ith_iter) + " (current cost = " + to_string(total_cost) + "):");
		
		// method 0: split and fill
		updatePeaks(peaks);
		if (_max_height == 0) break;

		// find the LSB with the second-high height
		int second_height = -1;
		for (int i = 0; i < _r; ++i) {
			if (_heights[i] != _max_height && _heights[i] > second_height) {
				second_height = _heights[i];
			}
		}
		int secnod_height_index = peaks[0];
		for (int i = 0; i < _r; ++i) {
			if (_heights[i] >= second_height) {
				secnod_height_index = i;
			}
		}

		remaining.clear();
		for (int i = 0; i < peaks.size(); i++) {
			bool is_successful = split(peaks[i], secnod_height_index * 2);
			if (!is_successful) {
				remaining.push_back(peaks[i]);
			}
		}
		if (remaining.empty()) {
			//if (to_print_info) { cout << "[system pause] >> "; string temp; cin >> temp; cout << endl; }
			continue;
		}

		// method 1 + 2 : counter + new (partial) adder
		int dealing_remaining_index = 0;
		vector<int> new_heights_counter = _heights;
		vector<int> new_n_carry_counter = _n_carry;
		vector<int> new_n_counter_counter = _n_counter;
		vector<vector<int>> new_counter_sizes_counter = _counter_sizes;
		int cost_counter = doCounter(new_heights_counter, new_n_carry_counter, new_n_counter_counter, new_counter_sizes_counter, remaining, dealing_remaining_index);
		
		// method 3 : single-gate
		unordered_set<int> new_excluded_single;
		int cost_single = doSingle(new_excluded_single, remaining);

		if (cost_counter <= cost_single) {
			total_cost += cost_counter;
			_heights = new_heights_counter;
			_n_carry = new_n_carry_counter;
			_n_counter = new_n_counter_counter;
			_counter_sizes = new_counter_sizes_counter;

			if (dealing_remaining_index < remaining.size()) {  // new adder is used
				total_cost -= countAdderCost(remaining[dealing_remaining_index]);	// avoid counting twice
				//if (to_print_info) { cout << "[system pause] >> "; string temp; cin >> temp; cout << endl; }
				break;
			}
		}
		else {
			total_cost += cost_single;
			for (int i = 0; i < _r; i++) {
				for (Bit& bit : _bit_table[i]) {
					if (new_excluded_single.count(bit.getGateId()) > 0 && _heights[i] > 0) {
						bit.setInactivate();
						_excluded[bit.getGateId()] += pow(2, (-1 - i)) * 2 * M_PI;
						_heights[i]--;
					}
				}
			}
			removeExcluded();
		}
		//if (to_print_info) { cout << "[system pause] >> "; string temp; cin >> temp; cout << endl; }
	}

	// calculate the remaining cost of adders
	for (int i = _r - 1; i >= 0; --i) {
		if (_heights[i] > n_adder) {
			total_cost += countAdderCost(i) * (_heights[i] - n_adder);
			n_adder = _heights[i];
		}
	}
	return make_pair(total_cost, n_adder);
}

// ==================================================

/* ===== Function Description:
	Find the gate index to split a qubit. 
	Preferentially choose gates that can compensate for lower bits.
	If no gate is found, return -1; if no preference, return _n;
*/
int Optimizer::findSplittedGate(int index, unordered_set<int>& pos_gates, unordered_set<int>& neg_gates, int index_bound) {		
	vector<int> n_needed_bits(_n + 1, 1);	// [_n] for general

	for (int end_index = index + 1; end_index < min(_r, index_bound); ++end_index) {
		if (_heights[end_index] == _max_height) return -1;

		for (int i = 0; i < _n + 1; ++i) {
			n_needed_bits[i] *= 2;
		}
		if (*min_element(n_needed_bits.begin(), n_needed_bits.end()) > _max_height) return -1;

		for (int i = 0; i < _bit_table[end_index].size(); ++i) {
			int gate_id = _bit_table[end_index][i].getGateId();
			bool can_discharge = false;
			if (_bit_table[end_index][i].isPos() && neg_gates.count(gate_id) != 0) {
				can_discharge = true;
			}
			else if (_bit_table[end_index][i].isNeg() && pos_gates.count(gate_id) != 0) {
				can_discharge = true;
			}

			if (can_discharge) {
				n_needed_bits[gate_id] -= 2;			// cancel out one and put the other in the same position
				if (n_needed_bits[gate_id] == 0) {
					return gate_id;
				}
			}
		}
		// discharged with a bit with oppisite sign (reverse Booth's encoding)

		for (int i = 0; i < _n + 1; ++i) {
			n_needed_bits[i] -= _max_height - 1 - _heights[end_index];
		}

		int mini = INT_MAX;
		int mini_index = -1;
		for (int i = 0; i < _n + 1; ++i) {
			if (pos_gates.count(i) == 0 && neg_gates.count(i) == 0 && i != _n)
				continue;

			if (n_needed_bits[i] <= mini) {
				mini = n_needed_bits[i];
				mini_index = i;
			}
		}
		if (mini <= 0) {
			return mini_index;
		}
	}
	return -1;
}

/* ===== Function Description:
	Split the bit with 'gate_id' at 'index' column into lower bits.
*/
void Optimizer::splitGate(int index, int gate_id) {
	assert(gate_id != -1);

	int n_needed_bit = 1;
	if (gate_id == _n) {	// any bit can be used; keep the flexibility
		_n_split_from[index]++;
		_heights[index]--;

		for (int i = index + 1; i < _r; ++i) {
			n_needed_bit *= 2;
			int capacity = _max_height - 1 - _heights[i];
			if (n_needed_bit <= capacity) {
				_n_split_to[i] += n_needed_bit;
				_heights[i] += n_needed_bit;
				return;
			}
			else {
				n_needed_bit -= capacity;
				_n_split_to[i] += capacity;
				_heights[i] += capacity;
			}
		}
	}
	else {
		// decide the bit type
		BITTYPE bit_type;
		BITTYPE bit_type_inv;
		bool flag = false;
		for (int sub_index = 0; sub_index < _bit_table[index].size(); ++sub_index) {
			if (_bit_table[index][sub_index].getGateId() == gate_id) {
				bit_type = _bit_table[index][sub_index].getType();
				bit_type_inv = _bit_table[index][sub_index].getInvType();
				_bit_table[index].erase(_bit_table[index].begin() + sub_index);
				_heights[index]--;
				flag = true;
				break;
			}
		}
		if (!flag) assert(false);

		for (int i = index + 1; i < _r; ++i) {
			n_needed_bit *= 2;

			for (int j = 0; j < _bit_table[i].size(); j++) {
				if (_bit_table[i][j].getGateId() == gate_id && _bit_table[i][j].getType() == bit_type_inv) {
					_bit_table[i][j].invType();
					n_needed_bit -= 2;
					if (n_needed_bit == 0) {
						return;
					}
				}
			}

			while (_heights[i] < _max_height - 1) {
				_bit_table[i].push_back(Bit(bit_type, _gate_list[gate_id]));
				_heights[i]++;
				n_needed_bit--;

				if (n_needed_bit == 0) {
					return;
				}
			}
		}
	}
	assert(false);	// should end in the for loop
}

/* ===== Function Description:
	Try to use the split method to reduce the heiginvTypeht at 'index' column.
*/
bool Optimizer::split(int index, int index_bound) {
	if (_heights[index] - _n_carry[index] - _counter_sizes[index].size() <= 0) return false;
	// notice that carry bits and counter bits cannot be splitted, but splitted-to bits can be further splitted

	unordered_set<int> pos_gates, neg_gates;
	for (Bit& bit : _bit_table[index]) {
		if (bit.isPos()) pos_gates.insert(bit.getGateId());
		else if (bit.isNeg()) neg_gates.insert(bit.getGateId());
	}

	int splitted_gate = findSplittedGate(index, pos_gates, neg_gates, index_bound);
	if (splitted_gate == -1) {
		return false;
	}

	splitGate(index, splitted_gate);
	return true;
}

// ==================================================

/* ===== Function Description:
	Find the minimum one and maintain 'counter_sizes' in decreasing order.
	E.g., 4 4 4 3 3 2 2 [1] 1 1 1
*/
int Optimizer::findTargetCounterMin(const vector<int>& counter_sizes) {
	for (int i = counter_sizes.size() - 1; i > 0; --i) {\
		if (counter_sizes[i] < counter_sizes[i - 1]) {
			return i;
		}
	}  
	return 0;
}

/* ===== Function Description:
	Try to merge two counters.
*/
int Optimizer::mergeCounter(vector<int>& new_heights, vector<int>& new_n_carry, vector<int>& new_n_counter, vector<vector<int>>& new_counter_sizes, const vector<int>& peaks, int& dealing_peak_index) {
	int index = peaks[dealing_peak_index];

	int adder_saved_cost = 0;
	int counter_saved_cost = 0;
	int counter_extra_cost = 0;
	int new_peak = peaks[dealing_peak_index];

	// remove the last counter
	int original_counter_size = new_counter_sizes[index].back();
	new_counter_sizes[index].pop_back();
	counter_saved_cost = countCounterCost(original_counter_size, index);
	int original_counter_bitlength = log2(original_counter_size) + 1;

	for (; dealing_peak_index < peaks.size(); ++dealing_peak_index) {
		if (peaks[dealing_peak_index] >= index - original_counter_bitlength + 1) {
			if (dealing_peak_index + 1 < peaks.size())	new_peak = peaks[dealing_peak_index + 1];
			else										new_peak = -1;
		} // this range is also saved
	}
	if (new_peak == -1) adder_saved_cost = countAdderCost(index);
	else			          adder_saved_cost = countAdderCost(index) - countAdderCost(new_peak);
	for (int i = 0; i < original_counter_bitlength && index - i >= 0; ++i) {
		new_heights[index - i]--;
		if (i > 0) new_n_carry[index - i]--;
	}

	// merge into other counters
	for (int i = 0; i < original_counter_size; ++i) {
		int target_counter = findTargetCounterMin(new_counter_sizes[index]);
		
		new_counter_sizes[index][target_counter]++;
		int new_bitlength = log2(new_counter_sizes[index][target_counter]) + 1;
		bool need_new_carry = (new_counter_sizes[index][target_counter] == pow(2, new_bitlength - 1));
		if (need_new_carry && index - new_bitlength + 1 >= 0) {
			new_n_carry[index - new_bitlength + 1]++;
			new_heights[index - new_bitlength + 1]++;
			if (new_heights[index - new_bitlength + 1] >= _max_height) {
				return INT_MAX;
			}  // cannot apply counter method	// this restriction may be losen
		}
		counter_extra_cost += countCounterCost(new_counter_sizes[index][target_counter], index) - countCounterCost(new_counter_sizes[index][target_counter] - 1, index);
	}

	if (adder_saved_cost + counter_saved_cost > counter_extra_cost) {
		return (counter_extra_cost - counter_saved_cost);
	}
	else {
		return INT_MAX;
	}  // do not apply counter
}

/* ===== Function Description:
	Try to use the counter method to reduce the height at each peak column.
	The remaining peak columns are dealt by adder method.
	Return the cost.
*/
int Optimizer::doCounter(vector<int>& new_heights, vector<int>& new_n_carry, vector<int>& new_n_counter, vector<vector<int>>& new_counter_sizes, const vector<int>& peaks, int& dealing_peak_index) {
	int cost_adder = 0;
	for (dealing_peak_index = 0; dealing_peak_index < peaks.size(); dealing_peak_index++) {
		int index = peaks[dealing_peak_index];
		if (new_heights[index] - new_n_carry[index] <= 0) {
			return INT_MAX;
		} // cannot apply counter

		if (new_heights[index] - new_n_carry[index] - new_counter_sizes[index].size() <= 0) {	// can only merge old counters
			if (new_counter_sizes[index].size() < 2) {	// nothing to merge
				break;
			}

			vector<int> temp_heights = new_heights;
			vector<int> temp_n_carry = new_n_carry;
			vector<int> temp_n_counter = new_n_counter;
			vector<vector<int>> temp_counter_sizes = new_counter_sizes;
			int temp_dealing_peak_index = dealing_peak_index;

			int cost = mergeCounter(temp_heights, temp_n_carry, temp_n_counter, temp_counter_sizes, peaks, temp_dealing_peak_index);
			if (cost == INT_MAX) break;	// failed

			cost_adder += cost;
			new_heights = temp_heights;
			new_n_carry = temp_n_carry;
			new_n_counter = temp_n_counter;
			new_counter_sizes = temp_counter_sizes;
			dealing_peak_index = temp_dealing_peak_index;
		}
		else {	// merge a new bit into existing counter 

			int adder_saved_cost = 0;
			int counter_extra_cost = 0;

			if (dealing_peak_index + 1 < peaks.size()) {
				adder_saved_cost = countAdderCost(peaks[dealing_peak_index]) - countAdderCost(peaks[dealing_peak_index + 1]);
			}
			else {
				adder_saved_cost = countAdderCost(peaks[dealing_peak_index]);
			}

			bool create_new_counter = true;
			if (new_heights[index] - new_n_carry[index] - (new_n_counter[index] - new_counter_sizes[index].size()) < 2) {
				create_new_counter = false;
			}
			if (index - 1 >= 0 && new_heights[index - 1] >= _max_height) {
				create_new_counter = false;
			}

			if (create_new_counter) { // merge two bits to form a counter
				counter_extra_cost = countCounterCost(2, index);
				if (counter_extra_cost >= adder_saved_cost) {
					break;
				}  // do not use counter method

				// update
				new_heights[index]--;
				new_n_counter[index] += 2;
				new_counter_sizes[index].emplace_back(2);
				if (index - 1 >= 0) {
					new_heights[index - 1]++;
					new_n_carry[index - 1]++;
				}
				cost_adder += counter_extra_cost;
			}
			else { // merge a new bit into an existing counter
				if (new_counter_sizes[index].empty()) {
					break;
				}  // no counter to merge

				int target_counter = findTargetCounterMin(new_counter_sizes[index]);

				int new_bitlength = log2(new_counter_sizes[index][target_counter] + 1) + 1;
				bool need_new_carry = (new_counter_sizes[index][target_counter] + 1 == pow(2, new_bitlength - 1));
				if (need_new_carry && index - new_bitlength + 1 >= 0) {
					new_heights[index - new_bitlength + 1]++;
					if (new_heights[index - new_bitlength + 1] >= _max_height) {
						break;
					}  // cannot apply counter method	// this restriction may be losen
				}
				counter_extra_cost = countCounterCost(new_counter_sizes[index][target_counter] + 1, index) - countCounterCost(new_counter_sizes[index][target_counter], index);
				if (counter_extra_cost >= adder_saved_cost) {
					break;
				}  // do not use counter method
				
				// update
				new_n_counter[index]++;
				new_counter_sizes[index][target_counter]++;
				new_heights[index]--;
				if (need_new_carry && index - new_bitlength >= 0) {
					new_heights[index - new_bitlength]++;
					new_n_carry[index - new_bitlength]++;
				}
				cost_adder += counter_extra_cost;
			}
		}
	}
 
	if (dealing_peak_index < peaks.size()) {
		cost_adder += countAdderCost(peaks[dealing_peak_index]);
	}	// new adder is still needed ('break' is triggered)

	return cost_adder;
}

// ==================================================

/* ===== Function Description: 
	Remove excluded bits.
*/
void Optimizer::removeExcluded() {
	for (int i = 0; i < _r; ++i) {
		_bit_table[i].erase(remove_if(_bit_table[i].begin(), _bit_table[i].end(), [this](Bit& b) {return (!b.isActivate()); }), _bit_table[i].end());
	}
}

/* ===== Function Description:
	Try the single-gate method to reduce the height at each peak column.
	Return the cost.
*/
int Optimizer::doSingle(unordered_set<int>& new_excluded, vector<int> peaks_remaining) {
	for (int index : peaks_remaining) {
		if (_heights[index] - _n_carry[index] <= 0) {
			return INT_MAX;
		}  // cannot use single-gate method
	}
	
	// find removed gates
	while (!peaks_remaining.empty()) {
		vector<unordered_set<int>> n_involved_peaks(_n);
		for (int index : peaks_remaining) {
			vector<bool> is_counted(_n, false);
			for (Bit& bit : _bit_table[index]) {
				if (!is_counted[bit.getGateId()]) {
					n_involved_peaks[bit.getGateId()].insert(index);
					is_counted[bit.getGateId()] = true;
				}
			}
		}

		int max_involved_gate = -1;
		int n_involved = 0;
		for (int i = 0; i < _n; i++) {
			if (new_excluded.count(i) == 0 && n_involved_peaks[i].size() > n_involved) {
				max_involved_gate = i;
				n_involved = n_involved_peaks[i].size();
			}
		}
		if (max_involved_gate == -1) {
			return INT_MAX;
		}

		new_excluded.insert(max_involved_gate);
		auto it = remove_if(peaks_remaining.begin(), peaks_remaining.end(), 
			[this, n_involved_peaks, max_involved_gate](int& index) {return (n_involved_peaks[max_involved_gate].count(index)); }
		);
		peaks_remaining.erase(it, peaks_remaining.end());
	}

	// calculate cost // some counters may be saved, but we ignore them for simplicity
	int extra_cost = _cost_single * new_excluded.size();
	return extra_cost;
}

// ==================================================

/* ===== Function Description: 
	Turn all flexibilities into concrete implementations.
*/
void Optimizer::concrete() {
	// bit-spliting
	for (int i = 0; i < _r; ++i) {
		while (_n_split_from[i] > 0) {
			splitGateAny(i);
		}
	}
	for (int i = 0; i < _r; i++) {
		assert(_n_split_to[i] == 0);
		assert(_n_split_from[i] == 0);
	}

	// carries
	for (int i = _r - 1; i >= 0; --i) {
		for (int counter_size : _counter_sizes[i]) {
			vector<Bit> carrys;
			for (int k = 0; k < counter_size; ++k) {
				carrys.push_back(_bit_table[i][k]);
			}
			_bit_table[i].erase(_bit_table[i].begin(), _bit_table[i].begin() + counter_size);

			_bit_table[i].emplace_back(Bit(carrys, 0));
			for (int k = 1; k < int(log2(counter_size) + 1) && i - k >= 0; ++k) {
				_bit_table[i - k].emplace_back(Bit(carrys, k));
				_n_carry[i - k] -= 1;
			}

			_n_counter[i] -= counter_size;
		}
	}

	for (int i = 0; i < _r; i++) {
		assert(_n_carry[i] == 0);
		assert(_n_counter[i] == 0);
	}

	for (int i = 0; i < _r; i++) {
		assert(_bit_table[i].size() == _heights[i]);
	}
}

/* ===== Function Description:
	Split any bit at 'index' column into lower bits.
	The last bit is selected.
	(used at the final concretion stage)
*/
void Optimizer::splitGateAny(int index) {
	int gate_id = _bit_table[index].back().getGateId();
	BITTYPE bit_type = _bit_table[index].back().getType();
	_bit_table[index].pop_back();
	_n_split_from[index]--;

	int n_needed_bits = 1;
	for (int i = index + 1; i < _r; ++i) {
		n_needed_bits *= 2;
		while (_n_split_to[i] > 0) {
			_bit_table[i].push_back(Bit(bit_type, _gate_list[gate_id]));
			_n_split_to[i]--;
			n_needed_bits--;

			if (n_needed_bits == 0) {
				return;
			}
		}
	}
	assert(false);	// should end in the FOR loop
}
