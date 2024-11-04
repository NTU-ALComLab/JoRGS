#include <iostream>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <string> 
#include <ctime>
#include <sstream>
#include <cmath>
#include <climits> // for INT_MAX
#include <cfenv> // for fmod

// for M_PI
#define _USE_MATH_DEFINES	
#include <math.h>

using namespace std;

extern float COST_TOFFOLI;

template <typename T>
void print(T t) {
	cout << t << endl;
}
template<typename T, typename... Args>
void print(T t, Args... args) { // recursive variadic function
	cout << t << ' ';
	print(args...);
}

//============================================================================

enum BITTYPE {
	POS,
	NEG,
	CAR		// carry bit
};

enum GATETYPE {
	RX, RY, RZ, RXX, RYY, RZZ, P, CP
};

class Gate {
public:
	Gate(int id, GATETYPE type, const vector<int>& qubits) : _id(id), _type(type), _qubits(qubits) {};
	int getId() { return _id; }
	string getTypeStr() {
		switch (_type) {
			case GATETYPE::RX:
				return "rx";
			case GATETYPE::RY:
				return "ry";
			case GATETYPE::RZ:
				return "rz";
			case GATETYPE::RXX:
				return "rxx";
			case GATETYPE::RYY:
				return "ryy";
			case GATETYPE::RZZ:
				return "rzz";
			case GATETYPE::P:
				return "p";
			case GATETYPE::CP:
				return "cp";
			default:
				return "N/A";
		}
	}
	int getQubit(int index) const { return _qubits[index]; }
	void setName(const string& name) { _name = name; }
	const string& getName() const { return _name; }
private:
	int _id;
	int _type;
	vector<int> _qubits;
	float _value = 0;  // only for single-gate
	string _name;
	//bool _activate = true;
};

class Bit {
public:
	Bit(BITTYPE type, Gate* gate) : _type(type), _gate(gate) {}
	Bit(const vector<Bit>& carry_ins, int power) : _type(BITTYPE::CAR), _gate(nullptr), _carry_ins(carry_ins), _power(power) {} //, _index(index) {}

	const BITTYPE getType() const { return _type; }
	void invType() {
		if (_type == BITTYPE::POS)		_type = BITTYPE::NEG;
		else if (_type == BITTYPE::NEG) _type = BITTYPE::POS;
		else							assert(false);	// should not be BITTYPE::CAR
	}
	BITTYPE getInvType() {
		if (_type == BITTYPE::POS) return BITTYPE::NEG;
		if (_type == BITTYPE::NEG) return BITTYPE::POS;
		assert(false);	// should not be BITTYPE::CAR
	}
	char getTypeChr() {
		if (_type == BITTYPE::POS) return '+';
		if (_type == BITTYPE::NEG) return '-';
		return 'c';
	}
	int getGateId() { return _gate->getId(); }
	const string& getName() const { return _gate->getName(); }
	bool isActivate() { return _is_activate; }// _gate->isActivate();
	void setInactivate() { _is_activate = false; }
	bool isPos() { return (_type == BITTYPE::POS); }
	bool isNeg() { return (_type == BITTYPE::NEG); }

	const vector<Bit>& getCarryIns() const { return _carry_ins; }	// for counter bits
	int getPower() const { return _power; }	// for counter bits

private:
	BITTYPE _type;
	Gate* _gate;
	bool _is_activate = true;
	vector<Bit> _carry_ins;		// for counter bits
	int _power = 0;					// for counter bits
};

class Optimizer {
public:
	// defined in 'optimize.cpp'
	Optimizer(int precision, float cost_single = INT_MAX, bool is_same = false); // : _r(precision), _is_same(is_same), _cost_single(cost_single) {}
	pair<float, int> optimize(bool to_print_info = false);
	void concrete();
	
	// defined in 'io.cpp'
	void importQasm(const string& file_name);
	
	float exportQasm(const string& file_name);
	float importBitList(const string& file_name);
	void printInfo(const string& header = "");
private:
	int _n;				// number of gates = _gate_list.size()
	int _r;				// number of bits (precision)
	bool _is_same;		// special mode for synthesize same angles
	double _last_angle;
	vector<Gate*> _gate_list;
	vector<vector<Bit>> _bit_table;		// _r * _n
	float _cost_single;
	set<int> _involved_qubits_x;
	set<int> _involved_qubits_y;
	set<int> _involved_qubits_z;
	vector<int> _heights;			// _heights[i] = _bit_table[i].size() - _n_split_from[i] + _n_split_to[i] + _n_carry[i] - _n_counter[i] + _counter_sizes[i].size() - #excluded
	int _max_height;
	vector<int> _n_split_from;
	vector<int> _n_split_to;
	vector<int> _n_carry;
	vector<int> _n_counter;					// _n_counter[i] = sum(_counter_sizes[i])
	vector<vector<int>> _counter_sizes;		// each in the decreasing order
	vector<pair<Gate*, int>> single_gates;
	unordered_map<int, float> _excluded;
	vector<string> _headers;
	float _cost = 0;
	
	// defined in 'optimize.cpp'
	void updatePeaks(vector<int>& peaks);

	bool split(int index, int index_bound);
	int findSplittedGate(int index, unordered_set<int>& pos_gates, unordered_set<int>& neg_gates, int index_bound);
	void splitGate(int index, int gate_id);

	int doCounter(vector<int>& new_height, vector<int>& new_n_carry, vector<int>& new_n_counter, vector<vector<int>>& new_counter_sizes, const vector<int>& peaks, int& dealing_peak_index);
	int mergeCounter(vector<int>& new_height, vector<int>& new_n_carry, vector<int>& new_n_counter, vector<vector<int>>& new_counter_sizes, const vector<int>& peaks, int& dealing_peak_index);
	int findTargetCounterMin(const vector<int>& counter_sizes);
	
	int doSingle(unordered_set<int>& new_excluded, vector<int> peaks_remaining);
	void removeExcluded();

	void splitGateAny(int index);

	// defined in 'io.cpp'
	void exportQasmFourierTrans(ofstream& ofs, bool is_reverted);
	void exportQasmRotTypeTrans(ofstream& ofs, bool is_reverted);
	void exportQasmSetAnc(ofstream& ofs, bool is_reverted);
	void exportQasmWriteAdder(ofstream& ofs);
	int exportQasmSetAdderBits(ofstream& ofs, int ith_adder, bool is_reverted);
	void exportCounter(ofstream& ofs, const vector<Bit>& carry_ins, vector<int>& selected, int k, string& target_name, bool is_reverted);
	void exportQasmWriteSingle(ofstream& ofs);
};

// defined in 'external.cpp'
int nCr(int n, int k);
int countAdderCost(int min_bit);
int countCounterCost(int counter_size, int dis_to_head);
void boothEncode(vector<int>& bit_string);