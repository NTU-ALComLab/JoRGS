#include "headers.h"

/* ===== Function Description:
	Read from a bit list file.	(for testing the program)
	Side feedback: return the cost of pure-adder method.
*/
float Optimizer::importBitList(const string& file_name) {
	float adder_cost = 0;

	ifstream in_file(file_name, ios::in);
	if (!in_file.good()) {
		cerr << "File \"" << file_name << "\" is not found\n";
		exit(-1);
	}

	string line;
	while (getline(in_file, line)) {
		Gate* new_gate = new Gate((int)_gate_list.size(), GATETYPE::RZ, {});
		_gate_list.emplace_back(new_gate);

		stringstream line_ss(line);
		string word;
		int i = 0;
		int lsb = 0;
		while (getline(line_ss, word, ' ')) {
			if (word == "1") {
				if (!_is_same) {
					_bit_table[i].emplace_back(Bit(BITTYPE::POS, new_gate));
					_heights[i]++;
				}
				lsb = i;
			}
			else if (word == "-1") {
				if (!_is_same) {
					_bit_table[i].emplace_back(Bit(BITTYPE::NEG, new_gate));
					_heights[i]++;
				}
				lsb = i;
			}
			i++;
		}

		if (_is_same) {
			_bit_table[lsb].emplace_back(Bit(BITTYPE::POS, new_gate));
			_heights[lsb]++;
		}

		adder_cost += countAdderCost(lsb);
	}

	_n = _gate_list.size();

	if (_is_same) {
		for (int i = 0; i < _r; ++i) {
			if (_bit_table[i].size() != 0) {
				_r = i + 1;
				while (_heights.size() > _r) {
					_heights.pop_back();
					_n_carry.pop_back();
					_n_counter.pop_back();
					_n_split_from.pop_back();
					_n_split_to.pop_back();
					_counter_sizes.pop_back();
					_bit_table.pop_back();
				}
				break;
			}
		}
	}

	return adder_cost;
}

/* ===== Function Description:
	Read from a openQASM file.
*/
void Optimizer::importQasm(const string& file_name) {		// support gate set: RX, RY, RZ, RXX, RYY, RZZ, P, CP // GATETYPE
	// parse qasm files (read bit table)
	bool first_angle = true;

	ifstream in_file(file_name, ios::in);
	if (!in_file.good()) {
		cerr << "File \"" << file_name << "\" is not found\n";
		exit(-1);
	}

	string line;
	vector<int> bit_string(_r);
	while (getline(in_file, line)) {
		line = line.substr(0, line.find("//"));
		replace(line.begin(), line.end(), '(', ' ');
		replace(line.begin(), line.end(), ')', ' ');
		if (line.find_first_not_of("\t\n ") == string::npos) continue;

		stringstream line_ss(line);
		string word;
		getline(line_ss, word, ' ');

		if (word == "rx" || word == "ry" || word == "rz" || word == "rxx" || word == "ryy" || word == "rzz" || word == "p" || word == "cp") {
			// gate type
			GATETYPE gate_type;
			if (word == "rx")		gate_type = GATETYPE::RX;
			else if (word == "ry")	gate_type = GATETYPE::RY;
			else if (word == "rz")	gate_type = GATETYPE::RZ;
			else if (word == "rxx") gate_type = GATETYPE::RXX;
			else if (word == "ryy") gate_type = GATETYPE::RYY;
			else if (word == "rzz") gate_type = GATETYPE::RZZ;
			else if (word == "p")	gate_type = GATETYPE::P;
			else if (word == "cp")	gate_type = GATETYPE::CP;

			// rotation angle
			getline(line_ss, word, ' ');
			double angle = stod(word);
			angle = angle / M_PI / 2 + 1 + pow(2, -1 - _r);  // for rounding
			if (angle >= 1) angle -= 1;		// between 0 and 1
			if (_is_same == true && first_angle == false && angle != _last_angle) {
				cerr << "All angles must be the same under the --all_same mode." << endl;
				exit(-1);
			}
			first_angle = false;
			_last_angle = angle;

			// qubits
			vector<int> qubits;
			getline(line_ss, word, '[');
			while (getline(line_ss, word, ']')) {
				int var = stoi(word);
				qubits.emplace_back(var);
				getline(line_ss, word, '[');

				if (gate_type == GATETYPE::RX || gate_type == GATETYPE::RXX) {
					if (_involved_qubits_y.count(var) > 0 || _involved_qubits_z.count(var) > 0) {
						cerr << "[Error]: Qubit " << var << " appears in gates with different rotation-axis type." << endl;
						exit(-1);
					}
					_involved_qubits_x.insert(var);
				}		
				else if (gate_type == GATETYPE::RY || gate_type == GATETYPE::RYY) {
					if (_involved_qubits_x.count(var) > 0 || _involved_qubits_z.count(var) > 0) {
						cerr << "[Error]: Qubit " << var << " appears in gates with different rotation-axis type." << endl;
						exit(-1);
					}
					_involved_qubits_y.insert(var);
				}
				else {
					if (_involved_qubits_x.count(var) > 0 || _involved_qubits_y.count(var) > 0) {
						cerr << "[Error]: Qubit " << var << " appears in gates with different rotation-axis type." << endl;
						exit(-1);
					}
					_involved_qubits_z.insert(var);
				}									
			}

			// process
			Gate* new_gate = new Gate((int)_gate_list.size(), gate_type, qubits);
			_gate_list.emplace_back(new_gate);

			for (int i = 0; i < _r; ++i) {
				angle = angle * 2;
				if (angle > 1) {
					bit_string[i] = 1;
					angle -= 1;
				}
				else {
					bit_string[i] = 0;
				}
			}

			if (_is_same) {
				int lsb;
				for (lsb = _r - 1; lsb >= 0; --lsb) {
					if (bit_string[lsb] == 1)
						break;
				}
				_bit_table[lsb].emplace_back(Bit(BITTYPE::POS, new_gate));
				_heights[lsb]++;
			}
			else {
				boothEncode(bit_string);
				for (int i = 0; i < _r; ++i) {
					if (bit_string[i] == 1) {
						_bit_table[i].emplace_back(Bit(BITTYPE::POS, new_gate));
						_heights[i]++;
					}
					else if (bit_string[i] == -1) {
						_bit_table[i].emplace_back(Bit(BITTYPE::NEG, new_gate));
						_heights[i]++;
					}
				}
			}
		}
		else if (word == "qreg" || word == "creg" || word == "OPENQASM" || word == "include") {
			_headers.emplace_back(line);
			continue;
		}
		else {
			cerr << "[Warning]: Syntax \'" << word << "\' is not supported in this simulator. The line is ignored ..." << endl;
		}
	}

	// initialization
	_n = _gate_list.size();

	// remove redundant bits for special case
	if (_is_same) {
		for (int i = 0; i < _r; ++i) {
			if (_bit_table[i].size() != 0) {
				_r = i + 1;
				while (_heights.size() > _r) {
					_heights.pop_back();
					_n_carry.pop_back();
					_n_counter.pop_back();
					_n_split_from.pop_back();
					_n_split_to.pop_back();
					_counter_sizes.pop_back();
					_bit_table.pop_back();
				}
				break;
			}
		}
	}
}

/* ===== Function Description:
	Do Fourier-state transformation for the special case.
*/
void Optimizer::exportQasmFourierTrans(ofstream& ofs, bool is_reverted) {
  double c = 1 - (int)(_last_angle * pow(2, _r));
  
	for (int i = 0; i < _r; ++i) {
    double angle = M_PI * c;//fmod(c, 2);
    
		if (is_reverted) {	// cost is not counted
			ofs << "p(" << -angle << ") frs[" << i << "];\n";
		}
		else {
			ofs << "p(" << angle << ") frs[" << i << "];\n";
			_cost += _cost_single;
		}
    c /= 2; 
	}
}


/* ===== Function Description:
	Do rotation-type transformation between x/y-type and z-type.
*/
void Optimizer::exportQasmRotTypeTrans(ofstream& ofs, bool is_reverted) {
	for (int qubit : _involved_qubits_x) {
		ofs << "h q[" << qubit << "];\n";
	}
	for (int qubit : _involved_qubits_y) {
		if (is_reverted) {
			ofs << "h q[" << qubit << "];\n";
			ofs << "s q[" << qubit << "];\n";
		}
		else {
			ofs << "sdg q[" << qubit << "];\n";
			ofs << "h q[" << qubit << "];\n";
		}
	}
}	

/* ===== Function Description:
	Set the representative ancilla qubits for two-qubit gates.
*/
void Optimizer::exportQasmSetAnc(ofstream& ofs, bool is_reverted) {
	int ith_anc = 0;
	for (Gate* gate : _gate_list) {
		if (gate->getTypeStr() == "rxx" || gate->getTypeStr() == "ryy" || gate->getTypeStr() == "rzz") {
			ofs << "cx q[" << gate->getQubit(0) << "], anc[" << ith_anc << "];\n";
			ofs << "cx q[" << gate->getQubit(1) << "], anc[" << ith_anc << "];\n";
			if (!is_reverted) gate->setName("anc[" + to_string(ith_anc) + "]");
			ith_anc++;
		}
		else if (gate->getTypeStr() == "cp") {
			ofs << "ccx q[" << gate->getQubit(0) << "], q[" << gate->getQubit(1) << "], anc[" << ith_anc << "];\n";
			if (!is_reverted) _cost += COST_TOFFOLI;
			if (!is_reverted) gate->setName("anc[" + to_string(ith_anc) + "]");
			ith_anc++;
		}
		else {
			if (!is_reverted) gate->setName("q[" + to_string(gate->getQubit(0)) + "]");
		}
	}
}

/* ===== Function Description:
	Write counter circuits.
*/
void Optimizer::exportCounter(ofstream& ofs, const vector<Bit>& carry_ins, vector<int>& selected, int k, string& target_name, bool is_reverted) {
	if (selected.size() == k) {
		set<string> pos_gates, neg_gates;
		for (int ith_bit : selected) {
			if (carry_ins[ith_bit].getType() == BITTYPE::POS) {
				pos_gates.insert(carry_ins[ith_bit].getName());
			}
			else {
				neg_gates.insert(carry_ins[ith_bit].getName());
			}
		}
		set<string> excluded;
		for (string s : pos_gates) {
			if (neg_gates.count(s) != 0) {
				excluded.insert(s);
			}
		}
		if (pos_gates.size() + neg_gates.size() - 2 * excluded.size() == 0) return;

		for (string s : neg_gates)
			if (excluded.count(s) == 0)
				ofs << "x " << s << ";\n";

		if (pos_gates.size() + neg_gates.size() - 2 * excluded.size() == 1)			ofs << "cx ";
		else if (pos_gates.size() + neg_gates.size() - 2 * excluded.size() == 2)	ofs << "ccx ";
		else 																		ofs << "mcx ";
			
		if (!is_reverted && pos_gates.size() + neg_gates.size() - 2 * excluded.size() > 1) _cost += COST_TOFFOLI;	// a bound
		// Note that by storing target bits of previous k/2-controlled Toffoli gates, 
		// k-controlled Toffoli gates with k > 2 can be obtained by a single 2-controlled Toffoli gate

		for (auto s : pos_gates)
			ofs << s << ", ";
		for (auto s : neg_gates)
			ofs << s << ", ";
		ofs << target_name << ";\n";

		for (string s : neg_gates)
			if (excluded.count(s) == 0)
				ofs << "x " << s << ";\n";
		return;
	}

	int start = selected.empty() ? 0 : (selected.back() + 1);
	for (int i = start; i < carry_ins.size(); ++i) {
		selected.emplace_back(i);
		exportCounter(ofs, carry_ins, selected, k, target_name, is_reverted);
		selected.pop_back();
	}
}

/* ===== Function Description:
	Set adder bits.
*/
int Optimizer::exportQasmSetAdderBits(ofstream& ofs, int ith_adder, bool is_reverted) {
	int last_bit = -1;
	for (int i = 0; i < _r; ++i) {
		if (_bit_table[i].size() > ith_adder) {
			last_bit = i;
			if (_bit_table[i][ith_adder].getType() == BITTYPE::POS) {
				ofs << "cx " << _bit_table[i][ith_adder].getName() << ", add[" << i << "];\n";
			}
			else if (_bit_table[i][ith_adder].getType() == BITTYPE::NEG) {
				ofs << "x add[" << i << "];\n";
				ofs << "cx " << _bit_table[i][ith_adder].getName() << ", add[" << i << "];\n";
			}
			else {	// BITTYPE::CAR
				vector<int> selected;
				string target = "add[" + to_string(i) + "]";
				exportCounter(ofs, _bit_table[i][ith_adder].getCarryIns(), selected, pow(2, _bit_table[i][ith_adder].getPower()), target, is_reverted);
			}
		}
	}
	return last_bit;
}

/* ===== Function Description:
	Write adders.
*/
void Optimizer::exportQasmWriteAdder(ofstream& ofs) {
	for (int ith_adder = 0; true; ++ith_adder) {
		int last_bit = exportQasmSetAdderBits(ofs, ith_adder, false);
		if (last_bit == -1) break;
		//ofs << "barrier;\n";

		// main adder
		for (int i = last_bit; i > 0; --i) { // MAJ
			ofs << "cx add[" << i << "], frs[" << i << "];\n";
			ofs << "cx add[" << i << "], add[" << i + 1 << "];\n";
			ofs << "ccx add[" << i + 1 << "], frs[" << i << "], add[" << i << "]; \n";
			_cost += COST_TOFFOLI;
		}
		ofs << "cx add[0], frs[0];\n";
		ofs << "cx add[1], frs[0];\n";
		for (int i = 1; i <= last_bit; ++i) { // UMS
			ofs << "ccx add[" << i + 1 << "], frs[" << i << "], add[" << i << "]; \n";
			ofs << "cx add[" << i << "], add[" << i + 1 << "];\n";
			ofs << "cx add[" << i + 1 << "], frs[" << i << "];\n";
		}

		//ofs << "barrier;\n";
		exportQasmSetAdderBits(ofs, ith_adder, true);	 // reverted
	}
}

/* ===== Function Description:
	Write excluded single rotation gates.
*/
void Optimizer::exportQasmWriteSingle(ofstream& ofs) {
	for (auto pair : _excluded) {
		Gate* gate = _gate_list[pair.first];
		float value = pair.second;
		ofs << "rz(" << value << ") " << gate->getName() << ";\n";
		_cost += _cost_single;
	}
}


/* ===== Function Description:
	Write the optimized circuit in openQASM format.
*/
float Optimizer::exportQasm(const string& file_name) {
	int n_ancilla = 0;
	for (Gate* gate : _gate_list) {
		if (gate->getTypeStr() == "rxx" || gate->getTypeStr() == "ryy" || gate->getTypeStr() == "rzz" || gate->getTypeStr() == "cp") {
			n_ancilla++;
		}
	}
	// use another ancilla qubit to represent the two-qubit gate

	ofstream ofs;
	ofs.open(file_name);
	for (string& line : _headers) {
		ofs << line << endl;
	}
	ofs << "qreg anc[" << n_ancilla << "];\n";
	ofs << "qreg add[" << _r + 1 << "];\n";
	ofs << "qreg frs[" << _r << "];\n";
	ofs << "// Notice: All Toffoli gates are recovered after the circuit,.\n";
	ofs << "//           and the method in [C. Gidney, 2018] can be applied.\n";
	ofs << "//         We use the method to calculate the T-count,\n";
	ofs << "//           but we keep the original circuit for clearity.\n";
	ofs << "//         Also, counter circuits can be easily simplified,\n";
	ofs << "//           but we keep the original circuit for clearity.\n";
	if (_is_same) ofs << "//         The cost of reverse Fourier state transform is not counted;\n";
	ofs << endl;

	if (_is_same) exportQasmFourierTrans(ofs, false);

	exportQasmRotTypeTrans(ofs, false);			// rotation type transformation
	exportQasmSetAnc(ofs, false);				// set representative ancilla qubits for two-qubit gates
	exportQasmWriteAdder(ofs);
	exportQasmSetAnc(ofs, true);				
	exportQasmRotTypeTrans(ofs, true);
	exportQasmWriteSingle(ofs);

	if (_is_same) exportQasmFourierTrans(ofs, true);

	return _cost;
}

/* ===== Function Description:
	Print current status.
*/
void Optimizer::printInfo(const string& header) {
	cout << header << endl;
	cout << "  excludedList: ";
	for (auto item : _excluded) {
		cout << endl << "    ID " << item.first << ", value = " << item.second;
	}
	if (_excluded.empty()) cout << "(empty)";
	cout << endl;

	cout << "  bitTable: " << endl;
	cout << "    MSB" << endl;
	for (int i = 0; i < _r; i++) {
		cout << "    ";
		if (_bit_table[i].empty()) {
			cout << "(empty)";
		}
		for (int j = 0; j < _bit_table[i].size(); j++) {
			cout << _bit_table[i][j].getTypeChr() << "";
			if (j % 10 == 9) cout << " ";
		}
		cout << endl;
	}
	cout << "    LSB" << endl;

	cout << "  Heights   : ";
	for (int h : _heights) {
		cout << h << " ";
	}
	cout << endl;
	cout << "  #SplitFrom: ";
	for (int n : _n_split_from) {
		cout << n << " ";
	} cout << endl;
	cout << "  #SplitTo  : ";
	for (int n : _n_split_to) {
		cout << n << " ";
	} cout << endl;
	cout << "  #Carry    : ";
	for (int n : _n_carry) {
		cout << n << " ";
	} cout << endl;
	cout << "  #Counter  : ";
	for (int n : _n_counter) {
		cout << n << " ";
	} cout << endl;

	cout << "  CounterSizeTable: " << endl;
	cout << "    MSB" << endl;
	for (int i = 0; i < _r; i++) {
		cout << "    ";
		if (_counter_sizes[i].empty()) {
			cout << "(empty)";
		}
		for (int j = 0; j < _counter_sizes[i].size(); j++) {
			cout << _counter_sizes[i][j] << " ";
		} cout << endl;
	}
	cout << "    LSB" << endl;

	cout << "==========================" << endl << endl;
}
