#include "headers.h"

float COST_TOFFOLI = 4;

/* ===== Function Description:  // O(k)
	Calculate the combination number of (n, k).
*/
int nCr(int n, int k) {
	int x = M_PI;
	if (n < k) {
		return 0;
	}
	if (k > n - k) {
		k = n - k;
	}

	int ans = 1;
	for (int i = 1; i <= k; i++) {
		ans *= n;
		ans /= i;
		n--;
	}
	return ans;
}

/* ===== Function Description:  // O(1)
	Calculate the cost of an adder circuit.
*/
int countAdderCost(int min_bit) {
	int n_toffoli = ((min_bit - 0 + 1) - 1);
	int t_count = n_toffoli * COST_TOFFOLI;
	return t_count;
}

/* ===== Function Description:  // O(min(log(counter_size), dis_to_head))
	Calculate the cost of a counter circuit.
	Note that by storing target bits of previous k/2-controlled Toffoli gates,
	k-controlled Toffoli gates with k > 2 can be obtained by a single 2-controlled Toffoli gate.
*/
int countCounterCost(int counter_size, int dis_to_head) {
	int n_toffoli = 0;

	int comb = 2;
	while (comb <= counter_size && dis_to_head > 0) {
		n_toffoli += nCr(counter_size, comb);
		comb *= 2;
		dis_to_head--;
	}
	// recover circuit has no cost by the previous method 

	return n_toffoli * COST_TOFFOLI;
}

/* ===== Function Description:
	Perform Booth encoding to a bit string.
*/
void boothEncode(vector<int>& bit_string) {
	bool flag = false;
	for (int i = bit_string.size() - 1; i >= 0; --i) {
		if (bit_string[i] == 0) {
			if (flag == true) {
				bit_string[i] = 1;
				flag = false;
			}
		}
		if (bit_string[i] == 1) {
			if (flag == true) {
				bit_string[i] = 0;
			}
			else {
				if (i > 0 && bit_string[i - 1] == 1) {
					flag = true;
					bit_string[i] = -1;
				}
			}
		}
	}
}
