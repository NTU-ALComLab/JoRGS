#include <boost/program_options.hpp>
#include "headers.h"

float getCostSingle(int precision) {
  	float HST = 3.0 * (precision + 1) + log2(precision + 1);
  	float RUS = 1.149 * (precision + 1) + 9.2;
  	float PQF = 1.0 * (precision + 1) + 4 * log2(precision + 1) + 1.187;
  	float cost_single = min(HST, min(RUS, PQF));
  	return cost_single;
}

int main(int argc, char** argv) {
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()
        ("help", "produce help message")
        ("in",  po::value<string>(), "qasm file string for synthesis")
        ("out", po::value<string>(), "qasm file string after synthesis")
        ("prec", po::value<unsigned int>()->default_value(30), "precision in bits (default: 30)")
        ("cost", po::value<double>()->default_value(1000), "T-count of applying an independent single-gate rotation. Set it to a large number to disable applying single-gate rotations (default: 1000)")
        ("same", "use Fourier state transformation for the same-angle special case")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);
    
    if (vm.count("help") || !vm.count("in") || !vm.count("out")) {
  	    std::cout << description << std::endl;
  	    return 1;
	  }
    
    string in_cir  = vm["in"].as<string>();
    string out_cir = vm["out"].as<string>();
    int prec = vm["prec"].as<unsigned int>();
    int cost = vm["cost"].as<double>();
    bool is_same = (bool)vm.count("same");

    Optimizer op(prec, cost, is_same);
		op.importQasm(in_cir);
		op.optimize();
		op.concrete();
		float t_count = op.exportQasm(out_cir);
		cout << "Finished. Final T-count = " << t_count << endl;
    
	  return 0;
}


