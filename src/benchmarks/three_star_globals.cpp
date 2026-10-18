// Upstream 3-star requires process-global state in exactly one translation
// unit.
#include <pixie/rank_select/three_star.h>

#include <chrono>
#include <string>
#include <vector>

uint64_t SAMPLEDIST_EXP = 0;
uint64_t SAMPLEDIST = 0;

bool show_overhead = false;
std::string BITGEN;
std::string BVNAME;
bool doBench = false;
bool doValid = false;
bool doValidw = false;
bool doMultibench = false;
uint64_t synth_bits = 0;
bool useSynth = false;
bool do_manual_bits = false;
bool do_fractals = false;
bool do_chunky = false;
bool do_smooth = false;
bool do_quicksmooth = false;
bool do_gidneysmooth = false;
bool do_alternate = false;
bool do_every_kth = false;
bool do_bimodal = false;
bool do_patches = false;
bool do_sinus = false;
int bimodal_factor = 0;
int every_kth_spacing = 0;
int synth_seqsize = 0;
int synth_01ratio = 0;
bool invert_all_bits = false;
bool do_random_queries = false;
uint64_t synth_accesses = 0;
uint64_t synth_ranks = 0;
uint64_t synth_selects = 0;
uint64_t synth_select1s = 0;
uint64_t synth_select0s = 0;
std::vector<uint64_t> query_type_counter;
uint64_t N_queries = 0;
uint64_t N_bits = 0;
uint64_t N_ones = 0;
uint64_t N_zeros = 0;
uint64_t N_words = 0;
uint64_t seed = 0;
double eff_01ratio = 0.0;
int parameter_cycles = 1;
int curr_parameter_cycle = 0;
int curr_instance = 0;
double sin_threshhold = 0.0;
std::vector<uint64_t, AlignedAllocator<uint64_t, ALIGNMENT>> bits;
std::vector<uint64_t, AlignedAllocator<uint64_t, ALIGNMENT>> queries;
std::vector<uint64_t, AlignedAllocator<uint64_t, ALIGNMENT>>
    bits_64bit_reversed;
uint64_t total_used_space_in_bits = 0;
tp tCopyFinished;
double dRead = 0.0;
double dCopy = 0.0;
double dBuild = 0.0;
double dQueries = 0.0;
double dPerQuery = 0.0;

std::string get_unixtimestamp() {
  return {};
}
