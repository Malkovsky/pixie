# Rank/select {#benchmarks_rank_select}

## Machine

AMD Ryzen 7 8845HS (8 cores, 16 threads), 15 GiB RAM, WSL2
(Linux 6.6.87.2-microsoft-standard-WSL2). GCC 13.3, C++20,
`-O3 -DNDEBUG -march=native`; benchmarks pinned to CPU 0.
Results captured October 13–16, 2026; SDSL was added in a separate run.

<details>
<summary>lscpu output</summary>

```text
Architecture:                         x86_64
CPU op-mode(s):                       32-bit, 64-bit
Address sizes:                        48 bits physical, 48 bits virtual
Byte Order:                           Little Endian
CPU(s):                               16
On-line CPU(s) list:                  0-15
Vendor ID:                            AuthenticAMD
Model name:                           AMD Ryzen 7 8845HS w/ Radeon 780M Graphics
CPU family:                           25
Model:                                117
Thread(s) per core:                   2
Core(s) per socket:                   8
Socket(s):                            1
Stepping:                             2
BogoMIPS:                             7585.20
Flags:                                fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 ht syscall nx mmxext fxsr_opt pdpe1gb rdtscp lm constant_tsc rep_good nopl tsc_reliable nonstop_tsc cpuid extd_apicid tsc_known_freq pni pclmulqdq ssse3 fma cx16 sse4_1 sse4_2 movbe popcnt aes xsave avx f16c rdrand hypervisor lahf_lm cmp_legacy cr8_legacy abm sse4a misalignsse 3dnowprefetch osvw topoext perfctr_core ssbd ibrs ibpb stibp vmmcall fsgsbase bmi1 avx2 smep bmi2 erms invpcid avx512f avx512dq rdseed adx smap avx512ifma clflushopt clwb avx512cd sha_ni avx512bw avx512vl xsaveopt xsavec xgetbv1 xsaves avx512_bf16 clzero xsaveerptr arat avx512vbmi umip avx512_vbmi2 gfni vaes vpclmulqdq avx512_vnni avx512_bitalg avx512_vpopcntdq rdpid fsrm
Hypervisor vendor:                    Microsoft
Virtualization type:                  full
L1d cache:                            256 KiB (8 instances)
L1i cache:                            256 KiB (8 instances)
L2 cache:                             8 MiB (8 instances)
L3 cache:                             16 MiB (1 instance)
NUMA node(s):                         1
NUMA node0 CPU(s):                    0-15
Vulnerability Gather data sampling:   Not affected
Vulnerability Itlb multihit:          Not affected
Vulnerability L1tf:                   Not affected
Vulnerability Mds:                    Not affected
Vulnerability Meltdown:               Not affected
Vulnerability Mmio stale data:        Not affected
Vulnerability Reg file data sampling: Not affected
Vulnerability Retbleed:               Not affected
Vulnerability Spec rstack overflow:   Vulnerable: Safe RET, no microcode
Vulnerability Spec store bypass:      Mitigation; Speculative Store Bypass disabled via prctl
Vulnerability Spectre v1:             Mitigation; usercopy/swapgs barriers and __user pointer sanitization
Vulnerability Spectre v2:             Mitigation; Retpolines; IBPB conditional; IBRS_FW; STIBP always-on; RSB filling; PBRSB-eIBRS Not affected; BHI Not affected
Vulnerability Srbds:                  Not affected
Vulnerability Tsx async abort:        Not affected
```

</details>

## Reproduction

Run from the repository root with CMake 3.25+ and Python with Matplotlib:

```bash
cmake --preset benchmark-all-backends
cmake --build --preset benchmark-all-backends --target rank_select_benchmarks -j 2
export PIXIE_RANK_SELECT_FIXED_SEED=1
bench=build/benchmark-all-backends/rank_select_benchmarks
out=build/rank-select-reproduction
mkdir -p "$out"

for scale in small large; do
  sizes='1024|16384|262144|4194304|67108864|1073741824'
  query_time=1s; build_time=1s
  if [ "$scale" = large ]; then
    sizes=17179869184; query_time=20000000x; build_time=2x
  fi
  taskset -c 0 "$bench" \
    --benchmark_filter="^rank_select.*_(rank1|rank0|select1|select0)_(12p5|50|87p5)/N:($sizes)(/|$)" \
    --benchmark_repetitions=3 --benchmark_min_time="$query_time" \
    --benchmark_out="$out/queries-$scale.json" --benchmark_out_format=json
  taskset -c 0 "$bench" \
    --benchmark_filter="^rank_select.*_build_both_(12p5|50|87p5)/N:($sizes)(/|$)" \
    --benchmark_repetitions=10 --benchmark_min_time="$build_time" \
    --benchmark_out="$out/build-$scale.json" --benchmark_out_format=json
done
python3 scripts/plot_rank_select.py "$out"/*.json -o "$out/graphs"
```

Median wall times; auxiliary space is index bytes / input bytes, excluding
owned source copies. Build time includes those copies; footprint inspection is
untimed in the current harness. SDSL uses rank-v and MCL select for both bit
values and reports serialized directory sizes. FastBitVec stops at 2^30;
3-star and Orzo require BMI/BMI2. Missing select0 curves reflect upstream APIs.

## Results

<img src="rank_select_rank1.svg" alt="rank1" />

<img src="rank_select_rank0.svg" alt="rank0" />

<img src="rank_select_select1.svg" alt="select1" />

<img src="rank_select_select0.svg" alt="select0" />

<img src="rank_select_build.svg" alt="Construction" />

<img src="rank_select_auxiliary.svg" alt="Auxiliary space" />
