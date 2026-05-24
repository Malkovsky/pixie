# setup-cpp-repo Examples

This file contains concrete examples for the C++ project scaffold skill. Keep
general setup workflow in `SKILL.md`.

## Succinct Library Scaffold

Create a header-only C++20 library scaffold for a succinct data structures
project:

```bash
python3 .kilo/skills/setup-cpp-repo/scripts/init_cpp_project.py \
    --name succinct-lib --namespace succinct --output-dir .
```

Then verify the generated project:

```bash
cd succinct-lib
cmake --preset release
cmake --build --preset release -j
./build/release/unittests
```
