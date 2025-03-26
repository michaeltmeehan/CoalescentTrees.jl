# CoalescentTrees.jl

**CoalescentTrees.jl** is a Julia package for simulating genealogical trees under the coalescent model with a constant effective population size.

It supports both the **standard (unbounded)** coalescent model and the **bounded coalescent** model, which conditions on a root time (or bound) and sampling configuration.

---

## ✨ Features

- Efficient simulation of binary trees under:
  - **Standard coalescent** model
  - **Bounded coalescent** model (with forward-backward sampling)
- Support for heterochronous sampling (sequences sampled at multiple times)
- Clean separation of public API and internal logic
- Fully tested implementation with high test coverage

---

## 📦 Installation

```julia
] add path/to/CoalescentTrees
```

(Or clone the repository and include via `Pkg.develop`.)

---

## 🚀 Usage

```julia
using CoalescentTrees

# Sample from the standard coalescent
sampled_sequences = [1, 2]
sequence_times = [1.0, 2.0]
Ne = 10.0

tree = sample_tree(sampled_sequences, sequence_times, Ne)

# Sample from the bounded coalescent
bound_time = 0.0
bounded_tree = sample_tree(sampled_sequences, sequence_times, bound_time, Ne)
```

---

## 🧬 Tree Structure

Each call to `sample_tree` returns a `Tree` object with:

- `tree.time::Vector{Float64}`: time of each node
- `tree.left::Vector{Int}`: left child (0 if leaf)
- `tree.right::Vector{Int}`: right child (0 if leaf)

Node 1 is the root by construction. The remaining nodes follow reverse chronological order.

---

## 📖 Public API

| Function | Description |
|----------|-------------|
| `sample_tree(sampled_sequences, sequence_times, Ne)` | Sample tree from the **standard** coalescent model |
| `sample_tree(sampled_sequences, sequence_times, bound_time, Ne)` | Sample tree from the **bounded** coalescent model |
| `Tree` | Struct representing the tree topology and node times |

---

## 🧪 Testing

Run the full test suite:

```julia
import Pkg
Pkg.test("CoalescentTrees")
```

Tests cover all core functions, edge cases, and internal utilities.

---

## 📊 Benchmarking

To reproduce benchmarking results (e.g. distribution of coalescent times), see:

```
benchmarks/carson_et_al.jl
```

This script simulates many trees under each model and compares coalescence time distributions.

---

## 📂 Project Structure

```
src/
  CoalescentTrees.jl      # Main module
test/
  runtests.jl             # Entry point for all tests
  test_*.jl               # Unit tests
benchmarks/
  compare_models.jl       # Benchmarking and comparison with published results
```

---

## 📚 References

This package implements the bounded coalescent sampling algorithm as described in:

Carson, J. et al. (2022). *The bounded coalescent model: Conditioning a genealogy on a minimum root date*. Journal of Theoretical Biology, 548, 111165.  
[https://www.sciencedirect.com/science/article/pii/S0022519322001849?via%3Dihub#b0235](https://www.sciencedirect.com/science/article/pii/S0022519322001849?via%3Dihub#b0235)

The standard (unbounded) coalescent model is discussed in more detail in:
Kingman, J. F. C (1982). *The coalescent*. Stochastic Processes and Their Applications, 13, 235.
[https://www.sciencedirect.com/science/article/pii/0304414982900114](https://www.sciencedirect.com/science/article/pii/0304414982900114)

Tavare, S. (1984). *Line-of-descent and genealogical processes, and their applications in population genetics models*. Theoretical Population Biology, 26, 119.
[https://www.sciencedirect.com/science/article/pii/0040580984900273?via%3Dihub](https://www.sciencedirect.com/science/article/pii/0040580984900273?via%3Dihub)

---

## 🙏 Acknowledgements

The implementation of the bounded coalescent model, including the forward-filtering backward-sampling approach, is based entirely on the algorithm presented in the above paper by Carson et al. (2022). This package would not exist without their elegant and rigorous treatment of the problem.

---

## 🧱 License

MIT License. See `LICENSE` file for details.
