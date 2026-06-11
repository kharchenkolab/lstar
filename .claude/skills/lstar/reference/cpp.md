# C++ core: libstar (reference)

Header-only, MIT-licensed core in `core/include/lstar/lstar.hpp` (vendors `nlohmann/json.hpp`). It
is the shared engine: the R package links it via cpp11, and the Python accelerator binds it via
pybind11. Build with C++17; OpenMP for the kernels; optional zlib (`-DLSTAR_HAVE_ZLIB`) for gzip.

## Model structs

```cpp
struct NdArray { std::string dtype; std::vector<int64_t> shape; std::vector<uint8_t> bytes;
                 int64_t nelem(); template<class T> const T* as() const; };
struct Axis  { std::string name, origin, role; std::vector<std::string> labels; json provenance; };
struct Field { std::string name, role, encoding, state, subtype; std::vector<std::string> span;
               NdArray dense; std::vector<std::string> strings;      // utf8
               NdArray data, indices, indptr; std::vector<int64_t> shape; };  // csr/csc
struct Dataset { std::string kind, spec_version; std::vector<std::string> profiles, dropped;
                 std::vector<Axis> axes; std::vector<Field> fields;
                 Axis* axis(name); Field* field(name); };
```

## Zarr IO

```cpp
lstar::Dataset ds = lstar::read(root);     // reads an arbitrary chunk grid, C-order, fill-padded
                                           // edge chunks, missing chunk -> fill 0; gzip/zlib chunks
lstar::write(ds, root);                    // single-chunk arrays + a consolidated .zmetadata
```
- `read_array(dir)` assembles a (possibly chunked, possibly compressed) zarr v2 array into a
  contiguous C-order buffer. `decode_chunk` / `inflate_stream` handle gzip/zlib (zlib windowBits 47,
  auto-detects either header) when built with `LSTAR_HAVE_ZLIB`.
- Strings via `read_strings`/`write_strings` (UTF-8 bytes + int64 offsets).

## dtype normalization

Index and value arrays come in multiple widths. Normalize before the kernels:
```cpp
std::vector<int64_t> as_i64(const NdArray&);   // <i4 / <i8 -> int64  (indptr/indices)
std::vector<double>  as_f64(const NdArray&);   // <f4 / <f8 -> double (values; or keep native)
```
Measures are commonly **float32** (`<f4`) in AnnData/Seurat — the kernels are templated on the value
dtype so you read float32 in place (no widening copy) and accumulate in double.

## Translation primitives (OpenMP)

```cpp
struct ColStats { std::vector<double> mean, var; std::vector<int64_t> nnz; };

template<class T>
ColStats csc_col_mean_var(const T* data, const int64_t* indptr,
                          int64_t ncols, int64_t nrows, int n_threads = 0, bool lognorm = false);
//   zero-aware per-column mean/variance of a CSC matrix; lognorm applies log1p per nonzero on the
//   fly (pagoda2 lazy-view pattern). n_threads: 1 serial, N, <=0 = OpenMP default. Thread-invariant.

template<class T>
struct CsxArrays { std::vector<T> data; std::vector<int64_t> indices, indptr; int64_t nrows, ncols; };
template<class T>
CsxArrays<T> csc_to_csr(const T* data, const int64_t* indices, const int64_t* indptr,
                        int64_t nrows, int64_t ncols);
//   O(nnz) storage transpose (orientation flip genes<->cells) preserving the value dtype.
```

Use these for the genes×cells reorientation the Seurat/SCE profiles need, and for per-gene
dispersion stats without densifying.

## Build & test

```
cmake -S core -B core/build -DCMAKE_BUILD_TYPE=Release && cmake --build core/build -j
core/build/test_crossimpl    # reads a Python-written store, validates, round-trips
core/build/test_chunked      # reads a chunked+gzip store; transpose + threaded reduction
core/build/bench_colstats    # OpenMP scaling of csc_col_mean_var
```
CMake finds OpenMP and zlib; defines `LSTAR_HAVE_ZLIB` when zlib is present.

## Still planned
blosc2 codec; Zarr v3 + sharding; multithreaded/out-of-core full transpose and collection gather;
a Python↔C++ zero-copy boundary for whole-store streaming.
