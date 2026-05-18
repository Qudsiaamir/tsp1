# Simplex TSP Solver

A C++ command-line Traveling Salesman Problem solver that uses simplex-based linear programming with an optional LaTeX/TikZ visualization demo.

The project implements the Dantzig-Fulkerson-Johnson style iterative approach for solving TSP instances. It also includes a USA map example that can render the current tour or relaxation solution into a PDF when `pdflatex` is available.

## Features

- Exact linear-programming based TSP workflow using a simplex implementation.
- Interactive commands for solving, adding subtour-elimination constraints, setting edge values, undoing constraints, and exiting.
- MST-based 2-approximation command for quickly generating a tour.
- Included USA demo with an adjacency matrix, TikZ map, and example PDF output.
- Simple Makefile build with no external C++ library dependencies.

## Tech Stack

- C++11
- Make
- Optional: LaTeX with `pdflatex` for rendering the visualization PDF

## Folder Structure

```text
.
├── Makefile
├── README.md
├── LICENSE
├── include/
│   ├── mst.h
│   └── simplex.h
├── src/
│   ├── mst.cpp
│   ├── simplex.cpp
│   └── tsp_solver.cpp
└── examples/
    └── usa-demo/
        ├── adjacency_matrix.txt
        ├── country.jpg
        ├── graph.tex
        ├── usa.tex
        └── usa.pdf
```

## Installation

Clone the repository and build the executable:

```sh
git clone <repository-url>
cd <repository-folder>
make
```

The compiled binary is written to:

```text
build/tsp_solver
```

By default the Makefile uses `g++`. To use another C++ compiler:

```sh
make CXX=clang++
```

## Run Locally

Start the solver:

```sh
make run
```

For the included USA demo, enter these values when prompted:

```text
examples/usa-demo/adjacency_matrix.txt
examples/usa-demo
graph.tex
usa.tex
```

The solver then accepts interactive commands:

```text
SOLVE
REM_LOOP N x_1 x_2 ... x_N
REM_LOOP_RNG lo hi
SET x_1 x_2 v
UNDO N
APPROX_MST
EXIT
```

## Example Usage

Run a quick startup check using the bundled demo data and immediately exit:

```sh
make smoke-test
```

Run the demo manually:

```sh
build/tsp_solver
```

Then provide the USA demo paths shown above. Try `APPROX_MST` for a fast approximate tour, or `SOLVE` to run the simplex solver on the current constraints.

When `SOLVE` runs, the program writes an LP dump to `lp.txt`, updates `examples/usa-demo/graph.tex`, and tries to rebuild `examples/usa-demo/usa.pdf` with `pdflatex`.

## Environment Variables

No environment variables are required.

## Testing and Verification

There is no formal test suite yet. Use these commands for local verification:

```sh
make
make smoke-test
```

Optional visualization check:

```sh
build/tsp_solver
```

Then enter the demo paths and run `APPROX_MST` or `SOLVE`. If `pdflatex` is installed, `examples/usa-demo/usa.pdf` should be regenerated.

## Troubleshooting

### `pdflatex` is not installed

The solver can still run, but PDF generation will fail with a warning. Install a LaTeX distribution such as TeX Live or MacTeX if you want visualization output.

### `make` cannot find a compiler

Install a C++ compiler such as `g++` or `clang++`, then rerun:

```sh
make
```

### Paths fail when running the example

Use the example folder path without quotes:

```text
examples/usa-demo
```

The program accepts this path with or without a trailing slash.

### Generated files appear after running

`lp.txt`, `build/`, and LaTeX auxiliary files are ignored by Git. The demo files `graph.tex` and `usa.pdf` are tracked because they are part of the included example.

## Roadmap

- Add small automated tests for the simplex and MST modules.
- Add non-interactive command-line arguments for easier scripting.
- Add a smaller example graph for quick full-solver demonstrations.
- Improve input validation and error messages for custom graph files.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).
