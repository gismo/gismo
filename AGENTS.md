# AGENTS.md — G+Smo

Instructions for AI coding agents (and humans) contributing to G+Smo, a C++
library for isogeometric analysis. This file holds the code conventions every
contributor should follow, whatever tooling they use.

> **Using the [gsAgents](https://github.com/gismo/gsAgents) Claude Code plugin?**
> Its skills and agents (`/gismo:build-target`, `/gismo:run-tests`,
> `/gismo:syntax-check`, `/gismo:plan`, `/gismo:implement`, …) take precedence
> for *how* to build, test, plan and review. This file only sets *what* the code
> should look like, and the plugin enforces the same rules.

## Repository layout

- `src/<module>/` — core library, e.g. `src/gsCore`, `src/gsNurbs`, `src/gsAssembler`.
- `optional/<module>/` — optional submodules (`gsKLShell`, `gsElasticity`, …).
  **Each is a separate git repository**: commit inside the submodule, never from
  the root repo. Enable them via `GISMO_OPTIONAL` or `submodules.txt`.
- `examples/`, `optional/<module>/examples/` — one `.cpp` file per executable.
- `unittests/`, `optional/<module>/unittests/` — UnitTest++ suites, all compiled
  into a single `unittests` binary.
- `filedata/` — XML input data (geometries, PDE definitions).

## C++ conventions

### Files and layout

- Class and file names start with `gs`, e.g. `gsBasis` in `gsBasis.h`. All library
  code lives in `namespace gismo`.
- Templates are split over three files:
  - `gsFoo.h` holds the declaration.
  - `gsFoo.hpp` holds the implementation.
  - `gsFoo_.cpp` holds the explicit instantiations, `CLASS_TEMPLATE_INST gsFoo<real_t>;`.
- Non-template free functions are declared with `GISMO_EXPORT` in the `.h` and
  defined in a `.cpp`. Symbols are hidden by default (`-fvisibility=hidden`).
- Use `#pragma once`. Headers must be self-contained, so they compile on their own.
- Every file starts with the standard header:

  ```cpp
  /** @file gsFoo.h

      @brief One-line summary.

      This file is part of the G+Smo library.

      This Source Code Form is subject to the terms of the Mozilla Public
      License, v. 2.0. If a copy of the MPL was not distributed with this
      file, You can obtain one at http://mozilla.org/MPL/2.0/.

      Author(s): ...
  */
  ```

### Formatting

- Follow `.clang-format`: LLVM base, 4-space indent, Allman braces, `T* p`
  (pointer to the type), no namespace indentation.
- When editing existing code, match the surrounding file instead of reformatting
  it. Keep diffs focused.

### Types and idioms

- Use `real_t` for scalars and `index_t` for sizes and indices, not hard-coded
  `double`/`int`. The library builds
  with float, double and multiprecision `real_t`, so do not assume `double`.
- Use `gsMatrix`/`gsVector` (Eigen-based). Prefer Eigen block and vectorised
  operations over element-wise loops in performance-critical code.
- Use `give(x)` (from `gsCore/gsMemory.h`), not `std::move(x)`.
- Use the class smart-pointer typedefs `typename gsFoo::uPtr` and `::Ptr`
  (`memory::unique_ptr` and `memory::shared_ptr`).
- For output, use `gsInfo`, `gsWarn` and `gsDebug`, never `std::cout`/`std::cerr`
  directly.
- Error macros:
  - `GISMO_ASSERT(cond, msg)` is debug-only. Use it for internal invariants and in
    hot loops.
  - `GISMO_ENSURE(cond, msg)` is always checked. Use it to validate user input at
    API boundaries, not inside inner loops.
  - `GISMO_ERROR(msg)` is for unreachable or unsupported paths.
- Avoid dynamic allocation and exceptions inside assembly or evaluation loops.

### Numerics

- Guard against degenerate input: empty or zero-size matrices, degree 0, a single
  knot span, coincident points, mismatched dimensions, and wrong orientation.
- Watch for cancellation, unguarded division, `index_t` overflow and silent
  narrowing between `real_t` and `index_t`.
- When complexity is not obvious, state it in a comment.

## Documentation and comments

- Doxygen on everything public: `\brief`, `\param`, `\return`, `\tparam`, `\sa`;
  maths as `\f$ ... \f$`. Keep the `@file` header intact.
- Link theory to code. Name the method, scheme or paper a solver or assembler
  implements, and never invent citations.
- Document contracts, not the obvious: units, index conventions, tensor or matrix
  shapes, ownership, valid ranges, I/O formats, and why a non-obvious formulation
  is the stable or correct one.
- **Comments explain the code. Commit messages and PRs explain the change.** Never
  write the following into source:
  - Narration of the change, e.g. "removed the old loop" or "previously this used `gsFoo`".
  - Task or review scaffolding, e.g. "added for step 3", "per the spec" or `TODO(review)`.
  - Comments that restate what the code says, or first-person hedging.
  - Commented-out code you replaced; it is in git.

  The test: *would this comment still be true and useful to someone who never saw
  the diff?* Real `TODO`/`FIXME` notes and warnings about real traps are welcome.

## Unit tests

- Use UnitTest++ via `#include "gismo_unittest.h"`. `unittests/gsTutorial.cpp` is
  the reference.
- Put `SUITE(gsFoo_test)` in `gsFoo_test.cpp`, so the suite name matches the file
  name. Use `TEST(descriptive_name)` with `CHECK`, `CHECK_EQUAL`, `CHECK_CLOSE`,
  `CHECK_ARRAY_CLOSE` and `CHECK_THROW`.
- Compare against reference solutions: analytic values, manufactured solutions or
  convergence rates. Never paste the code's own output back in as the expected
  value.
- Pick tolerances that are tight enough to fail when the result is wrong. Prefer
  tolerances scaled by `math::limits::epsilon()` over arbitrary constants.
- Keep suites fast: coarse meshes and few refinements, seconds rather than minutes.
- Bug fixes come with a regression test that fails on the unfixed code.

## Examples

- One `examples/foo.cpp` builds the make target `foo`. Start from the closest
  existing example.
- Parse arguments with `gsCmdLine`. Defaults must run in seconds without arguments,
  so heavy resolutions are opt-in.
- Read input from `filedata/` via `gsFileData` or `gsReadFile`, and print with
  `gsInfo`. ParaView output (`gsWriteParaview`) goes behind a `--plot` switch,
  off by default.
- Convergence studies print an error and rate table and state the expected order.

## Building and testing

These rules apply without the plugin. Plugin users should use its skills, which
enforce them.

- **Never run bare `make`**, since it builds every example. Build named targets:
  `make <target> -j4`, e.g. `make gismo`, `make unittests` or `make poisson_example`.
- Keep `-j` moderate (about `nproc/2` at most). Unbounded parallelism can exhaust
  RAM.
- After adding a new `.cpp` (example, test or instantiation file), run `cmake .`
  in the build directory so the target is picked up.
- Unit tests need `GISMO_BUILD_UNITTESTS=ON`. Run `./bin/unittests` for all tests,
  or `./bin/unittests <prefix>` to run suites, tests or files by name prefix.
  This is UnitTest++, not doctest, so there is no `-R` or `--test-case`.
- Build and run the affected unit tests or example before opening a PR.

## Pull requests

- PRs need at least one approving review from `@gismo/admins`.
- Write the commit and PR description as a list of prefixed lines:

  ```
  NEW: ...
  IMPROVED: ...
  FIXED: ...
  API: ...
  ```

- New code needs Doxygen documentation and tests or an example.
- Do not commit build directories, generated files or local agent state (`.claude/`).
