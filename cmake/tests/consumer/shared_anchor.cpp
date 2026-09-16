#include <gsdeptest/dep_lib.h>

// Exported in both the sources_export_symbols and sources_hidden_symbols arms,
// giving nm a positive control inside the same dynamic symbol table: an empty
// `nm -D` output must not make the hidden arm pass vacuously.
int gsdeptest_shared_anchor() { return gsdeptest_lib_answer(); }
