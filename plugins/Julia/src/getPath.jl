using TOML

parsed = TOML.parsefile("plugins/Julia/Gismo_path.toml")
path_to_lib = parsed["library_path"]