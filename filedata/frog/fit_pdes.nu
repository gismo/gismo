#!/usr/bin/nu

# A nushell script that run the frog_pde binary for all the different PDEs in the pde_input folder, fitting them for error analysis.
ls ../filedata/frog/pde-inputs/ | get name | each {|file|
  bin/frog_pde --steps 5 --frogdir "frog/quad/" --input $file
}
