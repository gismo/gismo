#!/usr/bin/nu

bin/frog_pde --steps 5 --frogdir "frog/lagrange/" --input ../filedata/frog/pde-inputs/trig_Val3_Scale35.xml --errors $"errors/pde_trig3.csv" --condition $"errors/cond_pde_trig3.csv"
bin/frog_pde --steps 5 --frogdir "frog/lagrange/" --input ../filedata/frog/pde-inputs/trig_Val5_Scale35.xml --errors $"errors/pde_trig5.csv" --condition $"errors/cond_pde_trig5.csv"
bin/frog_pde --steps 5 --frogdir "frog/lagrange/" --input ../filedata/frog/pde-inputs/trig_Val6_Scale35.xml --errors $"errors/pde_trig6.csv" --condition $"errors/cond_pde_trig6.csv"
