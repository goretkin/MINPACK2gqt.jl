#!/bin/bash

(julia test_print_estsv.jl F) | (tr -s " ") > out.F
(julia test_print_estsv.jl J) | (tr -s " ") > out.J
