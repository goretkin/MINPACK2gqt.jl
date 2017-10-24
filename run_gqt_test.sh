#!/bin/bash

(julia test_print_gqt.jl F) | (tr -s " ") > outF.txt
(julia test_print_gqt.jl J) | (tr -s " ") > outJ.txt
