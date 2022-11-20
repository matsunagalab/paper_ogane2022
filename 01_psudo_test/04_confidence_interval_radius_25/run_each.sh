#! /bin/bash

export PATH=$PATH:/bin:/usr/bin:/usr/local/bin
julia run.jl ${1} >arg_id_${1}.log

