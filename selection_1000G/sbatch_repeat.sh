#!/bin/bash

for i in $(seq 20); do
	sbatch "${@:1}" && break || sleep 15
done

