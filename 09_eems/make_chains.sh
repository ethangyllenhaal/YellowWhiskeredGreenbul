#!/bin/bash

for i in $(seq 5 20); do
    cp params_chain1 params_chain${i}
    sed -i "s/chain1/chain${i}/g" params_chain${i}
done
