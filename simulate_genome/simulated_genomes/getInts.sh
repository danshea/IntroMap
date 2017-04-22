#!/bin/bash

grep 'Introgression ' $1 | cut -d\  -f5,10 | sort -k1
