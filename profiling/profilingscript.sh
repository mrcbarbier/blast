echo $$
valgrind --tool=callgrind   ./simu 

python ./gprof2dot.py -f callgrind callgrind.out.$(( $$+1)) | dot -Tsvg -o output.svg
