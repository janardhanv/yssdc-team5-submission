awk '(NR%2==0)' | \
    python -c "import sys; print "11111";print ''.join(x.strip('\n') for x in sys.stdin.xreadlines())" 
