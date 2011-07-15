FILE=/data/gbk/ref_chr7_00.gbk 

g++ splitter.cpp -o splitter.out
hadoop jar /usr/lib/hadoop-0.20/contrib/streaming/hadoop-streaming-0.20.2-cdh3u0.jar \
    -files "splitter.out" \
    -mapper "./splitter.out" \
    -numReduceTasks 0 \
    -input $FILE \
    -output gbk-splitted \
