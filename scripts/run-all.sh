
DATA=$1
SPLITTED="gbk-splitted"
#SPLITTED="gbk-1"
MAPPER="genefinder"
HADOOP=hadoop/bin/hadoop

g++ -O2 splitter.cpp -o splitter
g++ -O2 genefinder.cpp -o genefinder

cat $DATA |./splitter > data-splitted.txt

$HADOOP dfs -rmr $SPLITTED
$HADOOP fs -put data-splitted.txt $SPLITTED

$HADOOP dfs -rmr gbk
$HADOOP jar hadoop/contrib/streaming/hadoop-streaming-0.20.203.0.jar \
    -numReduceTasks 0 \
    -file "$MAPPER" \
    -input "$SPLITTED" \
    -output "gbk" \
    -mapper "./$MAPPER" \
$HADOOP jar yssdc/team5.jar gbk/part-00000 tmp result.txt
