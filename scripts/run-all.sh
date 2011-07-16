
DATA=$1
SPLITTED="gbk-splitted"
MAPPER="genefinder"
HADOOP=hadoop/bin/hadoop

g++ -O2 splitter.cpp -o splitter
g++ -O2 genefinder.cpp -o genefinder

cat $DATA |./splitter | head -n 3 > data-splitted.txt




$HADOOP dfs -rmr $SPLITTED
$HADOOP fs -put data-splitted.txt $SPLITTED

$HADOOP dfs -rmr gbk
$HADOOP jar hadoop/contrib/streaming/hadoop-streaming-0.20.203.0.jar \
    -D mapred.map.tasks=300 \
    -file "$MAPPER" \
    -input "$SPLITTED" \
    -output "gbk" \
    -mapper "./$MAPPER" \
    -reducer org.apache.hadoop.mapred.lib.IdentityReducer \
    
$HADOOP jar yssdc/team5.jar gbk/part-00000 tmp result.txt
