export SPARK_CLASSPATH=target/bwa-spark-0.2.0.jar
export SPARK_JAVA_OPTS="-Dspark.executor.memory=32g"
java -Xmx32g $SPARK_JAVA_OPTS -jar target/bwa-spark-0.2.0.jar /home/ytchen/genomics/data/ERR013140_2.filt.fastq hdfs://Jc11:9000/user/ytchen/ERR013140_2.filt.fastq.test3
