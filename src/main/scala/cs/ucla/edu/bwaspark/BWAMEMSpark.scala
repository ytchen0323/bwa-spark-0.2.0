package cs.ucla.edu.bwaspark

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
//import org.apache.spark.storage.StorageLevel

import scala.io.{Source, BufferedSource}

import cs.ucla.edu.bwaspark.datatype._

object BWAMEMSpark {
   def main(args: Array[String]) {
      //val sc = new SparkContext("local[12]", "BWA-mem Spark",
         //"/home/hadoopmaster/spark/spark-0.9.0-incubating-bin-hadoop2-prebuilt/", List("/home/ytchen/incubator/bwa-spark-0.2.0/target/bwa-spark-0.2.0.jar"))
      //val sc = new SparkContext("spark://Jc11:7077", "BWA-mem Spark",
         //"/home/hadoopmaster/spark/spark-0.9.0-incubating-bin-hadoop2-prebuilt/", List("/home/ytchen/incubator/bwa-spark-0.2.0/target/bwa-spark-0.2.0.jar"))

      //val fastqLoader = new FASTQLocalFileLoader()
      //fastqLoader.storeFASTQInHDFS(sc, args(0), args(1))
      //val fastqRDDLoader = new FASTQRDDLoader(sc, "hdfs://Jc11:9000/user/ytchen/ERR013140_2.filt.fastq.test4/", 13)
      //val fastqRDD = fastqRDDLoader.RDDLoadAll()
      //val fastqRDD = fastqRDDLoader.RDDLoad("hdfs://Jc11:9000/user/ytchen/ERR013140_2.filt.fastq.test4/2/")

      val bwt = new BWTType()
      bwt.BWTLoad("/home/pengwei/genomics/ReferenceMetadata/human_g1k_v37.fasta.bwt")
   }
}
