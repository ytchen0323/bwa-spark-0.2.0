package cs.ucla.edu.bwaspark

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
//import org.apache.spark.storage.StorageLevel

import scala.collection.mutable.MutableList

import cs.ucla.edu.bwaspark.datatype._
import cs.ucla.edu.bwaspark.worker1.BWAMemWorker1._
import cs.ucla.edu.bwaspark.worker2.BWAMemWorker2._
import cs.ucla.edu.bwaspark.debug.DebugFlag._

import java.io.FileReader
import java.io.BufferedReader

object BWAMEMSpark {
  // load reads from the FASTQ file (for testing use)
  private def loadFASTQSeqs(fileName: String): MutableList[String] = {
    
    val reader = new BufferedReader(new FileReader(fileName))
    var line = reader.readLine
    var i = 0    
    var seqs = new MutableList[String]

    while(line != null) {
      if(i % 4 == 1) {
        seqs += line
      }
      i += 1
      line = reader.readLine
    } 

    seqs.foreach(println(_))
    seqs
  }

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


    //loading index files
    val bwaIdx = new BWAIdxType
    val prefix = "/home/pengwei/genomics/ReferenceMetadata/human_g1k_v37.fasta"
    bwaIdx.load(prefix, 0)

    //loading BWA MEM options
    val bwaMemOpt = new MemOptType
    bwaMemOpt.load

    //debugLevel = 1

    //loading reads
    var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_20reads.fq")

    val regsAllReads = seqs.map( seq => bwaMemWorker1(bwaMemOpt, bwaIdx.bwt, bwaIdx.bns, bwaIdx.pac, null, seq.length, seq) )

    // print regs for all reads
    var readNum = 0
    regsAllReads.foreach(read => {
      var i = 0
      println("Read " + readNum)
      read.foreach(r => {
        print("Reg " + i + "(")
        print(r.rBeg + ", " + r.rEnd + ", " + r.qBeg + ", " + r.qEnd + ", " + r.score + ", " + r.trueScore + ", ")
        println(r.sub + ", "  + r.csub + ", " + r.subNum + ", " + r.width + ", " + r.seedCov + ", " + r.secondary + ")")
        i += 1
      } )
      println("##############################################################")
      readNum += 1
    } )

  } 
}
