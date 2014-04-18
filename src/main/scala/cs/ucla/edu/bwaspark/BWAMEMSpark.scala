package cs.ucla.edu.bwaspark

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
//import org.apache.spark.storage.StorageLevel

import scala.io.{Source, BufferedSource}
import scala.collection.mutable.MutableList

import cs.ucla.edu.bwaspark.datatype._
import cs.ucla.edu.bwaspark.worker1.MemChainToAlign._

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

      val idx = new BWAIdxType()
      idx.load("/home/pengwei/genomics/ReferenceMetadata/human_g1k_v37.fasta", 1)

      //println("Hello!")
      val opt = new MemOptType()
      opt.load()
      //for ( i <- 0 to 24){
      //  println(opt.mat(i))
      //}

      readTestData("/home/ytchen/bwa/bwa-0.7.8/log_1read")
      //printAllReads

      // debugging message
      testReadChains.foreach( read => {
        print("Sequence ")
        for(i <- 0 to 100)
          print(read.seq(i))
        println
      } )


      //testReadChains.foreach( read => read.chains.map( chain => memChainToAln(opt, idx.bns.l_pac, idx.pac, 101, read.seq, chain) ) )

      
      testReadChains.foreach( read => {
        var regs = new MutableList[MemAlnRegType]

        read.chains.map( chain => {
          regs ++= memChainToAln(opt, idx.bns.l_pac, idx.pac, 101, read.seq, chain) 
          } ) 

        // print all regs for a single read
        var i = 0
        regs.foreach(r => {
          print("Reg " + i + "(")
          print(r.rBeg + ", " + r.rEnd + ", " + r.qBeg + ", " + r.qEnd + ", " + r.score + ", " + r.trueScore + ", ")
          println(r.sub + ", "  + r.csub + ", " + r.subNum + ", " + r.width + ", " + r.seedCov + ", " + r.secondary + ")")
          i += 1
          } )

        } )


            //bwt.load("/home/pengwei/genomics/ReferenceMetadata/human_g1k_v37.fasta")
   }
}
