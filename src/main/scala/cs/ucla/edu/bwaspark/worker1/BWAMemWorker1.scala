package cs.ucla.edu.bwaspark.worker1

import cs.ucla.edu.bwaspark.datatype._
import scala.collection.mutable.MutableList
import java.util.TreeSet
import java.util.Comparator
import cs.ucla.edu.bwaspark.worker1.MemChain._
import cs.ucla.edu.bwaspark.worker1.MemChainFilter._
import cs.ucla.edu.bwaspark.worker1.MemChainToAlign._
import cs.ucla.edu.bwaspark.worker1.MemSortAndDedup._

//this standalone object defines the main job of BWA MEM:
//1)for each read, generate all the possible seed chains
//2)using SW algorithm to extend each chain to all possible aligns
object BWAMemWorker1 {
  
  //the function which do the main task
  def bwaMemWorker1(opt: MemOptType, //BWA MEM options
                    bwt: BWTType, //BWT and Suffix Array
                    bns: BNTSeqType, //.ann, .amb files
                    pac: Array[Byte], //.pac file uint8_t
                    pes: Array[MemPeStat], //pes array
                    len: Int, //the length of the read
                    seq: Array[Int] //a read
                    ): MutableList[MemAlnRegType] = { //all possible alignment  

    //for paired alignment, to add
    //!!!to add!!!
    //for now, we only focus on single sequence alignment
    //if (!(opt.flag & MEM_F_PE)) {
    if (true) {

      //pre-process: transform A/C/G/T to 0,1,2,3

      def locusEncode(locus: Int): Byte = {
        //transforming from A/C/G/T to 0,1,2,3
        locus match {
          case 'A' => 0
          case 'a' => 0
          case 'C' => 1
          case 'c' => 1
          case 'G' => 2
          case 'g' => 2
          case 'T' => 3
          case 't' => 3
          case '-' => 5
          case _ => 4
        }
      }

      val read: Array[Byte] = seq.map(ele => locusEncode(ele))

      //first step: generate all possible MEM chains for this read
      val chains = generateChains(opt, bwt, bns.l_pac, len, read) 

      //second step: filter chains
      val chainsFiltered = memChainFilter(opt, chains)

      //third step: for each chain, from chain to aligns
      val regs = new MutableList[MemAlnRegType]
      val alignRegArray = chainsFiltered.map(ele => memChainToAln(opt, bns.l_pac, pac, len, read, ele, regs)).reduce( (a,b) => b)

      //last step: sorting and deduplication

      val pureRegArray = memSortAndDedup(alignRegArray, opt.maskLevelRedun)

      pureRegArray
    }
    else {
      assert (false)
      null
    }
  }
}
