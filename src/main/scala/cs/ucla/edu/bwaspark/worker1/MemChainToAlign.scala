package cs.ucla.edu.bwaspark.worker1

import scala.util.control.Breaks._
import scala.List

import cs.ucla.edu.bwaspark.datatype._

// Used for read test input data
import java.io.{FileReader, BufferedReader}

object MemChainToAlign {
  var testSeedChains: List[List[MemChainType]] = Nil
  
  /**
    *  Read the test chain data generated from bwa-0.7.8 (C version)
    *
    *  @param fileName the test data file name
    */
  def readTestData(fileName: String) {
    val reader = new BufferedReader(new FileReader(fileName))

    var line = reader.readLine
    var chains: List[MemChainType] = Nil
    var chainPos: Long = 0
    var seeds: List[MemSeedType] = Nil

    while(line != null) {
      val lineFields = line.split(" ")      

      // Find a sequence
      if(lineFields(0) == "Sequence") {
        chains = Nil
      }
      // Find a chain
      else if(lineFields(0) == "Chain") {
        seeds = Nil
        chainPos = lineFields(1).toLong
      }
      // Fina a seed
      else if(lineFields(0) == "Seed") {
        seeds = (new MemSeedType(lineFields(1).toLong, lineFields(2).toInt, lineFields(3).toInt)) :: seeds
      }
      // append the current list
      else if(lineFields(0) == "ChainEnd") {
        val cur_seeds = seeds.reverse
        chains = (new MemChainType(chainPos, cur_seeds)) :: chains
      }
      // append the current list
      else if(lineFields(0) == "SequenceEnd") {
        val cur_chains = chains.reverse
        testSeedChains = cur_chains :: testSeedChains
      }

      line = reader.readLine
    }

    testSeedChains = testSeedChains.reverse
  }


  /**
    *  Print all the chains (and seeds) from all input reads 
    *  Only for debugging use
    */
  def printAllReads() {
    def printChains(chains: List[MemChainType]) {
      println("Sequence");
      def printSeeds(seeds: List[MemSeedType]) {
        seeds.foreach(s => println("Seed " + s.rBeg + " " + s.qBeg + " " + s.len))
      }
    
      chains.map(p => {
        println("Chain " + p.pos + " " + p.seeds.length)
        printSeeds(p.seeds)
                      } )
    }

    testSeedChains.foreach(printChains(_))
  }


  // calculate the maximum possible span of this alignment
  private def calMaxGap(opt: MemOptType, qLen: Int): Int = {
    val lenDel = ((qLen * opt.a - opt.oDel).toDouble / opt.eDel.toDouble + 1.0).toInt
    val lenIns = ((qLen * opt.a - opt.oIns).toDouble / opt.eIns.toDouble + 1.0).toInt
    var len = -1

    if(lenDel > lenIns)
      len = lenDel
    else
      len = lenIns

    if(len <= 1) len = 1

    val tmp = opt.w << 1

    if(len < tmp) len
    else tmp
  }

  class SRTType(len_i: Int, index_i: Int) {
    var len: Int = len_i
    var index: Int = index_i
  }

  /**
    *  The main function of memChainToAlign class
    *
    *  @param opt the MemOptType object
    *  @param pacLen the pac length
    *  @param pac the pac array
    *  @param queryLen the query length (read length)
    *  @param chain one of the mem chains of the read
    */
  //def memChainToAln(opt: MemOptType, pacLen: Long, pac: Array[Byte], queryLen: Int, chain: MemChainType): RDD[MemAlnRegType] = {
  def memChainToAln(opt: MemOptType, pacLen: Long, pac: Array[Byte], queryLen: Int, chain: MemChainType) {
    var rmax: Array[Long] = new Array[Long](2)   
    var srt: Array[SRTType] = new Array[SRTType](chain.seeds.length) 
    var alnRegs: List[MemAlnRegType] = Nil

    // calculate the maximum possible span of this alignment
    rmax = getMaxSpan(opt, pacLen, queryLen, chain)
    println("rmax(0): " + rmax(0) + ", rmax(1): " + rmax(1))
    // retrieve the reference sequence
    //(rlen, rseq) = bnsGetSeq(pacLen, pac, rmax(0), rmax(1))
    //assert(rlen == rmax(1) - rmax(0))
   
    // Setup the value of srt array
    for(i <- 0 to (chain.seeds.length - 1)) {
      srt(i).len = chain.seeds(i).len 
      srt(i).index = i
    }

    srt = srt.sortBy(s => s.len)
    //srt.map(s => println("(" + s._1 + ", " + s._2 + ")") )  // debugging use

    // The main for loop
    for(k <- (chain.seeds.length - 1) to 0) {
      val seed = chain.seeds(k)
      var i = testExtension(opt, seed, alnRegs)

      if(i < alnRegs.length) i = checkOverlapping(k + 1, seed, chain, srt)
      
      // no overlapping seeds; then skip extension
      if(i == chain.seeds.length) {
        srt(k).index = 0  // mark that seed extension has not been performed
      }
      else {
        // push the current align reg into the output list
        // (need to write some codes here)
     
        if(seed.qBeg > 0) 
          leftExtension
        // else 
     
        if((seed.qBeg + seed.len) != queryLen) 
          rightExtension
  
        computeSeedCoverage
      }
    }
  }

  // retrieve the reference sequence
  // scala: l_pac, pac, rmax[0], rmax[1], return &rlen, return rseq
  // c: l_pac, pac, beg, end, len, return rseq
  private def bnsGetSeq(pacLen: Long, pac: Array[Byte], beg: Long, end: Long) : (Array[Byte], Long) = {
    var endVar: Long = 0//for swapping
    var begVar: Long = 0//for swapping 
    if(end < beg) {//if end is smaller, swap
      endVar = beg
      begVar = end
  }
  else {//else keep the value
    endVar = end
    begVar = beg
    }
    if(endVar > (pacLen<<1)) endVar = pacLen<<1
    if(begVar < 0) begVar = 0
    var rLen: Long = endVar - begVar// for return rlen
    var seq: Array[Byte] = new Array[Byte](rLen.toInt)//for return seq

    if(begVar >= pacLen || endVar <= pacLen) {
      var k: Long = 0
      var l: Long = 0
      if( begVar >= pacLen ) {//reverse strand
	var begF: Long = pacLen<<1 - 1 - endVar
	var endF: Long = pacLen<<1 - 1 - begVar
	for( k <- endF until begF) {
	  seq(l.toInt) = ( 3 - getPac(pac, k) ).toByte
	  l = l + 1
 	 }
	}
      else {
	for( k <- begVar until endVar ) {
	  seq(l.toInt) = ( getPac(pac, k) ).toByte
	  l = l + 1
	}
      }
	}
	else
	  rLen = 0

	(seq, rLen)//return two value
    }
  //#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
  private def getPac( pac: Array[Byte], l: Long) : Long = {
	var pacValue: Long = ( pac((l>>2).toInt) >> ((~(l)&3) <<1) ) & 3
	pacValue
  }
 	


  // get the max possible span
  private def getMaxSpan(opt: MemOptType, pacLen: Long, queryLen: Int, chain: MemChainType): Array[Long] = {
    var rmax: Array[Long] = new Array[Long](2)
    val doublePacLen = pacLen << 1
    rmax(0) = doublePacLen
    rmax(1) = 0

    val seedMinRBeg = chain.seeds.map(seed => 
      { seed.rBeg - ( seed.qBeg + calMaxGap(opt, seed.qBeg) ) } ).min
    val seedMaxREnd = chain.seeds.map(seed => 
      { seed.rBeg + seed.len + (queryLen - seed.qBeg - seed.len) + calMaxGap(opt, queryLen - seed.qBeg - seed.len) } ).max
   
    if(rmax(0) > seedMinRBeg) rmax(0) = seedMinRBeg
    if(rmax(1) < seedMaxREnd) rmax(1) = seedMaxREnd
      
    if(rmax(0) <= 0) rmax(0) = 0
    if(rmax(1) >= doublePacLen) rmax(1) = doublePacLen

    // crossing the forward-reverse boundary; then choose one side
    if(rmax(0) < pacLen && pacLen < rmax(1)) {
      // this works because all seeds are guaranteed to be on the same strand
      if(chain.seeds(0).rBeg < pacLen) rmax(1) = pacLen
      else rmax(0) = pacLen
    }

    rmax
  }
   
  // test whether extension has been made before
  // NOTE: need to be tested!!!
  private def testExtension(opt: MemOptType, seed: MemSeedType, regs: List[MemAlnRegType]): Int = {
    var rDist: Long = -1 
    var qDist: Int = -1
    var maxGap: Int = -1
    var minDist: Int = -1
    var w: Int = -1
    var breakIdx: Int = regs.length

    breakable {
      for(i <- 0 to (regs.length - 1)) {
        
        if(seed.rBeg >= regs(i).rBeg && (seed.rBeg + seed.len) <= regs(i).rEnd && seed.qBeg >= regs(i).qBeg && (seed.qBeg + seed.len) <= regs(i).qEnd) {
          // qDist: distance ahead of the seed on query; rDist: on reference
          qDist = seed.qBeg - regs(i).qBeg
          rDist = seed.rBeg - regs(i).rBeg

          if(qDist < rDist) minDist = qDist 
          else minDist = rDist.toInt

          // the maximal gap allowed in regions ahead of the seed
          maxGap = calMaxGap(opt, minDist)

          // bounded by the band width          
          if(maxGap < opt.w) w = maxGap
          else w = opt.w
          
          // the seed is "around" a previous hit
          if((qDist - rDist) < w && (rDist - qDist) < w) { 
            breakIdx = i 
            break
          }

          // the codes below are similar to the previous four lines, but this time we look at the region behind
          qDist = regs(i).qEnd - (seed.qBeg + seed.len)
          rDist = regs(i).rEnd - (seed.rBeg + seed.len)
          
          if(qDist < rDist) minDist = qDist
          else minDist = rDist.toInt

          maxGap = calMaxGap(opt, minDist)

          if(maxGap < opt.w) w = maxGap
          else w = opt.w

          if((qDist - rDist) < w && (rDist - qDist) < w) {
            breakIdx = i
            break
          }          
        }
      }
    }

    breakIdx
  }
    
  // check overlapping seeds in the same chain
  private def checkOverlapping(startIdx: Int, seed: MemSeedType, chain: MemChainType, srt: Array[SRTType]): Int = {
    var breakIdx = chain.seeds.length

    breakable {
      for(i <- startIdx to (chain.seeds.length - 1)) {
        if(srt(i).index != 0) {
          val targetSeed = chain.seeds(srt(i).index)

          // only check overlapping if t is long enough; TODO: more efficient by early stopping
          // NOTE: the original implementation may be not correct!!!
          if(targetSeed.len >= seed.len * 0.95) {
            if(seed.qBeg <= targetSeed.qBeg && (seed.qBeg + seed.len - targetSeed.qBeg) >= (seed.len>>2) && (targetSeed.qBeg - seed.qBeg) != (targetSeed.rBeg - seed.rBeg)) {
              breakIdx = i
              break
            }
            if(targetSeed.qBeg <= seed.qBeg && (targetSeed.qBeg + targetSeed.len - seed.qBeg) >= (seed.len>>2) && (seed.qBeg - targetSeed.qBeg) != (seed.rBeg - targetSeed.rBeg)) {
              breakIdx = i
              break
            }
          }
        }
      }
    }

    breakIdx
  }

  // left extension of the current seed
  private def leftExtension() {

  }

  // right extension of the current seed
  private def rightExtension() {

  }
    
  // compute seed coverage
  private def computeSeedCoverage() {

  }
 
  private def SWExtend(
    qLen: Int, query: Array[Byte], tLen: Int, target: Array[Byte], m: Int, mat: Array[Byte],
    oDel: Int, eDel: Int, oIns: Int, eInt: Int, w: Int, endBonus: Int, zdrop: Int, h0: Int): Array[Int] =  
  {
    var retArray: Array[Int] = new Array[Int](6)

    retArray
  }
}

