package cs.ucla.edu.bwaspark.worker1

import scala.util.control.Breaks._
import scala.List
import scala.math.abs
import scala.collection.mutable.MutableList

import cs.ucla.edu.bwaspark.datatype._

// Used for read test input data
import java.io.{FileReader, BufferedReader}

object MemChainToAlign {
  val MAX_BAND_TRY = 2    

  class ReadChain(chains_i: MutableList[MemChainType], seq_i: Array[Byte]) {
    var chains: MutableList[MemChainType] = chains_i
    var seq: Array[Byte] = seq_i
  }

  var testReadChains: MutableList[ReadChain] = new MutableList
  
  /**
    *  Read the test chain data generated from bwa-0.7.8 (C version)
    *
    *  @param fileName the test data file name
    */
  def readTestData(fileName: String) {
    val reader = new BufferedReader(new FileReader(fileName))

    var line = reader.readLine
    var chains: MutableList[MemChainType] = new MutableList
    var chainPos: Long = 0
    var seeds: MutableList[MemSeedType] = new MutableList
    var seq: Array[Byte] = new Array[Byte](101)
    //var seqLen: Int = 0

    while(line != null) {
      val lineFields = line.split(" ")      

      // Find a sequence
      if(lineFields(0) == "Sequence") {
        chains = new MutableList
        //seqLen = lineFields(1).toInt
        //seq = new Array[Byte](seqLen)
        seq = lineFields(2).getBytes
        seq = seq.map(s => (s - 48).toByte) // ASCII => Byte(Int)
      }
      // Find a chain
      else if(lineFields(0) == "Chain") {
        seeds = new MutableList
        chainPos = lineFields(1).toLong
      }
      // Fina a seed
      else if(lineFields(0) == "Seed") {
        seeds += (new MemSeedType(lineFields(1).toLong, lineFields(2).toInt, lineFields(3).toInt))
      }
      // append the current list
      else if(lineFields(0) == "ChainEnd") {
        val cur_seeds = seeds
        chains += (new MemChainType(chainPos, cur_seeds))
      }
      // append the current list
      else if(lineFields(0) == "SequenceEnd") {
        val cur_chains = chains
        val cur_seq = seq 
        testReadChains += (new ReadChain(cur_chains, seq))
      }

      line = reader.readLine
    }

  }


  /**
    *  Print all the chains (and seeds) from all input reads 
    *  Only for debugging use
    */
  def printAllReads() {
    def printChains(chains: MutableList[MemChainType]) {
      println("Sequence");
      def printSeeds(seeds: MutableList[MemSeedType]) {
        seeds.foreach(s => println("Seed " + s.rBeg + " " + s.qBeg + " " + s.len))
      }
    
      chains.map(p => {
        println("Chain " + p.pos + " " + p.seeds.length)
        printSeeds(p.seeds)
                      } )
    }

    testReadChains.foreach(r => printChains(r.chains))
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
  def memChainToAln(opt: MemOptType, pacLen: Long, pac: Array[Byte], queryLen: Int, query: Array[Byte], chain: MemChainType): MutableList[MemAlnRegType] = {
    var rmax: Array[Long] = new Array[Long](2)   
    var srt: Array[SRTType] = new Array[SRTType](chain.seeds.length) 
    var alnRegs: MutableList[MemAlnRegType] = new MutableList[MemAlnRegType]
    var aw: Array[Int] = new Array[Int](2)

 
    // calculate the maximum possible span of this alignment
    rmax = getMaxSpan(opt, pacLen, queryLen, chain)
    //println("rmax(0): " + rmax(0) + ", rmax(1): " + rmax(1))  // debugging message

    // retrieve the reference sequence
    val ret = bnsGetSeq(pacLen, pac, rmax(0), rmax(1))
    var rseq = ret._1
    val rlen = ret._2
    assert(rlen == rmax(1) - rmax(0))
    // debugging message
    //println(rlen)     
    //for(i <- 0 until rlen.toInt)
      //print(rseq(i).toInt)
    //println

    // Setup the value of srt array
    for(i <- 0 to (chain.seeds.length - 1)) 
      srt(i) = new SRTType(chain.seeds(i).len, i)

    srt = srt.sortBy(s => s.len)
    //srt.map(s => println("(" + s.len + ", " + s.index + ")") )  // debugging message


    // The main for loop
    for(k <- (chain.seeds.length - 1) to 0) {
      val seed = chain.seeds( srt(k).index )
      var i = testExtension(opt, seed, alnRegs)
     
      //println("Test1: " + i)

      if(i < alnRegs.length) i = checkOverlapping(k + 1, seed, chain, srt)
      
      //println("Test2: " + i + ", Chain seed length: " + chain.seeds.length)

      // no overlapping seeds; then skip extension
      if(i == chain.seeds.length) {
        srt(k).index = 0  // mark that seed extension has not been performed
      }
      else {
        // push the current align reg into the output list
        // initialize a new alnreg
        var reg = new MemAlnRegType
        reg.width = opt.w
        aw(0) = opt.w
        aw(1) = opt.w
        reg.score = -1
        reg.trueScore = -1
     
        // left extension
        println("s.qbeg: " + seed.qBeg)
        if(seed.qBeg > 0) {
          val ret = leftExtension(opt, seed, rmax, query, rseq, reg) 
          reg = ret._1
          aw(0) = ret._2
          println("LEX qbeg: " + reg.qBeg + ", rbeg: " + reg.rBeg)
        }
        else {
          reg.score = seed.len * opt.a
          reg.trueScore = seed.len * opt.a
          reg.qBeg = 0
          reg.rBeg = seed.rBeg
          println("NLEX qbeg: " + reg.qBeg + ", rbeg: " + reg.rBeg)
        }
            
        // right extension
        println("s.qbeg: " + seed.qBeg + ", s.len: " + seed.len + ", queryLen: " + queryLen)
        if((seed.qBeg + seed.len) != queryLen) {
          val ret = rightExtension(opt, seed, rmax, query, queryLen, rseq, reg)
          reg = ret._1
          aw(1) = ret._2
          println("REX qend: " + reg.qEnd + ", rend: " + reg.rEnd)
        }
        else {
          reg.qEnd = queryLen
          reg.rEnd = seed.rBeg + seed.len
          println("NREX qend: " + reg.qEnd + ", rend: " + reg.rEnd)
        }
  
        reg.seedCov = computeSeedCoverage(chain, reg)

        if(aw(0) > aw(1)) reg.width = aw(0)
        else reg.width = aw(1)

        // push the current align reg into the output list
        alnRegs += reg
      }
    }

    alnRegs
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
      var l: Int = 0
      if( begVar >= pacLen ) {//reverse strand
        var begF: Long = (pacLen<<1) - 1 - endVar
        var endF: Long = (pacLen<<1) - 1 - begVar
        for( k <- endF to (begF + 1) by -1) {
          seq(l) = (3 - getPac(pac, k)).toByte
          l = l + 1
        }
      }
      else {
        for( k <- begVar until endVar ) {
          seq(l) = getPac(pac, k).toByte
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
    var pacValue: Long = ( pac((l>>>2).toInt) >>> (((~l)&3) <<1) ) & 3
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
  private def testExtension(opt: MemOptType, seed: MemSeedType, regs: MutableList[MemAlnRegType]): Int = {
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
  private def leftExtension(opt: MemOptType, seed: MemSeedType, rmax: Array[Long], query: Array[Byte], rseq: Array[Byte], reg: MemAlnRegType): (MemAlnRegType, Int) = {
    var aw = 0
    val tmp = (seed.rBeg - rmax(0)).toInt
    var qs = new Array[Byte](seed.qBeg)
    var rs = new Array[Byte](tmp)
    var qle = -1
    var tle = -1
    var gtle = -1
    var gscore = -1
    var maxoff = -1

    var regResult = reg
    
    for(i <- 0 to (seed.qBeg - 1)) qs(i) = query(seed.qBeg - 1 - i)
    for(i <- 0 to (tmp - 1)) rs(i) = rseq(tmp - 1 - i)
    
    breakable {
      for(i <- 0 to (MAX_BAND_TRY - 1)) {
        var prev = regResult.score
        aw = opt.w << i
        val results = SWExtend(seed.qBeg, qs, tmp, rs, 5, opt.mat, opt.oDel, opt.eDel, opt.oIns, opt.eIns, aw, opt.penClip5, opt.zdrop, seed.len * opt.a)
        regResult.score = results(0)
        qle = results(1)
        tle = results(2)
        gtle = results(3)
        gscore = results(4)
        maxoff = results(5)

        if(regResult.score == prev || ( maxoff < (aw >> 1) + (aw >> 2) ) ) break
      }
    }

    // check whether we prefer to reach the end of the query
    // local extension
    if(gscore <= 0 || gscore <= (regResult.score - opt.penClip5)) {
      regResult.qBeg = seed.qBeg - qle
      regResult.rBeg = seed.rBeg - tle
      regResult.trueScore = regResult.score
    }
    // to-end extension
    else {
      regResult.qBeg = 0
      regResult.rBeg = seed.rBeg - gtle
      regResult.trueScore = gscore
    }

    (regResult, aw)
  }

  // right extension of the current seed
  private def rightExtension(opt: MemOptType, seed: MemSeedType, rmax: Array[Long], query: Array[Byte], queryLen: Int, rseq: Array[Byte], reg: MemAlnRegType): (MemAlnRegType, Int) = {
    var aw = 0
    var regResult = reg
    var qe = seed.qBeg + seed.len
    var re = seed.rBeg + seed.len - rmax(0)
    var sc0 = regResult.score
    var qle = -1
    var tle = -1
    var gtle = -1
    var gscore = -1
    var maxoff = -1

    assert(re >= 0)

    var qeArray = new Array[Byte](queryLen - qe)
    // fill qeArray
    for(i <- 0 to (queryLen - qe - 1)) qeArray(i) = query(qe + i)

    var reArray = new Array[Byte]((rmax(1) - rmax(0) - re).toInt)
    // fill reArray
    for(i <- 0 to (rmax(1) - rmax(0) - re - 1).toInt) reArray(i) = rseq(re.toInt + i)

    breakable {
      for(i <- 0 to (MAX_BAND_TRY - 1)) {
        var prev = regResult.score
        aw = opt.w << i
        val results = SWExtend(queryLen - qe, qeArray, (rmax(1) - rmax(0) - re).toInt, reArray, 5, opt.mat, opt.oDel, opt.eDel, opt.oIns, opt.eIns, aw, opt.penClip3, opt.zdrop, sc0)
        regResult.score = results(0)
        qle = results(1)
        tle = results(2)
        gtle = results(3)
        gscore = results(4)
        maxoff = results(5)

        if(regResult.score == prev || ( maxoff < (aw >> 1) + (aw >> 2) ) ) break
      }
    }
    // check whether we prefer to reach the end of the query
    // local extension
    if(gscore <= 0 || gscore <= (regResult.score - opt.penClip3)) {
      regResult.qEnd = qe + qle
      regResult.rEnd = rmax(0) + re + tle
      regResult.trueScore += regResult.score - sc0
    }
    else {
      regResult.qEnd = queryLen
      regResult.rEnd = rmax(0) + re + gtle
      regResult.trueScore += gscore - sc0
    }

    (regResult, aw)
  }
    
  // compute seed coverage
  private def computeSeedCoverage(chain: MemChainType, reg: MemAlnRegType): Int = {
    var seedcov = 0
    
    for(i <- 0 to (chain.seeds.length - 1)) {
      // seed fully contained
      if(chain.seeds(i).qBeg >= reg.qBeg && 
         chain.seeds(i).qBeg + chain.seeds(i).len <= reg.qEnd &&
         chain.seeds(i).rBeg >= reg.rBeg &&
         chain.seeds(i).rBeg + chain.seeds(i).len <= reg.rEnd)
        seedcov += chain.seeds(i).len   // this is not very accurate, but for approx. mapQ, this is good enough
    }

    seedcov
  }
 
  class EHType(e_i: Int, h_i: Int) {
    var e: Int = e_i
    var h: Int = h_i
  }

  private def SWExtend(
    qLen: Int, query: Array[Byte], tLen: Int, target: Array[Byte], m: Int, mat: Array[Byte],
    oDel: Int, eDel: Int, oIns: Int, eIns: Int, w_i: Int, endBonus: Int, zdrop: Int, h0: Int): Array[Int] =  
  {
    var retArray: Array[Int] = new Array[Int](6)
    var eh: Array[EHType] = new Array[EHType](qLen + 1) // score array
    var qp: Array[Byte] = new Array[Byte](qLen * m) // query profile
    var oeDel = oDel + eDel
    var oeIns = oIns + eIns
    var i = 0 
    var j = 0 
    var k = 0
    var w = w_i

    for(i <- 0 until (qLen + 1))
      eh(i) = new EHType(0, 0)

    // generate the query profile
    for(k <- 0 to (m - 1)) {
      val p = k * m
      
      for(j <- 0 to (qLen - 1)) {
        qp(i) = mat(p + query(j))
        i += 1
      }
    }

    // fill the first row
    eh(0).h = h0
    if(h0 > oeIns) eh(1).h = h0 - oeIns
    else eh(1).h = 0
    j = 2
    while(j <= qLen && eh(j-1).h > eIns) {
      eh(j).h = eh(j-1).h - eIns
      j += 1
    }

    // adjust $w if it is too large
    k = m * m
    var max = 0
    for(i <- 0 to (k - 1)) // get the max score
      if(max < mat(i))
        max = mat(i)
    var maxIns = ((qLen * max + endBonus - oIns).toDouble / eIns + 1.0).toInt
    if(maxIns < 1) maxIns = 1
    var maxDel = ((qLen * max + endBonus - oDel).toDouble / eDel + 1.0).toInt
    if(maxDel < 1) maxDel = 1
    if(w > maxDel) w = maxDel  // TODO: is this necessary? (in original C implementation)

    // DP loop
    max = h0
    var max_i = -1
    var max_j = -1
    var max_ie = -1
    var gscore = -1
    var max_off = 0
    var beg = 0
    var end = qLen

    breakable {
      for(i <- 0 to (tLen - 1)) {
        var t = -1
        var f = 0
        var h1 = -1
        var m = 0
        var mj = -1
        var qPtr = target(i) * qLen
        // compute the first column
        h1 = h0 - (oDel + eDel * (i + 1))
        if(h1 < 0) h1 = 0
        // apply the band and the constraint (if provided)
        if (beg < i - w) beg = i - w
        if (end > i + w + 1) end = i + w + 1
        if (end > qLen) end = qLen

      
        for(j <- beg to (end - 1)) {
          // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
          // Similar to SSE2-SW, cells are computed in the following order:
          //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
          //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
          //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
          var h = eh(j).h
          var e = eh(j).e   // get H(i-1,j-1) and E(i-1,j)
          eh(j).h = h1
          h += qp(qPtr + j)
          if(h < e) h = e
          if(h < f) h = f 
          h1 = h            // save H(i,j) to h1 for the next column
          if(m < h) { 
            mj = j          // record the position where max score is achieved
            m = h           // m is stored at eh[mj+1]
          }
          t = h - oeDel
          if(t < 0) t = 0
          e -= eDel
          if(e < t) e = t   // computed E(i+1,j)
          eh(j).e = e       // save E(i+1,j) for the next row
          t = h - oeIns
          if(t < 0) t = 0
          f -= eIns
          if(f < t) f = t
        }
      
        eh(end).h = h1
        eh(end).e = 0
        // end == j after the previous loop
        if(end == qLen) {
          if(gscore < h1) {
            max_ie = i
            gscore = h1
          }
        }

        if(m == 0) break

        if(m > max) {
          max = m
          max_i = i
          max_j = mj

          if(max_off < abs(mj - i)) max_off = abs(mj - i)
        }
        else if(zdrop > 0) {
          if((i - max_i) > (mj - max_j)) 
            if(max - m - ((i - max_i) - (mj - max_j)) * eDel > zdrop) break
          else
            if(max - m - ((mj - max_j) - (i - max_i)) * eIns > zdrop) break;
        }
        
        // update beg and end for the next round
        j = mj
        while(j >= beg && eh(j).h > 0) {
          beg = j + 1
          j -= 1
        }
        j = mj + 2
        while(j <= end && eh(j).h > 0) {
          end = j
          j += 1
        }
      }
    }

    retArray(0) = max
    retArray(1) = max_j + 1
    retArray(2) = max_i + 1
    retArray(3) = max_ie + 1
    retArray(4) = gscore
    retArray(5) = max_off
    
    retArray
  }
}

