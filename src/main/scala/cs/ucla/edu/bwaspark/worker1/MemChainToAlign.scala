package cs.ucla.edu.bwaspark.worker1

import cs.ucla.edu.bwaspark.datatype._

// Used for read test input data
import java.io.{FileReader, BufferedReader}

object MemChainToAlign {
  var testSeedChains: List[List[MemChainType]] = Nil

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
        //println("Sequence")
        chains = Nil
      }
      // Find a chain
      else if(lineFields(0) == "Chain") {
        //println("Chain")
        seeds = Nil
        chainPos = lineFields(1).toLong
      }
      // Fina a seed
      else if(lineFields(0) == "Seed") {
        //println("Seed")
        seeds = (new MemSeedType(lineFields(1).toLong, lineFields(2).toInt, lineFields(3).toInt)) :: seeds
      }
      // append the current list
      else if(lineFields(0) == "ChainEnd") {
        //println("ChainEnd")
        val cur_seeds = seeds.reverse
        chains = (new MemChainType(chainPos, cur_seeds)) :: chains
      }
      // append the current list
      else if(lineFields(0) == "SequenceEnd") {
        val cur_chains = chains.reverse
        testSeedChains = cur_chains :: testSeedChains
        //println("SequenceEnd")
      }

      line = reader.readLine
    }

    testSeedChains = testSeedChains.reverse
  }

  // Print the read results from the input file
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
}

