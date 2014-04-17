package cs.ucla.edu.bwaspark.worker1

import cs.ucla.edu.bwaspark.datatype._
import scala.math._
import scala.collection.mutable.MutableList
import java.util.TreeSet
import java.util.Comparator

//standalone object for generating all MEM chains for each read
object MemChain {

  //generate a chain tree for each read

  def generateChainTree(opt: MemOptType, l_pac: Long, smemItr: SMemItrType): TreeSet[MemChainType] = {

    //calculate splitLen
    val splitLen = min((opt.min_seed_len * opt.split_factor + 0.499).toInt, smemItr.len)

    //!!!be careful!!!
    //we temporarily assign startWidth as 1 other than 2
    //but it need be fixed
    //val startWidth = { if (opt.flag & MEM_F_NO_EXACT) 2 else 1 }
    val startWidth = 1

    //the mutable list is for storing all the bi-intervals generated from specific point
    var bwtIntvOnPoint: MutableList[BWTIntvType] = null

    //the return value: mutable TreeSet which maintains all the chains
    var chainTree: TreeSet[MemChainType] = null

    //if bi-intervals responding to specific point are generated
    //go through each seed, either
    //1) merge it to existing chain
    //2) generate new chain from it
    while ( (bwtIntvOnPoint = smemNext(smemItr, splitLen, opt.split_width, startWidth)) != null ) {
      //traverse all the seeds
      for (i <- 0 until bwtIntvOnPoint.length) {

        //end - start = length of the seed
        val seedLen = bwtIntvOnPoint(i).endPoint - bwtIntvOnPoint(i).startPoint

        //ignore the seed if it it too short or too repetitive
        if (seedLen < opt.min_seed_len || bwtIntvOnPoint(i).s > opt.max_occ) {
          //do nothing
        }
        //the statements in the else clause are the main part
        //it will traverse all the possible positions in the suffix array(reference)
        else {
          //traverse each aligned position
          for (j <- 0 until bwtIntvOnPoint(i).s) {

            //prepare for generating a new seed
            var rBeg = suffixArrayPos2ReferencePos()
            var qBeg = bwtIntvOnPoint(i).startPoint
            var len = seedLen

            //handle edge cases
            if (rBeg < l_pac && l_pac < rBeg + len) {
              //do nothing if a seed crossing the boundary
            }
            else {
              //generate a seed for this position
              val newSeed = new MemSeedType(rBeg, qBeg, len)

              //find the closest chain in the existing chain tree
              //the closest chain should satisfy
              //1)chain.pos <= seed.rbegin
              //2)the chain.pos is the largest one of all chains that satisfy 1)

              def findClosestChain(chainTree: TreeSet[MemChainType], refPoint: Long): MemChainType = {
                //if the tree is empty, return null
                if (chainTree == null) null
                else {
                  //create a temporary chain for finding the lower chain
                  //because lower.pos < tmpChain.pos
                  //having refPoint + 1 to handle if lower.pos == tmpChain.pos
                  val tmpChain = new MemChainType(refPoint + 1, null)
                  chainTree.lower(tmpChain)
                }
              }

              val targetChain = findClosestChain(chainTree, newSeed.rBeg)

              //test if the seed can be merged into some existing chain
              //return true/false
              //if return true, actually also DID the merging task

              //define tryMergeSeedToChain to test if a seed can be merged with some chain in the chain tree
              def tryMergeSeedToChain(opt: MemOptType, l_pac: Long, chain: MemChainType, seed: MemSeedType): = {

                //get query begin and end, reference begin and end
                //!!!to clarify!!!: the order of seeds in a chain
                //qBeg sorting? or rBeg sorting?
                val qBegChain = chain.seeds.head.qBeg
                val rBegChain = chain.seeds.head.rBeg
                val qEndChain = chain.seeds.last.qBeg + chain.seeds.last.len
                val rEndChain = chain.seeds.last.rBeg + chain.seeds.last.len

                //if the seed is fully contained by the chain, return true
                if (qBegChain <= seed.qBeg && qEndChain >= seed.qBeg + seed.len &&
                    rBegChain <= seed.rBeg && rEndChain >= seed.rBeg + seed.len)
                  true
                //if not in the same strand (crossing l_pac boundary), return false
                else if ( (rBegChain < l_pac || chain.seeds.last.rBeg < l_pac) &&
                          seed.rBeg >= l_pac)
                  false
                else {
                  //follow the conditions judged in original BWA test_and_merge function
                  val x = seed.qBeg - chain.seeds.last.qBeg // always non-negtive???
                  val y = seed.rBeg - chain.seeds.last.rBeg

                  if (y >= 0 &&
                      x - y <= opt.w &&
                      y - x <= opt.w &&
                      x - chain.seeds.last.len < opt.max_chain_gap &&
                      y - chain.seeds.last.len < opt.max_chain_gap) {
                    //all the conditions are satisfied? growing the chain
                    chain += seed
                    true //return true
                  }
                  false
                }
              }

              val isMergable = tryMergeSeedToChain(opt, l_pac, targetChain, newSeed)

              //add the seed as a new chain if not mergable
              if (!isMergable) {
                val newSeedList = new MutableList[MemSeedType](newSeed)
                val newChain = new MemChainType(rBeg, newSeedList)
                //Push the new chain to the chain tree
                //1)if the chainTree is empty
                if (chainTree == null) {
                  //using java style to new a TreeSet[MemChainType]
                  chainTree = new TreeSet[MemChainType](new Comparator[MemChainType]() {
                    def compare(a: MemChainType, b: MemChainType): Int = {
                      a.pos - b.pos
                    }
                  } )
                  //insert the chain to the tree
                  chainTree.add(newChain)
                }
                //2)if the chainTree is not empty, directly add it
                else chainTree.add(newChain)
              }
            }            
          }
        }
      }
    }

    //finally, return the tree
    chainTree

  }

  def traverseChainTree(chainTree: TreeSet[MemChainType]): Array[MemChainType] = {
    //if the tree is empty, return null
    if (chainTree == null) null
    
    //else, it's gonna be a simple map() task
    else {
      val itr = chainTree.iterator
      val chains = new Array[MemChainType](chainTree.size).map(ele => itr.next)
      chains
    }
  }

  //generate chains for each read
  def generateChains(opt: MemOptType, bwt: BWTType, l_pac: Long, len: Int, seq: Array[Int]): Array[MemChainType] = {
    
    //if the query is shorter than the seed length, no match, return null
    if (len < opt.min_seed_len) {
      null
    }
    //the else part the real meaty part for this function
    else {
      //generate a SMemItrType object for smemNext
      val smemItr = new SMemItrType(bwt, 
                                    query, 
                                    0, //the first startpoint for smemNext is 0
                                    len,
                                    null, //matches array
                                    null, //sub array
                                    null, //temporary array 0
                                    null) //temporary array 1

      //generate a tree for all chains
      val chainTree = generateChainTree(opt, l_pac, smemItr)

      //return value, the chains to be generated for a read
      val chains = traverseChainTree(chainTree)
      chains
    }
  }  
}
