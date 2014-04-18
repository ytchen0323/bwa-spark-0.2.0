package cs.ucla.edu.bwaspark.worker2

import scala.util.control.Breaks._
import scala.List
import scala.math.abs
import scala.collection.mutable.MutableList

import cs.ucla.edu.bwaspark.datatype._

//MemMarkPrimarySe
object MemMarkPrimarySe {
  
  
  /**
   * The main function of memMarkPrimarySe class

   * @param opt the MemOptType object
   * @param n the Int object
   * @param a the MemAlnRegType object
   * @param id the Long object
   */
  //??opt a list or a object??
  def memMarkPrimarySe(opt: MemOptType, n: Int, a: List[MemAlnRegType], id: Long) {
    var i: Int = 0
    var j: Int = 0
    var k: Int = 0
    var z: MutableList[Int] = new MutableList

    if( n != 0) {
      for( i <- 0 until n ) {
	a(i).sub = 0
	a(i).secondary = -1
	a(i).hash = hash64((id + i.toLong))
      }

      //ks_introsort(mem_ars_hash, n, a)
      //#define alnreg_hlt(a, b) ((a).score > (b).score || ((a).score == (b).score && (a).hash < (b).hash))
      a.sortWith( (x, y) => ((x.score > y.score) || ( x.score == y.score && x.hash < y.hash)  ) )

	

    }


  }

  private def hash64( key: Long ) : Long = {
    var keyVar: Long = key
    keyVar += ~(keyVar << 32)
    keyVar ^= (keyVar >>> 22)
    keyVar += ~(keyVar << 13)
    keyVar ^= (keyVar >>> 8)
    keyVar += (keyVar << 3)
    keyVar ^= (keyVar >>> 15)
    keyVar += ~(keyVar <<27)
    keyVar ^= (keyVar >>> 31)
    keyVar
  }
}
