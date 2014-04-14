package cs.ucla.edu.bwaspark.datatype

import java.io.{FileInputStream, IOException}
import java.nio.{ByteBuffer}
import java.nio.file.{Files, Path, Paths}
import java.nio.channels.FileChannel
import java.nio.ByteOrder

import scala.util.control.Breaks._
import scala.List


class BWTType {
  // Data Structure
  var primary: Long = _
  var L2: Array[Long] = new Array[Long](5)
  var seqLen: Long = _
  var bwtSize: Long = _
  var bwt: Array[Int] = _
  var cntTable: Array[Int] = new Array[Int](256)
  var saIntv: Int = _
  var numSa: Long = _
  var sa: Array[Long] = _

  val readBufSize = 0x80000

  def load(prefix: String) {
    // Load the .bwt file
    BWTLoad(prefix + ".bwt")
    // Load the .sa file
    SALoad(prefix + ".sa")
    // Generate CNT table
    genCNTTable
  }

  private def SALoad(saFile: String) {
    val fs = new FileInputStream(saFile)
    val fc = fs.getChannel
    val buf = ByteBuffer.allocate(readBufSize)
    buf.order(ByteOrder.nativeOrder)
    
    val smallBuf = ByteBuffer.allocate(8)
    smallBuf.order(ByteOrder.nativeOrder)

    // Read primary in .sa file
    fc.read(smallBuf)
    smallBuf.flip
    val primarySa = smallBuf.getLong
    smallBuf.rewind
    assert(primarySa == primary, "SA-BWT inconsistency: primary is not the same.")

    // Skipped 4 * sizeof(Long)
    val tmpBuf = ByteBuffer.allocate(32)
    tmpBuf.order(ByteOrder.nativeOrder)
    fc.read(tmpBuf)

    // Read saIntv member variable in BWTType
    fc.read(smallBuf)
    smallBuf.flip
    saIntv = smallBuf.getLong.asInstanceOf[Int]  // Long -> Int (done in the original c code)
    smallBuf.rewind

    // Read seqLen in .sa file
    fc.read(smallBuf)
    smallBuf.flip
    val seqLenSa = smallBuf.getLong
    smallBuf.rewind
    assert(seqLenSa == seqLen, "SA-BWT inconsistency: seq_len is not the same.")

    numSa = (seqLen + saIntv) / saIntv;
    sa = new Array[Long](numSa.asInstanceOf[Int])  // Long -> Int
    
    // Read the sa array member variable in BWTType
    sa(0) = -1
    var ret = 0
    var startIdx = 1

    while(ret >= 0) {
      ret = fc.read(buf)
      buf.flip

      // Fill the data from buf
      while(buf.hasRemaining) {
        val saPiece = buf.getLong
        sa(startIdx) = saPiece
        startIdx += 1
      }

      buf.rewind
    }
    
    fc.close

    // Debugging message
    println("saIntv: " + saIntv)
    println("numSa: " + numSa)
    for(i <- 0 to 19)
      println("sa(" + i + "): " + sa(i))
  }

  private def BWTLoad(bwtFile: String) {
    val fs = new FileInputStream(bwtFile)
    val fc = fs.getChannel
    val fileSize = fc.size
    val buf = ByteBuffer.allocate(readBufSize)
    buf.order(ByteOrder.nativeOrder)

    val smallBuf = ByteBuffer.allocate(8)
    smallBuf.order(ByteOrder.nativeOrder)

    // Read primary member variable in BWTType
    fc.read(smallBuf)
    smallBuf.flip
    primary = smallBuf.getLong
    smallBuf.rewind

    // Read L[1] - L[4] member variables in BWTType
    L2(0) = 0
    for(i <- 1 to 4) {
      fc.read(smallBuf)
      smallBuf.flip
      L2(i) = smallBuf.getLong
      smallBuf.rewind
    }
    seqLen = L2(4)

    // Read the bwt array in BWTType
    var ret = 0
    var bytesRead: Long = 0
    var numIntRead = 0
    bwtSize = (fileSize - 8 * 5) >> 2

    println("File size: " + fileSize)

    // Potential array boundary problem here
    // Enforce the type of numElem to be Int
    // This is because the size of a scala array cannot exceed the Int upper bound
    bwt = new Array[Int](bwtSize.asInstanceOf[Int])

    while(ret >= 0) {
      ret = fc.read(buf)
      buf.flip

      // Fill the data from buf
      while(buf.hasRemaining) {
        val bwtSeg = buf.getInt
        bwt(numIntRead) = bwtSeg
        numIntRead += 1
      }

      buf.rewind
    } 

    fc.close

    // Debugging messages
    println("Num Integers Read: " + numIntRead)
    println("primary: " + primary)
    println("bwtSize: " + bwtSize)
    println("seqLen: " + seqLen)
    for(i <- 0 to 4) 
      println("L2(" + i + "): " + L2(i))
    for(i <- 0 to 19)
      println("bwt(" + i + "): " + bwt(i))
  }

  implicit def bool2int(b:Boolean) = if (b) 1 else 0

  private def genCNTTable() {
    for(i <- 0 to 255) {
      var x = 0
      for(j <- 0 to 3) {
        x = x | ((((i&3) == j) + ((i>>>2&3) == j) + ((i>>>4&3) == j) + (i>>>6 == j)) << (j<<3))
      }   
      cntTable(i) = x
    }

    for(i <- 0 to 255) 
      println("cntTable(" + i + "): " + cntTable(i))
  }
}

