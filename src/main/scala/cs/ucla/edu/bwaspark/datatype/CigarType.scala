package cs.ucla.edu.bwaspark.datatype

import scala.collection.immutable.Vector

class CigarType {
  var cigarSegs: Vector[CigarSegType] = scala.collection.immutable.Vector.empty
  var cigarStr: String = new String
}

