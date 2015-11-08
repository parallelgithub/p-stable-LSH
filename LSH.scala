import scala.io.Source
import scala.util.Random
import scala.math


object LSH {

	val prime = (1L << 29) - 33

	val binWidthW = 4.0
	val radiusR = 1.0
	val successProbability = 0.9

	//from http://picomath.org/scala/Erf.scala.html
	//   & http://www.johndcook.com/blog/cpp_erf/
	//It is correct, can be better...
	def ERF(x: Double): Double =  {
		val a1: Double =  0.254829592;
		val a2: Double = -0.284496736;
		val a3: Double =  1.421413741;
		val a4: Double = -1.453152027;
		val a5: Double =  1.061405429;
		val p: Double  =  0.3275911;
		// Save the sign of x
		val sign = if (x < 0) -1 else 1
		val absx =  math.abs(x)

		// A&S formula 7.1.26, rational approximation of error function
		val t = 1.0/(1.0 + p*absx);
		val y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x);
		sign*y
	}

	//from E2LSH code in SelfTuring.cpp
	//w: bin width
	//c: | vector_x - vector_y | , i.e. distance of 2 point
	//it is OK
	def computeFunctionP(w: Double, c: Double):Double = {
		  val x = w / c;
			val M_SQRT2 = math.sqrt(2.0)
			val M_2_SQRTPI = 2.0 / math.sqrt(math.Pi)
			ERF(x / M_SQRT2) - M_2_SQRTPI / M_SQRT2 / x * (1 - math.exp(-(x*x) / 2));
	}

	//from E2LSH & manual
	def computeLfromKP(k: Int, successProb: Double):Int = {
			val p1 = computeFunctionP(binWidthW, radiusR)
		  math.ceil(math.log(1 - successProb) / math.log(1 - math.pow(p1, k))).toInt;
	}

	//from the original paper
	//   & http://cseweb.ucsd.edu/~dasgupta/254-embeddings/lawrence.pdf
	def computeK(p2: Double,n: Int): Int = {
		(math.log(n)/math.log(1.0/p2)).toInt
	}

	class HashFunctions(numberOfFunctionG: Int,
	                    dimensionOfVectorG: Int,
	                    dimensionOfVectorA: Int, 
	                    tableSize: Int, 
	                    inWidthW: Double ){

		private case class lshFunctionParameter(vectorA: Vector[Double], parameterB: Double)
		

		private val functionG = List.fill(numberOfFunctionG, dimensionOfVectorG)(
				new lshFunctionParameter(
					Vector.fill(dimensionOfVectorA)(Random.nextGaussian()),
					Random.nextDouble() * binWidthW)
		)

		def dotHash(vectorV: Vector[Double]) = {
			functionG.map( 
				_.map( x => {
					val dot = (vectorV zip x.vectorA).map{case (a,b) => a * b}.foldLeft(0.0)(_+_)
					math.floor((dot + x.parameterB)/binWidthW).toInt
				}).toVector
			)
			
		}

		private val functionH1 = IndexedSeq.fill(dimensionOfVectorG)((Random.nextDouble()*prime).toLong)
		private val functionH2 = IndexedSeq.fill(dimensionOfVectorG)((Random.nextDouble()*prime).toLong)

		private def primaryHash(v: IndexedSeq[Int], tableSize: Int) = {
			val dot = (v zip functionH1).map{case (a,b) => a * b}.foldLeft(0L)(_+_)
			val x = dot % prime
			if (x < 0)
				((x + prime) % tableSize).toInt
			else 
				(x % tableSize).toInt
			//val hashValue: Int = ((dot % prime)%tableSize).toInt
			//hashValue
		}

		//dot = the value of inner product
		//the return value is dot%prime (or dot%prime + prime when negative)
		//  where prime = (1L << 29) - 33
		private def secondaryHash(v: IndexedSeq[Int]) = {
			val zipVector = (v zip functionH2).map{case (a,b) => a * b}

			val bitFilter =  (1L << 29) - 1
			val dot = zipVector.foldLeft(0L) ( (x,y) => {
				val z = x + y
				(z & bitFilter) +	33 * (z >> 29)
				}	
			)
			val hashValue: Int = if (dot>=prime) (dot-prime).toInt else dot.toInt
			hashValue
		}

		//vectorG means g_i(v) for given v and for all i belonging to L
		def tableHash(reducedVectors: List[Vector[Int]], 
		                 nameOfVec: Int,
										 table: Vector[Vector[List[Any]]]): Vector[Vector[List[Any]]] = {

			//For a given original vector v, assume g_i is the i-th lsh function
			//v is reduced to reducedVectors(0) by g_0
			//v is reduced to reducedVectors(1) by g_1 ...

			/*
			reducedVectors.zipWithIndex.foreach { case(reducedVec, indexOflshFunction) => {
				val hashValue1 = primaryHash(reducedVec, tableSize)
				val hashValue2 = secondaryHash(reducedVec)
				
				val oldList = table(indexOflshFunction)(hashValue1) 
				val newElement = (hashValue2, indexOfVec)
				 table(indexOflshFunction).updated(hashValue1, newElement )
			}}
			*/

			def computeHashValueToTable(reducedVecs: List[(Vector[Int], Int)],
								   table: Vector[Vector[List[Any]]]): Vector[Vector[List[Any]]]= {
				reducedVecs match {
					case Nil => table
					case (head, indexOflshFunction) :: tail => {
						val locationOfVec= primaryHash(head, tableSize)
						val fingerprintOfVec= secondaryHash(head)
						
						val oldList = table(indexOflshFunction)(locationOfVec) 
						val newElement = (fingerprintOfVec, nameOfVec)
						val newtable = table.updated(indexOflshFunction,table(indexOflshFunction).updated(locationOfVec, newElement :: oldList))

						computeHashValueToTable(tail, newtable)
					}
					case _ => println("Error");table
				}
			}

			computeHashValueToTable(reducedVectors.zipWithIndex, table)

		} //end of def tablingHas
	}

	def main(args: Array[String]){
		//val baseVectors = Source.fromFile("../dataset/sift/siftsmall/base.input").getLines.toList.map{_.split(" ").map{_.toDouble}.toVector}
		//val queryVectors = Source.fromFile("../dataset/sift/siftsmall/query.input").getLines.toList.map{_.split(" ").toVector}
		val learnVectors = Source.fromFile("../dataset/learn.input").getLines.toList.map{_.split(" ").map{_.toDouble}.toVector}

		val numberOfVectors = learnVectors.length
		val dimensionD = learnVectors(1).length

		//from the original paper & we let p2 = p(radiusR=1.0)
		//val dimensionK = computeK(computeFunctionP(binWidthW, radiusR), numberOfVectors)
		val dimensionK = 10 
		//val numberOflshL = computeLfromKP(dimensionK,successProbability)
		val numberOflshL = 5
		val tableSize = numberOfVectors


		//?The value will be minus because radius is less than bin width
		//println("p(R=1.0) = computeFunctionP(4.0,1.0) = " + computeFunctionP(4.0,1.0))


		//learnVectors: 25000 * 128 List

		
		val hashFn = new HashFunctions(numberOflshL, dimensionK, dimensionD, tableSize, binWidthW)
		
/*
		//改成讀一行即掉
		learnVectors.zipWithIndex.foreach( {case (vector, indexOfVector) => {

			val reducedVectors = hashFn.dotHash(vector)

			hashFn.tableingHash(reducedVectors, indexOfVector)
				
		}}) //end of foreach
*/
		def computeHashValue(learnVec: List[(Vector[Double], Int)],
		                     table: Vector[Vector[List[Any]]]): Vector[Vector[List[Any]]]= {
			learnVec match {
				case Nil => table
				case (vec,nameOfVec) :: tail => {
					val reducedVectors = hashFn.dotHash(vec)
					val newTable = hashFn.tableHash(reducedVectors, nameOfVec, table)
					computeHashValue(tail, newTable)
				}
				case _ => println("Error");table
			}
		}
		val emptyTable = 
		      Vector.fill(numberOflshL, tableSize)(Nil)
		val table = computeHashValue(learnVectors.zipWithIndex , emptyTable)

		println(table)

	} //end of def main()
}
