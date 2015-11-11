import scala.io.Source
import scala.util.Random
import scala.math


object LSH {

	val prime = (1L << 29) - 33

	//the case of w=4 and r=1 is not good
	val binWidthW = 400.0
	val radiusR = 100.0
	val successProbability = 0.9
	println("Bin width w = " + binWidthW)
	println("Success Probability = " + successProbability)
	println("Radius R = " + radiusR)

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
		val x = w / c
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

	def innderProduct(v1: Vector[Double], v2: Vector[Double]): Double = {
		var sum = 0.0
		val len = v1.length
		for(i <- 0 until len) {
			sum = sum + v1(i) * v2(i)
		}	
		sum
	}

	def innderProductInt(v1: Vector[Int], v2: Vector[Long]): Long = {
		var sum: Long = 0
		val len = v1.length
		for(i <- 0 until len) {
			sum = sum + v1(i) * v2(i)
		}	
		sum
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
					//!!! the functional coding style is bottleneck
					//val dot = (vectorV zip x.vectorA).map{case (a,b) => a * b}.foldLeft(0.0)(_+_)

					//MUCH faster than code above
					val dot = innderProduct(vectorV, x.vectorA)

					math.floor((dot + x.parameterB)/binWidthW).toInt
				}).toVector
			)
			
		}

		private val functionH1 = Vector.fill(dimensionOfVectorG)((Random.nextDouble()*prime).toLong)
		private val functionH2 = Vector.fill(dimensionOfVectorG)((Random.nextDouble()*prime).toLong)

		private def primaryHash(v: Vector[Int]) = {
			val dot = (v zip functionH1).map{case (a,b) => a * b}.foldLeft(0L)(_+_)
			//val dot = innderProductInt(v, functionH1)
			//test println(dot == dou)

			val x = dot % prime
			if (x < 0)
				((x + prime) % tableSize).toInt
			else 
				(x % tableSize).toInt
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

/*		//version of two dimension table (i.e. table: Vector[Vector[List]])

		//vectorG means g_i(v) for given v and for all i belonging to L
		def tableHash(reducedVectors: List[Vector[Int]], 
		                 nameOfVec: Int,
										 table: Vector[Vector[List[Any]]]): Vector[Vector[List[Any]]] = {

			//For a given original vector v, assume g_i is the i-th lsh function
			//v is reduced to reducedVectors(0) by g_0
			//v is reduced to reducedVectors(1) by g_1 ...

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
*/

		//vectorG means g_i(v) for given v and for all i belonging to L
		def tableHash(reducedVectors: List[Vector[Int]], 
		                nameOfVec: Int,
						table: Vector[List[Any]]): Vector[List[Any]] = {

			//For a given original vector v, assume g_i is the i-th lsh function
			//v is reduced to reducedVectors(0) by g_0
			//v is reduced to reducedVectors(1) by g_1 ...

			def computeHashValueToTable(reducedVecs: List[Vector[Int]],
								   table: Vector[List[Any]]): Vector[List[Any]] = {
				reducedVecs match {
					case Nil => table
					case head :: tail => {
						val locationOfVec= primaryHash(head)
						val fingerprintOfVec= secondaryHash(head)
						
						val oldList = table(locationOfVec) 
						val newElement = (fingerprintOfVec, nameOfVec)
						val newtable = table.updated(locationOfVec, newElement :: oldList)

						computeHashValueToTable(tail, newtable)
					}
					case _ => println("Error");table
				}
			}

			computeHashValueToTable(reducedVectors, table)

		} //end of def tablingHas


		def computeOneQuery (reducedQuery: List[Vector[Int]],
			           table: Vector[List[Any]]): Set[Int] = {

			def computeEachReducedQuery(queries: List[Vector[Int]]): Set[Int] ={
				queries match {
					case head :: tail => {

						val location = primaryHash(head)
						val fingerprintOfQuery = secondaryHash(head)

						def findFingerprint(pointList: List[Any]): Set[Int]= 
							pointList match {
								case (fingerprint, nameOfVec: Int) :: tail => {
									if (fingerprint == fingerprintOfQuery){
										findFingerprint(tail) + nameOfVec
									}else
										findFingerprint(tail)
									}
								case _ => Set.empty
							}
						
						computeEachReducedQuery(tail) ++ findFingerprint(table(location))
					}
					case _ => Set.empty
				}
			}

			computeEachReducedQuery(reducedQuery)

		} //end of def findQuery

	} //end of class HashFunctions

	def main(args: Array[String]){
		//val baseVectors = Source.fromFile("../dataset/sift/siftsmall/base.input").getLines.toList.map{_.split(" ").map{_.toDouble}.toVector}
		val queryVectors = Source.fromFile("../dataset/my_query.input").getLines.toList.map{_.split(" ").map{_.toDouble}.toVector}
		val learnVectors = Source.fromFile("../dataset/siftsmall_learn.input").getLines.toList.map{_.split(" ").map{_.toDouble}.toVector}

		val numberOfVectors = learnVectors.length
		val dimensionD = learnVectors(1).length
		println("Number of vectors from learn_file: " + numberOfVectors)
		println("Dimension of the vector from learn_file : " + dimensionD)

		//need to find a good way to dertermine parameter k
		//val dimensionK = computeK(computeFunctionP(binWidthW, radiusR), numberOfVectors)
		val dimensionK = 10 
		//val numberOflshL = computeLfromKP(dimensionK,successProbability)
		val numberOflshL = computeLfromKP(dimensionK,successProbability).min(20)
		println("Reduced dimension k = " + dimensionK)
		println("Number of tables L = " + numberOflshL)

		val tableSize = numberOfVectors

		//?The value will be minus because radius is less than bin width
		//println("p(R=1.0) = computeFunctionP(4.0,1.0) = " + computeFunctionP(4.0,1.0))

		
		val hashFn = new HashFunctions(numberOflshL, dimensionK, dimensionD, tableSize, binWidthW)
		
		def computeHashValue(learnVec: List[(Vector[Double], Int)],
		                     table: Vector[List[Any]]): Vector[List[Any]]= {
			learnVec match {
				case Nil => table
				case (vec,nameOfVec) :: tail => {
					//println("Vector " + nameOfVec + " of " + numberOfVectors)

					val reducedVectors = hashFn.dotHash(vec)
					//test val reducedVectors = List.fill(numberOflshL)(Vector.fill(dimensionK)(1))

					//a little slow
					val newTable = hashFn.tableHash(reducedVectors, nameOfVec, table)

					computeHashValue(tail, newTable)
				}
				case _ => println("Error");table
			}
		}
		val emptyTable = 
		      Vector.fill(tableSize)(Nil)
		val table = computeHashValue(learnVectors.zipWithIndex , emptyTable)

		println
		println("Process queries...")
		println("Number of vectors from query_file: " + queryVectors.length)
		println("Dimension of the vector from query_file : " + queryVectors(1).length)
		for(query <- queryVectors){
			print("Query " + (queryVectors.indexOf(query)+1) + " : ")
			val reducedQuery = hashFn.dotHash(query)	
			val resultOfQuery = hashFn.computeOneQuery(reducedQuery, table)
			println(resultOfQuery)
			
		}


	} //end of def main()
}
