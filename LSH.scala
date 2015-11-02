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
	                    binWidthW: Double ){

		case class lshFunctionParameter(vectorA: Vector[Double], parameterB: Double)
		

		val functionG = List.fill(numberOfFunctionG)(
			List.fill(dimensionOfVectorG)(
				new lshFunctionParameter(Vector.fill(
															 dimensionOfVectorA)(Random.nextGaussian()),
														     Random.nextDouble() * binWidthW)
			)
		)

		def dotHash(vectorV: Vector[Double]) = {
			functionG.map( 
				_.map( x => {
					val dot = (vectorV zip x.vectorA).map{case (a,b) => a * b}.foldLeft(0.0)(_+_)
					math.floor((dot + x.parameterB)/binWidthW).toInt
				}).toVector
			)
			
		}

		val functionH1 = IndexedSeq.fill(dimensionOfVectorG)((Random.nextDouble()*prime).toLong)
		val functionH2 = IndexedSeq.fill(dimensionOfVectorG)((Random.nextDouble()*prime).toLong)

			def primaryHash(v: IndexedSeq[Int], tableSize: Int) = {
				val dot = (v zip functionH1).map{case (a,b) => a * b}.foldLeft(0L)(_+_)
				val hashValue: Int = ((dot % prime)%tableSize).toInt
				hashValue
			}

			//dot = the value of inner product
			//the return value is dot%prime (or dot%prime + prime when negative)
			//  where prime = (1L << 29) - 33
			def secondaryHash(v: IndexedSeq[Int]) = {
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

			def tableHash(vectorG: List[Vector[Int]], tableSize: Int) = {

				vectorG.foreach(x => {
					primaryHash(x, tableSize)
					secondaryHash(x)
				})
			}
	}

	def main(args: Array[String]){
		//val baseVectors = Source.fromFile("../dataset/sift/siftsmall/base.input").getLines.toList.map{_.split(" ").map{_.toDouble}.toVector}
		//val queryVectors = Source.fromFile("../dataset/sift/siftsmall/query.input").getLines.toList.map{_.split(" ").toVector}
		val learnVectors = Source.fromFile("../dataset/sift/siftsmall/learn.input").getLines.toList.map{_.split(" ").map{_.toDouble}.toVector}

		val numberOfVectors = learnVectors.length
		val dimensionD = learnVectors(1).length

		//from the original paper & we let p2 = p(radiusR=1.0)
		//val dimensionK = computeK(computeFunctionP(binWidthW, radiusR), numberOfVectors)
		val dimensionK = 20
		//val numberOflshL = computeLfromKP(dimensionK,successProbability)
		val numberOflshL = 10
		val tableSize = numberOfVectors


		//?The value will be minus because radius is less than bin width
		//println("p(R=1.0) = computeFunctionP(4.0,1.0) = " + computeFunctionP(4.0,1.0))


		//learnVectors: 25000 * 128 List

		
		val hashFn = new HashFunctions(numberOflshL, dimensionK, dimensionD, binWidthW)
		
		var a = 0

		//改成讀一行即掉
		for (vector <- learnVectors){

/*
			def lshFunctionH(): Int = {
				val vectorA = IndexedSeq.fill(dimensionD)(Random.nextGaussian())
				val parameterB = Random.nextDouble() * binWidthW
				val dot = (vector zip vectorA).map{case (a,b) => a * b}.foldLeft(0.0)(_+_)
				math.floor((dot + parameterB)/binWidthW).toInt
			}
			val vectorG = List.fill(dimensionK)(lshFunctionH)
*/
			val vectorG = hashFn.dotHash(vector)


			//if (a < 10){
				//hashFn.computePrimaryHash(hashFn.functionH1)
				//hashFn.computeSecondaryHash(hashFn.functionH2)

				hashFn.tableHash(vectorG, tableSize)
				//a += 1
			//}
				
		} //end of for( <- learnvector)


	} //end of def main()
}
